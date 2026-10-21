"""PyTorch neural network and ordinary linear model training, plus self-contained ONNX export."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import onnx
import onnxruntime as ort
import torch
from onnx import TensorProto, helper
from torch import nn

from preface import __version__
from preface.hardware import get_onnx_providers


class PrefacePredictor(nn.Module):
    """Core PREFACE regression network (single-hidden-layer Sigmoid MLP or linear OLM)."""

    def __init__(self, n_feat: int, hidden: int = 2, is_olm: bool = False) -> None:
        super().__init__()
        self.n_feat = n_feat
        self.hidden = hidden
        self.is_olm = is_olm

        if is_olm or hidden <= 0:
            self.hidden_layer: nn.Linear | None = None
            self.activation: nn.Sigmoid | None = None
            self.output_layer = nn.Linear(n_feat, 1)
        else:
            self.hidden_layer = nn.Linear(n_feat, hidden)
            self.activation = nn.Sigmoid()
            self.output_layer = nn.Linear(hidden, 1)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        if self.hidden_layer is not None and self.activation is not None:
            x = self.activation(self.hidden_layer(x))
        return self.output_layer(x)


class PrefaceFullPipeline(nn.Module):
    """End-to-end ONNX-exportable PREFACE inference module."""

    def __init__(
        self,
        predictor: PrefacePredictor,
        mean_features: np.ndarray,
        pca_center: np.ndarray,
        pca_rotation: np.ndarray,
        the_intercept: float,
        the_slope: float,
    ) -> None:
        super().__init__()
        self.predictor = predictor
        self.register_buffer(
            "mean_features",
            torch.as_tensor(mean_features, dtype=torch.float32).unsqueeze(0),
        )
        self.register_buffer(
            "pca_center",
            torch.as_tensor(pca_center, dtype=torch.float32).unsqueeze(0),
        )
        self.register_buffer(
            "pca_rotation",
            torch.as_tensor(pca_rotation, dtype=torch.float32),
        )
        self.register_buffer(
            "the_intercept",
            torch.tensor(float(the_intercept), dtype=torch.float32),
        )
        self.register_buffer(
            "the_slope",
            torch.tensor(float(the_slope), dtype=torch.float32),
        )

    def forward(self, ratio: torch.Tensor) -> torch.Tensor:
        imputed = torch.where(torch.isnan(ratio), self.mean_features, ratio)
        projected = (imputed - self.pca_center) @ self.pca_rotation
        raw_pred = self.predictor(projected).squeeze(-1)
        return self.the_intercept + self.the_slope * raw_pred


@dataclass(frozen=True)
class ONNXModelMetadata:
    """Metadata embedded in model.onnx metadata_props."""

    possible_features: list[str]
    the_intercept_x: float
    the_slope_x: float
    n_feat: int
    is_olm: bool
    version: str


@dataclass(frozen=True)
class FFYModelMetadata:
    """Metadata embedded in a saved FFY ONNX model."""

    bin_loci: list[str]
    informative_indices: list[int]
    female_baseline: float
    male_100_baseline: float
    yfrac_threshold: float
    binsize: int
    version: str


def train_olm(
    x_train: np.ndarray,
    y_train: np.ndarray,
    device: torch.device,
) -> PrefacePredictor:
    """Fit unregularized Gaussian linear model matching glmnet(..., family='gaussian', lambda=0)."""
    n_samples, n_feat = x_train.shape
    design = np.column_stack([np.ones(n_samples, dtype=np.float64), x_train])
    coeffs, _, _, _ = np.linalg.lstsq(design, y_train.ravel(), rcond=None)

    model = PrefacePredictor(n_feat=n_feat, hidden=0, is_olm=True)
    with torch.no_grad():
        model.output_layer.bias.fill_(float(coeffs[0]))
        model.output_layer.weight.copy_(
            torch.as_tensor(coeffs[1:], dtype=torch.float32).unsqueeze(0)
        )
    return model.to(device)


def train_neural_net(
    x_train: np.ndarray,
    y_train: np.ndarray,
    hidden: int,
    device: torch.device,
    stepmax: int = 100_000,
    threshold: float = 0.01,
) -> PrefacePredictor:
    """Train single-hidden-layer Sigmoid neural network using full-batch Rprop matching R's neuralnet."""
    n_samples, n_feat = x_train.shape
    train_device = (
        torch.device("cpu")
        if (device.type == "mps" and n_samples * n_feat < 100_000)
        else device
    )

    model = PrefacePredictor(n_feat=n_feat, hidden=hidden, is_olm=False).to(
        train_device
    )
    for param in model.parameters():
        nn.init.normal_(param, mean=0.0, std=1.0)

    x_t = torch.as_tensor(x_train, dtype=torch.float32, device=train_device)
    y_t = torch.as_tensor(y_train, dtype=torch.float32, device=train_device).reshape(
        n_samples, 1
    )

    optimizer = torch.optim.Rprop(
        model.parameters(),
        lr=0.01,
        etas=(0.5, 1.2),
        step_sizes=(1e-6, 50.0),
    )

    converged = False
    params = list(model.parameters())

    for _ in range(stepmax):
        optimizer.zero_grad(set_to_none=True)
        pred = model(x_t)
        loss = 0.5 * torch.sum((pred - y_t) ** 2)
        loss.backward()

        if all(
            p.grad is not None
            and bool((p.grad.detach().abs().max() < threshold).item())
            for p in params
        ):
            converged = True
            break

        optimizer.step()

    if not converged:
        raise RuntimeError(
            "Neural network did not converge. Re-run and decrease --hidden or optimize --nfeat. "
            "Alternatively, use --olm."
        )

    return model.to(device)


def fit_predictor(
    x_train: np.ndarray,
    y_train: np.ndarray,
    hidden: int,
    is_olm: bool,
    device: torch.device,
    log_prefix: str = "",
) -> PrefacePredictor:
    """Dispatch training to either OLM or Rprop neural network."""
    if is_olm:
        print(f"{log_prefix}Training ordinary linear model ...")
        return train_olm(x_train, y_train, device=device)
    print(f"{log_prefix}Training neural network ...")
    return train_neural_net(x_train, y_train, hidden=hidden, device=device)


@torch.no_grad()
def predict_model(
    model: PrefacePredictor,
    x_eval: np.ndarray,
    device: torch.device,
) -> np.ndarray:
    """Evaluate a trained PrefacePredictor on PCA-projected features."""
    model.eval()
    eval_device = next(model.parameters()).device
    x_t = torch.as_tensor(x_eval, dtype=torch.float32, device=eval_device)
    return model(x_t).squeeze(-1).cpu().numpy().astype(np.float64)


def _save_onnx_with_metadata(
    onnx_model: onnx.ModelProto,
    meta_entries: dict[str, str],
    onnx_path: Path,
) -> None:
    for key, val in meta_entries.items():
        entry = onnx_model.metadata_props.add()
        entry.key = key
        entry.value = val
    onnx.checker.check_model(onnx_model)
    onnx_path.parent.mkdir(parents=True, exist_ok=True)
    onnx.save(onnx_model, str(onnx_path))


def _read_onnx_props(onnx_path: Path, required_keys: set[str]) -> dict[str, str]:
    onnx_model = onnx.load(str(onnx_path), load_external_data=False)
    props = {entry.key: entry.value for entry in onnx_model.metadata_props}
    missing = required_keys - props.keys()
    if missing:
        raise ValueError(
            f"ONNX model '{onnx_path}' is missing required metadata: {missing}"
        )
    return props


def _run_single_onnx_tensor(
    onnx_path: Path,
    input_name: str,
    output_name: str,
    values: np.ndarray,
    requested_device: str = "auto",
) -> float:
    session = ort.InferenceSession(
        str(onnx_path), providers=get_onnx_providers(requested_device)
    )
    input_tensor = np.asarray(values, dtype=np.float32).reshape(1, -1)
    outputs = session.run([output_name], {input_name: input_tensor})
    return float(np.asarray(outputs[0]).ravel()[0])


def export_pipeline_to_onnx(
    predictor: PrefacePredictor,
    mean_features: np.ndarray,
    pca_center: np.ndarray,
    pca_rotation: np.ndarray,
    the_intercept: float,
    the_slope: float,
    the_intercept_x: float,
    the_slope_x: float,
    possible_features: list[str],
    onnx_path: Path,
) -> None:
    """Export the full PREFACE preprocessing + model + skew-correction pipeline to a single ONNX file."""
    cpu_predictor = predictor.cpu().eval()
    pipeline = PrefaceFullPipeline(
        predictor=cpu_predictor,
        mean_features=mean_features,
        pca_center=pca_center,
        pca_rotation=pca_rotation,
        the_intercept=the_intercept,
        the_slope=the_slope,
    ).eval()

    dummy_input = torch.zeros((1, len(possible_features)), dtype=torch.float32)
    batch_dim = torch.export.Dim("batch_size", min=1)

    onnx_program = torch.onnx.export(
        pipeline,
        (dummy_input,),
        input_names=["ratio"],
        output_names=["preface_prediction"],
        dynamic_shapes={"ratio": {0: batch_dim}},
        opset_version=18,
        dynamo=True,
    )
    if onnx_program is None:
        raise RuntimeError("ONNX export failed to produce an ONNXProgram.")

    _save_onnx_with_metadata(
        onnx_program.model_proto,
        {
            "possible_features": json.dumps(possible_features),
            "the_intercept_X": repr(float(the_intercept_x)),
            "the_slope_X": repr(float(the_slope_x)),
            "n_feat": str(cpu_predictor.n_feat),
            "is_olm": "1" if cpu_predictor.is_olm else "0",
            "version": __version__,
        },
        onnx_path,
    )


def load_onnx_metadata(onnx_path: Path) -> ONNXModelMetadata:
    """Read PREFACE metadata properties from a saved model.onnx file."""
    props = _read_onnx_props(
        onnx_path,
        {"possible_features", "the_intercept_X", "the_slope_X", "n_feat", "is_olm"},
    )
    return ONNXModelMetadata(
        possible_features=json.loads(props["possible_features"]),
        the_intercept_x=float(props["the_intercept_X"]),
        the_slope_x=float(props["the_slope_X"]),
        n_feat=int(props["n_feat"]),
        is_olm=props["is_olm"] == "1",
        version=props.get("version", __version__),
    )


def run_onnx_inference(
    onnx_path: Path,
    feature_ratios: np.ndarray,
    requested_device: str = "auto",
) -> float:
    """Execute ONNX Runtime inference on raw feature ratios (including NaNs)."""
    return _run_single_onnx_tensor(
        onnx_path, "ratio", "preface_prediction", feature_ratios, requested_device
    )


def export_ffy_model_to_onnx(
    intercept: float,
    weights: np.ndarray,
    bin_loci: list[str],
    informative_indices: list[int],
    female_baseline: float,
    male_100_baseline: float,
    yfrac_threshold: float,
    binsize: int,
    onnx_path: Path,
) -> None:
    """Export a fitted FFY Lasso model and calibration baselines to a self-contained ONNX file."""
    n_bins = len(bin_loci)
    w_arr = np.asarray(weights, dtype=np.float32).reshape(n_bins, 1)

    graph = helper.make_graph(
        [
            helper.make_node("MatMul", ["y_norm", "lasso_weights"], ["matmul_out"]),
            helper.make_node("Add", ["matmul_out", "lasso_intercept"], ["raw_out"]),
            helper.make_node("Squeeze", ["raw_out", "squeeze_axes"], ["squeezed_out"]),
            helper.make_node(
                "Clip", ["squeezed_out", "clip_min", "clip_max"], ["ffy_prediction"]
            ),
        ],
        "PREFACE_FFY_Graph",
        [
            helper.make_tensor_value_info(
                "y_norm", TensorProto.FLOAT, ["batch_size", n_bins]
            )
        ],
        [
            helper.make_tensor_value_info(
                "ffy_prediction", TensorProto.FLOAT, ["batch_size"]
            )
        ],
        [
            helper.make_tensor(
                "lasso_weights",
                TensorProto.FLOAT,
                [n_bins, 1],
                w_arr.flatten().tolist(),
            ),
            helper.make_tensor(
                "lasso_intercept", TensorProto.FLOAT, [1], [float(intercept)]
            ),
            helper.make_tensor("clip_min", TensorProto.FLOAT, [], [0.0]),
            helper.make_tensor("clip_max", TensorProto.FLOAT, [], [100.0]),
            helper.make_tensor("squeeze_axes", TensorProto.INT64, [1], [1]),
        ],
    )
    onnx_model = helper.make_model(
        graph, producer_name="preface", opset_imports=[helper.make_opsetid("", 18)]
    )
    onnx_model.ir_version = 9

    _save_onnx_with_metadata(
        onnx_model,
        {
            "model_type": "PREFACE_FFY",
            "bin_loci": json.dumps(bin_loci),
            "informative_indices": json.dumps(informative_indices),
            "female_baseline": repr(float(female_baseline)),
            "male_100_baseline": repr(float(male_100_baseline)),
            "yfrac_threshold": repr(float(yfrac_threshold)),
            "binsize": str(binsize),
            "version": __version__,
        },
        onnx_path,
    )


def load_ffy_onnx_metadata(onnx_path: Path) -> FFYModelMetadata:
    """Read FFY model metadata properties from a saved ONNX file."""
    props = _read_onnx_props(
        onnx_path,
        {
            "bin_loci",
            "informative_indices",
            "female_baseline",
            "male_100_baseline",
            "yfrac_threshold",
            "binsize",
        },
    )
    return FFYModelMetadata(
        bin_loci=json.loads(props["bin_loci"]),
        informative_indices=json.loads(props["informative_indices"]),
        female_baseline=float(props["female_baseline"]),
        male_100_baseline=float(props["male_100_baseline"]),
        yfrac_threshold=float(props["yfrac_threshold"]),
        binsize=int(props["binsize"]),
        version=props.get("version", __version__),
    )


def run_ffy_onnx_inference(
    onnx_path: Path,
    aligned_y: np.ndarray,
    requested_device: str = "auto",
) -> float:
    """Execute ONNX Runtime inference on normalized chromosome Y bin values."""
    pred = _run_single_onnx_tensor(
        onnx_path, "y_norm", "ffy_prediction", aligned_y, requested_device
    )
    return float(np.clip(pred, 0.0, 100.0))
