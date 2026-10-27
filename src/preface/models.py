"""PyTorch neural network and ordinary linear model training, plus self-contained ONNX export."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import onnx
import onnxruntime as ort
import torch
from torch import nn

from preface import __version__
from preface.hardware import get_onnx_providers
from preface.logging_utils import get_logger
from preface.math_utils import format_percentage

logger = get_logger(__name__)

MODEL_TYPE_ALL = "PREFACE_ALL"
MODEL_TYPE_FFY = "PREFACE_FFY"
MODEL_TYPE_PREFACE = "PREFACE_ONLY"


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
    """End-to-end ONNX-exportable PREFACE + FFX inference module."""

    def __init__(
        self,
        predictor: PrefacePredictor,
        mean_features: np.ndarray,
        pca_center: np.ndarray,
        pca_rotation: np.ndarray,
        the_intercept: float,
        the_slope: float,
        the_intercept_x: float = 0.0,
        the_slope_x: float = 1.0,
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
        self.register_buffer(
            "the_intercept_x",
            torch.tensor(float(the_intercept_x), dtype=torch.float32),
        )
        safe_slope_x = float(the_slope_x) if abs(float(the_slope_x)) > 1e-12 else 1.0
        self.register_buffer(
            "the_slope_x",
            torch.tensor(safe_slope_x, dtype=torch.float32),
        )

    def forward(
        self, ratio: torch.Tensor, x_ratio: torch.Tensor
    ) -> tuple[torch.Tensor, torch.Tensor]:
        imputed = torch.where(torch.isnan(ratio), self.mean_features, ratio)
        projected = (imputed - self.pca_center) @ self.pca_rotation
        raw_pred = self.predictor(projected).squeeze(-1)
        preface_pred = self.the_intercept + self.the_slope * raw_pred
        ffx_pred = (x_ratio - self.the_intercept_x) / self.the_slope_x
        return ffx_pred, preface_pred


class FFYPipeline(nn.Module):
    """End-to-end ONNX-exportable FFY sex classification and Lasso regression module."""

    def __init__(
        self,
        ffy_intercept: float,
        ffy_weights: np.ndarray,
        informative_indices: list[int],
        yfrac_threshold: float,
    ) -> None:
        super().__init__()
        n_y_bins = int(np.asarray(ffy_weights).size)
        info_idx = informative_indices if informative_indices else list(range(n_y_bins))
        score_mask = np.zeros((n_y_bins, 1), dtype=np.float32)
        score_mask[info_idx, 0] = 1.0 / max(1, len(info_idx))

        self.register_buffer(
            "ffy_weights",
            torch.as_tensor(ffy_weights, dtype=torch.float32).reshape(n_y_bins, 1),
        )
        self.register_buffer(
            "ffy_intercept",
            torch.tensor(float(ffy_intercept), dtype=torch.float32),
        )
        self.register_buffer(
            "score_weights",
            torch.as_tensor(score_mask, dtype=torch.float32),
        )
        self.register_buffer(
            "yfrac_threshold",
            torch.tensor(float(yfrac_threshold), dtype=torch.float32),
        )

    def forward(
        self, y_norm: torch.Tensor, no_female_gate: torch.Tensor
    ) -> tuple[torch.Tensor, torch.Tensor]:
        y_clean = torch.where(torch.isnan(y_norm), torch.zeros_like(y_norm), y_norm)
        y_score = (y_clean @ self.score_weights).squeeze(-1)
        is_male_bool = y_score >= self.yfrac_threshold
        sex_is_male = is_male_bool.to(torch.int64)

        raw_ffy = (y_clean @ self.ffy_weights).squeeze(-1) + self.ffy_intercept
        clipped_ffy = torch.clamp(raw_ffy, min=0.0, max=1.0)
        keep_ffy = is_male_bool | no_female_gate
        ffy_pred = torch.where(keep_ffy, clipped_ffy, torch.zeros_like(clipped_ffy))
        return sex_is_male, ffy_pred


class UnifiedPrefacePipeline(nn.Module):
    """Unified ONNX-exportable FFY + FFX + PREFACE inference module."""

    def __init__(
        self,
        predictor: PrefacePredictor,
        mean_features: np.ndarray,
        pca_center: np.ndarray,
        pca_rotation: np.ndarray,
        the_intercept: float,
        the_slope: float,
        the_intercept_x: float,
        the_slope_x: float,
        ffy_intercept: float,
        ffy_weights: np.ndarray,
        informative_indices: list[int],
        yfrac_threshold: float,
    ) -> None:
        super().__init__()
        self.ffy_pipe = FFYPipeline(
            ffy_intercept=ffy_intercept,
            ffy_weights=ffy_weights,
            informative_indices=informative_indices,
            yfrac_threshold=yfrac_threshold,
        )
        self.preface_pipe = PrefaceFullPipeline(
            predictor=predictor,
            mean_features=mean_features,
            pca_center=pca_center,
            pca_rotation=pca_rotation,
            the_intercept=the_intercept,
            the_slope=the_slope,
            the_intercept_x=the_intercept_x,
            the_slope_x=the_slope_x,
        )

    def forward(
        self,
        ratio: torch.Tensor,
        x_ratio: torch.Tensor,
        y_norm: torch.Tensor,
        no_female_gate: torch.Tensor,
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
        sex_is_male, ffy_pred = self.ffy_pipe(y_norm, no_female_gate)
        ffx_pred, preface_pred = self.preface_pipe(ratio, x_ratio)
        return sex_is_male, ffy_pred, ffx_pred, preface_pred


@dataclass(frozen=True)
class ONNXModelMetadata:
    """Metadata embedded in a PREFACE or unified ONNX model."""

    model_type: str
    possible_features: list[str]
    the_intercept_x: float
    the_slope_x: float
    n_feat: int
    is_olm: bool
    binsize: int
    ffy_binsize: int
    bin_loci: list[str]
    informative_indices: list[int]
    female_baseline: float
    male_100_baseline: float
    yfrac_threshold: float
    no_female_gate: bool
    version: str


@dataclass(frozen=True)
class UnifiedPredictionResult:
    """Result from running inference on an integrated, FFY-only, or PREFACE-only ONNX model."""

    sex: str | None = None
    ffy: float | None = None
    ffx: float | None = None
    preface: float | None = None

    def to_json_dict(self) -> dict[str, str | float]:
        out: dict[str, str | float] = {}
        if self.sex is not None:
            out["sex"] = self.sex
        if self.ffy is not None:
            out["FFY"] = self.ffy
        if self.ffx is not None:
            out["FFX"] = self.ffx
        if self.preface is not None:
            out["PREFACE"] = self.preface
        return out

    def to_text_lines(self) -> list[str]:
        lines: list[str] = []
        if self.sex is not None:
            lines.append(f"sex = {self.sex}")
        if self.ffy is not None:
            lines.append(f"FFY = {format_percentage(self.ffy)}")
        if self.ffx is not None:
            lines.append(f"FFX = {format_percentage(self.ffx)}")
        if self.preface is not None:
            lines.append(f"PREFACE = {format_percentage(self.preface)}")
        return lines


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
    threshold: float = 1e-4,
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
        logger.info("%sTraining ordinary linear model ...", log_prefix)
        return train_olm(x_train, y_train, device=device)
    logger.info("%sTraining neural network ...", log_prefix)
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
    logger.info("Exported ONNX model to '%s'", onnx_path)


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
    binsize: int = 100_000,
) -> None:
    """Export the PREFACE + FFX preprocessing, model, and skew-correction pipeline to ONNX."""
    cpu_predictor = predictor.cpu().eval()
    pipeline = PrefaceFullPipeline(
        predictor=cpu_predictor,
        mean_features=mean_features,
        pca_center=pca_center,
        pca_rotation=pca_rotation,
        the_intercept=the_intercept,
        the_slope=the_slope,
        the_intercept_x=the_intercept_x,
        the_slope_x=the_slope_x,
    ).eval()

    dummy_ratio = torch.zeros((1, len(possible_features)), dtype=torch.float32)
    dummy_x_ratio = torch.zeros((1,), dtype=torch.float32)
    batch_r = torch.export.Dim("batch_size", min=1)
    batch_x = torch.export.Dim("batch_size_x", min=1)

    onnx_program = torch.onnx.export(
        pipeline,
        (dummy_ratio, dummy_x_ratio),
        input_names=["ratio", "x_ratio"],
        output_names=["ffx_prediction", "preface_prediction"],
        dynamic_shapes={"ratio": {0: batch_r}, "x_ratio": {0: batch_x}},
        opset_version=18,
        dynamo=True,
    )
    if onnx_program is None:
        raise RuntimeError("ONNX export failed to produce an ONNXProgram.")

    _save_onnx_with_metadata(
        onnx_program.model_proto,
        {
            "model_type": MODEL_TYPE_PREFACE,
            "possible_features": json.dumps(possible_features),
            "the_intercept_X": repr(float(the_intercept_x)),
            "the_slope_X": repr(float(the_slope_x)),
            "n_feat": str(cpu_predictor.n_feat),
            "is_olm": "1" if cpu_predictor.is_olm else "0",
            "binsize": str(binsize),
            "version": __version__,
        },
        onnx_path,
    )


def export_unified_pipeline_to_onnx(
    predictor: PrefacePredictor,
    mean_features: np.ndarray,
    pca_center: np.ndarray,
    pca_rotation: np.ndarray,
    the_intercept: float,
    the_slope: float,
    the_intercept_x: float,
    the_slope_x: float,
    possible_features: list[str],
    ffy_intercept: float,
    ffy_weights: np.ndarray,
    bin_loci: list[str],
    informative_indices: list[int],
    female_baseline: float,
    male_100_baseline: float,
    yfrac_threshold: float,
    onnx_path: Path,
    binsize: int = 100_000,
    ffy_binsize: int = 5_000,
    no_female_gate: bool = False,
) -> None:
    """Export the integrated FFY + FFX + PREFACE pipeline to a single ONNX file."""
    cpu_predictor = predictor.cpu().eval()
    pipeline = UnifiedPrefacePipeline(
        predictor=cpu_predictor,
        mean_features=mean_features,
        pca_center=pca_center,
        pca_rotation=pca_rotation,
        the_intercept=the_intercept,
        the_slope=the_slope,
        the_intercept_x=the_intercept_x,
        the_slope_x=the_slope_x,
        ffy_intercept=ffy_intercept,
        ffy_weights=ffy_weights,
        informative_indices=informative_indices,
        yfrac_threshold=yfrac_threshold,
    ).eval()

    dummy_ratio = torch.zeros((1, len(possible_features)), dtype=torch.float32)
    dummy_x_ratio = torch.zeros((1,), dtype=torch.float32)
    dummy_y_norm = torch.zeros((1, len(bin_loci)), dtype=torch.float32)
    dummy_no_gate = torch.tensor([bool(no_female_gate)], dtype=torch.bool)
    batch_r = torch.export.Dim("batch_size", min=1)
    batch_x = torch.export.Dim("batch_size_x", min=1)
    batch_y = torch.export.Dim("batch_size_y", min=1)

    onnx_program = torch.onnx.export(
        pipeline,
        (dummy_ratio, dummy_x_ratio, dummy_y_norm, dummy_no_gate),
        input_names=["ratio", "x_ratio", "y_norm", "no_female_gate"],
        output_names=[
            "sex_is_male",
            "ffy_prediction",
            "ffx_prediction",
            "preface_prediction",
        ],
        dynamic_shapes={
            "ratio": {0: batch_r},
            "x_ratio": {0: batch_x},
            "y_norm": {0: batch_y},
            "no_female_gate": {},
        },
        opset_version=18,
        dynamo=True,
    )
    if onnx_program is None:
        raise RuntimeError("ONNX export failed to produce an ONNXProgram.")

    _save_onnx_with_metadata(
        onnx_program.model_proto,
        {
            "model_type": MODEL_TYPE_ALL,
            "possible_features": json.dumps(possible_features),
            "the_intercept_X": repr(float(the_intercept_x)),
            "the_slope_X": repr(float(the_slope_x)),
            "n_feat": str(cpu_predictor.n_feat),
            "is_olm": "1" if cpu_predictor.is_olm else "0",
            "binsize": str(binsize),
            "ffy_binsize": str(ffy_binsize),
            "bin_loci": json.dumps(bin_loci),
            "informative_indices": json.dumps(informative_indices),
            "female_baseline": repr(float(female_baseline)),
            "male_100_baseline": repr(float(male_100_baseline)),
            "yfrac_threshold": repr(float(yfrac_threshold)),
            "no_female_gate": "1" if no_female_gate else "0",
            "version": __version__,
        },
        onnx_path,
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
    no_female_gate: bool = False,
) -> None:
    """Export a fitted FFY Lasso model, sex classifier, and baselines to a self-contained ONNX file."""
    pipeline = FFYPipeline(
        ffy_intercept=intercept,
        ffy_weights=weights,
        informative_indices=informative_indices,
        yfrac_threshold=yfrac_threshold,
    ).eval()

    dummy_y_norm = torch.zeros((1, len(bin_loci)), dtype=torch.float32)
    dummy_no_gate = torch.tensor([bool(no_female_gate)], dtype=torch.bool)
    batch_y = torch.export.Dim("batch_size", min=1)

    onnx_program = torch.onnx.export(
        pipeline,
        (dummy_y_norm, dummy_no_gate),
        input_names=["y_norm", "no_female_gate"],
        output_names=["sex_is_male", "ffy_prediction"],
        dynamic_shapes={"y_norm": {0: batch_y}, "no_female_gate": {}},
        opset_version=18,
        dynamo=True,
    )
    if onnx_program is None:
        raise RuntimeError("ONNX export failed to produce an ONNXProgram.")

    _save_onnx_with_metadata(
        onnx_program.model_proto,
        {
            "model_type": MODEL_TYPE_FFY,
            "bin_loci": json.dumps(bin_loci),
            "informative_indices": json.dumps(informative_indices),
            "female_baseline": repr(float(female_baseline)),
            "male_100_baseline": repr(float(male_100_baseline)),
            "yfrac_threshold": repr(float(yfrac_threshold)),
            "binsize": str(binsize),
            "ffy_binsize": str(binsize),
            "no_female_gate": "1" if no_female_gate else "0",
            "version": __version__,
        },
        onnx_path,
    )


def load_onnx_metadata(onnx_path: Path) -> ONNXModelMetadata:
    """Read metadata properties from any PREFACE ONNX model (integrated, FFY-only, or PREFACE-only)."""
    onnx_model = onnx.load(str(onnx_path), load_external_data=False)
    props = {entry.key: entry.value for entry in onnx_model.metadata_props}

    model_type = props.get("model_type", "")
    if not model_type:
        if "possible_features" in props and "bin_loci" in props:
            model_type = MODEL_TYPE_ALL
        elif "possible_features" in props:
            model_type = MODEL_TYPE_PREFACE
        elif "bin_loci" in props:
            model_type = MODEL_TYPE_FFY
        else:
            raise ValueError(
                f"ONNX model '{onnx_path}' is missing required PREFACE metadata."
            )

    if model_type in {MODEL_TYPE_ALL, MODEL_TYPE_PREFACE}:
        req = {
            "possible_features",
            "the_intercept_X",
            "the_slope_X",
            "n_feat",
            "is_olm",
        }
        if missing := req - props.keys():
            raise ValueError(
                f"ONNX model '{onnx_path}' is missing required metadata: {missing}"
            )
    if model_type in {MODEL_TYPE_ALL, MODEL_TYPE_FFY}:
        req_ffy = {
            "bin_loci",
            "informative_indices",
            "female_baseline",
            "male_100_baseline",
            "yfrac_threshold",
        }
        if missing_ffy := req_ffy - props.keys():
            raise ValueError(
                f"ONNX model '{onnx_path}' is missing required FFY metadata: {missing_ffy}"
            )

    binsize = int(props.get("binsize", "100000"))
    ffy_binsize = int(
        props.get(
            "ffy_binsize",
            props.get("binsize", "5000") if model_type == MODEL_TYPE_FFY else "5000",
        )
    )

    return ONNXModelMetadata(
        model_type=model_type,
        possible_features=json.loads(props["possible_features"])
        if "possible_features" in props
        else [],
        the_intercept_x=float(props.get("the_intercept_X", "0.0")),
        the_slope_x=float(props.get("the_slope_X", "1.0")),
        n_feat=int(props.get("n_feat", "0")),
        is_olm=props.get("is_olm", "0") == "1",
        binsize=binsize,
        ffy_binsize=ffy_binsize,
        bin_loci=json.loads(props["bin_loci"]) if "bin_loci" in props else [],
        informative_indices=json.loads(props["informative_indices"])
        if "informative_indices" in props
        else [],
        female_baseline=float(props.get("female_baseline", "0.0")),
        male_100_baseline=float(props.get("male_100_baseline", "0.5")),
        yfrac_threshold=float(props.get("yfrac_threshold", "0.01")),
        no_female_gate=props.get("no_female_gate", "0") == "1",
        version=props.get("version", __version__),
    )


def run_unified_onnx_inference(
    onnx_path: Path,
    meta: ONNXModelMetadata,
    feature_ratios: np.ndarray | None = None,
    x_ratio: float | None = None,
    aligned_y: np.ndarray | None = None,
    requested_device: str = "auto",
    no_female_gate_override: bool | None = None,
) -> UnifiedPredictionResult:
    """Execute a single unified ONNX Runtime inference call for any PREFACE model type."""
    session = ort.InferenceSession(
        str(onnx_path), providers=get_onnx_providers(requested_device)
    )
    input_names = {inp.name for inp in session.get_inputs()}
    output_names = [out.name for out in session.get_outputs()]

    feed: dict[str, np.ndarray] = {}
    if "ratio" in input_names:
        if feature_ratios is None:
            raise ValueError("Model requires autosomal feature ratios ('ratio').")
        feed["ratio"] = np.asarray(feature_ratios, dtype=np.float32).reshape(1, -1)
    if "x_ratio" in input_names:
        if x_ratio is None:
            raise ValueError("Model requires chrX ratio ('x_ratio').")
        feed["x_ratio"] = np.array([float(x_ratio)], dtype=np.float32)
    if "y_norm" in input_names:
        if aligned_y is None:
            raise ValueError("Model requires normalized chrY bins ('y_norm').")
        feed["y_norm"] = np.asarray(aligned_y, dtype=np.float32).reshape(1, -1)
    if "no_female_gate" in input_names:
        use_no_gate = (
            bool(no_female_gate_override)
            if no_female_gate_override is not None
            else bool(meta.no_female_gate)
        )
        feed["no_female_gate"] = np.array([use_no_gate], dtype=bool)

    raw_outputs = session.run(output_names, feed)
    out_map = dict(zip(output_names, raw_outputs))

    sex = (
        ("M" if int(np.asarray(out_map["sex_is_male"]).ravel()[0]) == 1 else "F")
        if "sex_is_male" in out_map
        else None
    )
    ffy = (
        float(
            np.clip(float(np.asarray(out_map["ffy_prediction"]).ravel()[0]), 0.0, 1.0)
        )
        if "ffy_prediction" in out_map
        else None
    )
    ffx = (
        float(np.asarray(out_map["ffx_prediction"]).ravel()[0])
        if "ffx_prediction" in out_map
        else None
    )
    preface = (
        float(np.asarray(out_map["preface_prediction"]).ravel()[0])
        if "preface_prediction" in out_map
        else None
    )
    return UnifiedPredictionResult(sex=sex, ffy=ffy, ffx=ffx, preface=preface)


def run_onnx_inference(
    onnx_path: Path,
    feature_ratios: np.ndarray,
    x_ratio: float = 1.0,
    requested_device: str = "auto",
) -> float:
    """Execute ONNX Runtime inference on raw feature ratios (including NaNs) and return PREFACE prediction."""
    meta = load_onnx_metadata(onnx_path)
    dummy_y = (
        np.zeros(max(1, len(meta.bin_loci)), dtype=np.float64)
        if meta.bin_loci
        else None
    )
    res = run_unified_onnx_inference(
        onnx_path,
        meta,
        feature_ratios=feature_ratios,
        x_ratio=x_ratio,
        aligned_y=dummy_y,
        requested_device=requested_device,
    )
    assert res.preface is not None
    return res.preface
