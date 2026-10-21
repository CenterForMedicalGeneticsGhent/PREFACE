"""End-to-end training and prediction pipelines for PREFACE."""

from __future__ import annotations

import time
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import torch

from preface.hardware import configure_hardware, resolve_device
from preface.io import (
    emit_prediction_output,
    load_config,
    load_inference_bed,
    load_training_dataset,
)
from preface.math_utils import ols_fit, pca_irlba_torch, rlm_huber
from preface.models import (
    export_pipeline_to_onnx,
    fit_predictor,
    load_onnx_metadata,
    predict_model,
    run_onnx_inference,
)

if TYPE_CHECKING:
    from preface.plotting import PerformanceMetrics


def _format_stat(val: float) -> str:
    if not np.isfinite(val):
        return "NA"
    return f"{val:.15g}"


def _run_cross_validation(
    training_frame: np.ndarray,
    ff_all: np.ndarray,
    labeled_mask: np.ndarray,
    labeled_indices: np.ndarray,
    n_feat: int,
    hidden: int,
    olm: bool,
    device: torch.device,
    repeats_dir: Path,
    repeats: int = 10,
) -> np.ndarray:
    """Execute k-fold cross-validation across labeled samples and save per-repeat diagnostic plots."""
    from preface.plotting import plot_performance

    n_total = training_frame.shape[0]
    n_labeled = labeled_indices.size
    fold_predictions: list[np.ndarray] = []

    for rep in range(1, repeats + 1):
        print(f"Model training | Repeat {rep}/{repeats} ...")

        start_idx = (rep - 1) * n_labeled // repeats
        end_idx = rep * n_labeled // repeats
        test_index_overall = labeled_indices[start_idx:end_idx]

        is_test_overall = np.zeros(n_total, dtype=bool)
        is_test_overall[test_index_overall] = True
        non_test_indices = np.flatnonzero(~is_test_overall)
        train_index_subset = np.flatnonzero(labeled_mask[non_test_indices])

        print("\tExecuting principal component analysis ...")
        pca_fold = pca_irlba_torch(
            training_frame[non_test_indices, :],
            n_components=min(n_feat * 10, non_test_indices.size - 1),
            device=device,
        )

        x_train = pca_fold.x[train_index_subset, :n_feat]
        y_train = ff_all[non_test_indices[train_index_subset]]
        x_test = pca_fold.project(training_frame[test_index_overall, :], n_feat)
        y_test = ff_all[test_index_overall]

        fold_model = fit_predictor(
            x_train,
            y_train,
            hidden=hidden,
            is_olm=olm,
            device=device,
            log_prefix="\t",
        )
        pred_fold = predict_model(fold_model, x_test, device=device)
        plot_performance(
            pred_fold,
            y_test,
            pca_fold.importance_prop,
            n_feat,
            "PREFACE (%)",
            "FF (%)",
            repeats_dir / f"repeat_{rep}.png",
        )
        fold_predictions.append(pred_fold)

    return np.concatenate(fold_predictions)


def _write_training_statistics(
    stats_path: Path,
    outdir_str: str,
    predictions: np.ndarray,
    labeled_ff: np.ndarray,
    labeled_ids: list[str],
    metrics: PerformanceMetrics,
    the_intercept: float,
    the_slope: float,
    elapsed_seconds: int,
) -> None:
    """Format and write training_statistics.txt matching PREFACE.R."""
    index_10 = labeled_ff < (20.0 - predictions)
    deviations_10 = np.abs(predictions[index_10] - labeled_ff[index_10])
    mae_10 = float(np.mean(deviations_10)) if deviations_10.size > 0 else float("nan")
    sd_10 = (
        float(np.std(deviations_10, ddof=1)) if deviations_10.size > 1 else float("nan")
    )

    deviations = np.abs(predictions - labeled_ff)
    outlier_threshold = metrics.mae + 3.0 * metrics.sd_diff
    outlier_indices = np.flatnonzero(deviations > outlier_threshold)

    lines: list[str] = ["PREFACE - PREdict FetAl ComponEnt", ""]
    if outlier_indices.size > 0:
        lines.extend(
            [
                "Below, some of the top candidates for outlier removal are listed.",
                "If you know some of these are low quality/have sex aberrations (when using FFY as response variable), remove them from the config file and re-run.",
                "Avoid removing other cases, as this will result in inaccurate performance statistics and possible overfitting towards irrelevant models.",
                "",
                "_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_",
                "ID\tFF (%) - PREFACE (%)",
            ]
        )
        outlier_abs_vals = deviations[outlier_indices]
        outlier_signed_vals = (predictions - labeled_ff)[outlier_indices]
        for idx in np.argsort(outlier_abs_vals)[::-1]:
            sample_id = labeled_ids[int(outlier_indices[idx])]
            lines.append(f"{sample_id}\t{_format_stat(outlier_signed_vals[idx])}")
        lines.extend(["_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_", ""])

    lines.extend(
        [
            f"Training time: {elapsed_seconds} seconds",
            f"Overall correlation (r): {_format_stat(metrics.cor)}",
            f"Overall mean absolute error (MAE): {_format_stat(metrics.mae)} ± {_format_stat(metrics.sd_diff)}",
            f"FF < 10% mean absolute error (MAE): {_format_stat(mae_10)} ± {_format_stat(sd_10)}",
            "Correction for skew: ",
            f"\tIntercept: {_format_stat(the_intercept)}",
            f"\tSlope: {_format_stat(the_slope)}",
            "",
            (
                f"Do not forget to verify whether the '--nfeat' parameter captures the first 'random' phase "
                f"(and not too much of the 'non-random' phase) at '{outdir_str}overall_performance.png'.\n"
                "If you believe this parameter is not located in an optimal position, decrease/increase '--nfeat' and re-run."
            ),
        ]
    )
    stats_path.write_text("\n".join(lines), encoding="utf-8")


def run_train(
    config_path: Path,
    outdir: Path,
    n_feat: int = 50,
    hidden: int = 2,
    cpus: int = 1,
    femprop: bool = False,
    olm: bool = False,
    noskewcorrect: bool = False,
    device_name: str = "auto",
    seed: int = 1,
) -> Path:
    """Execute the full PREFACE training workflow and export model.onnx."""
    from preface.plotting import plot_ffx, plot_performance

    start_time = time.perf_counter()

    if not config_path.is_file():
        raise FileNotFoundError(f"The file '{config_path}' does not exist.")

    outdir.mkdir(parents=True, exist_ok=True)
    if not outdir.is_dir():
        raise OSError(f"Could not create directory '{outdir}'.")

    resolved_cpus = configure_hardware(cpus=cpus, seed=seed)
    device = resolve_device(device_name)
    rng = np.random.default_rng(seed)

    train_genders = {"M", "F"} if femprop else {"M"}
    config_df = load_config(config_path, rng=rng)
    genders = config_df["gender"].to_list()
    labeled_mask = np.array([g in train_genders for g in genders], dtype=bool)
    n_labeled = int(np.sum(labeled_mask))

    if n_labeled < n_feat:
        raise ValueError(f"Please provide at least {n_feat} labeled samples.")

    dataset = load_training_dataset(config_df=config_df, cpus=resolved_cpus)
    training_frame = dataset.training_frame
    ff_all = config_df["FF"].to_numpy().astype(np.float64, copy=False)
    ids_all = config_df["ID"].to_list()

    repeats_dir = outdir / "training_repeats"
    repeats_dir.mkdir(parents=True, exist_ok=True)

    repeats = 10
    max_feat = n_labeled - int(n_labeled / repeats) - 1
    if n_feat > max_feat:
        print(
            f"Too few samples were provided for --nfeat {n_feat}, using --nfeat {max_feat}"
        )
        n_feat = max_feat

    labeled_indices = np.flatnonzero(labeled_mask)
    predictions = _run_cross_validation(
        training_frame=training_frame,
        ff_all=ff_all,
        labeled_mask=labeled_mask,
        labeled_indices=labeled_indices,
        n_feat=n_feat,
        hidden=hidden,
        olm=olm,
        device=device,
        repeats_dir=repeats_dir,
        repeats=repeats,
    )
    labeled_ff = ff_all[labeled_indices]

    if not noskewcorrect:
        sub_n = max(2, predictions.size // 4)
        sub_idx = rng.permutation(predictions.size)[:sub_n]
        the_intercept, the_slope = ols_fit(predictions[sub_idx], labeled_ff[sub_idx])
    else:
        the_intercept, the_slope = 0.0, 1.0

    male_mask = np.array([g == "M" for g in genders], dtype=bool)
    v1_male_ff = ff_all[male_mask]
    v2_male_x = dataset.x_ratios[male_mask]
    the_intercept_x, the_slope_x, wr = rlm_huber(v1_male_ff, v2_male_x)
    plot_ffx(
        v1_male_ff, v2_male_x, the_intercept_x, the_slope_x, wr, outdir / "FFX.png"
    )

    predictions = the_intercept + the_slope * predictions

    print("Executing final principal component analysis ...")
    pca_final = pca_irlba_torch(
        training_frame,
        n_components=min(n_feat * 10, training_frame.shape[0] - 1),
        device=device,
    )
    final_model = fit_predictor(
        pca_final.x[labeled_indices, :n_feat],
        labeled_ff,
        hidden=hidden,
        is_olm=olm,
        device=device,
    )

    metrics = plot_performance(
        predictions,
        labeled_ff,
        pca_final.importance_prop,
        n_feat,
        "PREFACE (%)",
        "FF (%)",
        outdir / "overall_performance.png",
    )

    outdir_str = f"{outdir.as_posix().rstrip('/')}/"
    _write_training_statistics(
        stats_path=outdir / "training_statistics.txt",
        outdir_str=outdir_str,
        predictions=predictions,
        labeled_ff=labeled_ff,
        labeled_ids=[ids_all[i] for i in labeled_indices],
        metrics=metrics,
        the_intercept=the_intercept,
        the_slope=the_slope,
        elapsed_seconds=round(time.perf_counter() - start_time),
    )

    onnx_path = outdir / "model.onnx"
    export_pipeline_to_onnx(
        predictor=final_model,
        mean_features=dataset.mean_features,
        pca_center=pca_final.center,
        pca_rotation=pca_final.rotation[:, :n_feat],
        the_intercept=the_intercept,
        the_slope=the_slope,
        the_intercept_x=the_intercept_x,
        the_slope_x=the_slope_x,
        possible_features=dataset.possible_features,
        onnx_path=onnx_path,
    )

    print(
        f"Finished! Consult '{outdir_str}training_statistics.txt' to analyse your model's performance."
    )
    return onnx_path


def run_predict(
    infile: Path,
    model_path: Path,
    json_arg: str | None = None,
    device_name: str = "auto",
) -> tuple[float, float]:
    """Run PREFACE prediction on a single sample BED file using a trained ONNX model."""
    if not infile.is_file():
        raise FileNotFoundError(f"The file '{infile}' does not exist.")
    if not model_path.is_file():
        raise FileNotFoundError(f"The file '{model_path}' does not exist.")

    meta = load_onnx_metadata(model_path)
    sample_data = load_inference_bed(infile, meta.possible_features)

    ffx = (sample_data.x_ratio - meta.the_intercept_x) / meta.the_slope_x
    prediction = run_onnx_inference(
        onnx_path=model_path,
        feature_ratios=sample_data.feature_ratios,
        requested_device=device_name,
    )

    emit_prediction_output(
        [f"FFX = {ffx:.4g}%", f"PREFACE = {prediction:.4g}%"],
        {"FFX": ffx / 100.0, "PREFACE": prediction / 100.0},
        json_arg,
    )
    return ffx, prediction
