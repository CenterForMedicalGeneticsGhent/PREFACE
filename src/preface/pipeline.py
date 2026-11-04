"""End-to-end training and prediction pipelines for PREFACE."""

from __future__ import annotations

import json
import time
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import torch

from preface.ffy import (
    FFYTrainingResult,
    fit_cohort_ffy_from_profiles,
    format_ffy_statistics_lines,
)
from preface.hardware import configure_hardware, resolve_device
from preface.io import (
    SAMPLESHEET_SUFFIXES,
    align_to_loci,
    build_training_dataset_from_unified,
    extract_unified_sample,
    load_cohort_unified_data,
    load_samplesheet,
    write_samplesheet,
)
from preface.logging_utils import get_logger
from preface.math_utils import format_stat, ols_fit, pca_irlba_torch, rlm_huber
from preface.models import (
    MODEL_TYPE_ALL,
    MODEL_TYPE_FFY,
    MODEL_TYPE_PREFACE,
    UnifiedPredictionResult,
    export_ffy_model_to_onnx,
    export_pipeline_to_onnx,
    export_unified_pipeline_to_onnx,
    fit_predictor,
    load_onnx_metadata,
    predict_model,
    run_unified_onnx_inference,
)

if TYPE_CHECKING:
    from preface.plotting import PerformanceMetrics

logger = get_logger(__name__)


class TrainTarget(StrEnum):
    """Supported training targets for `preface train`."""

    ALL = "all"
    FFY = "ffy"
    PREFACE = "preface"
    FFX = "ffx"


@dataclass(frozen=True)
class PrefixPaths:
    """Resolved output paths derived from a CLI --prefix argument."""

    onnx_path: Path
    samplesheet_path: Path
    training_stats_path: Path
    ffy_stats_path: Path
    overall_plot_path: Path
    ffx_plot_path: Path
    repeats_dir: Path
    predict_json_path: Path


def resolve_prefix_paths(
    prefix: Path | str, samplesheet_suffix: str = ".yaml"
) -> PrefixPaths:
    """Resolve output artifact paths from a file prefix or directory path."""
    raw = str(prefix)
    p = Path(raw)
    ext = (
        samplesheet_suffix
        if samplesheet_suffix.lower() in SAMPLESHEET_SUFFIXES
        else ".yaml"
    )
    if raw.endswith(("/", "\\")) or (p.exists() and p.is_dir()):
        p.mkdir(parents=True, exist_ok=True)
        return PrefixPaths(
            p / "model.onnx",
            p / f"samplesheet{ext}",
            p / "training_statistics.txt",
            p / "ffy_statistics.txt",
            p / "overall_performance.png",
            p / "FFX.png",
            p / "training_repeats",
            p / "predictions.json",
        )
    base = next(
        (
            raw[: -len(s)]
            for s in (".onnx", ".json", ".yaml", ".yml")
            if raw.lower().endswith(s) and len(raw) > len(s)
        ),
        raw,
    )
    Path(base).parent.mkdir(parents=True, exist_ok=True)
    return PrefixPaths(
        p if raw.lower().endswith(".onnx") else Path(f"{base}.onnx"),
        Path(f"{base}.samplesheet{ext}"),
        Path(f"{base}.training_statistics.txt"),
        Path(f"{base}.ffy_statistics.txt"),
        Path(f"{base}.overall_performance.png"),
        Path(f"{base}.FFX.png"),
        Path(f"{base}.training_repeats"),
        p if raw.lower().endswith(".json") else Path(f"{base}.json"),
    )


def emit_prediction_output(
    text_lines: list[str],
    json_dict: dict[str, str | float],
    prefix: Path | str | None = None,
) -> None:
    """Write prediction results as JSON to <prefix>.json or plain-text percentages to stdout."""
    if prefix is None:
        for line in text_lines:
            print(line)
        return
    paths = resolve_prefix_paths(prefix)
    payload = json.dumps(json_dict, indent=2) + "\n"
    paths.predict_json_path.parent.mkdir(parents=True, exist_ok=True)
    paths.predict_json_path.write_text(payload, encoding="utf-8")
    logger.info("Wrote predictions to '%s'", paths.predict_json_path)


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
        logger.info("Model training | Repeat %d/%d ...", rep, repeats)

        start_idx = (rep - 1) * n_labeled // repeats
        end_idx = rep * n_labeled // repeats
        test_index_overall = labeled_indices[start_idx:end_idx]

        is_test_overall = np.zeros(n_total, dtype=bool)
        is_test_overall[test_index_overall] = True
        non_test_indices = np.flatnonzero(~is_test_overall)
        train_index_subset = np.flatnonzero(labeled_mask[non_test_indices])

        logger.info("\tExecuting principal component analysis ...")
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
            "PREFACE",
            "FF",
            repeats_dir / f"repeat_{rep}.png",
        )
        fold_predictions.append(pred_fold)

    return np.concatenate(fold_predictions)


def _write_training_statistics(
    stats_path: Path,
    overall_plot_path: Path,
    predictions: np.ndarray,
    labeled_ff: np.ndarray,
    labeled_ids: list[str],
    metrics: PerformanceMetrics,
    the_intercept: float,
    the_slope: float,
    elapsed_seconds: int,
    ffy_result: FFYTrainingResult | None = None,
) -> None:
    """Format and write training_statistics.txt (including FFY stats when trained together)."""
    index_10 = labeled_ff < (0.20 - predictions)
    deviations_10 = np.abs(predictions[index_10] - labeled_ff[index_10])
    mae_10 = float(np.mean(deviations_10)) if deviations_10.size > 0 else float("nan")
    sd_10 = (
        float(np.std(deviations_10, ddof=1)) if deviations_10.size > 1 else float("nan")
    )

    deviations = np.abs(predictions - labeled_ff)
    outlier_threshold = metrics.mae + 3.0 * metrics.sd_diff
    outlier_indices = np.flatnonzero(deviations > outlier_threshold)

    lines: list[str] = ["PREFACE - PREdict FetAl ComponEnt", ""]
    if ffy_result is not None:
        lines.extend(format_ffy_statistics_lines(ffy_result))
        lines.append("")

    if outlier_indices.size > 0:
        lines.extend(
            [
                "Below, some of the top candidates for outlier removal are listed.",
                "If you know some of these are low quality/have sex aberrations (when using FFY as response variable), remove them from the samplesheet and re-run.",
                "Avoid removing other cases, as this will result in inaccurate performance statistics and possible overfitting towards irrelevant models.",
                "",
                "_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_",
                "ID\tFF - PREFACE",
            ]
        )
        outlier_abs_vals = deviations[outlier_indices]
        outlier_signed_vals = (predictions - labeled_ff)[outlier_indices]
        for idx in np.argsort(outlier_abs_vals)[::-1]:
            sample_id = labeled_ids[int(outlier_indices[idx])]
            lines.append(f"{sample_id}\t{format_stat(outlier_signed_vals[idx])}")
        lines.extend(["_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_", ""])

    lines.extend(
        [
            f"Training time: {elapsed_seconds} seconds",
            f"Overall correlation (r): {format_stat(metrics.cor)}",
            f"Overall mean absolute error (MAE): {format_stat(metrics.mae)} ± {format_stat(metrics.sd_diff)}",
            f"FF < 0.1 mean absolute error (MAE): {format_stat(mae_10)} ± {format_stat(sd_10)}",
            "Correction for skew: ",
            f"\tIntercept: {format_stat(the_intercept)}",
            f"\tSlope: {format_stat(the_slope)}",
            "",
            (
                f"Do not forget to verify whether the '--nfeat' parameter captures the first 'random' phase "
                f"(and not too much of the 'non-random' phase) at '{overall_plot_path.as_posix()}'.\n"
                "If you believe this parameter is not located in an optimal position, decrease/increase '--nfeat' and re-run."
            ),
        ]
    )
    stats_path.parent.mkdir(parents=True, exist_ok=True)
    stats_path.write_text("\n".join(lines), encoding="utf-8")


def run_train(
    samplesheet_path: Path,
    prefix: Path | str,
    target: TrainTarget | str = TrainTarget.ALL,
    n_feat: int = 50,
    hidden: int = 2,
    cpus: int = 1,
    femprop: bool = False,
    olm: bool = False,
    noskewcorrect: bool = False,
    binsize: int = 100_000,
    ffy_binsize: int = 5_000,
    mapq: int = 1,
    yfrac: float | None = None,
    lasso_lambda: float = 1e-4,
    no_female_gate: bool = False,
    reference: Path | None = None,
    male_ref_paths: list[Path] | None = None,
    cache_counts: Path | None = None,
    device_name: str = "auto",
    seed: int = 1,
) -> Path:
    """Execute the single-pass FFY and/or PREFACE training workflow and export <prefix>.onnx."""
    from preface.plotting import plot_ffx, plot_performance

    start_time = time.perf_counter()
    resolved_target = TrainTarget(str(target).lower())
    run_ffy = resolved_target in {TrainTarget.ALL, TrainTarget.FFY}
    run_preface = resolved_target in {
        TrainTarget.ALL,
        TrainTarget.PREFACE,
        TrainTarget.FFX,
    }

    raw_entries = load_samplesheet(samplesheet_path)
    paths = resolve_prefix_paths(
        prefix, samplesheet_suffix=samplesheet_path.suffix.lower()
    )

    resolved_cpus = configure_hardware(cpus=cpus, seed=seed)
    device = resolve_device(device_name)
    rng = np.random.default_rng(seed)

    train_sexes = {"M", "F"} if femprop else {"M"}
    if not run_ffy:
        missing_labels = any(
            e.sex is None or (e.sex in train_sexes and e.ff is None)
            for e in raw_entries
        )
        if missing_labels:
            raise ValueError(
                "Samplesheet has missing 'sex' or 'ff' labels required for 'preface train preface'. "
                "Run 'preface train all' or 'preface train ffy' first to populate baseline FFY labels."
            )

    logger.info(
        "Reading %d cohort sample(s) in a single pass (include_ffy=%s, include_cn=%s) ...",
        len(raw_entries),
        run_ffy,
        run_preface,
    )
    unified_samples, lb_profiles = load_cohort_unified_data(
        entries=raw_entries,
        binsize=binsize,
        ffy_binsize=ffy_binsize,
        mapq=mapq,
        reference=reference,
        male_ref_paths=male_ref_paths,
        cpus=resolved_cpus,
        include_cn=run_preface,
        include_ffy=run_ffy,
        cache_counts=cache_counts,
    )

    ffy_result: FFYTrainingResult | None = None
    if run_ffy:
        cohort_y_profiles = [
            s.y_profile for s in unified_samples if s.y_profile is not None
        ]
        ffy_result = fit_cohort_ffy_from_profiles(
            entries=raw_entries,
            cohort_profiles=cohort_y_profiles,
            lb_profiles=lb_profiles,
            binsize=ffy_binsize,
            yfrac=yfrac,
            lasso_lambda=lasso_lambda,
            no_female_gate=no_female_gate,
            seed=seed,
        )
        raw_entries = ffy_result.entries
        default_fmt = "json" if samplesheet_path.suffix.lower() == ".json" else "yaml"
        write_samplesheet(
            raw_entries,
            outfile=paths.samplesheet_path,
            default_format=default_fmt,
        )

        if not run_preface:
            elapsed_ffy = round(time.perf_counter() - start_time)
            ffy_lines = format_ffy_statistics_lines(ffy_result) + [
                f"Training time: {elapsed_ffy} seconds",
            ]
            paths.ffy_stats_path.parent.mkdir(parents=True, exist_ok=True)
            paths.ffy_stats_path.write_text("\n".join(ffy_lines), encoding="utf-8")
            export_ffy_model_to_onnx(
                intercept=ffy_result.intercept,
                weights=ffy_result.weights,
                bin_loci=ffy_result.bin_loci,
                informative_indices=ffy_result.informative_indices,
                female_baseline=ffy_result.female_baseline,
                male_100_baseline=ffy_result.male_100_baseline,
                yfrac_threshold=ffy_result.yfrac_threshold,
                binsize=ffy_result.binsize,
                onnx_path=paths.onnx_path,
                no_female_gate=no_female_gate,
            )
            logger.info(
                "Finished FFY training! Consult '%s' and '%s'.",
                paths.ffy_stats_path,
                paths.samplesheet_path,
            )
            return paths.onnx_path

    perm = rng.permutation(len(raw_entries))
    entries = [raw_entries[int(i)] for i in perm]

    sexes = [e.sex for e in entries]
    ff_all = np.array(
        [float(e.ff) if e.ff is not None else np.nan for e in entries],
        dtype=np.float64,
    )
    labeled_mask = np.array(
        [s in train_sexes and np.isfinite(f) for s, f in zip(sexes, ff_all)],
        dtype=bool,
    )
    n_labeled = int(np.sum(labeled_mask))

    if n_labeled < n_feat:
        raise ValueError(f"Please provide at least {n_feat} labeled samples.")

    dataset = build_training_dataset_from_unified(
        entries=entries, unified_samples=unified_samples
    )
    training_frame = dataset.training_frame
    ids_all = [e.id for e in entries]

    paths.repeats_dir.mkdir(parents=True, exist_ok=True)

    repeats = 10
    max_feat = n_labeled - int(n_labeled / repeats) - 1
    if n_feat > max_feat:
        logger.warning(
            "Too few samples were provided for --nfeat %d, using --nfeat %d",
            n_feat,
            max_feat,
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
        repeats_dir=paths.repeats_dir,
        repeats=repeats,
    )
    labeled_ff = ff_all[labeled_indices]

    if not noskewcorrect:
        sub_n = max(2, predictions.size // 4)
        sub_idx = rng.permutation(predictions.size)[:sub_n]
        the_intercept, the_slope = ols_fit(predictions[sub_idx], labeled_ff[sub_idx])
    else:
        the_intercept, the_slope = 0.0, 1.0

    male_mask = np.array(
        [s == "M" and np.isfinite(f) for s, f in zip(sexes, ff_all)], dtype=bool
    )
    v1_male_ff = ff_all[male_mask]
    v2_male_x = dataset.x_ratios[male_mask]
    the_intercept_x, the_slope_x, wr = rlm_huber(v1_male_ff, v2_male_x)
    plot_ffx(
        v1_male_ff,
        v2_male_x,
        the_intercept_x,
        the_slope_x,
        wr,
        paths.ffx_plot_path,
    )

    predictions = the_intercept + the_slope * predictions

    logger.info("Executing final principal component analysis ...")
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
        "PREFACE",
        "FF",
        paths.overall_plot_path,
    )

    _write_training_statistics(
        stats_path=paths.training_stats_path,
        overall_plot_path=paths.overall_plot_path,
        predictions=predictions,
        labeled_ff=labeled_ff,
        labeled_ids=[ids_all[i] for i in labeled_indices],
        metrics=metrics,
        the_intercept=the_intercept,
        the_slope=the_slope,
        elapsed_seconds=round(time.perf_counter() - start_time),
        ffy_result=ffy_result,
    )

    if ffy_result is not None:
        export_unified_pipeline_to_onnx(
            predictor=final_model,
            mean_features=dataset.mean_features,
            pca_center=pca_final.center,
            pca_rotation=pca_final.rotation[:, :n_feat],
            the_intercept=the_intercept,
            the_slope=the_slope,
            the_intercept_x=the_intercept_x,
            the_slope_x=the_slope_x,
            possible_features=dataset.possible_features,
            ffy_intercept=ffy_result.intercept,
            ffy_weights=ffy_result.weights,
            bin_loci=ffy_result.bin_loci,
            informative_indices=ffy_result.informative_indices,
            female_baseline=ffy_result.female_baseline,
            male_100_baseline=ffy_result.male_100_baseline,
            yfrac_threshold=ffy_result.yfrac_threshold,
            onnx_path=paths.onnx_path,
            binsize=binsize,
            ffy_binsize=ffy_binsize,
            no_female_gate=no_female_gate,
        )
    else:
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
            onnx_path=paths.onnx_path,
            binsize=binsize,
        )

    logger.info(
        "Finished! Consult '%s' to analyse your model's performance.",
        paths.training_stats_path,
    )
    return paths.onnx_path


def run_predict(
    infile: Path,
    model_path: Path,
    prefix: Path | str | None = None,
    mapq: int = 1,
    reference: Path | None = None,
    no_female_gate: bool | None = None,
    device_name: str = "auto",
) -> UnifiedPredictionResult:
    """Run unified FFY / FFX / PREFACE prediction on a single alignment file using an ONNX model."""
    if not infile.is_file():
        raise FileNotFoundError(f"The file '{infile}' does not exist.")
    if not model_path.is_file():
        raise FileNotFoundError(f"The file '{model_path}' does not exist.")

    logger.debug("Loading ONNX model from '%s'", model_path)
    meta = load_onnx_metadata(model_path)

    need_cn = meta.model_type in {MODEL_TYPE_ALL, MODEL_TYPE_PREFACE}
    need_ffy = meta.model_type in {MODEL_TYPE_ALL, MODEL_TYPE_FFY}

    unified_sample = extract_unified_sample(
        infile,
        sample_id=infile.name,
        binsize=meta.binsize,
        ffy_binsize=meta.ffy_binsize,
        mapq=mapq,
        reference=reference,
        include_cn=need_cn,
        include_ffy=need_ffy,
    )

    selected_ratios: np.ndarray | None = None
    if need_cn:
        assert unified_sample.feature_loci is not None
        assert unified_sample.feature_ratios is not None
        selected_ratios = align_to_loci(
            unified_sample.feature_loci,
            unified_sample.feature_ratios,
            meta.possible_features,
        )

    aligned_y: np.ndarray | None = None
    if need_ffy:
        assert unified_sample.y_profile is not None
        aligned_y = align_to_loci(
            unified_sample.y_profile.bin_loci,
            unified_sample.y_profile.y_norm_bins,
            meta.bin_loci,
            fill_nan=0.0,
        )

    result = run_unified_onnx_inference(
        onnx_path=model_path,
        meta=meta,
        feature_ratios=selected_ratios,
        x_ratio=unified_sample.x_ratio,
        aligned_y=aligned_y,
        requested_device=device_name,
        no_female_gate_override=no_female_gate,
    )

    logger.debug(
        "Predicted sex=%s, FFY=%s, FFX=%s, PREFACE=%s for '%s'",
        result.sex,
        result.ffy,
        result.ffx,
        result.preface,
        infile,
    )

    emit_prediction_output(
        result.to_text_lines(),
        result.to_json_dict(),
        prefix=prefix,
    )
    return result
