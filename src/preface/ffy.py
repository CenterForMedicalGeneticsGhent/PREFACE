"""Baseline FFY (chromosome Y fetal fraction) calculation per PMC6771918 Section 3.4."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from preface.io import (
    FFY_TARGET_CHRS,
    SampleEntry,
    compute_autosomal_density,
    emit_prediction_output,
    load_samplesheet,
    parallel_load_samples,
    read_alignment_bins,
    write_samplesheet,
)
from preface.logging_utils import get_logger

logger = get_logger(__name__)

CHRY_PAR1_END: int = 2_781_479
CHRY_HETEROCHROMATIN_START: int = 28_800_000


@dataclass(frozen=True)
class SampleYProfile:
    """Normalized chromosome Y bin depths and coordinates for a single sample."""

    bin_loci: list[str]
    bin_starts: np.ndarray
    bin_ends: np.ndarray
    y_norm_bins: np.ndarray


def _coordinate_informative_mask(starts: np.ndarray, ends: np.ndarray) -> np.ndarray:
    """Mask out PAR1, PAR2, and Yq12 heterochromatin when coordinates span human chrY."""
    if ends.size == 0:
        return np.zeros(0, dtype=bool)
    if int(np.max(ends)) > CHRY_PAR1_END:
        mask = (starts >= CHRY_PAR1_END) & (ends <= CHRY_HETEROCHROMATIN_START)
        if np.any(mask):
            return mask
    return np.ones(starts.size, dtype=bool)


def extract_sample_y_profile(
    input_path: Path,
    binsize: int = 5000,
    mapq: int = 1,
    reference: Path | None = None,
    index: Path | None = None,
) -> SampleYProfile:
    """Extract diploid-normalized chromosome Y bin profile from an alignment file."""
    ref_lengths, bin_counts = read_alignment_bins(
        input_path,
        binsize=binsize,
        mapq=mapq,
        reference=reference,
        target_chrs=FFY_TARGET_CHRS,
        index=index,
    )
    if "Y" not in ref_lengths or "Y" not in bin_counts:
        raise ValueError(
            "Alignment header does not contain chromosome Y ('Y' or 'chrY')."
        )

    r_auto = compute_autosomal_density(ref_lengths, bin_counts, binsize)

    y_len = ref_lengths["Y"]
    y_starts0 = np.arange(0, y_len, binsize, dtype=np.int64)
    y_starts = y_starts0 + 1
    y_ends = np.minimum(y_starts0 + binsize, y_len)
    y_widths = np.maximum(1, y_ends - y_starts0).astype(np.float64)
    y_counts = bin_counts["Y"].astype(np.float64)

    y_norm = np.where(
        np.isfinite(y_counts),
        y_counts / (r_auto * y_widths),
        np.nan,
    )
    bin_loci = [f"Y:{int(s)}-{int(e)}" for s, e in zip(y_starts, y_ends)]

    logger.debug(
        "Extracted chrY profile from '%s': %d chrY bins",
        input_path,
        len(bin_loci),
    )
    return SampleYProfile(
        bin_loci=bin_loci,
        bin_starts=y_starts,
        bin_ends=y_ends,
        y_norm_bins=y_norm,
    )


def _align_profile_to_loci(
    profile: SampleYProfile, target_loci: list[str]
) -> np.ndarray:
    if profile.bin_loci == target_loci:
        return np.nan_to_num(profile.y_norm_bins, nan=0.0)
    locus_map = dict(zip(profile.bin_loci, profile.y_norm_bins))
    arr = np.array([locus_map.get(loc, 0.0) for loc in target_loci], dtype=np.float64)
    return np.nan_to_num(arr, nan=0.0)


def _extract_profile_worker(
    args: tuple[Path, int, int, Path | None, Path | None],
) -> SampleYProfile:
    path, binsize, mapq, reference, index = args
    return extract_sample_y_profile(
        path, binsize=binsize, mapq=mapq, reference=reference, index=index
    )


def run_ffy_single(
    infile: Path,
    model_path: Path | None = None,
    binsize: int = 5000,
    mapq: int = 1,
    yfrac: float | None = None,
    reference: Path | None = None,
    json_arg: str | None = None,
    device_name: str = "auto",
) -> tuple[str, float]:
    """Calculate sex and baseline FFY (0-1 fraction) for a single BAM/CRAM/SAM sample."""
    from preface.models import load_ffy_onnx_metadata, run_ffy_onnx_inference

    if model_path is not None:
        if not model_path.is_file():
            raise FileNotFoundError(f"The file '{model_path}' does not exist.")
        meta = load_ffy_onnx_metadata(model_path)
        profile = extract_sample_y_profile(
            infile, binsize=meta.binsize, mapq=mapq, reference=reference
        )
        aligned_y = _align_profile_to_loci(profile, meta.bin_loci)
        info_idx = np.asarray(meta.informative_indices, dtype=np.int64)
        y_score = (
            float(np.mean(aligned_y[info_idx]))
            if info_idx.size > 0
            else float(np.mean(aligned_y))
        )
        thresh = yfrac if yfrac is not None else meta.yfrac_threshold
        sex = "M" if y_score >= thresh else "F"
        ffy = (
            run_ffy_onnx_inference(model_path, aligned_y, device_name)
            if sex == "M"
            else 0.0
        )
    else:
        profile = extract_sample_y_profile(
            infile, binsize=binsize, mapq=mapq, reference=reference
        )
        coord_mask = _coordinate_informative_mask(profile.bin_starts, profile.bin_ends)
        vals = profile.y_norm_bins[coord_mask]
        finite_vals = vals[np.isfinite(vals)]
        if finite_vals.size == 0:
            raise ValueError(f"Sample '{infile}' has no finite chromosome Y bins.")
        y_score = float(np.mean(finite_vals))
        thresh = yfrac if yfrac is not None else 0.01
        sex = "M" if y_score >= thresh else "F"
        ffy = float(np.clip(2.0 * y_score, 0.0, 1.0)) if sex == "M" else 0.0

    logger.debug("Single-sample FFY for '%s': sex=%s, FFY=%.6f", infile, sex, ffy)
    emit_prediction_output(
        [f"sex = {sex}", f"FFY = {ffy:.4g}"],
        {"FFY": ffy, "sex": sex},
        json_arg,
    )
    return sex, ffy


def estimate_cohort_ffy(
    entries: Sequence[SampleEntry],
    male_ref_paths: list[Path] | None = None,
    save_model_path: Path | None = None,
    binsize: int = 5000,
    mapq: int = 1,
    yfrac: float | None = None,
    lasso_lambda: float = 1e-4,
    reference: Path | None = None,
    cpus: int = 1,
    seed: int = 1,
) -> list[SampleEntry]:
    """Estimate baseline FFY (0-1 fraction) and sex across a cohort of SampleEntry objects."""
    from preface.math_utils import fit_lasso_fista, infer_sex_gmm_1d
    from preface.models import export_ffy_model_to_onnx

    entry_list = list(entries)
    n_samples = len(entry_list)
    if n_samples == 0:
        raise ValueError("Samplesheet contains no samples.")

    rng = np.random.default_rng(seed)

    combined_items: list[
        tuple[str, tuple[Path, int, int, Path | None, Path | None]]
    ] = [
        (entry.id, (entry.filepath, binsize, mapq, reference, entry.index))
        for entry in entry_list
    ]
    if male_ref_paths:
        combined_items.extend(
            (f"male_ref:{p.name}", (p, binsize, mapq, reference, None))
            for p in male_ref_paths
        )

    logger.info("Extracting chrY profiles for %d cohort samples ...", n_samples)
    all_profiles = parallel_load_samples(
        combined_items, _extract_profile_worker, cpus=cpus
    )
    cohort_profiles = all_profiles[:n_samples]
    lb_profiles = all_profiles[n_samples:]

    first_prof = cohort_profiles[0]
    ref_loci = first_prof.bin_loci
    y_matrix = np.vstack([_align_profile_to_loci(p, ref_loci) for p in cohort_profiles])

    coord_mask = _coordinate_informative_mask(
        first_prof.bin_starts, first_prof.bin_ends
    )
    non_empty_mask = np.max(y_matrix, axis=0) > 0.0
    init_info_mask = (
        coord_mask & non_empty_mask
        if np.any(coord_mask & non_empty_mask)
        else coord_mask
    )
    raw_scores = np.mean(y_matrix[:, init_info_mask], axis=1)

    inferred_all, learned_thresh = infer_sex_gmm_1d(raw_scores, yfrac=yfrac)
    sexes = [
        entry.sex.strip().upper()
        if entry.sex is not None and entry.sex.strip().upper() in {"M", "F"}
        else inf
        for entry, inf in zip(entry_list, inferred_all)
    ]

    male_mask = np.array([s == "M" for s in sexes], dtype=bool)
    female_mask = ~male_mask
    n_males = int(np.sum(male_mask))
    n_females = int(np.sum(female_mask))
    logger.info(
        "Classified cohort fetal sex: %d male, %d female (threshold=%.6f)",
        n_males,
        n_females,
        learned_thresh,
    )

    info_mask = init_info_mask.copy()
    if np.any(male_mask) and np.any(female_mask):
        med_f = np.median(y_matrix[female_mask, :], axis=0)
        mean_m = np.mean(y_matrix[male_mask, :], axis=0)
        refined = coord_mask & (mean_m > med_f + 1e-4) & (med_f < 0.15)
        if np.any(refined):
            info_mask = refined

    info_indices = np.flatnonzero(info_mask).tolist()
    logger.info(
        "Selected %d informative chrY bins out of %d total bins",
        len(info_indices),
        len(ref_loci),
    )
    y_scores = np.mean(y_matrix[:, info_mask], axis=1)

    if np.any(male_mask) and np.any(female_mask) and yfrac is None:
        learned_thresh = float(
            0.5 * (np.median(y_scores[female_mask]) + np.min(y_scores[male_mask]))
        )

    female_baseline = (
        float(np.median(y_scores[female_mask])) if np.any(female_mask) else 0.0
    )

    if lb_profiles:
        logger.info(
            "Using %d male reference (LB) sample(s) for 100%% male baseline ...",
            len(lb_profiles),
        )
        lb_matrix = np.vstack(
            [_align_profile_to_loci(p, ref_loci) for p in lb_profiles]
        )
        male_100_baseline = float(np.median(np.mean(lb_matrix[:, info_mask], axis=1)))
    else:
        male_100_baseline = female_baseline + 0.5

    logger.debug(
        "Calibration baselines: female=%.6f, male_100=%.6f",
        female_baseline,
        male_100_baseline,
    )

    denom = max(male_100_baseline - female_baseline, 1e-6)
    prior_ffy = np.zeros(n_samples, dtype=np.float64)
    prior_ffy[male_mask] = np.clip(
        (y_scores[male_mask] - female_baseline) / denom, 0.0, 1.0
    )

    x_lasso = np.zeros_like(y_matrix)
    x_lasso[:, info_mask] = y_matrix[:, info_mask]

    if n_males >= 4 and n_samples >= 4:
        k_folds = min(10, n_samples)
        logger.info(
            "Running %d-fold cross-validated Lasso (lambda=%g) across %d samples ...",
            k_folds,
            lasso_lambda,
            n_samples,
        )
        perm = rng.permutation(n_samples)
        cv_preds = np.zeros(n_samples, dtype=np.float64)
        for fold in range(k_folds):
            test_idx = perm[
                fold * n_samples // k_folds : (fold + 1) * n_samples // k_folds
            ]
            train_mask = np.ones(n_samples, dtype=bool)
            train_mask[test_idx] = False
            b0, w_fold = fit_lasso_fista(
                x_lasso[train_mask, :], prior_ffy[train_mask], lam=lasso_lambda
            )
            cv_preds[test_idx] = b0 + x_lasso[test_idx, :] @ w_fold
        final_ffy = np.where(male_mask, np.clip(cv_preds, 0.0, 1.0), 0.0)
    else:
        logger.warning(
            "Fewer than 4 male/total samples (%d male, %d total); skipping CV Lasso and using prior FFY.",
            n_males,
            n_samples,
        )
        final_ffy = prior_ffy.copy()

    if save_model_path is not None:
        full_b0, full_w = fit_lasso_fista(x_lasso, prior_ffy, lam=lasso_lambda)
        export_ffy_model_to_onnx(
            intercept=full_b0,
            weights=full_w,
            bin_loci=ref_loci,
            informative_indices=info_indices,
            female_baseline=female_baseline,
            male_100_baseline=male_100_baseline,
            yfrac_threshold=learned_thresh,
            binsize=binsize,
            onnx_path=save_model_path,
        )

    return [
        SampleEntry(
            id=entry.id,
            filepath=entry.filepath,
            index=entry.index,
            sex=sex,
            ff=float(entry.ff) if entry.ff is not None else float(ffy_val),
        )
        for entry, sex, ffy_val in zip(entry_list, sexes, final_ffy)
    ]


def run_ffy_cohort(
    samplesheet_path: Path,
    outfile: Path | None = None,
    male_ref_paths: list[Path] | None = None,
    save_model_path: Path | None = None,
    binsize: int = 5000,
    mapq: int = 1,
    yfrac: float | None = None,
    lasso_lambda: float = 1e-4,
    reference: Path | None = None,
    cpus: int = 1,
    seed: int = 1,
) -> list[SampleEntry]:
    """Calculate baseline FFY and sex across a cohort samplesheet per PMC6771918 Section 3.4."""
    from preface.hardware import configure_hardware

    resolved_cpus = configure_hardware(cpus=cpus, seed=seed)
    entries = load_samplesheet(samplesheet_path)
    updated_entries = estimate_cohort_ffy(
        entries=entries,
        male_ref_paths=male_ref_paths,
        save_model_path=save_model_path,
        binsize=binsize,
        mapq=mapq,
        yfrac=yfrac,
        lasso_lambda=lasso_lambda,
        reference=reference,
        cpus=resolved_cpus,
        seed=seed,
    )
    default_fmt = "json" if samplesheet_path.suffix.lower() == ".json" else "yaml"
    write_samplesheet(updated_entries, outfile=outfile, default_format=default_fmt)
    return updated_entries
