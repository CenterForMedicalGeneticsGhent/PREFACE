"""Baseline FFY (chromosome Y fetal fraction) calculation per PMC6771918 Section 3.4."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from preface.io import (
    SampleEntry,
    SampleYProfile,
    align_to_loci,
    extract_unified_sample,
)
from preface.logging_utils import get_logger
from preface.math_utils import (
    fit_lasso_fista,
    format_stat,
    get_m_diff,
    get_sd_diff,
    infer_sex_gmm_1d,
    pearson_cor,
)

logger = get_logger(__name__)

CHRY_PAR1_END: int = 2_781_479
CHRY_HETEROCHROMATIN_START: int = 28_800_000


@dataclass(frozen=True)
class FFYTrainingResult:
    """Calibrated FFY model parameters, updated samplesheet entries, and summary statistics."""

    entries: list[SampleEntry]
    intercept: float
    weights: np.ndarray
    bin_loci: list[str]
    informative_indices: list[int]
    female_baseline: float
    male_100_baseline: float
    yfrac_threshold: float
    binsize: int
    no_female_gate: bool
    n_males: int
    n_females: int
    cv_cor: float
    cv_mae: float
    cv_sd: float


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
    data = extract_unified_sample(
        input_path,
        ffy_binsize=binsize,
        mapq=mapq,
        reference=reference,
        index=index,
        include_cn=False,
        include_ffy=True,
    )
    assert data.y_profile is not None
    logger.debug(
        "Extracted chrY profile from '%s': %d chrY bins",
        input_path,
        len(data.y_profile.bin_loci),
    )
    return data.y_profile


def fit_cohort_ffy_from_profiles(
    entries: Sequence[SampleEntry],
    cohort_profiles: Sequence[SampleYProfile],
    lb_profiles: Sequence[SampleYProfile] | None = None,
    binsize: int = 5000,
    yfrac: float | None = None,
    lasso_lambda: float = 1e-4,
    no_female_gate: bool = False,
    seed: int = 1,
) -> FFYTrainingResult:
    """Fit 2-GMM sex classification and non-negative Lasso FFY calibration on extracted chrY profiles."""
    entry_list = list(entries)
    n_samples = len(entry_list)
    if n_samples == 0:
        raise ValueError("Samplesheet contains no samples.")
    if len(cohort_profiles) != n_samples:
        raise ValueError("Mismatch between cohort entries and extracted chrY profiles.")

    rng = np.random.default_rng(seed)

    first_prof = cohort_profiles[0]
    ref_loci = first_prof.bin_loci
    y_matrix = np.vstack(
        [
            align_to_loci(p.bin_loci, p.y_norm_bins, ref_loci, fill_nan=0.0)
            for p in cohort_profiles
        ]
    )

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
            [
                align_to_loci(p.bin_loci, p.y_norm_bins, ref_loci, fill_nan=0.0)
                for p in lb_profiles
            ]
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
    unclamped_prior = np.clip((y_scores - female_baseline) / denom, 0.0, 1.0)
    if no_female_gate:
        prior_ffy = unclamped_prior
    else:
        prior_ffy = np.where(male_mask, unclamped_prior, 0.0)

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
        clipped_cv = np.clip(cv_preds, 0.0, 1.0)
        final_ffy = (
            clipped_cv if no_female_gate else np.where(male_mask, clipped_cv, 0.0)
        )
    else:
        logger.warning(
            "Fewer than 4 male/total samples (%d male, %d total); skipping CV Lasso and using prior FFY.",
            n_males,
            n_samples,
        )
        final_ffy = prior_ffy.copy()

    full_b0, full_w = fit_lasso_fista(x_lasso, prior_ffy, lam=lasso_lambda)

    eval_target = np.array(
        [
            float(e.ff) if e.ff is not None else float(p)
            for e, p in zip(entry_list, prior_ffy)
        ],
        dtype=np.float64,
    )
    eval_mask = male_mask if np.any(male_mask) else np.ones(n_samples, dtype=bool)
    cv_cor = pearson_cor(final_ffy[eval_mask], eval_target[eval_mask])
    cv_mae = get_m_diff(final_ffy[eval_mask], eval_target[eval_mask], use_abs=True)
    cv_sd = get_sd_diff(final_ffy[eval_mask], eval_target[eval_mask])

    updated_entries = [
        SampleEntry(
            id=entry.id,
            filepath=entry.filepath,
            index=entry.index,
            sex=sex,
            ff=float(entry.ff) if entry.ff is not None else float(ffy_val),
        )
        for entry, sex, ffy_val in zip(entry_list, sexes, final_ffy)
    ]

    return FFYTrainingResult(
        entries=updated_entries,
        intercept=full_b0,
        weights=full_w,
        bin_loci=ref_loci,
        informative_indices=info_indices,
        female_baseline=female_baseline,
        male_100_baseline=male_100_baseline,
        yfrac_threshold=learned_thresh,
        binsize=binsize,
        no_female_gate=no_female_gate,
        n_males=n_males,
        n_females=n_females,
        cv_cor=cv_cor,
        cv_mae=cv_mae,
        cv_sd=cv_sd,
    )


def format_ffy_statistics_lines(result: FFYTrainingResult) -> list[str]:
    """Format human-readable FFY model statistics lines."""
    return [
        "FFY - Baseline Chromosome Y Fetal Fraction Model",
        f"Cohort sex classification: {result.n_males} male, {result.n_females} female",
        f"Sex classification threshold (yfrac): {format_stat(result.yfrac_threshold)}",
        f"Female baseline depth: {format_stat(result.female_baseline)}",
        f"100% male baseline depth: {format_stat(result.male_100_baseline)}",
        f"Informative chrY bins: {len(result.informative_indices)} / {len(result.bin_loci)} (binsize={result.binsize} bp)",
        f"Female FFY gating disabled (--no-female-gate): {result.no_female_gate}",
        f"FFY Lasso CV correlation (r): {format_stat(result.cv_cor)}",
        f"FFY Lasso CV MAE: {format_stat(result.cv_mae)} ± {format_stat(result.cv_sd)}",
    ]
