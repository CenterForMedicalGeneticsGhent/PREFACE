"""Diagnostic plot generation matching PREFACE.R (overall_performance.png, repeat_i.png, FFX.png)."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from preface import COLOR_A, COLOR_B, COLOR_C
from preface.math_utils import get_m_diff, get_sd_diff, ols_fit, pearson_cor


def _signif(value: float, digits: int = 3) -> str:
    """Format floating-point numbers to significant digits matching R's signif()."""
    if not np.isfinite(value):
        return "NA"
    return f"{value:.{digits}g}"


def plot_performance(
    v1: np.ndarray,
    v2: np.ndarray,
    importance_prop: np.ndarray,
    n_feat: int,
    xlab: str,
    ylab: str,
    path: Path,
) -> tuple[float, float, float, float, float]:
    """Generate 3-panel performance diagnostic plot matching plot.performance in PREFACE.R.

    Returns (intercept, slope, mae, sd_diff, cor).
    """
    v1_arr = np.asarray(v1, dtype=np.float64).ravel()
    v2_arr = np.asarray(v2, dtype=np.float64).ravel()
    imp = np.asarray(importance_prop, dtype=np.float64).ravel()

    fig, axes = plt.subplots(1, 3, figsize=(7.6, 2.65), dpi=1024)
    fig.subplots_adjust(left=0.09, right=0.98, bottom=0.20, top=0.86, wspace=0.30)

    # Panel 1: PCA Proportion of variance (log-log)
    ax_pca = axes[0]
    nonzero_imp = imp[imp > 0]
    y_min = float(np.min(nonzero_imp)) if nonzero_imp.size > 0 else 1e-6
    y_max = float(np.max(nonzero_imp)) if nonzero_imp.size > 0 else 1.0
    n_pcs = max(1, imp.size)
    x_idx = np.arange(1, n_pcs + 1, dtype=np.float64)

    log_x = np.log(x_idx)
    log_y = np.log(np.clip(imp, y_min, None))
    ax_pca.plot(log_x, log_y, color=COLOR_A, linewidth=1.8)

    log_nfeat = float(np.log(max(1, n_feat)))
    ax_pca.plot(
        [log_nfeat, log_nfeat],
        [np.log(y_min), np.log(y_max * 0.99)],
        color=COLOR_C,
        linewidth=2.0,
        linestyle=":",
    )
    ax_pca.text(
        log_nfeat,
        np.log(y_max),
        "Number of features",
        color=COLOR_C,
        ha="center",
        va="bottom",
        fontsize=7.5,
    )

    ax_pca.set_xlim(0.0, float(np.log(max(2, n_pcs))))
    ax_pca.set_ylim(float(np.log(y_min)), float(np.log(y_max * 1.05)))
    ax_pca.set_title("PCA", fontsize=9.5, fontweight="bold")
    ax_pca.set_xlabel("Principal components", fontsize=8)
    ax_pca.set_ylabel("Proportion of variance", fontsize=8)

    x_ticks_raw = np.unique(np.round(np.linspace(1, n_pcs, num=min(6, n_pcs))).astype(int))
    ax_pca.set_xticks(np.log(x_ticks_raw))
    ax_pca.set_xticklabels([str(t) for t in x_ticks_raw], fontsize=7)
    ax_pca.set_yticks([np.log(y_min), np.log(y_max)])
    ax_pca.set_yticklabels([f"{y_min:.3g}", f"{y_max:.3g}"], fontsize=7)
    ax_pca.spines["top"].set_visible(False)
    ax_pca.spines["right"].set_visible(False)

    # Panel 2: Scatter plot
    ax_scat = axes[1]
    mx = float(max(np.max(v1_arr, initial=1.0), np.max(v2_arr, initial=1.0)))
    ax_scat.scatter(v1_arr, v2_arr, s=10, color="black", edgecolors="none")

    intercept, slope = ols_fit(v1_arr, v2_arr)
    cor_val = pearson_cor(v1_arr, v2_arr)

    ax_scat.plot([0.0, mx], [0.0, mx], color=COLOR_B, linewidth=2.0, linestyle=":", label="f(x)=x")
    ax_scat.plot(
        [0.0, mx],
        [intercept, intercept + mx * slope],
        color=COLOR_A,
        linewidth=2.0,
        linestyle="--",
        label="OLS fit",
    )
    ax_scat.set_xlim(0.0, mx)
    ax_scat.set_ylim(0.0, mx)
    ax_scat.set_title("Scatter plot", fontsize=9.5, fontweight="bold")
    ax_scat.set_xlabel(xlab, fontsize=8)
    ax_scat.set_ylabel(ylab, fontsize=8)
    ax_scat.tick_params(labelsize=7)
    ax_scat.text(0.0, mx * 1.03, f"(r = {_signif(cor_val, 3)})", fontsize=8, ha="left")
    ax_scat.legend(loc="upper left", frameon=False, fontsize=7.5)
    ax_scat.spines["top"].set_visible(False)
    ax_scat.spines["right"].set_visible(False)

    # Panel 3: Histogram of (v1 - v2)
    ax_hist = axes[2]
    diffs = v1_arr - v2_arr
    n_bins = max(20, int(v1_arr.size / 10))
    counts, breaks, _ = ax_hist.hist(diffs, bins=n_bins, color="black", edgecolor="white", linewidth=0.4)
    hist_mx = float(np.max(counts)) if counts.size > 0 else 1.0

    mean_signed_diff = get_m_diff(v1_arr, v2_arr, use_abs=False)
    mae_val = get_m_diff(v1_arr, v2_arr, use_abs=True)
    sd_val = get_sd_diff(v1_arr, v2_arr)

    ax_hist.plot([0.0, 0.0], [0.0, hist_mx], color=COLOR_B, linewidth=2.0, linestyle=":", label="x=0")
    ax_hist.plot(
        [mean_signed_diff, mean_signed_diff],
        [0.0, hist_mx],
        color=COLOR_A,
        linewidth=2.0,
        linestyle="--",
        label="mean error",
    )
    ax_hist.set_title("Histogram", fontsize=9.5, fontweight="bold")
    ax_hist.set_xlabel(f"{xlab} - {ylab}", fontsize=8)
    ax_hist.set_ylabel("Density", fontsize=8)
    ax_hist.tick_params(labelsize=7)
    ax_hist.text(
        float(np.min(breaks)),
        hist_mx * 1.03,
        f"(MAE = {_signif(mae_val, 3)} ± {_signif(sd_val, 3)})",
        fontsize=8,
        ha="left",
    )
    ax_hist.legend(loc="upper left", frameon=False, fontsize=7.5)
    ax_hist.spines["top"].set_visible(False)
    ax_hist.spines["right"].set_visible(False)

    fig.savefig(path, dpi=1024)
    plt.close(fig)

    return intercept, slope, mae_val, sd_val, cor_val


def plot_ffx(
    v1_ff: np.ndarray,
    v2_x_ratios: np.ndarray,
    fit_intercept: float,
    fit_slope: float,
    weighted_r: float,
    path: Path,
) -> None:
    """Generate 2-panel FFX.png plot matching PREFACE.R lines 301-342."""
    v1 = np.asarray(v1_ff, dtype=np.float64).ravel()
    v2 = np.asarray(v2_x_ratios, dtype=np.float64).ravel()

    fig, axes = plt.subplots(1, 2, figsize=(5.0, 2.8), dpi=1024)
    fig.subplots_adjust(left=0.13, right=0.97, bottom=0.18, top=0.90, wspace=0.34)

    min_v1 = float(np.min(v1)) if v1.size > 0 else 0.0
    max_v1 = float(np.max(v1)) if v1.size > 0 else 1.0

    # Left panel: FF (%) vs mean(ratio X)
    ax_left = axes[0]
    ax_left.scatter(v1, v2, s=8, color="black", edgecolors="none")
    ax_left.plot(
        [min_v1, max_v1],
        [fit_intercept + min_v1 * fit_slope, fit_intercept + max_v1 * fit_slope],
        color=COLOR_A,
        linewidth=2.0,
        linestyle="--",
        label="RLS fit",
    )
    ax_left.plot([], [], color=COLOR_B, linewidth=2.0, linestyle=":", label="f(x)=x")
    ax_left.plot([], [], color="none", label=f"(wr = {_signif(weighted_r, 4)})")
    ax_left.set_xlim(0.0, max_v1)
    ax_left.set_xlabel("FF (%)", fontsize=8)
    ax_left.set_ylabel("μ(ratio X)", fontsize=8)
    ax_left.tick_params(labelsize=7)
    ax_left.legend(loc="upper right", frameon=False, fontsize=6.5)
    ax_left.spines["top"].set_visible(False)
    ax_left.spines["right"].set_visible(False)

    # Right panel: FF (%) vs FFX (%)
    ax_right = axes[1]
    ffx_vals = (v2 - fit_intercept) / fit_slope if abs(fit_slope) > 1e-12 else np.zeros_like(v2)
    ax_right.scatter(v1, ffx_vals, s=8, color="black", edgecolors="none")
    ax_right.plot(
        [min_v1, max_v1],
        [min_v1, max_v1],
        color=COLOR_B,
        linewidth=2.0,
        linestyle=":",
    )
    ax_right.set_xlim(0.0, max_v1)
    ax_right.set_xlabel("FF (%)", fontsize=8)
    ax_right.set_ylabel("FFX (%)", fontsize=8)
    ax_right.tick_params(labelsize=7)
    ax_right.spines["top"].set_visible(False)
    ax_right.spines["right"].set_visible(False)

    fig.savefig(path, dpi=1024)
    plt.close(fig)
