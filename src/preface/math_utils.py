"""Statistical and linear algebra routines matching PREFACE.R (PCA, OLS, Huber RLM)."""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
import torch


@dataclass(frozen=True)
class PCAResult:
    """Principal Component Analysis result matching R's prcomp_irlba + summary()."""

    x: np.ndarray
    center: np.ndarray
    rotation: np.ndarray
    sdev: np.ndarray
    importance_prop: np.ndarray

    def project(self, matrix: np.ndarray, n_feat: int) -> np.ndarray:
        """Center and project samples onto the top n_feat principal components."""
        return (matrix - self.center) @ self.rotation[:, :n_feat]


def pca_irlba_torch(
    matrix: np.ndarray,
    n_components: int,
    device: torch.device,
) -> PCAResult:
    """Compute unscaled, centered truncated PCA matching R's prcomp_irlba(..., scale.=FALSE)."""
    n_samples, n_features = matrix.shape
    k = max(1, min(n_components, n_samples - 1, n_features))

    center = np.mean(matrix, axis=0, dtype=np.float64)
    centered_np = matrix - center
    total_var = float(np.sum(np.var(centered_np, axis=0, ddof=1)))

    use_dtype = torch.float32 if device.type == "mps" else torch.float64
    linalg_device = torch.device("cpu") if device.type == "mps" else device
    x_t = torch.as_tensor(centered_np, dtype=use_dtype, device=device)

    if n_samples <= n_features:
        gram = (x_t @ x_t.T).to(device=linalg_device, dtype=torch.float64)
        evals, evecs = torch.linalg.eigh(gram)
        evals = evals.flip(0)[:k].clamp(min=0.0)
        u_k = evecs.flip(1)[:, :k]
        s_k = torch.sqrt(evals)
        safe_s = torch.where(s_k > 1e-12, s_k, torch.ones_like(s_k))
        x_linalg = torch.as_tensor(
            centered_np, dtype=torch.float64, device=linalg_device
        )
        v_k = ((x_linalg.T @ u_k) / safe_s.unsqueeze(0)).cpu()
        u_k = u_k.cpu()
        s_k = s_k.cpu()
    else:
        svd_input = x_t.to(device=linalg_device, dtype=torch.float64)
        u_full, s_full, vh_full = torch.linalg.svd(svd_input, full_matrices=False)
        u_k = u_full[:, :k].cpu()
        s_k = s_full[:k].cpu()
        v_k = vh_full[:k, :].T.cpu()

    rotation = v_k.numpy().astype(np.float64, copy=False)
    s_np = s_k.numpy().astype(np.float64, copy=False)
    scores = (u_k.numpy() * s_np[np.newaxis, :]).astype(np.float64, copy=False)

    max_abs_idx = np.argmax(np.abs(rotation), axis=0)
    signs = np.sign(rotation[max_abs_idx, np.arange(k)])
    signs[signs == 0] = 1.0
    rotation = rotation * signs[np.newaxis, :]
    scores = scores * signs[np.newaxis, :]

    sdev = s_np / np.sqrt(max(1, n_samples - 1))
    importance_prop = (sdev**2) / total_var if total_var > 0 else np.zeros_like(sdev)

    return PCAResult(
        x=scores,
        center=center,
        rotation=rotation,
        sdev=sdev,
        importance_prop=importance_prop,
    )


def get_m_diff(v1: np.ndarray, v2: np.ndarray, use_abs: bool = True) -> float:
    """Compute mean difference or mean absolute difference matching get.m.diff in PREFACE.R."""
    diff = np.asarray(v1, dtype=np.float64) - np.asarray(v2, dtype=np.float64)
    return float(np.mean(np.abs(diff) if use_abs else diff))


def get_sd_diff(v1: np.ndarray, v2: np.ndarray) -> float:
    """Compute sample standard deviation (ddof=1) of absolute errors matching get.sd.diff."""
    abs_diff = np.abs(
        np.asarray(v1, dtype=np.float64) - np.asarray(v2, dtype=np.float64)
    )
    return float(np.std(abs_diff, ddof=1)) if abs_diff.size > 1 else 0.0


def pearson_cor(v1: np.ndarray, v2: np.ndarray) -> float:
    """Compute Pearson correlation coefficient matching R's cor(v1, v2)."""
    x = np.asarray(v1, dtype=np.float64).ravel()
    y = np.asarray(v2, dtype=np.float64).ravel()
    if x.size <= 1 or np.std(x) == 0.0 or np.std(y) == 0.0:
        return 0.0
    return float(np.corrcoef(x, y)[0, 1])


def ols_fit(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    """Ordinary least squares fit matching R's coef(lsfit(x, y)) -> (intercept, slope)."""
    xv = np.asarray(x, dtype=np.float64).ravel()
    yv = np.asarray(y, dtype=np.float64).ravel()
    design = np.column_stack([np.ones_like(xv), xv])
    coeffs, _, _, _ = np.linalg.lstsq(design, yv, rcond=None)
    return float(coeffs[0]), float(coeffs[1])


def rlm_huber(
    x: np.ndarray,
    y: np.ndarray,
    k: float = 1.345,
    max_iter: int = 20,
    acc: float = 1e-4,
) -> tuple[float, float, float]:
    """Robust linear model using Huber M-estimation matching MASS::rlm(y ~ x) and r2wls."""
    xv = np.asarray(x, dtype=np.float64).ravel()
    yv = np.asarray(y, dtype=np.float64).ravel()
    design = np.column_stack([np.ones_like(xv), xv])

    coeffs, _, _, _ = np.linalg.lstsq(design, yv, rcond=None)
    fitted = design @ coeffs
    resid = yv - fitted
    weights = np.ones_like(yv)

    for _ in range(max_iter):
        prev_resid = resid.copy()
        scale = float(np.median(np.abs(resid)) / 0.6745)
        if scale < 1e-12:
            break
        u = np.abs(resid / scale)
        weights = np.where(u <= k, 1.0, k / np.maximum(u, 1e-12))

        sqrt_w = np.sqrt(weights)
        coeffs, _, _, _ = np.linalg.lstsq(
            design * sqrt_w[:, np.newaxis],
            yv * sqrt_w,
            rcond=None,
        )
        fitted = design @ coeffs
        resid = yv - fitted

        denom = np.sum(prev_resid**2)
        if denom > 0 and np.sqrt(np.sum((resid - prev_resid) ** 2) / denom) < acc:
            break

    sse = float(np.sum(weights * (resid**2)))
    observed = resid + fitted
    w_mean = float(np.average(observed, weights=weights))
    sst = float(np.sum(weights * ((observed - w_mean) ** 2)))
    r2 = max(0.0, 1.0 - sse / sst) if sst > 0 else 0.0
    return float(coeffs[0]), float(coeffs[1]), float(np.sqrt(r2))


def infer_sex_gmm_1d(
    y_scores: np.ndarray,
    yfrac: float | None = None,
    fallback_threshold: float = 0.01,
) -> tuple[list[str], float]:
    """Classify samples into female ('F') and male ('M') via explicit threshold or 1D 2-GMM."""
    scores = np.asarray(y_scores, dtype=np.float64).ravel()
    if yfrac is not None:
        thresh = float(yfrac)
        return ["M" if float(s) >= thresh else "F" for s in scores], thresh
    if scores.size < 2 or float(np.max(scores) - np.min(scores)) < 1e-4:
        return [
            "M" if float(s) >= fallback_threshold else "F" for s in scores
        ], fallback_threshold

    mu0 = float(np.min(scores))
    above = scores[scores > mu0 + 0.005]
    mu1 = float(np.mean(above)) if above.size > 0 else float(np.max(scores))
    var0 = 1e-6
    var1 = max(
        float(np.var(above)) if above.size > 1 else float(np.var(scores)) * 0.50,
        1e-5,
    )
    pi0 = 0.5

    for _ in range(100):
        prev_mu0, prev_mu1 = mu0, mu1
        p0 = pi0 * np.exp(-0.5 * ((scores - mu0) ** 2) / var0) / math.sqrt(var0)
        p1 = (1.0 - pi0) * np.exp(-0.5 * ((scores - mu1) ** 2) / var1) / math.sqrt(var1)
        denom = np.maximum(p0 + p1, 1e-300)
        gamma1 = p1 / denom
        gamma0 = 1.0 - gamma1

        n0 = float(np.sum(gamma0))
        n1 = float(np.sum(gamma1))
        if n0 < 1e-3 or n1 < 1e-3:
            break

        pi0 = float(np.clip(n0 / scores.size, 0.05, 0.95))
        mu0 = float(np.sum(gamma0 * scores) / n0)
        mu1 = float(np.sum(gamma1 * scores) / n1)
        var0 = max(float(np.sum(gamma0 * ((scores - mu0) ** 2)) / n0), 1e-6)
        var1 = max(float(np.sum(gamma1 * ((scores - mu1) ** 2)) / n1), 1e-6)

        if abs(mu0 - prev_mu0) < 1e-7 and abs(mu1 - prev_mu1) < 1e-7:
            break

    if mu0 > mu1:
        mu0, mu1 = mu1, mu0
        var0, var1 = var1, var0
        pi0 = 1.0 - pi0

    if mu1 - mu0 < 0.005:
        return [
            "M" if float(s) >= fallback_threshold else "F" for s in scores
        ], fallback_threshold

    grid = np.linspace(mu0, mu1, 500)
    d0 = pi0 * np.exp(-0.5 * ((grid - mu0) ** 2) / var0) / math.sqrt(var0)
    d1 = (1.0 - pi0) * np.exp(-0.5 * ((grid - mu1) ** 2) / var1) / math.sqrt(var1)
    threshold = float(grid[int(np.argmin(np.abs(d0 - d1)))])
    return ["M" if float(s) >= threshold else "F" for s in scores], threshold


def fit_lasso_fista(
    x_mat: np.ndarray,
    y_prior: np.ndarray,
    lam: float = 1e-4,
    max_iter: int = 500,
    tol: float = 1e-6,
) -> tuple[float, np.ndarray]:
    """Fit non-negative Lasso regression on standardized chrY bins matching PMC6771918 Eq. 3."""
    n_samples, n_bins = x_mat.shape
    y_frac = np.asarray(y_prior, dtype=np.float64).ravel()
    y_mean = float(np.mean(y_frac))
    y_centered = y_frac - y_mean

    x_clean = np.nan_to_num(x_mat, nan=0.0)
    x_mean = np.mean(x_clean, axis=0)
    x_std_dev = np.std(x_clean, axis=0)
    active = x_std_dev > 1e-8

    if not np.any(active) or float(np.std(y_frac)) < 1e-8:
        return y_mean, np.zeros(n_bins, dtype=np.float64)

    x_act = (x_clean[:, active] - x_mean[active]) / x_std_dev[active]
    m_active = int(np.sum(active))

    v = np.ones(m_active, dtype=np.float64) / math.sqrt(m_active)
    for _ in range(15):
        v_next = x_act.T @ (x_act @ v)
        norm_v = float(np.linalg.norm(v_next))
        if norm_v <= 1e-12:
            break
        v = v_next / norm_v
    lipschitz = max(float(np.sum((x_act @ v) ** 2)) / max(1, n_samples), 1e-6)
    step = 1.0 / lipschitz
    effective_lam = lam * max(1.0, math.sqrt(m_active / max(1, n_samples)))

    beta = np.zeros(m_active, dtype=np.float64)
    z = np.zeros(m_active, dtype=np.float64)
    t = 1.0

    for _ in range(max_iter):
        grad = (x_act.T @ (x_act @ z - y_centered)) / n_samples
        beta_next = np.maximum(0.0, z - step * grad - step * effective_lam)
        max_diff = float(np.max(np.abs(beta_next - beta)))
        t_next = 0.5 * (1.0 + math.sqrt(1.0 + 4.0 * t * t))
        z = beta_next + ((t - 1.0) / t_next) * (beta_next - beta)
        beta = beta_next
        t = t_next
        if max_diff < tol:
            break

    scaled_beta = beta / x_std_dev[active]
    intercept = y_mean - float(np.dot(x_mean[active], scaled_beta))
    full_weights = np.zeros(n_bins, dtype=np.float64)
    full_weights[active] = scaled_beta
    return intercept, full_weights
