"""Statistical and linear algebra routines matching PREFACE.R (PCA, OLS, Huber RLM)."""

from __future__ import annotations

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
