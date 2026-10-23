"""Hardware discovery, RNG seeding, and execution provider selection."""

from __future__ import annotations

import os
import random

import numpy as np
import onnxruntime as ort
import torch

from preface.logging_utils import get_logger

logger = get_logger(__name__)


def resolve_device(requested: str = "auto") -> torch.device:
    """Resolve the best available PyTorch compute device."""
    req = requested.strip().lower()
    if req == "auto":
        if torch.cuda.is_available():
            dev = torch.device("cuda")
        elif torch.backends.mps.is_available():
            dev = torch.device("mps")
        else:
            dev = torch.device("cpu")
    elif req == "cuda":
        if not torch.cuda.is_available():
            raise ValueError("Requested device 'cuda' is not available on this system.")
        dev = torch.device("cuda")
    elif req == "mps":
        if not torch.backends.mps.is_available():
            raise ValueError("Requested device 'mps' is not available on this system.")
        dev = torch.device("mps")
    elif req == "cpu":
        dev = torch.device("cpu")
    else:
        raise ValueError(
            f"Unsupported device '{requested}'. Choose from: auto, cpu, cuda, mps."
        )
    logger.debug("Resolved compute device '%s' -> %s", requested, dev)
    return dev


def configure_hardware(cpus: int, seed: int = 1) -> int:
    """Configure thread pools, matmul precision, and deterministic seeds."""
    resolved_cpus = max(1, cpus if cpus > 0 else (os.cpu_count() or 1))
    os.environ["POLARS_MAX_THREADS"] = str(resolved_cpus)
    os.environ["OMP_NUM_THREADS"] = str(resolved_cpus)
    torch.set_num_threads(resolved_cpus)
    torch.set_float32_matmul_precision("high")

    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)

    logger.debug("Configured hardware: cpus=%d, seed=%d", resolved_cpus, seed)
    return resolved_cpus


def get_onnx_providers(requested_device: str = "auto") -> list[str]:
    """Return ordered ONNX Runtime execution providers matching the requested device."""
    available = set(ort.get_available_providers())
    req = requested_device.strip().lower()

    preferred: list[str] = []
    if req in ("auto", "cuda") and "CUDAExecutionProvider" in available:
        preferred.append("CUDAExecutionProvider")
    if req in ("auto", "mps") and "CoreMLExecutionProvider" in available:
        preferred.append("CoreMLExecutionProvider")
    preferred.append("CPUExecutionProvider")
    logger.debug(
        "Selected ONNX Runtime providers for '%s': %s", requested_device, preferred
    )
    return preferred
