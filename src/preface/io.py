"""Data loading and preprocessing for PREFACE config and copy number alteration .bed files."""

from __future__ import annotations

import json
import os
from collections.abc import Callable, Sequence
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import polars as pl

from preface import EXCLUDE_CHRS
from preface.logging_utils import get_logger

logger = get_logger(__name__)

REQUIRED_CONFIG_COLS: tuple[str, ...] = ("ID", "filepath", "gender", "FF")
REQUIRED_BED_COLS: tuple[str, ...] = ("chr", "start", "end", "ratio")
NA_VALUES: list[str] = ["NA", "NaN", "nan", ""]


@dataclass(frozen=True)
class TrainingDataset:
    """Preprocessed training matrix and sample metadata."""

    config_df: pl.DataFrame
    training_frame: np.ndarray
    x_ratios: np.ndarray
    possible_features: list[str]
    mean_features: np.ndarray


@dataclass(frozen=True)
class SampleBedData:
    """Parsed single-sample BED data for prediction."""

    x_ratio: float
    feature_ratios: np.ndarray


def _locus_id_expr() -> pl.Expr:
    return (
        pl.col("chr")
        + ":"
        + pl.col("start").cast(pl.String)
        + "-"
        + pl.col("end").cast(pl.String)
    )


def load_config(config_path: Path, rng: np.random.Generator) -> pl.DataFrame:
    """Load and deterministically shuffle the training config TSV file."""
    df = pl.read_csv(
        config_path,
        separator="\t",
        null_values=NA_VALUES,
        schema_overrides={
            "ID": pl.String,
            "filepath": pl.String,
            "gender": pl.String,
            "FF": pl.Float64,
        },
    )
    missing = [c for c in REQUIRED_CONFIG_COLS if c not in df.columns]
    if missing:
        raise ValueError(
            f"Config file '{config_path}' is missing required columns: {missing}"
        )

    return df[rng.permutation(df.height)]


def _read_bed_table(bed_path: Path) -> pl.DataFrame:
    if not bed_path.is_file():
        raise FileNotFoundError(f"The file '{bed_path}' does not exist.")
    return pl.read_csv(
        bed_path,
        separator="\t",
        null_values=NA_VALUES,
        schema_overrides={
            "chr": pl.String,
            "start": pl.Int64,
            "end": pl.Int64,
            "ratio": pl.Float64,
        },
        columns=list(REQUIRED_BED_COLS),
    )


def _init_loader_worker() -> None:
    """Restrict per-worker thread pools to prevent oversubscription across processes."""
    os.environ["POLARS_MAX_THREADS"] = "1"
    os.environ["OMP_NUM_THREADS"] = "1"


def parallel_load_samples[T, R](
    items: Sequence[tuple[str, T]],
    worker_fn: Callable[[T], R],
    cpus: int,
) -> list[R]:
    """Load sample files in parallel across worker processes while preserving input order."""
    total = len(items)
    if total == 0:
        return []

    max_workers = min(max(1, cpus), total)
    if max_workers == 1:
        results: list[R] = []
        for done_count, (sample_id, payload) in enumerate(items, start=1):
            logger.info(
                "Loading sample %s | %d/%d remaining ...",
                sample_id,
                total - done_count,
                total,
            )
            results.append(worker_fn(payload))
        return results

    ordered_results: dict[int, R] = {}
    with ProcessPoolExecutor(
        max_workers=max_workers,
        initializer=_init_loader_worker,
    ) as pool:
        future_to_meta = {
            pool.submit(worker_fn, payload): (idx, sample_id)
            for idx, (sample_id, payload) in enumerate(items)
        }
        done_count = 0
        try:
            for future in as_completed(future_to_meta):
                idx, sample_id = future_to_meta[future]
                ordered_results[idx] = future.result()
                done_count += 1
                logger.info(
                    "Loading sample %s | %d/%d remaining ...",
                    sample_id,
                    total - done_count,
                    total,
                )
        except Exception:
            for pending in future_to_meta:
                pending.cancel()
            raise

    return [ordered_results[i] for i in range(total)]


def _load_single_sample_bed(
    args: tuple[str, bool],
) -> tuple[np.ndarray, list[str] | None, list[str] | None]:
    filepath, include_meta = args
    df = _read_bed_table(Path(filepath)).filter(pl.col("chr") != "Y")
    ratios = df["ratio"].to_numpy().astype(np.float64, copy=False)
    if not include_meta:
        return ratios, None, None

    excluded = list(EXCLUDE_CHRS)
    chrs = df["chr"].to_list()
    feature_loci = (
        df.filter(~pl.col("chr").is_in(excluded))
        .select(_locus_id_expr().alias("locus"))["locus"]
        .to_list()
    )
    return ratios, chrs, feature_loci


def load_training_dataset(
    config_df: pl.DataFrame,
    cpus: int,
) -> TrainingDataset:
    """Load all training BED files concurrently, compute X ratios, filter NAs, and mean-impute."""
    n_samples = config_df.height
    items = [
        (
            str(config_df["ID"][i]),
            (str(config_df["filepath"][i]), i == 0),
        )
        for i in range(n_samples)
    ]
    loaded = parallel_load_samples(items, _load_single_sample_bed, cpus=cpus)

    sample_columns = [r for r, _, _ in loaded]
    ref_chrs = loaded[0][1]
    feature_names_all = loaded[0][2]
    assert ref_chrs is not None and feature_names_all is not None

    sub_matrix = np.column_stack(sample_columns)

    x_mask = np.array([c == "X" for c in ref_chrs], dtype=bool)
    with np.errstate(invalid="ignore"):
        x_ratios = np.power(2.0, np.nanmean(sub_matrix[x_mask, :], axis=0))

    logger.info("Creating training frame ...")

    keep_chr_mask = np.array([c not in EXCLUDE_CHRS for c in ref_chrs], dtype=bool)

    training_frame = np.ascontiguousarray(
        sub_matrix[keep_chr_mask, :].T, dtype=np.float64
    )
    valid_feature_mask = np.isnan(training_frame).sum(axis=0) < (n_samples * 0.01)

    training_frame = training_frame[:, valid_feature_mask]
    possible_features = [
        name for name, keep in zip(feature_names_all, valid_feature_mask) if keep
    ]
    mean_features = np.nanmean(training_frame, axis=0)

    na_rows, na_cols = np.where(np.isnan(training_frame))
    if na_cols.size > 0:
        training_frame[na_rows, na_cols] = mean_features[na_cols]

    logger.info(
        "Remaining training features after 'NA' filtering: %d", len(possible_features)
    )

    return TrainingDataset(
        config_df=config_df,
        training_frame=training_frame,
        x_ratios=x_ratios,
        possible_features=possible_features,
        mean_features=mean_features,
    )


def load_inference_bed(bed_path: Path, possible_features: list[str]) -> SampleBedData:
    """Parse a single BED file for inference, returning X.ratio and ordered feature ratios."""
    df = _read_bed_table(bed_path)
    x_vals = (
        df.filter(pl.col("chr") == "X")["ratio"]
        .to_numpy()
        .astype(np.float64, copy=False)
    )
    x_ratio = float(2.0 ** np.nanmean(x_vals))

    feat_ids = df.select(_locus_id_expr().alias("locus"))["locus"].to_list()
    ratios = df["ratio"].to_numpy().astype(np.float64, copy=False)

    locus_map = dict(zip(feat_ids, ratios))
    selected_ratios = np.array(
        [locus_map.get(feat, np.nan) for feat in possible_features],
        dtype=np.float64,
    )

    return SampleBedData(x_ratio=x_ratio, feature_ratios=selected_ratios)


def emit_prediction_output(
    text_lines: list[str],
    json_dict: dict[str, str | float],
    json_arg: str | None,
) -> None:
    """Write single-sample prediction results as human-readable lines or JSON to stdout/file."""
    if json_arg is None:
        for line in text_lines:
            print(line)
        return
    payload = json.dumps(json_dict)
    if json_arg == "":
        print(payload)
    else:
        Path(json_arg).write_text(payload, encoding="utf-8")
