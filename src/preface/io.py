"""Data loading and preprocessing for PREFACE config and copy number alteration .bed files."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import polars as pl

from preface import EXCLUDE_CHRS

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
        raise ValueError(f"Config file '{config_path}' is missing required columns: {missing}")

    perm = rng.permutation(df.height)
    return df[perm]


def _read_bed_table(bed_path: Path) -> pl.DataFrame:
    if not bed_path.is_file():
        raise FileNotFoundError(f"The file '{bed_path}' does not exist.")
    df = pl.read_csv(
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
    return df


def _load_single_sample_ratios(args: tuple[int, int, str, str]) -> np.ndarray:
    idx, total, sample_id, filepath = args
    remaining = total - (idx + 1)
    print(f"Loading sample {sample_id} | {remaining}/{total} remaining ...")
    df = _read_bed_table(Path(filepath))
    filtered = df.filter(pl.col("chr") != "Y")
    return filtered["ratio"].to_numpy().astype(np.float64, copy=False)


def load_training_dataset(
    config_df: pl.DataFrame,
    cpus: int,
) -> TrainingDataset:
    """Load all training BED files concurrently, compute X ratios, filter NAs, and mean-impute."""
    n_samples = config_df.height
    first_path = Path(config_df["filepath"][0])
    ref_df = _read_bed_table(first_path).filter(pl.col("chr") != "Y")

    chrs = ref_df["chr"].to_list()
    starts = ref_df["start"].to_list()
    ends = ref_df["end"].to_list()

    tasks = [
        (i, n_samples, str(config_df["ID"][i]), str(config_df["filepath"][i]))
        for i in range(n_samples)
    ]

    with ThreadPoolExecutor(max_workers=max(1, cpus)) as pool:
        sample_columns = list(pool.map(_load_single_sample_ratios, tasks))

    # Shape: (n_loci_non_Y, n_samples)
    sub_matrix = np.column_stack(sample_columns)

    x_mask = np.array([c == "X" for c in chrs], dtype=bool)
    with np.errstate(invalid="ignore"):
        x_means = np.nanmean(sub_matrix[x_mask, :], axis=0)
    x_ratios = np.power(2.0, x_means)

    print("Creating training frame ...")

    keep_chr_mask = np.array([c not in EXCLUDE_CHRS for c in chrs], dtype=bool)
    feature_names_all = [
        f"{chrs[i]}:{starts[i]}-{ends[i]}"
        for i in range(len(chrs))
        if keep_chr_mask[i]
    ]

    # Transpose to shape: (n_samples, n_autosomal_features)
    training_frame = np.ascontiguousarray(sub_matrix[keep_chr_mask, :].T, dtype=np.float64)

    na_counts = np.isnan(training_frame).sum(axis=0)
    valid_feature_mask = na_counts < (n_samples * 0.01)

    training_frame = training_frame[:, valid_feature_mask]
    possible_features = [
        name for name, keep in zip(feature_names_all, valid_feature_mask) if keep
    ]
    mean_features = np.nanmean(training_frame, axis=0)

    na_rows, na_cols = np.where(np.isnan(training_frame))
    if na_cols.size > 0:
        training_frame[na_rows, na_cols] = mean_features[na_cols]

    print(f"Remaining training features after 'NA' filtering: {len(possible_features)}")

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
    x_df = df.filter(pl.col("chr") == "X")
    x_vals = x_df["ratio"].to_numpy().astype(np.float64, copy=False)
    x_ratio = float(2.0 ** np.nanmean(x_vals))

    feat_ids = (
        df["chr"] + ":" + df["start"].cast(pl.String) + "-" + df["end"].cast(pl.String)
    ).to_list()
    ratios = df["ratio"].to_numpy().astype(np.float64, copy=False)

    locus_map = dict(zip(feat_ids, ratios))
    selected_ratios = np.array(
        [locus_map.get(feat, np.nan) for feat in possible_features],
        dtype=np.float64,
    )

    return SampleBedData(x_ratio=x_ratio, feature_ratios=selected_ratios)
