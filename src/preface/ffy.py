"""Baseline FFY (chromosome Y fetal fraction) calculation per PMC6771918 Section 3.4."""

from __future__ import annotations

import importlib.util
import shutil
import struct
import subprocess
import zlib
from dataclasses import dataclass
from pathlib import Path
from typing import BinaryIO

import numpy as np
import polars as pl

from preface import EXCLUDE_CHRS
from preface.io import (
    NA_VALUES,
    emit_prediction_output,
    open_text_file,
    parallel_load_samples,
)
from preface.logging_utils import get_logger

logger = get_logger(__name__)

AUTOSOME_CHRS: frozenset[str] = frozenset(str(i) for i in range(1, 23)) - EXCLUDE_CHRS
CHRY_PAR1_END: int = 2_781_479
CHRY_HETEROCHROMATIN_START: int = 28_800_000
SKIP_BAM_FLAGS: int = 0x4 | 0x100 | 0x200 | 0x400 | 0x800
ALIGNMENT_SUFFIXES: frozenset[str] = frozenset({".bam", ".cram", ".sam"})


@dataclass(frozen=True)
class SampleYProfile:
    """Normalized chromosome Y bin depths and coordinates for a single sample."""

    bin_loci: list[str]
    bin_starts: np.ndarray
    bin_ends: np.ndarray
    y_norm_bins: np.ndarray


def _normalize_chr(raw: str) -> str:
    cleaned = raw.strip()
    if cleaned.lower().startswith("chr"):
        cleaned = cleaned[3:]
    upper = cleaned.upper()
    if upper in {"X", "Y", "M", "MT"}:
        return "M" if upper == "MT" else upper
    return cleaned


def _coordinate_informative_mask(starts: np.ndarray, ends: np.ndarray) -> np.ndarray:
    """Mask out PAR1, PAR2, and Yq12 heterochromatin when coordinates span human chrY."""
    if ends.size == 0:
        return np.zeros(0, dtype=bool)
    if int(np.max(ends)) > CHRY_PAR1_END:
        mask = (starts >= CHRY_PAR1_END) & (ends <= CHRY_HETEROCHROMATIN_START)
        if np.any(mask):
            return mask
    return np.ones(starts.size, dtype=bool)


def _compute_y_profile_from_bins(
    auto_starts: np.ndarray,
    auto_ends: np.ndarray,
    auto_counts: np.ndarray,
    y_starts: np.ndarray,
    y_ends: np.ndarray,
    y_counts: np.ndarray,
) -> SampleYProfile:
    """Normalize chromosome Y bin counts by mean mappable diploid autosomal read density."""
    if y_counts.size == 0:
        raise ValueError("Input contains no chromosome Y bins.")
    if np.any(auto_counts < 0) or np.any(y_counts < 0):
        raise ValueError(
            "Input contains negative bin values. Provide raw/binned read counts rather than log2 ratios."
        )

    auto_widths = np.maximum(1, auto_ends - auto_starts).astype(np.float64)
    valid_auto = np.isfinite(auto_counts) & (auto_counts > 0)
    if not np.any(valid_auto):
        raise ValueError(
            "Input contains no positive autosomal read counts for normalization."
        )

    r_auto = float(np.mean(auto_counts[valid_auto] / auto_widths[valid_auto]))
    if r_auto <= 0.0 or not np.isfinite(r_auto):
        raise ValueError("Computed autosomal read density is non-positive.")

    order = np.argsort(y_starts, kind="mergesort")
    y_starts_sorted = y_starts[order].astype(np.int64, copy=False)
    y_ends_sorted = y_ends[order].astype(np.int64, copy=False)
    y_counts_sorted = y_counts[order].astype(np.float64, copy=False)

    y_widths = np.maximum(1, y_ends_sorted - y_starts_sorted).astype(np.float64)
    y_norm = np.where(
        np.isfinite(y_counts_sorted),
        y_counts_sorted / (r_auto * y_widths),
        np.nan,
    )
    bin_loci = [f"Y:{int(s)}-{int(e)}" for s, e in zip(y_starts_sorted, y_ends_sorted)]

    return SampleYProfile(
        bin_loci=bin_loci,
        bin_starts=y_starts_sorted,
        bin_ends=y_ends_sorted,
        y_norm_bins=y_norm,
    )


def _pick_col(
    col_map: dict[str, str], candidates: tuple[str, ...], fallback: str
) -> str:
    return next((col_map[k] for k in candidates if k in col_map), fallback)


def _read_binned_counts_bed(bed_path: Path) -> SampleYProfile:
    """Parse a headered or headerless 4-column binned read count BED/TSV file."""
    first_data_line = ""
    skip_lines = 0
    with open_text_file(bed_path) as fh:
        for line in fh:
            stripped = line.strip()
            if not stripped:
                skip_lines += 1
                continue
            if stripped.startswith("#"):
                clean = stripped.lstrip("#").strip().lower()
                if clean.startswith(("chr", "chrom", "chromosome")):
                    first_data_line = stripped
                    break
                skip_lines += 1
                continue
            first_data_line = stripped
            break

    if not first_data_line:
        raise ValueError(f"Binned read count file '{bed_path}' is empty.")

    parts = first_data_line.split("\t")
    if len(parts) < 4:
        parts = first_data_line.split()
    if len(parts) < 4:
        raise ValueError(
            f"File '{bed_path}' must have at least 4 tab-separated columns (chr, start, end, count)."
        )

    try:
        int(parts[1])
        int(parts[2])
        has_header = False
    except ValueError:
        has_header = True

    df = pl.read_csv(
        bed_path,
        separator="\t",
        skip_lines=skip_lines,
        has_header=has_header,
        null_values=NA_VALUES,
        infer_schema_length=1000,
    )
    if has_header:
        col_map = {c.strip().lstrip("#").lower(): c for c in df.columns}
        cols = [
            _pick_col(col_map, ("chr", "chrom", "chromosome"), df.columns[0]),
            _pick_col(col_map, ("start", "chromstart"), df.columns[1]),
            _pick_col(col_map, ("end", "chromend"), df.columns[2]),
            _pick_col(
                col_map,
                (
                    "count",
                    "counts",
                    "reads",
                    "read_count",
                    "depth",
                    "coverage",
                    "ratio",
                ),
                df.columns[3],
            ),
        ]
    else:
        cols = df.columns[:4]

    df = df.select(
        pl.col(cols[0]).cast(pl.String).alias("chr"),
        pl.col(cols[1]).cast(pl.Int64).alias("start"),
        pl.col(cols[2]).cast(pl.Int64).alias("end"),
        pl.col(cols[3]).cast(pl.Float64).alias("count"),
    ).with_columns(pl.Series("chr", [_normalize_chr(c) for c in df[cols[0]].to_list()]))

    auto_df = df.filter(pl.col("chr").is_in(list(AUTOSOME_CHRS)))
    y_df = df.filter(pl.col("chr") == "Y")

    return _compute_y_profile_from_bins(
        auto_starts=auto_df["start"].to_numpy(),
        auto_ends=auto_df["end"].to_numpy(),
        auto_counts=auto_df["count"].to_numpy(),
        y_starts=y_df["start"].to_numpy(),
        y_ends=y_df["end"].to_numpy(),
        y_counts=y_df["count"].to_numpy(),
    )


def _iter_bgzf_blocks(fh: BinaryIO) -> bytes:
    """Decompress all BGZF/gzip blocks from a BAM file stream."""
    chunks: list[bytes] = []
    while True:
        header = fh.read(12)
        if not header:
            break
        if len(header) < 12:
            raise ValueError("Truncated BGZF block header in BAM file.")
        id1, id2, cm, flg, _mtime, _xfl, _os, xlen = struct.unpack("<BBBBIBBH", header)
        if id1 != 0x1F or id2 != 0x8B or cm != 8:
            raise ValueError("Invalid gzip/BGZF magic bytes in BAM file.")
        extra = fh.read(xlen) if (flg & 0x04) else b""
        bsize: int | None = None
        idx = 0
        while idx + 4 <= len(extra):
            si1, si2, slen = struct.unpack_from("<BBH", extra, idx)
            idx += 4
            if si1 == 66 and si2 == 67 and slen == 2 and idx + 2 <= len(extra):
                bsize = struct.unpack_from("<H", extra, idx)[0]
            idx += slen

        if bsize is None:
            raise ValueError("BAM file is missing BGZF BC subfield.")
        cdata_len = (bsize + 1) - 12 - xlen - 8
        if cdata_len < 0:
            raise ValueError("Corrupt BGZF block size in BAM file.")
        cdata = fh.read(cdata_len)
        trailer = fh.read(8)
        if len(cdata) < cdata_len or len(trailer) < 8:
            raise ValueError("Truncated BGZF block payload in BAM file.")
        _crc, isize = struct.unpack("<II", trailer)
        if isize > 0:
            chunks.append(zlib.decompress(cdata, -zlib.MAX_WBITS))
    return b"".join(chunks)


def _init_ref_bins(
    refs: list[tuple[str, int]], binsize: int
) -> tuple[dict[int, str], dict[str, int], dict[str, np.ndarray]]:
    idx_to_chr: dict[int, str] = {}
    ref_lengths: dict[str, int] = {}
    bin_counts: dict[str, np.ndarray] = {}
    for idx, (raw_name, length) in enumerate(refs):
        norm_chr = _normalize_chr(raw_name)
        if (norm_chr in AUTOSOME_CHRS or norm_chr == "Y") and length > 0:
            idx_to_chr[idx] = norm_chr
            ref_lengths[norm_chr] = length
            bin_counts[norm_chr] = np.zeros(
                (length + binsize - 1) // binsize, dtype=np.int64
            )
    return idx_to_chr, ref_lengths, bin_counts


def _record_hit(
    bin_counts: dict[str, np.ndarray],
    chrom: str | None,
    pos0: int,
    binsize: int,
) -> None:
    if chrom is not None and pos0 >= 0:
        arr = bin_counts[chrom]
        bin_idx = pos0 // binsize
        if bin_idx < arr.size:
            arr[bin_idx] += 1


def _finalize_bam_bin_counts(
    ref_lengths: dict[str, int],
    bin_counts: dict[str, np.ndarray],
    binsize: int,
) -> SampleYProfile:
    """Convert per-chromosome bin count arrays from aligned reads into a SampleYProfile."""
    if "Y" not in ref_lengths or "Y" not in bin_counts:
        raise ValueError(
            "Alignment header does not contain chromosome Y ('Y' or 'chrY')."
        )

    auto_starts_list: list[np.ndarray] = []
    auto_ends_list: list[np.ndarray] = []
    auto_counts_list: list[np.ndarray] = []

    for chrom in sorted(AUTOSOME_CHRS):
        if chrom not in ref_lengths:
            continue
        rlen = ref_lengths[chrom]
        starts = np.arange(0, rlen, binsize, dtype=np.int64)
        auto_starts_list.append(starts + 1)
        auto_ends_list.append(np.minimum(starts + binsize, rlen))
        auto_counts_list.append(bin_counts[chrom].astype(np.float64))

    if not auto_counts_list:
        raise ValueError(
            "Alignment file contains no reference autosomes for normalization."
        )

    y_len = ref_lengths["Y"]
    y_starts = np.arange(0, y_len, binsize, dtype=np.int64)
    return _compute_y_profile_from_bins(
        auto_starts=np.concatenate(auto_starts_list),
        auto_ends=np.concatenate(auto_ends_list),
        auto_counts=np.concatenate(auto_counts_list),
        y_starts=y_starts + 1,
        y_ends=np.minimum(y_starts + binsize, y_len),
        y_counts=bin_counts["Y"].astype(np.float64),
    )


def _read_bam_native(bam_path: Path, binsize: int, mapq: int) -> SampleYProfile:
    """Read and bin a BAM file using Python's built-in zlib and struct modules."""
    with bam_path.open("rb") as fh:
        raw = _iter_bgzf_blocks(fh)

    if len(raw) < 12 or raw[:4] != b"BAM\x01":
        raise ValueError(f"File '{bam_path}' is not a valid BAM file.")

    l_text = struct.unpack_from("<i", raw, 4)[0]
    offset = 8 + l_text
    n_ref = struct.unpack_from("<i", raw, offset)[0]
    offset += 4

    refs: list[tuple[str, int]] = []
    for _ in range(n_ref):
        l_name = struct.unpack_from("<i", raw, offset)[0]
        offset += 4
        name = raw[offset : offset + l_name].rstrip(b"\x00").decode("ascii")
        offset += l_name
        l_ref = struct.unpack_from("<i", raw, offset)[0]
        offset += 4
        refs.append((name, l_ref))

    ref_id_to_chr, ref_lengths, bin_counts = _init_ref_bins(refs, binsize)
    total_len = len(raw)
    while offset + 4 <= total_len:
        block_size = struct.unpack_from("<i", raw, offset)[0]
        rec_start = offset + 4
        offset = rec_start + block_size
        if offset > total_len or block_size < 16:
            break
        ref_id, pos0, _lq, rec_mapq, _bin, _nc, flag = struct.unpack_from(
            "<iiBBHHH", raw, rec_start
        )
        if (flag & SKIP_BAM_FLAGS) == 0 and rec_mapq >= mapq:
            _record_hit(bin_counts, ref_id_to_chr.get(ref_id), pos0, binsize)

    return _finalize_bam_bin_counts(ref_lengths, bin_counts, binsize)


def _read_alignment_external(
    aln_path: Path, binsize: int, mapq: int, reference: Path | None
) -> SampleYProfile:
    """Read BAM/CRAM/SAM via pysam when available, falling back to samtools CLI."""
    if importlib.util.find_spec("pysam") is not None:
        import pysam

        ref_arg = str(reference) if reference is not None else None
        with pysam.AlignmentFile(str(aln_path), reference_filename=ref_arg) as samfile:
            refs = [
                (samfile.get_reference_name(i), int(samfile.lengths[i]))
                for i in range(samfile.nreferences)
            ]
            tid_to_chr, ref_lengths, bin_counts = _init_ref_bins(refs, binsize)
            for read in samfile.fetch(until_eof=True):
                if (read.flag & SKIP_BAM_FLAGS) == 0 and read.mapping_quality >= mapq:
                    _record_hit(
                        bin_counts,
                        tid_to_chr.get(read.reference_id),
                        read.reference_start,
                        binsize,
                    )
        return _finalize_bam_bin_counts(ref_lengths, bin_counts, binsize)

    samtools_bin = shutil.which("samtools")
    if samtools_bin is None:
        raise RuntimeError(
            f"Cannot read '{aln_path}': install 'pysam' or 'samtools', or provide a '.bam' or binned '.bed' file."
        )
    ref_flags = ["-T", str(reference)] if reference is not None else []
    hdr_out = subprocess.run(
        [samtools_bin, "view", "-H", *ref_flags, str(aln_path)],
        capture_output=True,
        text=True,
        check=True,
    ).stdout
    refs_sam: list[tuple[str, int]] = []
    for line in hdr_out.splitlines():
        if line.startswith("@SQ"):
            fields = {
                f[:2]: f[3:] for f in line.split("\t")[1:] if f[:3] in {"SN:", "LN:"}
            }
            if "SN" in fields and "LN" in fields:
                refs_sam.append((fields["SN"], int(fields["LN"])))

    idx_to_chr, ref_lengths, bin_counts = _init_ref_bins(refs_sam, binsize)
    raw_to_chr = {
        name: idx_to_chr[i] for i, (name, _) in enumerate(refs_sam) if i in idx_to_chr
    }
    view_cmd = [
        samtools_bin,
        "view",
        "-F",
        str(SKIP_BAM_FLAGS),
        "-q",
        str(mapq),
        *ref_flags,
        str(aln_path),
    ]
    with subprocess.Popen(view_cmd, stdout=subprocess.PIPE, text=True) as proc:
        assert proc.stdout is not None
        for line in proc.stdout:
            parts = line.split("\t", 5)
            if len(parts) >= 4:
                _record_hit(
                    bin_counts, raw_to_chr.get(parts[2]), int(parts[3]) - 1, binsize
                )
        if proc.wait() != 0:
            raise RuntimeError(f"samtools view failed on '{aln_path}'.")

    return _finalize_bam_bin_counts(ref_lengths, bin_counts, binsize)


def extract_sample_y_profile(
    input_path: Path,
    binsize: int = 5000,
    mapq: int = 1,
    reference: Path | None = None,
) -> SampleYProfile:
    """Extract diploid-normalized chromosome Y bin profile from an alignment or binned BED file."""
    if not input_path.is_file():
        raise FileNotFoundError(f"The file '{input_path}' does not exist.")
    if binsize <= 0:
        raise ValueError("--binsize must be a positive integer.")

    suffix = input_path.suffix.lower()
    if suffix in ALIGNMENT_SUFFIXES:
        if suffix == ".bam" and importlib.util.find_spec("pysam") is None:
            profile = _read_bam_native(input_path, binsize=binsize, mapq=mapq)
        else:
            profile = _read_alignment_external(
                input_path, binsize=binsize, mapq=mapq, reference=reference
            )
    else:
        profile = _read_binned_counts_bed(input_path)

    logger.debug(
        "Extracted chrY profile from '%s': %d chrY bins",
        input_path,
        len(profile.bin_loci),
    )
    return profile


def _align_profile_to_loci(
    profile: SampleYProfile, target_loci: list[str]
) -> np.ndarray:
    if profile.bin_loci == target_loci:
        return np.nan_to_num(profile.y_norm_bins, nan=0.0)
    locus_map = dict(zip(profile.bin_loci, profile.y_norm_bins))
    arr = np.array([locus_map.get(loc, 0.0) for loc in target_loci], dtype=np.float64)
    return np.nan_to_num(arr, nan=0.0)


def _extract_profile_worker(
    args: tuple[Path, int, int, Path | None],
) -> SampleYProfile:
    path, binsize, mapq, reference = args
    return extract_sample_y_profile(
        path, binsize=binsize, mapq=mapq, reference=reference
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
    """Calculate gender and baseline FFY (%) for a single BAM/CRAM/BED sample."""
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
        gender = "M" if y_score >= thresh else "F"
        ffy_pct = (
            run_ffy_onnx_inference(model_path, aligned_y, device_name)
            if gender == "M"
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
        gender = "M" if y_score >= thresh else "F"
        ffy_pct = (
            float(np.clip(2.0 * y_score * 100.0, 0.0, 100.0)) if gender == "M" else 0.0
        )

    logger.debug(
        "Single-sample FFY for '%s': gender=%s, FFY=%.4f%%", infile, gender, ffy_pct
    )
    emit_prediction_output(
        [f"gender = {gender}", f"FFY = {ffy_pct:.4g}%"],
        {"FFY": ffy_pct / 100.0, "gender": gender},
        json_arg,
    )
    return gender, ffy_pct


def run_ffy_cohort(
    config_path: Path,
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
) -> pl.DataFrame:
    """Calculate baseline FFY and gender across a cohort per PMC6771918 Section 3.4."""
    from preface.hardware import configure_hardware
    from preface.math_utils import fit_lasso_fista, infer_sex_gmm_1d
    from preface.models import export_ffy_model_to_onnx

    if not config_path.is_file():
        raise FileNotFoundError(f"The file '{config_path}' does not exist.")

    resolved_cpus = configure_hardware(cpus=cpus, seed=seed)
    rng = np.random.default_rng(seed)

    df = pl.read_csv(
        config_path, separator="\t", null_values=NA_VALUES, infer_schema_length=1000
    )
    if "ID" not in df.columns or (
        "filepath" not in df.columns and "reads_filepath" not in df.columns
    ):
        raise ValueError(
            f"Cohort file '{config_path}' is missing required column 'ID' or 'filepath/reads_filepath'."
        )

    n_samples = df.height
    if n_samples == 0:
        raise ValueError(f"Cohort file '{config_path}' contains no samples.")

    reads_col = "reads_filepath" if "reads_filepath" in df.columns else "filepath"
    sample_ids = [str(s) for s in df["ID"].to_list()]
    input_paths = [Path(str(p)) for p in df[reads_col].to_list()]

    combined_items = [
        (sid, (path, binsize, mapq, reference))
        for sid, path in zip(sample_ids, input_paths)
    ]
    if male_ref_paths:
        combined_items.extend(
            (f"male_ref:{p.name}", (p, binsize, mapq, reference))
            for p in male_ref_paths
        )

    logger.info("Extracting chrY profiles for %d cohort samples ...", n_samples)
    all_profiles = parallel_load_samples(
        combined_items, _extract_profile_worker, cpus=resolved_cpus
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
    if "gender" in df.columns:
        genders = [
            str(g).strip().upper()
            if g is not None and str(g).strip().upper() in {"M", "F"}
            else inf
            for g, inf in zip(df["gender"].to_list(), inferred_all)
        ]
    else:
        genders = inferred_all

    male_mask = np.array([g == "M" for g in genders], dtype=bool)
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
        ((y_scores[male_mask] - female_baseline) / denom) * 100.0, 0.0, 100.0
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
        final_ffy = np.where(male_mask, np.clip(cv_preds, 0.0, 100.0), 0.0)
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

    out_path_col = "filepath" if "filepath" in df.columns else reads_col
    out_df = pl.DataFrame(
        {
            "ID": df["ID"].cast(pl.String),
            "filepath": df[out_path_col].cast(pl.String),
            "gender": genders,
            "FF": final_ffy,
        }
    )

    tsv_text = out_df.write_csv(separator="\t")
    if outfile is not None:
        outfile.parent.mkdir(parents=True, exist_ok=True)
        outfile.write_text(tsv_text, encoding="utf-8")
        logger.info("Wrote cohort FFY config (%d samples) to '%s'", n_samples, outfile)
    else:
        print(tsv_text, end="")

    return out_df
