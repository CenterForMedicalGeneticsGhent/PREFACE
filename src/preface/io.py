"""Samplesheet validation, alignment binning, and feature extraction for PREFACE."""

from __future__ import annotations

import importlib.util
import json
import os
import shutil
import struct
import subprocess
import zlib
from collections.abc import Callable, Sequence
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Any, BinaryIO

import jsonschema
import numpy as np
import yaml

from preface import EXCLUDE_CHRS
from preface.logging_utils import get_logger

logger = get_logger(__name__)

SAMPLESHEET_SCHEMA_PATH: Path = Path(__file__).with_name("samplesheet.schema.json")
SAMPLESHEET_SUFFIXES: frozenset[str] = frozenset({".yaml", ".yml", ".json"})
ALIGNMENT_SUFFIXES: frozenset[str] = frozenset({".bam", ".cram", ".sam"})
SKIP_BAM_FLAGS: int = 0x4 | 0x100 | 0x200 | 0x400 | 0x800

ALL_AUTOSOME_CHRS: frozenset[str] = frozenset(str(i) for i in range(1, 23))
AUTOSOME_CHRS: frozenset[str] = ALL_AUTOSOME_CHRS - EXCLUDE_CHRS
CN_TARGET_CHRS: frozenset[str] = ALL_AUTOSOME_CHRS | {"X"}
FFY_TARGET_CHRS: frozenset[str] = AUTOSOME_CHRS | {"Y"}
CANONICAL_FEATURE_CHRS: tuple[str, ...] = tuple(
    str(i) for i in range(1, 23) if str(i) not in EXCLUDE_CHRS
)


@dataclass(frozen=True)
class SampleEntry:
    """Single sample entry from a PREFACE samplesheet."""

    id: str
    filepath: Path
    sex: str | None = None
    ff: float | None = None

    def to_dict(self) -> dict[str, str | float | None]:
        return {
            "id": self.id,
            "filepath": self.filepath.as_posix(),
            "sex": self.sex,
            "ff": self.ff,
        }


@dataclass(frozen=True)
class TrainingDataset:
    """Preprocessed training matrix and sample metadata."""

    entries: list[SampleEntry]
    training_frame: np.ndarray
    x_ratios: np.ndarray
    possible_features: list[str]
    mean_features: np.ndarray


@dataclass(frozen=True)
class SampleAlignmentData:
    """Parsed single-sample alignment features for prediction."""

    x_ratio: float
    feature_ratios: np.ndarray


def load_samplesheet_schema() -> dict[str, Any]:
    """Load the bundled JSON Schema for PREFACE samplesheets."""
    return json.loads(SAMPLESHEET_SCHEMA_PATH.read_text(encoding="utf-8"))


def validate_samplesheet_data(
    data: object, source: Path | str = "<samplesheet>"
) -> list[SampleEntry]:
    """Validate parsed samplesheet data against the JSON Schema and return SampleEntry objects."""
    schema = load_samplesheet_schema()
    validator = jsonschema.Draft202012Validator(schema)
    errors = sorted(
        validator.iter_errors(data), key=lambda err: list(err.absolute_path)
    )
    if errors:
        formatted: list[str] = []
        for err in errors:
            loc = "$" + "".join(
                f"[{p}]" if isinstance(p, int) else f".{p}" for p in err.absolute_path
            )
            formatted.append(f"{loc}: {err.message}")
        detail = "; ".join(formatted)
        raise ValueError(f"Samplesheet '{source}' failed schema validation: {detail}")

    assert isinstance(data, list)
    entries: list[SampleEntry] = []
    seen_ids: set[str] = set()
    for idx, item in enumerate(data):
        sample_id = str(item["id"])
        if sample_id in seen_ids:
            raise ValueError(
                f"Samplesheet '{source}' contains duplicate sample id '{sample_id}' at $[{idx}]."
            )
        seen_ids.add(sample_id)
        sex_val = item.get("sex")
        ff_val = item.get("ff")
        entries.append(
            SampleEntry(
                id=sample_id,
                filepath=Path(str(item["filepath"])),
                sex=str(sex_val) if sex_val is not None else None,
                ff=float(ff_val) if ff_val is not None else None,
            )
        )
    return entries


def load_samplesheet(
    samplesheet_path: Path, rng: np.random.Generator | None = None
) -> list[SampleEntry]:
    """Load, lint, and optionally shuffle a YAML or JSON samplesheet file."""
    if not samplesheet_path.is_file():
        raise FileNotFoundError(f"The file '{samplesheet_path}' does not exist.")

    suffix = samplesheet_path.suffix.lower()
    if suffix not in SAMPLESHEET_SUFFIXES:
        raise ValueError(
            f"Unsupported samplesheet format '{samplesheet_path}'. Expected .yaml, .yml, or .json."
        )

    raw_text = samplesheet_path.read_text(encoding="utf-8")
    try:
        if suffix == ".json":
            parsed = json.loads(raw_text)
        else:
            parsed = yaml.safe_load(raw_text)
    except (json.JSONDecodeError, yaml.YAMLError) as exc:
        raise ValueError(
            f"Failed to parse samplesheet '{samplesheet_path}': {exc}"
        ) from exc

    entries = validate_samplesheet_data(parsed, source=samplesheet_path)
    if rng is not None:
        perm = rng.permutation(len(entries))
        return [entries[int(i)] for i in perm]
    return entries


def write_samplesheet(
    entries: Sequence[SampleEntry],
    outfile: Path | None = None,
    default_format: str = "yaml",
) -> str:
    """Serialize validated SampleEntry items to YAML or JSON."""
    payload = [entry.to_dict() for entry in entries]
    validate_samplesheet_data(
        payload, source=outfile if outfile is not None else "<output>"
    )

    fmt = default_format.lower()
    if outfile is not None:
        out_suffix = outfile.suffix.lower()
        if out_suffix == ".json":
            fmt = "json"
        elif out_suffix in {".yaml", ".yml"}:
            fmt = "yaml"
        else:
            raise ValueError(
                f"Unsupported output samplesheet extension '{outfile}'. Expected .yaml, .yml, or .json."
            )

    if fmt == "json":
        serialized = json.dumps(payload, indent=2) + "\n"
    else:
        serialized = yaml.safe_dump(payload, sort_keys=False)

    if outfile is not None:
        outfile.parent.mkdir(parents=True, exist_ok=True)
        outfile.write_text(serialized, encoding="utf-8")
        logger.info("Wrote samplesheet (%d samples) to '%s'", len(entries), outfile)
    else:
        print(serialized, end="")
    return serialized


def normalize_chr(raw: str) -> str:
    """Normalize reference chromosome names by stripping any 'chr' prefix."""
    cleaned = raw.strip()
    if cleaned.lower().startswith("chr"):
        cleaned = cleaned[3:]
    upper = cleaned.upper()
    if upper in {"X", "Y", "M", "MT"}:
        return "M" if upper == "MT" else upper
    return cleaned


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
    refs: list[tuple[str, int]],
    binsize: int,
    target_chrs: frozenset[str],
) -> tuple[dict[int, str], dict[str, int], dict[str, np.ndarray]]:
    idx_to_chr: dict[int, str] = {}
    ref_lengths: dict[str, int] = {}
    bin_counts: dict[str, np.ndarray] = {}
    for idx, (raw_name, length) in enumerate(refs):
        norm_chr = normalize_chr(raw_name)
        if norm_chr in target_chrs and length > 0:
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


def _read_bam_native(
    bam_path: Path,
    binsize: int,
    mapq: int,
    target_chrs: frozenset[str],
) -> tuple[dict[str, int], dict[str, np.ndarray]]:
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

    ref_id_to_chr, ref_lengths, bin_counts = _init_ref_bins(refs, binsize, target_chrs)
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

    return ref_lengths, bin_counts


def _read_alignment_external(
    aln_path: Path,
    binsize: int,
    mapq: int,
    reference: Path | None,
    target_chrs: frozenset[str],
) -> tuple[dict[str, int], dict[str, np.ndarray]]:
    """Read BAM/CRAM/SAM via pysam when available, falling back to samtools CLI."""
    if importlib.util.find_spec("pysam") is not None:
        import pysam

        ref_arg = str(reference) if reference is not None else None
        with pysam.AlignmentFile(str(aln_path), reference_filename=ref_arg) as samfile:
            refs = [
                (samfile.get_reference_name(i), int(samfile.lengths[i]))
                for i in range(samfile.nreferences)
            ]
            tid_to_chr, ref_lengths, bin_counts = _init_ref_bins(
                refs, binsize, target_chrs
            )
            for read in samfile.fetch(until_eof=True):
                if (read.flag & SKIP_BAM_FLAGS) == 0 and read.mapping_quality >= mapq:
                    _record_hit(
                        bin_counts,
                        tid_to_chr.get(read.reference_id),
                        read.reference_start,
                        binsize,
                    )
        return ref_lengths, bin_counts

    samtools_bin = shutil.which("samtools")
    if samtools_bin is None:
        raise RuntimeError(
            f"Cannot read '{aln_path}': install 'pysam' or 'samtools', or provide a '.bam' file."
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

    idx_to_chr, ref_lengths, bin_counts = _init_ref_bins(refs_sam, binsize, target_chrs)
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

    return ref_lengths, bin_counts


def read_alignment_bins(
    aln_path: Path,
    binsize: int,
    mapq: int,
    reference: Path | None,
    target_chrs: frozenset[str],
) -> tuple[dict[str, int], dict[str, np.ndarray]]:
    """Validate and bin aligned reads from a BAM, CRAM, or SAM file."""
    if not aln_path.is_file():
        raise FileNotFoundError(f"The file '{aln_path}' does not exist.")
    if binsize <= 0:
        raise ValueError("--binsize must be a positive integer.")

    suffix = aln_path.suffix.lower()
    if suffix not in ALIGNMENT_SUFFIXES:
        raise ValueError(
            f"Unsupported input file '{aln_path}'. Only aligned read files (.bam, .cram, .sam) are supported."
        )

    if suffix == ".bam" and importlib.util.find_spec("pysam") is None:
        return _read_bam_native(
            aln_path, binsize=binsize, mapq=mapq, target_chrs=target_chrs
        )
    return _read_alignment_external(
        aln_path,
        binsize=binsize,
        mapq=mapq,
        reference=reference,
        target_chrs=target_chrs,
    )


def compute_autosomal_density(
    ref_lengths: dict[str, int],
    bin_counts: dict[str, np.ndarray],
    binsize: int,
) -> float:
    """Compute mean mappable diploid autosomal read density (reads per base pair)."""
    auto_widths_list: list[np.ndarray] = []
    auto_counts_list: list[np.ndarray] = []

    for chrom in sorted(AUTOSOME_CHRS):
        if chrom not in ref_lengths:
            continue
        rlen = ref_lengths[chrom]
        starts = np.arange(0, rlen, binsize, dtype=np.int64)
        ends = np.minimum(starts + binsize, rlen)
        auto_widths_list.append(np.maximum(1, ends - starts).astype(np.float64))
        auto_counts_list.append(bin_counts[chrom].astype(np.float64))

    if not auto_counts_list:
        raise ValueError(
            "Alignment file contains no reference autosomes for normalization."
        )

    auto_widths = np.concatenate(auto_widths_list)
    auto_counts = np.concatenate(auto_counts_list)
    valid_auto = np.isfinite(auto_counts) & (auto_counts > 0.0)
    if not np.any(valid_auto):
        raise ValueError(
            "Alignment file contains no positive autosomal read counts for normalization."
        )

    r_auto = float(np.mean(auto_counts[valid_auto] / auto_widths[valid_auto]))
    if r_auto <= 0.0 or not np.isfinite(r_auto):
        raise ValueError("Computed autosomal read density is non-positive.")
    return r_auto


def _chrom_bin_log2_ratios(
    rlen: int, counts: np.ndarray, binsize: int, r_auto: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    starts0 = np.arange(0, rlen, binsize, dtype=np.int64)
    ends = np.minimum(starts0 + binsize, rlen)
    widths = np.maximum(1, ends - starts0).astype(np.float64)
    counts_f = counts.astype(np.float64)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratios = np.where(counts_f > 0.0, np.log2(counts_f / (r_auto * widths)), np.nan)
    return starts0 + 1, ends, ratios


def _align_ratios_to_loci(
    loci: list[str], ratios: np.ndarray, target_loci: list[str]
) -> np.ndarray:
    if loci == target_loci:
        return ratios
    locus_map = dict(zip(loci, ratios))
    return np.array(
        [locus_map.get(loc, np.nan) for loc in target_loci], dtype=np.float64
    )


def _extract_sample_cn_profile(
    aln_path: Path,
    binsize: int,
    mapq: int,
    reference: Path | None,
) -> tuple[float, list[str], np.ndarray]:
    """Extract chrX copy-number ratio and autosomal log2 ratio features from an alignment file."""
    ref_lengths, bin_counts = read_alignment_bins(
        aln_path,
        binsize=binsize,
        mapq=mapq,
        reference=reference,
        target_chrs=CN_TARGET_CHRS,
    )
    if "X" not in ref_lengths or "X" not in bin_counts:
        raise ValueError(
            f"Alignment file '{aln_path}' does not contain chromosome X ('X' or 'chrX')."
        )

    r_auto = compute_autosomal_density(ref_lengths, bin_counts, binsize)

    _, _, x_log2 = _chrom_bin_log2_ratios(
        ref_lengths["X"], bin_counts["X"], binsize, r_auto
    )
    with np.errstate(invalid="ignore"):
        x_ratio = float(np.power(2.0, np.nanmean(x_log2)))

    feature_loci: list[str] = []
    feature_ratios_list: list[np.ndarray] = []
    for chrom in CANONICAL_FEATURE_CHRS:
        if chrom not in ref_lengths:
            continue
        starts, ends, ratios = _chrom_bin_log2_ratios(
            ref_lengths[chrom], bin_counts[chrom], binsize, r_auto
        )
        feature_loci.extend(f"{chrom}:{int(s)}-{int(e)}" for s, e in zip(starts, ends))
        feature_ratios_list.append(ratios)

    if not feature_ratios_list:
        raise ValueError(
            f"Alignment file '{aln_path}' contains no non-excluded autosomal bins."
        )

    return x_ratio, feature_loci, np.concatenate(feature_ratios_list)


def _init_loader_worker() -> None:
    """Restrict per-worker thread pools to prevent oversubscription across processes."""
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


def _load_single_sample_alignment(
    args: tuple[Path, int, int, Path | None],
) -> tuple[float, list[str], np.ndarray]:
    filepath, binsize, mapq, reference = args
    return _extract_sample_cn_profile(
        filepath, binsize=binsize, mapq=mapq, reference=reference
    )


def load_training_dataset(
    entries: Sequence[SampleEntry],
    binsize: int = 100_000,
    mapq: int = 1,
    reference: Path | None = None,
    cpus: int = 1,
) -> TrainingDataset:
    """Load all training alignment files concurrently, compute X ratios, filter NAs, and mean-impute."""
    entry_list = list(entries)
    n_samples = len(entry_list)
    if n_samples == 0:
        raise ValueError("Samplesheet contains no samples.")

    items = [
        (entry.id, (entry.filepath, binsize, mapq, reference)) for entry in entry_list
    ]
    loaded = parallel_load_samples(items, _load_single_sample_alignment, cpus=cpus)

    x_ratios = np.array([x_r for x_r, _, _ in loaded], dtype=np.float64)
    feature_names_all = loaded[0][1]
    aligned_rows = [
        _align_ratios_to_loci(loci, ratios, feature_names_all)
        for _, loci, ratios in loaded
    ]

    logger.info("Creating training frame ...")
    training_frame = np.ascontiguousarray(np.vstack(aligned_rows), dtype=np.float64)
    valid_feature_mask = np.isnan(training_frame).sum(axis=0) < (n_samples * 0.01)

    training_frame = training_frame[:, valid_feature_mask]
    possible_features = [
        name for name, keep in zip(feature_names_all, valid_feature_mask) if keep
    ]
    if not possible_features:
        raise ValueError(
            "No autosomal genomic bins remained after filtering missing values across the cohort."
        )

    mean_features = np.nanmean(training_frame, axis=0)
    na_rows, na_cols = np.where(np.isnan(training_frame))
    if na_cols.size > 0:
        training_frame[na_rows, na_cols] = mean_features[na_cols]

    logger.info(
        "Remaining training features after 'NA' filtering: %d", len(possible_features)
    )

    return TrainingDataset(
        entries=entry_list,
        training_frame=training_frame,
        x_ratios=x_ratios,
        possible_features=possible_features,
        mean_features=mean_features,
    )


def load_inference_alignment(
    aln_path: Path,
    possible_features: list[str],
    binsize: int = 100_000,
    mapq: int = 1,
    reference: Path | None = None,
) -> SampleAlignmentData:
    """Extract X ratio and ordered autosomal bin ratios from a single alignment file for inference."""
    x_ratio, loci, ratios = _extract_sample_cn_profile(
        aln_path, binsize=binsize, mapq=mapq, reference=reference
    )
    selected_ratios = _align_ratios_to_loci(loci, ratios, possible_features)
    return SampleAlignmentData(x_ratio=x_ratio, feature_ratios=selected_ratios)


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
