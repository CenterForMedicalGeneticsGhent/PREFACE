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
INDEX_SUFFIXES: frozenset[str] = frozenset({".gzi", ".bai", ".csi", ".crai"})
SKIP_BAM_FLAGS: int = 0x4 | 0x100 | 0x200 | 0x400 | 0x800

ALL_AUTOSOME_CHRS: frozenset[str] = frozenset(str(i) for i in range(1, 23))
AUTOSOME_CHRS: frozenset[str] = ALL_AUTOSOME_CHRS - EXCLUDE_CHRS
CN_TARGET_CHRS: frozenset[str] = ALL_AUTOSOME_CHRS | {"X"}
FFY_TARGET_CHRS: frozenset[str] = AUTOSOME_CHRS | {"Y"}
ALL_TARGET_CHRS: frozenset[str] = CN_TARGET_CHRS | FFY_TARGET_CHRS
CANONICAL_FEATURE_CHRS: tuple[str, ...] = tuple(
    str(i) for i in range(1, 23) if str(i) not in EXCLUDE_CHRS
)


@dataclass(frozen=True)
class SampleEntry:
    """Single sample entry from a PREFACE samplesheet."""

    id: str
    filepath: Path
    index: Path | None = None
    sex: str | None = None
    ff: float | None = None

    def to_dict(self) -> dict[str, str | float | None]:
        data: dict[str, str | float | None] = {
            "id": self.id,
            "filepath": self.filepath.as_posix(),
        }
        if self.index is not None:
            data["index"] = self.index.as_posix()
        data["sex"] = self.sex
        data["ff"] = self.ff
        return data


@dataclass(frozen=True)
class SampleYProfile:
    """Normalized chromosome Y bin depths and coordinates for a single sample."""

    bin_loci: list[str]
    bin_starts: np.ndarray
    bin_ends: np.ndarray
    y_norm_bins: np.ndarray


@dataclass(frozen=True)
class UnifiedSampleData:
    """Single-pass extracted features for both FFY and PREFACE/FFX."""

    sample_id: str
    x_ratio: float | None
    feature_loci: list[str] | None
    feature_ratios: np.ndarray | None
    y_profile: SampleYProfile | None


@dataclass(frozen=True)
class TrainingDataset:
    """Preprocessed training matrix and sample metadata."""

    entries: list[SampleEntry]
    training_frame: np.ndarray
    x_ratios: np.ndarray
    possible_features: list[str]
    mean_features: np.ndarray


def load_samplesheet_schema() -> dict[str, Any]:
    """Load the bundled JSON Schema for PREFACE samplesheets."""
    return json.loads(SAMPLESHEET_SCHEMA_PATH.read_text(encoding="utf-8"))


def _resolve_and_check_readable_file(
    raw_path: str,
    field_name: str,
    sample_id: str,
    source: Path | str,
    base_dir: Path | None = None,
) -> Path:
    path = Path(raw_path)
    if (
        not path.is_absolute()
        and not path.exists()
        and base_dir is not None
        and (cand := base_dir / path).exists()
    ):
        path = cand
    if not path.is_file():
        raise FileNotFoundError(
            f"Samplesheet '{source}' sample '{sample_id}' {field_name} '{raw_path}' does not exist."
        )
    if not os.access(path, os.R_OK):
        raise PermissionError(
            f"Samplesheet '{source}' sample '{sample_id}' {field_name} '{raw_path}' is not readable."
        )
    return path


def validate_samplesheet_data(
    data: object, source: Path | str = "<samplesheet>", base_dir: Path | None = None
) -> list[SampleEntry]:
    """Validate parsed samplesheet data against the JSON Schema and return SampleEntry objects."""
    validator = jsonschema.Draft202012Validator(load_samplesheet_schema())
    errors = sorted(
        validator.iter_errors(data), key=lambda err: list(err.absolute_path)
    )
    if errors:
        formatted = [
            "$"
            + "".join(
                f"[{p}]" if isinstance(p, int) else f".{p}" for p in err.absolute_path
            )
            + f": {err.message}"
            for err in errors
        ]
        raise ValueError(
            f"Samplesheet '{source}' failed schema validation: {'; '.join(formatted)}"
        )

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
        filepath = _resolve_and_check_readable_file(
            str(item["filepath"]), "filepath", sample_id, source, base_dir=base_dir
        )
        index_val = item.get("index")
        index_path: Path | None = None
        if index_val is not None:
            raw_idx = str(index_val)
            if Path(raw_idx).suffix.lower() not in INDEX_SUFFIXES:
                raise ValueError(
                    f"Samplesheet '{source}' sample '{sample_id}' index '{raw_idx}' must end with .gzi, .bai, .csi, or .crai."
                )
            index_path = _resolve_and_check_readable_file(
                raw_idx, "index", sample_id, source, base_dir=base_dir
            )
        sex_val, ff_val = item.get("sex"), item.get("ff")
        entries.append(
            SampleEntry(
                sample_id,
                filepath,
                index_path,
                str(sex_val) if sex_val is not None else None,
                float(ff_val) if ff_val is not None else None,
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
        parsed = json.loads(raw_text) if suffix == ".json" else yaml.safe_load(raw_text)
    except (json.JSONDecodeError, yaml.YAMLError) as exc:
        raise ValueError(
            f"Failed to parse samplesheet '{samplesheet_path}': {exc}"
        ) from exc

    entries = validate_samplesheet_data(
        parsed, source=samplesheet_path, base_dir=samplesheet_path.parent
    )
    return (
        [entries[int(i)] for i in rng.permutation(len(entries))]
        if rng is not None
        else entries
    )


def write_samplesheet(
    entries: Sequence[SampleEntry],
    outfile: Path | None = None,
    default_format: str = "yaml",
) -> str:
    """Serialize validated SampleEntry items to YAML or JSON."""
    payload = [entry.to_dict() for entry in entries]
    validate_samplesheet_data(
        payload,
        source=outfile if outfile is not None else "<output>",
        base_dir=outfile.parent if outfile is not None else None,
    )
    fmt = default_format.lower()
    if outfile is not None:
        out_suffix = outfile.suffix.lower()
        if out_suffix not in SAMPLESHEET_SUFFIXES:
            raise ValueError(
                f"Unsupported output samplesheet extension '{outfile}'. Expected .yaml, .yml, or .json."
            )
        fmt = "json" if out_suffix == ".json" else "yaml"

    serialized = (
        json.dumps(payload, indent=2) + "\n"
        if fmt == "json"
        else yaml.safe_dump(payload, sort_keys=False)
    )
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
    return (
        ("M" if upper == "MT" else upper) if upper in {"X", "Y", "M", "MT"} else cleaned
    )


def _iter_bgzf_blocks(fh: BinaryIO) -> bytes:
    """Decompress all BGZF/gzip blocks from a BAM file stream."""
    chunks: list[bytes] = []
    while header := fh.read(12):
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
        cdata, trailer = fh.read(cdata_len), fh.read(8)
        if len(cdata) < cdata_len or len(trailer) < 8:
            raise ValueError("Truncated BGZF block payload in BAM file.")
        _crc, isize = struct.unpack("<II", trailer)
        if isize > 0:
            chunks.append(zlib.decompress(cdata, -zlib.MAX_WBITS))
    return b"".join(chunks)


def _init_multi_ref_bins(
    refs: list[tuple[str, int]], bin_specs: Sequence[tuple[int, frozenset[str]]]
) -> tuple[dict[int, str], dict[str, int], list[dict[str, np.ndarray]]]:
    all_targets: frozenset[str] = frozenset().union(*(chrs for _, chrs in bin_specs))
    idx_to_chr: dict[int, str] = {}
    ref_lengths: dict[str, int] = {}
    spec_counts: list[dict[str, np.ndarray]] = [{} for _ in bin_specs]
    for idx, (raw_name, length) in enumerate(refs):
        norm_chr = normalize_chr(raw_name)
        if norm_chr in all_targets and length > 0:
            idx_to_chr[idx] = norm_chr
            ref_lengths[norm_chr] = length
            for s_idx, (bsize, target_chrs) in enumerate(bin_specs):
                if norm_chr in target_chrs:
                    spec_counts[s_idx][norm_chr] = np.zeros(
                        (length + bsize - 1) // bsize, dtype=np.int64
                    )
    return idx_to_chr, ref_lengths, spec_counts


def _record_multi_hit(
    spec_counts: list[dict[str, np.ndarray]],
    bin_sizes: list[int],
    chrom: str | None,
    pos0: int,
) -> None:
    if chrom is not None and pos0 >= 0:
        for counts_map, bsize in zip(spec_counts, bin_sizes):
            if (arr := counts_map.get(chrom)) is not None and (
                bin_idx := pos0 // bsize
            ) < arr.size:
                arr[bin_idx] += 1


def _read_bam_native_multi(
    bam_path: Path, bin_specs: Sequence[tuple[int, frozenset[str]]], mapq: int
) -> tuple[dict[str, int], list[dict[str, np.ndarray]]]:
    """Read and bin a BAM file across one or more bin sizes in a single pass."""
    with bam_path.open("rb") as fh:
        raw = _iter_bgzf_blocks(fh)
    if len(raw) < 12 or raw[:4] != b"BAM\x01":
        raise ValueError(f"File '{bam_path}' is not a valid BAM file.")
    offset = 8 + struct.unpack_from("<i", raw, 4)[0]
    n_ref = struct.unpack_from("<i", raw, offset)[0]
    offset += 4
    refs: list[tuple[str, int]] = []
    for _ in range(n_ref):
        l_name = struct.unpack_from("<i", raw, offset)[0]
        name = raw[offset + 4 : offset + 4 + l_name].rstrip(b"\x00").decode("ascii")
        l_ref = struct.unpack_from("<i", raw, offset + 4 + l_name)[0]
        offset += 8 + l_name
        refs.append((name, l_ref))

    ref_id_to_chr, ref_lengths, spec_counts = _init_multi_ref_bins(refs, bin_specs)
    bin_sizes = [bsize for bsize, _ in bin_specs]
    total_len = len(raw)
    while offset + 4 <= total_len:
        block_size = struct.unpack_from("<i", raw, offset)[0]
        rec_start, offset = offset + 4, offset + 4 + block_size
        if offset > total_len or block_size < 16:
            break
        ref_id, pos0, _lq, rec_mapq, _bin, _nc, flag = struct.unpack_from(
            "<iiBBHHH", raw, rec_start
        )
        if (flag & SKIP_BAM_FLAGS) == 0 and rec_mapq >= mapq:
            _record_multi_hit(spec_counts, bin_sizes, ref_id_to_chr.get(ref_id), pos0)
    return ref_lengths, spec_counts


def _read_alignment_external_multi(
    aln_path: Path,
    bin_specs: Sequence[tuple[int, frozenset[str]]],
    mapq: int,
    reference: Path | None,
    index: Path | None = None,
) -> tuple[dict[str, int], list[dict[str, np.ndarray]]]:
    """Read BAM/CRAM/SAM via pysam or samtools across multiple bin specs in a single pass."""
    bin_sizes = [bsize for bsize, _ in bin_specs]
    if importlib.util.find_spec("pysam") is not None:
        import pysam

        with pysam.AlignmentFile(
            str(aln_path),
            index_filename=str(index) if index else None,
            reference_filename=str(reference) if reference else None,
        ) as samfile:
            refs = [
                (samfile.get_reference_name(i), int(samfile.lengths[i]))
                for i in range(samfile.nreferences)
            ]
            tid_to_chr, ref_lengths, spec_counts = _init_multi_ref_bins(refs, bin_specs)
            for read in samfile.fetch(until_eof=True):
                if (read.flag & SKIP_BAM_FLAGS) == 0 and read.mapping_quality >= mapq:
                    _record_multi_hit(
                        spec_counts,
                        bin_sizes,
                        tid_to_chr.get(read.reference_id),
                        read.reference_start,
                    )
        return ref_lengths, spec_counts

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

    idx_to_chr, ref_lengths, spec_counts = _init_multi_ref_bins(refs_sam, bin_specs)
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
            if len(parts := line.split("\t", 5)) >= 4:
                _record_multi_hit(
                    spec_counts, bin_sizes, raw_to_chr.get(parts[2]), int(parts[3]) - 1
                )
        if proc.wait() != 0:
            raise RuntimeError(f"samtools view failed on '{aln_path}'.")
    return ref_lengths, spec_counts


def read_alignment_multi_bins(
    aln_path: Path,
    bin_specs: Sequence[tuple[int, frozenset[str]]],
    mapq: int,
    reference: Path | None,
    index: Path | None = None,
) -> tuple[dict[str, int], list[dict[str, np.ndarray]]]:
    """Validate and bin aligned reads across multiple bin sizes in a single pass."""
    if not aln_path.is_file():
        raise FileNotFoundError(f"The file '{aln_path}' does not exist.")
    if not os.access(aln_path, os.R_OK):
        raise PermissionError(f"The file '{aln_path}' is not readable.")
    if any(bsize <= 0 for bsize, _ in bin_specs):
        raise ValueError("--binsize must be a positive integer.")
    if aln_path.suffix.lower() not in ALIGNMENT_SUFFIXES:
        raise ValueError(
            f"Unsupported input file '{aln_path}'. Only aligned read files (.bam, .cram, .sam) are supported."
        )
    if index is not None:
        if index.suffix.lower() not in INDEX_SUFFIXES:
            raise ValueError(
                f"Unsupported index file '{index}'. Expected .gzi, .bai, .csi, or .crai."
            )
        if not index.is_file():
            raise FileNotFoundError(f"The index file '{index}' does not exist.")
        if not os.access(index, os.R_OK):
            raise PermissionError(f"The index file '{index}' is not readable.")
    if aln_path.suffix.lower() == ".bam" and importlib.util.find_spec("pysam") is None:
        return _read_bam_native_multi(aln_path, bin_specs=bin_specs, mapq=mapq)
    return _read_alignment_external_multi(
        aln_path, bin_specs=bin_specs, mapq=mapq, reference=reference, index=index
    )


def _read_bam_native(
    bam_path: Path, binsize: int, mapq: int, target_chrs: frozenset[str]
) -> tuple[dict[str, int], dict[str, np.ndarray]]:
    """Read and bin a BAM file using Python's built-in zlib and struct modules."""
    ref_lengths, spec_counts = _read_bam_native_multi(
        bam_path, [(binsize, target_chrs)], mapq
    )
    return ref_lengths, spec_counts[0]


def _bin_bounds_and_widths(
    rlen: int, binsize: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    starts0 = np.arange(0, rlen, binsize, dtype=np.int64)
    ends = np.minimum(starts0 + binsize, rlen)
    return starts0 + 1, ends, np.maximum(1, ends - starts0).astype(np.float64)


def compute_autosomal_density(
    ref_lengths: dict[str, int], bin_counts: dict[str, np.ndarray], binsize: int
) -> float:
    """Compute mean mappable diploid autosomal read density (reads per base pair)."""
    widths_parts = [
        _bin_bounds_and_widths(ref_lengths[c], binsize)[2]
        for c in sorted(AUTOSOME_CHRS)
        if c in ref_lengths and c in bin_counts
    ]
    counts_parts = [
        bin_counts[c].astype(np.float64)
        for c in sorted(AUTOSOME_CHRS)
        if c in ref_lengths and c in bin_counts
    ]
    if not counts_parts:
        raise ValueError(
            "Alignment file contains no reference autosomes for normalization."
        )
    auto_widths, auto_counts = (
        np.concatenate(widths_parts),
        np.concatenate(counts_parts),
    )
    valid = np.isfinite(auto_counts) & (auto_counts > 0.0)
    if not np.any(valid):
        raise ValueError(
            "Alignment file contains no positive autosomal read counts for normalization."
        )
    r_auto = float(np.mean(auto_counts[valid] / auto_widths[valid]))
    if r_auto <= 0.0 or not np.isfinite(r_auto):
        raise ValueError("Computed autosomal read density is non-positive.")
    return r_auto


def _chrom_bin_log2_ratios(
    rlen: int, counts: np.ndarray, binsize: int, r_auto: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    starts, ends, widths = _bin_bounds_and_widths(rlen, binsize)
    counts_f = counts.astype(np.float64)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratios = np.where(counts_f > 0.0, np.log2(counts_f / (r_auto * widths)), np.nan)
    return starts, ends, ratios


def align_to_loci(
    loci: list[str],
    values: np.ndarray,
    target_loci: list[str],
    fill_nan: float | None = None,
) -> np.ndarray:
    """Align bin values to target genomic loci, optionally replacing missing/NaN entries."""
    if loci == target_loci:
        arr = values.astype(np.float64, copy=False)
    else:
        default_val = np.nan if fill_nan is None else fill_nan
        locus_map = dict(zip(loci, values))
        arr = np.array(
            [locus_map.get(loc, default_val) for loc in target_loci], dtype=np.float64
        )
    return arr if fill_nan is None else np.nan_to_num(arr, nan=fill_nan)


def _build_y_profile_from_counts(
    ref_lengths: dict[str, int], bin_counts: dict[str, np.ndarray], binsize: int
) -> SampleYProfile:
    if "Y" not in ref_lengths or "Y" not in bin_counts:
        raise ValueError(
            "Alignment header does not contain chromosome Y ('Y' or 'chrY')."
        )
    r_auto = compute_autosomal_density(ref_lengths, bin_counts, binsize)
    y_starts, y_ends, y_widths = _bin_bounds_and_widths(ref_lengths["Y"], binsize)
    y_counts = bin_counts["Y"].astype(np.float64)
    y_norm = np.where(np.isfinite(y_counts), y_counts / (r_auto * y_widths), np.nan)
    return SampleYProfile(
        [f"Y:{int(s)}-{int(e)}" for s, e in zip(y_starts, y_ends)],
        y_starts,
        y_ends,
        y_norm,
    )


def _build_cn_profile_from_counts(
    aln_path: Path,
    ref_lengths: dict[str, int],
    bin_counts: dict[str, np.ndarray],
    binsize: int,
) -> tuple[float, list[str], np.ndarray]:
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
        if chrom in ref_lengths and chrom in bin_counts:
            starts, ends, ratios = _chrom_bin_log2_ratios(
                ref_lengths[chrom], bin_counts[chrom], binsize, r_auto
            )
            feature_loci.extend(
                f"{chrom}:{int(s)}-{int(e)}" for s, e in zip(starts, ends)
            )
            feature_ratios_list.append(ratios)
    if not feature_ratios_list:
        raise ValueError(
            f"Alignment file '{aln_path}' contains no non-excluded autosomal bins."
        )
    return x_ratio, feature_loci, np.concatenate(feature_ratios_list)


def extract_unified_sample(
    aln_path: Path,
    sample_id: str = "sample",
    binsize: int = 100_000,
    ffy_binsize: int = 5_000,
    mapq: int = 1,
    reference: Path | None = None,
    index: Path | None = None,
    include_cn: bool = True,
    include_ffy: bool = True,
) -> UnifiedSampleData:
    """Read an alignment file in a single pass and extract both FFY and CN/FFX features."""
    if include_cn and include_ffy and binsize == ffy_binsize:
        specs: list[tuple[int, frozenset[str]]] = [(binsize, ALL_TARGET_CHRS)]
    else:
        specs = [
            *([(binsize, CN_TARGET_CHRS)] if include_cn else []),
            *([(ffy_binsize, FFY_TARGET_CHRS)] if include_ffy else []),
        ]
    ref_lengths, spec_counts = read_alignment_multi_bins(
        aln_path, specs, mapq=mapq, reference=reference, index=index
    )
    x_ratio, feature_loci, feature_ratios = (
        _build_cn_profile_from_counts(aln_path, ref_lengths, spec_counts[0], binsize)
        if include_cn
        else (None, None, None)
    )
    y_profile = (
        _build_y_profile_from_counts(ref_lengths, spec_counts[-1], ffy_binsize)
        if include_ffy
        else None
    )
    return UnifiedSampleData(
        sample_id, x_ratio, feature_loci, feature_ratios, y_profile
    )


def _init_loader_worker() -> None:
    """Restrict per-worker thread pools to prevent oversubscription across processes."""
    os.environ["OMP_NUM_THREADS"] = "1"


def parallel_load_samples[T, R](
    items: Sequence[tuple[str, T]], worker_fn: Callable[[T], R], cpus: int
) -> list[R]:
    """Load sample files in parallel across worker processes while preserving input order."""
    total = len(items)
    if total == 0:
        return []
    max_workers = min(max(1, cpus), total)
    if max_workers == 1:
        results: list[R] = []
        for done, (sid, payload) in enumerate(items, start=1):
            logger.info(
                "Loading sample %s | %d/%d remaining ...", sid, total - done, total
            )
            results.append(worker_fn(payload))
        return results

    ordered: dict[int, R] = {}
    with ProcessPoolExecutor(
        max_workers=max_workers, initializer=_init_loader_worker
    ) as pool:
        futs = {
            pool.submit(worker_fn, payload): (idx, sid)
            for idx, (sid, payload) in enumerate(items)
        }
        try:
            for done, fut in enumerate(as_completed(futs), start=1):
                idx, sid = futs[fut]
                ordered[idx] = fut.result()
                logger.info(
                    "Loading sample %s | %d/%d remaining ...", sid, total - done, total
                )
        except Exception:
            for pending in futs:
                pending.cancel()
            raise
    return [ordered[i] for i in range(total)]


def _extract_unified_worker(
    args: tuple[str, Path, int, int, int, Path | None, Path | None, bool, bool],
) -> UnifiedSampleData:
    sid, fpath, bsize, ffy_bsize, mapq, ref, idx, inc_cn, inc_ffy = args
    return extract_unified_sample(
        fpath, sid, bsize, ffy_bsize, mapq, ref, idx, inc_cn, inc_ffy
    )


def save_cohort_counts_cache(
    cache_path: Path,
    samples: Sequence[UnifiedSampleData],
    lb_profiles: Sequence[SampleYProfile],
    binsize: int,
    ffy_binsize: int,
) -> None:
    """Persist single-pass extracted cohort features and chrY profiles to a compressed .npz file."""
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    has_cn = bool(samples and all(s.feature_ratios is not None for s in samples))
    has_ffy = bool(samples and all(s.y_profile is not None for s in samples))
    payload: dict[str, Any] = {
        "sample_ids": np.array([s.sample_id for s in samples], dtype=str),
        "binsize": np.array(binsize, dtype=np.int64),
        "ffy_binsize": np.array(ffy_binsize, dtype=np.int64),
        "has_cn": np.array(has_cn, dtype=bool),
        "has_ffy": np.array(has_ffy, dtype=bool),
    }
    if has_cn and (ref_cn := samples[0].feature_loci) is not None:
        payload["cn_loci"] = np.array(ref_cn, dtype=str)
        payload["x_ratios"] = np.array(
            [float(s.x_ratio) for s in samples if s.x_ratio is not None],
            dtype=np.float64,
        )
        payload["feature_ratios"] = np.vstack(
            [
                align_to_loci(s.feature_loci, s.feature_ratios, ref_cn)
                for s in samples
                if s.feature_loci is not None and s.feature_ratios is not None
            ]
        )
    if has_ffy and (first_y := samples[0].y_profile) is not None:
        ref_y = first_y.bin_loci
        payload["y_loci"] = np.array(ref_y, dtype=str)
        payload["y_starts"], payload["y_ends"] = first_y.bin_starts, first_y.bin_ends
        payload["y_norm_matrix"] = np.vstack(
            [
                align_to_loci(s.y_profile.bin_loci, s.y_profile.y_norm_bins, ref_y)
                for s in samples
                if s.y_profile is not None
            ]
        )
        if lb_profiles:
            payload["lb_y_norm_matrix"] = np.vstack(
                [align_to_loci(p.bin_loci, p.y_norm_bins, ref_y) for p in lb_profiles]
            )
    np.savez_compressed(cache_path, **payload)
    logger.info("Saved single-pass read count cache to '%s'", cache_path)


def load_cohort_counts_cache(
    cache_path: Path,
    expected_ids: Sequence[str],
    binsize: int,
    ffy_binsize: int,
    need_cn: bool,
    need_ffy: bool,
) -> tuple[list[UnifiedSampleData], list[SampleYProfile]] | None:
    """Load cached single-pass cohort features if compatible with the requested run."""
    if not cache_path.is_file():
        return None
    with np.load(cache_path, allow_pickle=False) as d:
        cached_ids = [str(x) for x in d["sample_ids"].tolist()]
        has_cn, has_ffy = bool(d["has_cn"]), bool(d["has_ffy"])
        if (
            (need_cn and (not has_cn or int(d["binsize"]) != binsize))
            or (need_ffy and (not has_ffy or int(d["ffy_binsize"]) != ffy_binsize))
            or set(cached_ids) != set(expected_ids)
            or len(cached_ids) != len(expected_ids)
        ):
            return None
        id_map = {sid: i for i, sid in enumerate(cached_ids)}
        cn_loci = [str(x) for x in d["cn_loci"].tolist()] if has_cn else None
        y_loci = [str(x) for x in d["y_loci"].tolist()] if has_ffy else None
        samples = [
            UnifiedSampleData(
                sid,
                float(d["x_ratios"][id_map[sid]]) if has_cn else None,
                cn_loci,
                d["feature_ratios"][id_map[sid]] if has_cn else None,
                SampleYProfile(
                    y_loci, d["y_starts"], d["y_ends"], d["y_norm_matrix"][id_map[sid]]
                )
                if (has_ffy and y_loci is not None)
                else None,
            )
            for sid in expected_ids
        ]
        lb_profiles = (
            [
                SampleYProfile(y_loci, d["y_starts"], d["y_ends"], row)
                for row in d["lb_y_norm_matrix"]
            ]
            if (has_ffy and "lb_y_norm_matrix" in d and y_loci is not None)
            else []
        )
    logger.info("Loaded single-pass read count cache from '%s'", cache_path)
    return samples, lb_profiles


def load_cohort_unified_data(
    entries: Sequence[SampleEntry],
    binsize: int = 100_000,
    ffy_binsize: int = 5_000,
    mapq: int = 1,
    reference: Path | None = None,
    male_ref_paths: list[Path] | None = None,
    cpus: int = 1,
    include_cn: bool = True,
    include_ffy: bool = True,
    cache_counts: Path | None = None,
) -> tuple[list[UnifiedSampleData], list[SampleYProfile]]:
    """Read all cohort alignment files in a single pass per sample (with optional disk caching)."""
    entry_list = list(entries)
    if not entry_list:
        raise ValueError("Samplesheet contains no samples.")
    if cache_counts is not None and cache_counts.is_file():
        cached = load_cohort_counts_cache(
            cache_counts,
            [e.id for e in entry_list],
            binsize,
            ffy_binsize,
            include_cn,
            include_ffy,
        )
        if cached is not None:
            return cached

    work_items = [
        (
            e.id,
            (
                e.id,
                e.filepath,
                binsize,
                ffy_binsize,
                mapq,
                reference,
                e.index,
                include_cn,
                include_ffy,
            ),
        )
        for e in entry_list
    ]
    if include_ffy and male_ref_paths:
        work_items.extend(
            (
                f"male_ref:{p.name}",
                (
                    f"male_ref:{p.name}",
                    p,
                    binsize,
                    ffy_binsize,
                    mapq,
                    reference,
                    None,
                    False,
                    True,
                ),
            )
            for p in male_ref_paths
        )
    all_loaded = parallel_load_samples(work_items, _extract_unified_worker, cpus=cpus)
    cohort_samples = all_loaded[: len(entry_list)]
    lb_profiles = [
        s.y_profile for s in all_loaded[len(entry_list) :] if s.y_profile is not None
    ]
    if cache_counts is not None:
        save_cohort_counts_cache(
            cache_counts, cohort_samples, lb_profiles, binsize, ffy_binsize
        )
    return cohort_samples, lb_profiles


def build_training_dataset_from_unified(
    entries: Sequence[SampleEntry], unified_samples: Sequence[UnifiedSampleData]
) -> TrainingDataset:
    """Construct TrainingDataset from pre-extracted UnifiedSampleData in the given entry order."""
    entry_list = list(entries)
    if not entry_list:
        raise ValueError("Samplesheet contains no samples.")
    sample_by_id = {s.sample_id: s for s in unified_samples}
    ordered = [sample_by_id[e.id] for e in entry_list]
    x_ratios = np.array(
        [float(s.x_ratio) if s.x_ratio is not None else np.nan for s in ordered],
        dtype=np.float64,
    )
    first_loci = ordered[0].feature_loci
    if first_loci is None:
        raise ValueError("Missing autosomal copy-number features for training.")

    aligned_rows = [
        align_to_loci(s.feature_loci, s.feature_ratios, first_loci)
        for s in ordered
        if s.feature_loci is not None and s.feature_ratios is not None
    ]
    logger.info("Creating training frame ...")
    training_frame = np.ascontiguousarray(np.vstack(aligned_rows), dtype=np.float64)
    valid_mask = np.isnan(training_frame).sum(axis=0) < (len(entry_list) * 0.01)
    training_frame = training_frame[:, valid_mask]
    possible_features = [name for name, keep in zip(first_loci, valid_mask) if keep]
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
        entry_list, training_frame, x_ratios, possible_features, mean_features
    )


def load_training_dataset(
    entries: Sequence[SampleEntry],
    binsize: int = 100_000,
    mapq: int = 1,
    reference: Path | None = None,
    cpus: int = 1,
) -> TrainingDataset:
    """Load all training alignment files concurrently, compute X ratios, filter NAs, and mean-impute."""
    unified_samples, _ = load_cohort_unified_data(
        entries,
        binsize=binsize,
        mapq=mapq,
        reference=reference,
        cpus=cpus,
        include_cn=True,
        include_ffy=False,
    )
    return build_training_dataset_from_unified(entries, unified_samples)
