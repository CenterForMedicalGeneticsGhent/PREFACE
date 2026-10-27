"""Unit and integration tests for PREFACE."""

from __future__ import annotations

import json
import struct
import zlib
from pathlib import Path

import numpy as np
import pytest
import torch
import yaml
from typer.testing import CliRunner

from preface.cli import _normalize_argv, app
from preface.ffy import (
    extract_sample_y_profile,
    run_ffy_single,
)
from preface.hardware import configure_hardware, get_onnx_providers, resolve_device
from preface.io import (
    FFY_TARGET_CHRS,
    _read_bam_native,
    load_samplesheet,
    load_training_dataset,
    parallel_load_samples,
    validate_samplesheet_data,
    write_samplesheet,
)
from preface.logging_utils import LogLevel, configure_logging
from preface.math_utils import (
    get_m_diff,
    get_sd_diff,
    ols_fit,
    pca_irlba_torch,
    pearson_cor,
    rlm_huber,
)
from preface.models import (
    export_pipeline_to_onnx,
    load_ffy_onnx_metadata,
    load_onnx_metadata,
    predict_model,
    run_onnx_inference,
    train_neural_net,
    train_olm,
)
from preface.pipeline import run_predict, run_train

runner = CliRunner()


def _make_bgzf_block(payload: bytes) -> bytes:
    """Wrap uncompressed bytes into a valid BGZF block."""
    comp_obj = zlib.compressobj(level=6, method=zlib.DEFLATED, wbits=-zlib.MAX_WBITS)
    cdata = comp_obj.compress(payload) + comp_obj.flush()
    bsize = 12 + 6 + len(cdata) + 8 - 1
    header = struct.pack(
        "<BBBBIBBHBBHH",
        0x1F,
        0x8B,
        8,
        0x04,
        0,
        0,
        0xFF,
        6,
        66,
        67,
        2,
        bsize,
    )
    crc = zlib.crc32(payload) & 0xFFFFFFFF
    trailer = struct.pack("<II", crc, len(payload))
    return header + cdata + trailer


def _write_synthetic_bam_from_bin_counts(
    bam_path: Path,
    refs: list[tuple[str, int]],
    bin_reads: list[tuple[int, int, int, int, int]],
) -> None:
    """Write a valid BGZF-compressed BAM file from (ref_id, pos0, mapq, flag, count) tuples."""
    parts: list[bytes] = [b"BAM\x01", struct.pack("<i", 0)]
    parts.append(struct.pack("<i", len(refs)))
    for name, length in refs:
        name_b = name.encode("ascii") + b"\x00"
        parts.append(
            struct.pack("<i", len(name_b)) + name_b + struct.pack("<i", length)
        )

    qname = b"r\x00"
    for ref_id, pos0, mapq, flag, count in bin_reads:
        if count <= 0:
            continue
        rec_body = (
            struct.pack(
                "<iiBBHHHiiii",
                ref_id,
                pos0,
                len(qname),
                mapq,
                4680,
                0,
                flag,
                0,
                -1,
                -1,
                0,
            )
            + qname
        )
        rec = struct.pack("<i", len(rec_body)) + rec_body
        parts.append(rec * count)

    uncompressed = b"".join(parts)
    blocks: list[bytes] = []
    chunk_size = 60_000
    for idx in range(0, len(uncompressed), chunk_size):
        blocks.append(_make_bgzf_block(uncompressed[idx : idx + chunk_size]))
    blocks.append(_make_bgzf_block(b""))
    bam_path.write_bytes(b"".join(blocks))


def test_hardware_resolution() -> None:
    dev = resolve_device("auto")
    assert dev.type in {"cpu", "cuda", "mps"}
    assert resolve_device("cpu").type == "cpu"
    with pytest.raises(ValueError, match="Unsupported device"):
        resolve_device("tpu")

    cpus = configure_hardware(cpus=2, seed=42)
    assert cpus == 2
    providers = get_onnx_providers("cpu")
    assert providers == ["CPUExecutionProvider"]


def test_math_utils_pca_ols_rlm() -> None:
    rng = np.random.default_rng(7)
    x = rng.normal(size=(30, 80))
    pca_res = pca_irlba_torch(x, n_components=10, device=torch.device("cpu"))
    assert pca_res.x.shape == (30, 10)
    assert pca_res.rotation.shape == (80, 10)
    assert pca_res.center.shape == (80,)
    assert pca_res.importance_prop.shape == (10,)
    assert np.all(np.diff(pca_res.importance_prop) <= 1e-12)
    assert 0.0 < float(np.sum(pca_res.importance_prop)) <= 1.0 + 1e-6

    t = np.linspace(0.01, 0.20, 40)
    y = 0.035 + 2.2 * t
    intercept, slope = ols_fit(t, y)
    assert np.isclose(intercept, 0.035)
    assert np.isclose(slope, 2.2)

    y_outlier = y.copy()
    y_outlier[0] = 1.5
    r_int, r_slope, wr = rlm_huber(t, y_outlier)
    assert np.isclose(r_int, 0.035, atol=0.01)
    assert np.isclose(r_slope, 2.2, atol=0.05)
    assert 0.9 < wr <= 1.0

    assert np.isclose(get_m_diff(t, t + 0.02, use_abs=True), 0.02)
    assert np.isclose(get_m_diff(t, t + 0.02, use_abs=False), -0.02)
    assert np.isclose(get_sd_diff(t, t + 0.02), 0.0)
    assert np.isclose(pearson_cor(t, y), 1.0)


def test_models_and_onnx_roundtrip(tmp_path: Path) -> None:
    rng = np.random.default_rng(1)
    n_samples, n_raw_features, n_feat = 40, 25, 4
    raw_x = rng.normal(size=(n_samples, n_raw_features))
    mean_features = np.mean(raw_x, axis=0)
    pca_res = pca_irlba_torch(raw_x, n_components=n_feat, device=torch.device("cpu"))

    x_pc = pca_res.x[:, :n_feat]
    true_y = np.clip(0.10 + 0.02 * x_pc[:, 0] - 0.01 * x_pc[:, 1], 0.01, 0.30)

    olm_model = train_olm(x_pc, true_y, device=torch.device("cpu"))
    olm_preds = predict_model(olm_model, x_pc, device=torch.device("cpu"))
    assert np.allclose(olm_preds, true_y, atol=1e-4)

    torch.manual_seed(1)
    nn_model = train_neural_net(x_pc, true_y, hidden=2, device=torch.device("cpu"))
    nn_preds = predict_model(nn_model, x_pc, device=torch.device("cpu"))
    assert pearson_cor(nn_preds, true_y) > 0.95

    onnx_path = tmp_path / "model.onnx"
    feature_names = [
        f"1:{i * 5000 + 1}-{(i + 1) * 5000}" for i in range(n_raw_features)
    ]
    export_pipeline_to_onnx(
        predictor=nn_model,
        mean_features=mean_features,
        pca_center=pca_res.center,
        pca_rotation=pca_res.rotation[:, :n_feat],
        the_intercept=0.005,
        the_slope=1.05,
        the_intercept_x=0.99,
        the_slope_x=-0.5,
        possible_features=feature_names,
        onnx_path=onnx_path,
        binsize=5000,
    )

    meta = load_onnx_metadata(onnx_path)
    assert meta.possible_features == feature_names
    assert np.isclose(meta.the_intercept_x, 0.99)
    assert np.isclose(meta.the_slope_x, -0.5)
    assert meta.n_feat == n_feat
    assert meta.binsize == 5000
    assert not meta.is_olm

    sample_with_nan = raw_x[0].copy()
    sample_with_nan[3] = np.nan
    expected_imputed = raw_x[0].copy()
    expected_imputed[3] = mean_features[3]
    expected_pc = (expected_imputed - pca_res.center) @ pca_res.rotation[:, :n_feat]
    expected_raw = float(
        predict_model(nn_model, expected_pc.reshape(1, -1), device=torch.device("cpu"))[
            0
        ]
    )
    expected_final = 0.005 + 1.05 * expected_raw

    onnx_pred = run_onnx_inference(onnx_path, sample_with_nan, requested_device="cpu")
    assert np.isclose(onnx_pred, expected_final, atol=1e-4)


def test_samplesheet_lint_and_validation(tmp_path: Path) -> None:
    example_sheet = (
        Path(__file__).resolve().parents[1] / "examples" / "samplesheet.yaml"
    )
    entries = load_samplesheet(example_sheet)
    assert len(entries) == 10
    assert entries[0].id == "ID_0"
    assert entries[0].index is not None and entries[0].index.name == "sample.bam.bai"
    assert entries[0].sex == "M"
    assert entries[0].ff is not None and 0.0 <= entries[0].ff <= 1.0

    res_lint = runner.invoke(app, ["lint", str(example_sheet)])
    assert res_lint.exit_code == 0

    json_sheet = tmp_path / "samplesheet.json"
    write_samplesheet(entries, outfile=json_sheet)
    loaded_json = load_samplesheet(json_sheet)
    assert loaded_json == entries
    assert runner.invoke(app, ["lint", str(json_sheet)]).exit_code == 0

    valid_bam = tmp_path / "valid.bam"
    valid_bam.write_bytes(b"BAM\x01")
    for ext in (".gzi", ".bai", ".csi", ".crai"):
        idx_file = tmp_path / f"valid{ext}"
        idx_file.write_bytes(b"idx")
        parsed_idx = validate_samplesheet_data(
            [{"id": "S1", "filepath": str(valid_bam), "index": str(idx_file)}]
        )
        assert parsed_idx[0].index == idx_file

    # Test rejection of invalid schema inputs
    invalid_cases = [
        ([], "non-empty"),
        ([{"id": "S1"}], "'filepath' is a required property"),
        ([{"id": "S1", "filepath": str(valid_bam), "sex": "X"}], "is not one of"),
        (
            [{"id": "S1", "filepath": str(valid_bam), "ff": 1.5}],
            "greater than the maximum of 1.0",
        ),
        (
            [{"id": "S1", "filepath": str(valid_bam), "ff": -0.1}],
            "less than the minimum of 0.0",
        ),
        (
            [{"id": "S1", "filepath": str(valid_bam), "index": "bad.idx"}],
            "does not match",
        ),
        (
            [{"id": "S1", "filepath": str(valid_bam), "extra": 1}],
            "Additional properties are not allowed",
        ),
        (
            [
                {"id": "S1", "filepath": str(valid_bam)},
                {"id": "S1", "filepath": str(valid_bam)},
            ],
            "duplicate sample id 'S1'",
        ),
    ]
    for payload, expected_err in invalid_cases:
        with pytest.raises(ValueError, match=expected_err):
            validate_samplesheet_data(payload)

    # Test file existence and readability checks for filepath and index
    with pytest.raises(FileNotFoundError, match="filepath .* does not exist"):
        validate_samplesheet_data(
            [{"id": "S1", "filepath": str(tmp_path / "missing.bam")}]
        )

    with pytest.raises(FileNotFoundError, match="index .* does not exist"):
        validate_samplesheet_data(
            [
                {
                    "id": "S1",
                    "filepath": str(valid_bam),
                    "index": str(tmp_path / "missing.bai"),
                }
            ]
        )

    unreadable_bam = tmp_path / "unreadable.bam"
    unreadable_bam.write_bytes(b"BAM\x01")
    unreadable_bam.chmod(0o000)
    try:
        with pytest.raises(PermissionError, match="not readable"):
            validate_samplesheet_data([{"id": "S1", "filepath": str(unreadable_bam)}])
    finally:
        unreadable_bam.chmod(0o644)

    bad_yaml = tmp_path / "bad.yaml"
    bad_yaml.write_text("- id: S1\n  filepath: a.bam\n  ff: 12.5\n", encoding="utf-8")
    res_bad = runner.invoke(app, ["lint", str(bad_yaml)])
    assert res_bad.exit_code == 1

    txt_file = tmp_path / "legacy.txt"
    txt_file.write_text("ID\tfilepath\n", encoding="utf-8")
    with pytest.raises(ValueError, match="Unsupported samplesheet format"):
        load_samplesheet(txt_file)


def test_alignment_binning_and_bed_rejection(tmp_path: Path) -> None:
    bed_file = tmp_path / "legacy.bed"
    bed_file.write_text("chr1\t0\t5000\t100\n", encoding="utf-8")
    with pytest.raises(ValueError, match="Only aligned read files"):
        extract_sample_y_profile(bed_file)

    bam_path = tmp_path / "sample.bam"
    refs = [("chr1", 10000), ("chrY", 10000)]
    bin_reads = [
        (0, 100, 30, 0, 200),
        (0, 5100, 30, 0, 200),
        (1, 200, 30, 0, 16),
        (1, 5200, 30, 0, 16),
        (1, 200, 30, 0x400, 50),  # duplicate flag -> skipped
        (1, 200, 0, 0, 50),  # mapq 0 -> skipped when mapq=1
    ]
    _write_synthetic_bam_from_bin_counts(bam_path, refs, bin_reads)

    # Test native BAM reader directly
    ref_lengths, bin_counts = _read_bam_native(
        bam_path, binsize=5000, mapq=1, target_chrs=FFY_TARGET_CHRS
    )
    assert ref_lengths == {"1": 10000, "Y": 10000}
    assert np.array_equal(bin_counts["1"], [200, 200])
    assert np.array_equal(bin_counts["Y"], [16, 16])

    # Test extract_sample_y_profile & run_ffy_single (0-1 fraction)
    prof_bam = extract_sample_y_profile(bam_path, binsize=5000, mapq=1)
    assert np.allclose(prof_bam.y_norm_bins, [0.08, 0.08])
    bam_sex, bam_ffy = run_ffy_single(bam_path, binsize=5000, mapq=1)
    assert bam_sex == "M"
    assert np.isclose(bam_ffy, 0.16)


def _create_synthetic_bam_cohort(
    base_dir: Path, n_samples: int = 24, binsize: int = 5000
) -> tuple[Path, Path, list[str], list[float]]:
    """Create synthetic BAM files and a YAML samplesheet for end-to-end FFY + train + predict testing."""
    bams_dir = base_dir / "bams"
    bams_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(123)

    n_chr1_bins = 20
    n_chr2_bins = 20
    n_chr13_bins = 4
    n_chrx_bins = 8
    n_chry_bins = 10
    refs = [
        ("chr1", n_chr1_bins * binsize),
        ("chr2", n_chr2_bins * binsize),
        ("chr13", n_chr13_bins * binsize),
        ("chrX", n_chrx_bins * binsize),
        ("chrY", n_chry_bins * binsize),
    ]

    n_auto_feat = n_chr1_bins + n_chr2_bins
    latent_signal = rng.uniform(-1.2, 1.2, size=n_auto_feat)
    base_depth = 400.0

    true_sexes: list[str] = []
    true_ffs: list[float] = []
    sheet_rows: list[dict[str, str | float]] = []

    for i in range(n_samples):
        sex = "M" if i < 18 else "F"
        ff = float(rng.uniform(0.04, 0.22)) if sex == "M" else 0.0
        true_sexes.append(sex)
        true_ffs.append(ff)

        bin_reads: list[tuple[int, int, int, int, int]] = []
        # chr1 & chr2 (autosomal features encoding latent fetal signal)
        for b in range(n_chr1_bins):
            ratio_log2 = ff * latent_signal[b] + float(rng.normal(0.0, 0.005))
            cnt = max(1, round(base_depth * (2.0**ratio_log2)))
            bin_reads.append((0, b * binsize + 100, 30, 0, cnt))
        for b in range(n_chr2_bins):
            if b == 0:
                # Always 0 reads in chr2 bin 0 -> tests NaN filtering/handling
                continue
            ratio_log2 = ff * latent_signal[n_chr1_bins + b] + float(
                rng.normal(0.0, 0.005)
            )
            cnt = max(1, round(base_depth * (2.0**ratio_log2)))
            bin_reads.append((1, b * binsize + 100, 30, 0, cnt))
        # chr13 (excluded autosome)
        for b in range(n_chr13_bins):
            bin_reads.append((2, b * binsize + 100, 30, 0, int(base_depth)))
        # chrX (copy-number ratio decreases linearly with male FF)
        x_scale = (1.0 - 0.5 * ff) if sex == "M" else 1.0
        for b in range(n_chrx_bins):
            cnt = max(
                1,
                round(base_depth * (x_scale + float(rng.normal(0.0, 0.002)))),
            )
            bin_reads.append((3, b * binsize + 100, 30, 0, cnt))
        # chrY (bins 2..7 informative for male FF, bins 0..1 baseline noise)
        for b in range(n_chry_bins):
            if 2 <= b < 8:
                y_frac = (0.002 + 0.5 * ff) if sex == "M" else 0.002
            else:
                y_frac = 0.001
            cnt = max(0, round(base_depth * y_frac))
            bin_reads.append((4, b * binsize + 100, 30, 0, cnt))

        bam_path = bams_dir / f"sample_{i}.bam"
        _write_synthetic_bam_from_bin_counts(bam_path, refs, bin_reads)
        sheet_rows.append(
            {
                "id": f"ID_{i}",
                "filepath": bam_path.as_posix(),
                "sex": sex,
                "ff": round(ff, 6),
            }
        )

    samplesheet_path = base_dir / "samplesheet.yaml"
    samplesheet_path.write_text(
        yaml.safe_dump(sheet_rows, sort_keys=False), encoding="utf-8"
    )
    return samplesheet_path, bams_dir / "sample_0.bam", true_sexes, true_ffs


def test_ffy_cohort_lasso_onnx_and_cli(tmp_path: Path) -> None:
    samplesheet_path, sample0_bam, true_sexes, true_ffs = _create_synthetic_bam_cohort(
        tmp_path, n_samples=24, binsize=5000
    )
    entries = load_samplesheet(samplesheet_path)
    unlabeled_sheet = tmp_path / "unlabeled.yaml"
    write_samplesheet(
        [type(e)(id=e.id, filepath=e.filepath, sex=None, ff=None) for e in entries],
        outfile=unlabeled_sheet,
    )

    out_sheet = tmp_path / "samplesheet_out.json"
    ffy_onnx = tmp_path / "ffy_model.onnx"

    res_cli = runner.invoke(
        app,
        [
            "ffy",
            "-s",
            str(unlabeled_sheet),
            "--outfile",
            str(out_sheet),
            "--save-model",
            str(ffy_onnx),
            "--binsize",
            "5000",
            "--cpus",
            "2",
            "--seed",
            "1",
        ],
    )
    assert res_cli.exit_code == 0, res_cli.stdout
    assert out_sheet.is_file()
    assert ffy_onnx.is_file()

    out_entries = load_samplesheet(out_sheet)
    assert [e.sex for e in out_entries] == true_sexes
    est_ffs = np.array([float(e.ff or 0.0) for e in out_entries])
    assert np.all((est_ffs >= 0.0) & (est_ffs <= 1.0))
    assert np.allclose(est_ffs[18:], 0.0)
    assert pearson_cor(est_ffs[:18], np.array(true_ffs[:18])) > 0.98

    meta = load_ffy_onnx_metadata(ffy_onnx)
    assert len(meta.bin_loci) == 10
    assert meta.binsize == 5000

    norm_args = _normalize_argv(
        [
            "ffy",
            "--infile",
            str(sample0_bam),
            "--model",
            str(ffy_onnx),
            "--device",
            "cpu",
            "--json",
        ]
    )
    res_single = runner.invoke(app, norm_args)
    assert res_single.exit_code == 0
    payload = json.loads(res_single.stdout.strip())
    assert payload["sex"] == "M"
    assert 0.0 < payload["FFY"] <= 1.0
    assert np.isclose(payload["FFY"], true_ffs[0], atol=0.02)

    fem_bam = entries[20].filepath
    fem_sex, fem_ffy = run_ffy_single(fem_bam, model_path=ffy_onnx, device_name="cpu")
    assert fem_sex == "F"
    assert fem_ffy == 0.0


def test_end_to_end_pipeline_and_cli(tmp_path: Path) -> None:
    samplesheet_path, sample_bam, _, _ = _create_synthetic_bam_cohort(
        tmp_path, n_samples=24, binsize=5000
    )
    outdir_nn = tmp_path / "out_nn"

    onnx_path = run_train(
        samplesheet_path=samplesheet_path,
        outdir=outdir_nn,
        n_feat=3,
        hidden=2,
        cpus=2,
        olm=False,
        binsize=5000,
        device_name="cpu",
        seed=1,
    )
    assert onnx_path.is_file()
    assert (outdir_nn / "overall_performance.png").is_file()
    assert (outdir_nn / "FFX.png").is_file()
    assert (outdir_nn / "training_statistics.txt").is_file()
    assert (outdir_nn / "training_repeats" / "repeat_1.png").is_file()
    assert (outdir_nn / "training_repeats" / "repeat_10.png").is_file()

    stats_text = (outdir_nn / "training_statistics.txt").read_text(encoding="utf-8")
    assert "PREFACE - PREdict FetAl ComponEnt" in stats_text
    assert "Overall correlation (r):" in stats_text
    assert "FF < 0.1 mean absolute error (MAE):" in stats_text

    # Test CLI predict in default, bare --json, and --json <file> modes
    res_std = runner.invoke(
        app,
        [
            "predict",
            "--infile",
            str(sample_bam),
            "--model",
            str(onnx_path),
            "--device",
            "cpu",
        ],
    )
    assert res_std.exit_code == 0
    assert "FFX =" in res_std.stdout
    assert "PREFACE =" in res_std.stdout
    assert "%" not in res_std.stdout

    norm_args = _normalize_argv(
        ["predict", "--infile", str(sample_bam), "--model", str(onnx_path), "--json"]
    )
    res_json_stdout = runner.invoke(app, norm_args)
    assert res_json_stdout.exit_code == 0
    parsed = json.loads(res_json_stdout.stdout.strip())
    assert "FFX" in parsed and "PREFACE" in parsed
    assert 0.0 < parsed["PREFACE"] < 1.0

    json_file = tmp_path / "pred.json"
    res_json_file = runner.invoke(
        app,
        [
            "predict",
            "--infile",
            str(sample_bam),
            "--model",
            str(onnx_path),
            "--json",
            str(json_file),
        ],
    )
    assert res_json_file.exit_code == 0
    assert json_file.is_file()
    assert json.loads(json_file.read_text(encoding="utf-8")) == parsed

    # Test CLI train with --olm and auto-FFY estimation on unlabeled samplesheet
    entries = load_samplesheet(samplesheet_path)
    unlabeled_sheet = tmp_path / "unlabeled_train.yaml"
    write_samplesheet(
        [type(e)(id=e.id, filepath=e.filepath, sex=None, ff=None) for e in entries],
        outfile=unlabeled_sheet,
    )
    outdir_olm = tmp_path / "out_olm"
    res_train_olm = runner.invoke(
        app,
        [
            "train",
            "-s",
            str(unlabeled_sheet),
            "--outdir",
            str(outdir_olm),
            "--nfeat",
            "3",
            "--olm",
            "--binsize",
            "5000",
            "--ffy-binsize",
            "5000",
            "--cpus",
            "2",
            "--device",
            "cpu",
        ],
    )
    assert res_train_olm.exit_code == 0, res_train_olm.stdout
    assert (outdir_olm / "model.onnx").is_file()
    ffx_val, pref_val = run_predict(
        sample_bam, outdir_olm / "model.onnx", device_name="cpu"
    )
    assert np.isfinite(ffx_val)
    assert np.isfinite(pref_val)


def test_logging_and_cli_loglevel(tmp_path: Path) -> None:
    logger = configure_logging("debug")
    assert logger.level == 10
    logger = configure_logging(LogLevel.WARNING)
    assert logger.level == 30
    with pytest.raises(ValueError, match="Unsupported log level"):
        configure_logging("TRACE")

    bam_path = tmp_path / "log_sample.bam"
    refs = [("chr1", 10000), ("chrY", 10000)]
    bin_reads = [
        (0, 100, 30, 0, 200),
        (0, 5100, 30, 0, 200),
        (1, 200, 30, 0, 16),
        (1, 5200, 30, 0, 16),
    ]
    _write_synthetic_bam_from_bin_counts(bam_path, refs, bin_reads)

    res_debug = runner.invoke(
        app,
        [
            "ffy",
            "--infile",
            str(bam_path),
            "--loglevel",
            "DEBUG",
        ],
    )
    assert res_debug.exit_code == 0
    assert "sex = M" in res_debug.stdout
    assert "Single-sample FFY" in res_debug.output

    res_root_err = runner.invoke(
        app,
        [
            "--loglevel",
            "ERROR",
            "ffy",
            "--infile",
            str(bam_path),
        ],
    )
    assert res_root_err.exit_code == 0
    assert "sex = M" in res_root_err.stdout
    assert "Single-sample FFY" not in res_root_err.output

    res_invalid = runner.invoke(
        app,
        [
            "ffy",
            "--infile",
            str(bam_path),
            "--loglevel",
            "INVALID",
        ],
    )
    assert res_invalid.exit_code != 0


def _square_or_fail_worker(val: int) -> int:
    if val < 0:
        raise FileNotFoundError(f"Missing sample {val}")
    return val * val


def test_parallel_load_samples(tmp_path: Path) -> None:
    assert parallel_load_samples([], _square_or_fail_worker, cpus=2) == []

    items = [(f"S_{i}", i) for i in range(8)]
    res_seq = parallel_load_samples(items, _square_or_fail_worker, cpus=1)
    res_par = parallel_load_samples(items, _square_or_fail_worker, cpus=2)
    assert res_seq == [i * i for i in range(8)]
    assert res_par == res_seq

    bad_items = [("S_0", 2), ("S_bad", -1), ("S_2", 4)]
    with pytest.raises(FileNotFoundError, match="Missing sample -1"):
        parallel_load_samples(bad_items, _square_or_fail_worker, cpus=2)

    samplesheet_path, _, _, _ = _create_synthetic_bam_cohort(
        tmp_path, n_samples=8, binsize=5000
    )
    entries = load_samplesheet(samplesheet_path, rng=np.random.default_rng(1))
    ds_1 = load_training_dataset(entries, binsize=5000, cpus=1)
    ds_2 = load_training_dataset(entries, binsize=5000, cpus=2)
    assert ds_1.possible_features == ds_2.possible_features
    assert np.allclose(ds_1.training_frame, ds_2.training_frame)
    assert np.allclose(ds_1.x_ratios, ds_2.x_ratios)
