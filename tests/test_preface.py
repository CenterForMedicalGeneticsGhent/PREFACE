"""Unit and integration tests for the PREFACE Python rewrite."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import polars as pl
import pytest
import torch
from typer.testing import CliRunner

from preface.cli import _normalize_argv, app
from preface.hardware import configure_hardware, get_onnx_providers, resolve_device
from preface.math_utils import get_m_diff, get_sd_diff, ols_fit, pca_irlba_torch, pearson_cor, rlm_huber
from preface.models import (
    export_pipeline_to_onnx,
    load_onnx_metadata,
    predict_model,
    run_onnx_inference,
    train_neural_net,
    train_olm,
)
from preface.pipeline import run_predict, run_train

runner = CliRunner()


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

    t = np.linspace(1.0, 20.0, 40)
    y = 3.5 + 2.2 * t
    intercept, slope = ols_fit(t, y)
    assert np.isclose(intercept, 3.5)
    assert np.isclose(slope, 2.2)

    # Inject an outlier for robust linear regression
    y_Outlier = y.copy()
    y_Outlier[0] = 150.0
    r_int, r_slope, wr = rlm_huber(t, y_Outlier)
    assert np.isclose(r_int, 3.5, atol=0.3)
    assert np.isclose(r_slope, 2.2, atol=0.05)
    assert 0.9 < wr <= 1.0

    assert np.isclose(get_m_diff(t, t + 2.0, use_abs=True), 2.0)
    assert np.isclose(get_m_diff(t, t + 2.0, use_abs=False), -2.0)
    assert np.isclose(get_sd_diff(t, t + 2.0), 0.0)
    assert np.isclose(pearson_cor(t, y), 1.0)


def test_models_and_onnx_roundtrip(tmp_path: Path) -> None:
    rng = np.random.default_rng(1)
    n_samples, n_raw_features, n_feat = 40, 25, 4
    raw_x = rng.normal(size=(n_samples, n_raw_features))
    mean_features = np.mean(raw_x, axis=0)
    pca_res = pca_irlba_torch(raw_x, n_components=n_feat, device=torch.device("cpu"))

    x_pc = pca_res.x[:, :n_feat]
    true_y = 8.0 + 1.5 * x_pc[:, 0] - 0.8 * x_pc[:, 1]

    olm_model = train_olm(x_pc, true_y, device=torch.device("cpu"))
    olm_preds = predict_model(olm_model, x_pc, device=torch.device("cpu"))
    assert np.allclose(olm_preds, true_y, atol=1e-4)

    torch.manual_seed(1)
    nn_model = train_neural_net(x_pc, true_y, hidden=2, device=torch.device("cpu"))
    nn_preds = predict_model(nn_model, x_pc, device=torch.device("cpu"))
    assert pearson_cor(nn_preds, true_y) > 0.95

    onnx_path = tmp_path / "model.onnx"
    feature_names = [f"1:{i * 100 + 1}-{(i + 1) * 100}" for i in range(n_raw_features)]
    export_pipeline_to_onnx(
        predictor=nn_model,
        mean_features=mean_features,
        pca_center=pca_res.center,
        pca_rotation=pca_res.rotation[:, :n_feat],
        the_intercept=0.5,
        the_slope=1.1,
        the_intercept_x=0.95,
        the_slope_x=0.005,
        possible_features=feature_names,
        onnx_path=onnx_path,
    )

    meta = load_onnx_metadata(onnx_path)
    assert meta.possible_features == feature_names
    assert np.isclose(meta.the_intercept_x, 0.95)
    assert np.isclose(meta.the_slope_x, 0.005)
    assert meta.n_feat == n_feat
    assert not meta.is_olm

    sample_with_nan = raw_x[0].copy()
    sample_with_nan[3] = np.nan
    expected_imputed = raw_x[0].copy()
    expected_imputed[3] = mean_features[3]
    expected_pc = (expected_imputed - pca_res.center) @ pca_res.rotation[:, :n_feat]
    expected_raw = float(predict_model(nn_model, expected_pc.reshape(1, -1), device=torch.device("cpu"))[0])
    expected_final = 0.5 + 1.1 * expected_raw

    onnx_pred = run_onnx_inference(onnx_path, sample_with_nan, requested_device="cpu")
    assert np.isclose(onnx_pred, expected_final, atol=1e-4)


def _create_synthetic_cohort(base_dir: Path, n_samples: int = 24) -> tuple[Path, Path]:
    """Create synthetic BED files and config.txt for end-to-end testing."""
    beds_dir = base_dir / "beds"
    beds_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(123)

    chrs = (["1"] * 40) + (["2"] * 40) + (["13"] * 5) + (["X"] * 10) + (["Y"] * 5)
    starts = [i * 100000 + 1 for i in range(len(chrs))]
    ends = [(i + 1) * 100000 for i in range(len(chrs))]

    latent_signal = rng.normal(0.0, 1.0, size=len(chrs))
    config_rows: list[str] = ["ID\tfilepath\tgender\tFF"]

    for i in range(n_samples):
        gender = "M" if i < 20 else "F"
        ff = float(rng.uniform(4.0, 18.0)) if gender == "M" else 0.0
        ratios = (ff / 100.0) * latent_signal + rng.normal(0.0, 0.01, size=len(chrs))

        # Encode X copy number ratio relationship for males
        for idx, c in enumerate(chrs):
            if c == "X":
                ratios[idx] = np.log2(max(0.5, 1.0 - 0.004 * ff + rng.normal(0.0, 0.001)))

        # Always NaN on locus 0 (simulating centromere/telomere)
        ratios[0] = np.nan

        bed_path = beds_dir / f"sample_{i}.bed"
        df = pl.DataFrame(
            {
                "chr": chrs,
                "start": starts,
                "end": ends,
                "ratio": [
                    "NaN" if np.isnan(val) else f"{val:.8f}"
                    for val in ratios
                ],
            }
        )
        df.write_csv(bed_path, separator="\t")
        config_rows.append(f"ID_{i}\t{bed_path.as_posix()}\t{gender}\t{ff}")

    config_path = base_dir / "config.txt"
    config_path.write_text("\n".join(config_rows) + "\n", encoding="utf-8")
    return config_path, beds_dir / "sample_0.bed"


def test_end_to_end_pipeline_and_cli(tmp_path: Path) -> None:
    config_path, sample_bed = _create_synthetic_cohort(tmp_path, n_samples=24)
    outdir_nn = tmp_path / "out_nn"

    onnx_path = run_train(
        config_path=config_path,
        outdir=outdir_nn,
        n_feat=3,
        hidden=2,
        cpus=2,
        olm=False,
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

    # Test CLI predict in default, bare --json, and --json <file> modes
    res_std = runner.invoke(
        app,
        ["predict", "--infile", str(sample_bed), "--model", str(onnx_path), "--device", "cpu"],
    )
    assert res_std.exit_code == 0
    assert "FFX =" in res_std.stdout
    assert "PREFACE =" in res_std.stdout

    norm_args = _normalize_argv(
        ["predict", "--infile", str(sample_bed), "--model", str(onnx_path), "--json"]
    )
    res_json_stdout = runner.invoke(app, norm_args)
    assert res_json_stdout.exit_code == 0
    parsed = json.loads(res_json_stdout.stdout.strip())
    assert "FFX" in parsed and "PREFACE" in parsed

    json_file = tmp_path / "pred.json"
    res_json_file = runner.invoke(
        app,
        [
            "predict",
            "--infile",
            str(sample_bed),
            "--model",
            str(onnx_path),
            "--json",
            str(json_file),
        ],
    )
    assert res_json_file.exit_code == 0
    assert json_file.is_file()
    assert json.loads(json_file.read_text(encoding="utf-8")) == parsed

    # Test CLI train with --olm
    outdir_olm = tmp_path / "out_olm"
    res_train_olm = runner.invoke(
        app,
        [
            "train",
            "--config",
            str(config_path),
            "--outdir",
            str(outdir_olm),
            "--nfeat",
            "3",
            "--olm",
            "--cpus",
            "2",
            "--device",
            "cpu",
        ],
    )
    assert res_train_olm.exit_code == 0
    assert (outdir_olm / "model.onnx").is_file()
    ffx_val, pref_val = run_predict(sample_bed, outdir_olm / "model.onnx", device_name="cpu")
    assert np.isfinite(ffx_val)
    assert np.isfinite(pref_val)
