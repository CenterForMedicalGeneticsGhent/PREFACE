"""Typer CLI entrypoint for PREFACE."""

from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import Annotated

import typer

from preface import __version__

app = typer.Typer(
    name="preface",
    help="PREFACE (PREdict FetAl ComponEnt) - Predict fetal fraction from shallow WGS NIPT data.",
    add_completion=False,
    no_args_is_help=True,
)


def _version_callback(value: bool) -> None:
    if value:
        typer.echo(__version__)
        raise typer.Exit()


@app.callback()
def common_options(
    version: Annotated[
        bool,
        typer.Option(
            "--version",
            "--v",
            "-v",
            help="Show version and exit.",
            callback=_version_callback,
            is_eager=True,
        ),
    ] = False,
) -> None:
    """PREFACE command line interface."""


@app.command("train")
def train_cmd(
    config: Annotated[
        Path,
        typer.Option(
            "--config", help="Path to tab-separated training config.txt file."
        ),
    ],
    outdir: Annotated[
        Path,
        typer.Option(
            "--outdir", help="Output directory for trained model.onnx and plots."
        ),
    ],
    nfeat: Annotated[
        int,
        typer.Option(
            "--nfeat", help="Number of principal components to use during modeling."
        ),
    ] = 50,
    hidden: Annotated[
        int,
        typer.Option("--hidden", help="Number of hidden units used in neural network."),
    ] = 2,
    cpus: Annotated[
        int,
        typer.Option("--cpus", help="Number of CPU threads for parallel processing."),
    ] = os.cpu_count() or 1,
    femprop: Annotated[
        bool,
        typer.Option(
            "--femprop",
            help="Include female fetuses in supervised learning phase (requires proportional FF labels).",
        ),
    ] = False,
    olm: Annotated[
        bool,
        typer.Option(
            "--olm", help="Use an ordinary linear model instead of a neural network."
        ),
    ] = False,
    noskewcorrect: Annotated[
        bool,
        typer.Option(
            "--noskewcorrect",
            help="Disable linear skew correction on cross-validated predictions.",
        ),
    ] = False,
    device: Annotated[
        str,
        typer.Option("--device", help="Compute device to use (auto, cpu, cuda, mps)."),
    ] = "auto",
    seed: Annotated[
        int,
        typer.Option("--seed", help="Random seed for reproducibility."),
    ] = 1,
) -> None:
    """Train a PREFACE model and export it to ONNX."""
    from preface.pipeline import run_train

    try:
        run_train(
            config_path=config,
            outdir=outdir,
            n_feat=nfeat,
            hidden=hidden,
            cpus=cpus,
            femprop=femprop,
            olm=olm,
            noskewcorrect=noskewcorrect,
            device_name=device,
            seed=seed,
        )
    except (FileNotFoundError, OSError, ValueError, RuntimeError) as exc:
        typer.echo(str(exc))
        raise typer.Exit(code=1) from exc


@app.command("predict")
def predict_cmd(
    infile: Annotated[
        Path,
        typer.Option(
            "--infile", help="Path to sample copy number alteration .bed file."
        ),
    ],
    model: Annotated[
        Path,
        typer.Option("--model", help="Path to trained model.onnx file."),
    ],
    json_out: Annotated[
        str | None,
        typer.Option(
            "--json",
            help="Output predictions in JSON format. Pass empty string '' (or bare --json) for stdout, or a file path.",
        ),
    ] = None,
    device: Annotated[
        str,
        typer.Option(
            "--device", help="Execution device for ONNX Runtime (auto, cpu, cuda, mps)."
        ),
    ] = "auto",
) -> None:
    """Predict fetal fraction (FFX and PREFACE) for an input BED file using a trained ONNX model."""
    from preface.pipeline import run_predict

    try:
        run_predict(
            infile=infile,
            model_path=model,
            json_arg=json_out,
            device_name=device,
        )
    except (FileNotFoundError, OSError, ValueError, RuntimeError) as exc:
        typer.echo(str(exc))
        raise typer.Exit(code=1) from exc


@app.command("ffy")
def ffy_cmd(
    config: Annotated[
        Path | None,
        typer.Option(
            "--config",
            help="Path to cohort TSV file (ID, filepath, optional reads_filepath/gender).",
        ),
    ] = None,
    infile: Annotated[
        Path | None,
        typer.Option(
            "--infile",
            help="Path to single-sample BAM/CRAM or binned read count BED/TSV file.",
        ),
    ] = None,
    outfile: Annotated[
        Path | None,
        typer.Option(
            "--outfile",
            help="Output path for generated cohort config.txt (defaults to stdout).",
        ),
    ] = None,
    model: Annotated[
        Path | None,
        typer.Option(
            "--model",
            help="Path to calibrated FFY ONNX model for single-sample inference.",
        ),
    ] = None,
    save_model: Annotated[
        Path | None,
        typer.Option(
            "--save-model",
            help="Output path to save calibrated FFY Lasso ONNX model in cohort mode.",
        ),
    ] = None,
    male_ref: Annotated[
        list[Path] | None,
        typer.Option(
            "--male-ref",
            help="Paths to 100% male reference BAM/CRAM/BED samples (LBs) per PMC6771918 Eq. 1.",
        ),
    ] = None,
    binsize: Annotated[
        int,
        typer.Option(
            "--binsize",
            help="Bin width in base pairs when binning aligned BAM/CRAM reads.",
        ),
    ] = 5000,
    mapq: Annotated[
        int,
        typer.Option(
            "--mapq",
            help="Minimum mapping quality for aligned BAM/CRAM reads.",
        ),
    ] = 1,
    yfrac: Annotated[
        float | None,
        typer.Option(
            "--yfrac",
            help="Explicit normalized chrY read depth threshold for male/female classification.",
        ),
    ] = None,
    lasso_lambda: Annotated[
        float,
        typer.Option(
            "--lambda",
            help="Lasso regularization penalty lambda for chrY bin weighting (PMC6771918 Eq. 3).",
        ),
    ] = 1e-4,
    reference: Annotated[
        Path | None,
        typer.Option(
            "--reference",
            help="Reference FASTA file path for CRAM alignment decoding.",
        ),
    ] = None,
    json_out: Annotated[
        str | None,
        typer.Option(
            "--json",
            help="Output single-sample FFY in JSON format. Pass bare --json for stdout, or a file path.",
        ),
    ] = None,
    cpus: Annotated[
        int,
        typer.Option("--cpus", help="Number of CPU threads for parallel processing."),
    ] = os.cpu_count() or 1,
    device: Annotated[
        str,
        typer.Option(
            "--device",
            help="Execution device for ONNX Runtime (auto, cpu, cuda, mps).",
        ),
    ] = "auto",
    seed: Annotated[
        int,
        typer.Option("--seed", help="Random seed for reproducibility."),
    ] = 1,
) -> None:
    """Calculate baseline chromosome Y fetal fraction (FFY) from binned read counts or aligned reads."""
    from preface.ffy import run_ffy_cohort, run_ffy_single

    try:
        if (config is None) == (infile is None):
            raise ValueError(
                "Provide exactly one of --config (cohort mode) or --infile (single-sample mode)."
            )
        if config is not None:
            run_ffy_cohort(
                config_path=config,
                outfile=outfile,
                male_ref_paths=male_ref,
                save_model_path=save_model,
                binsize=binsize,
                mapq=mapq,
                yfrac=yfrac,
                lasso_lambda=lasso_lambda,
                reference=reference,
                cpus=cpus,
                seed=seed,
            )
        elif infile is not None:
            run_ffy_single(
                infile=infile,
                model_path=model,
                binsize=binsize,
                mapq=mapq,
                yfrac=yfrac,
                reference=reference,
                json_arg=json_out,
                device_name=device,
            )
    except (FileNotFoundError, OSError, ValueError, RuntimeError) as exc:
        typer.echo(str(exc))
        raise typer.Exit(code=1) from exc


def _normalize_argv(argv: list[str]) -> list[str]:
    """Support PREFACE.R's bare `--json` flag (without value) alongside `--json path/to/file.json`."""
    normalized: list[str] = []
    i = 0
    while i < len(argv):
        arg = argv[i]
        if arg == "--json" and (i + 1 >= len(argv) or argv[i + 1].startswith("-")):
            normalized.extend(["--json", ""])
            i += 1
            continue
        normalized.append(arg)
        i += 1
    return normalized


def main() -> None:
    """CLI entrypoint normalizing legacy R-style optional flag arguments."""
    sys.argv[1:] = _normalize_argv(sys.argv[1:])
    app()


if __name__ == "__main__":
    main()
