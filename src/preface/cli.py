"""Typer CLI entrypoint for PREFACE."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Annotated

import typer

from preface import __version__
from preface.logging_utils import LogLevel, configure_logging, get_logger
from preface.pipeline import TrainTarget

logger = get_logger(__name__)

LogLevelOption = Annotated[
    LogLevel | None,
    typer.Option(
        "--loglevel",
        case_sensitive=False,
        help="Logging verbosity level (DEBUG, INFO, WARNING, ERROR, CRITICAL).",
    ),
]

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
    loglevel: Annotated[
        LogLevel,
        typer.Option(
            "--loglevel",
            case_sensitive=False,
            help="Logging verbosity level (DEBUG, INFO, WARNING, ERROR, CRITICAL).",
        ),
    ] = LogLevel.INFO,
) -> None:
    """PREFACE command line interface."""
    configure_logging(loglevel)


@app.command("lint")
def lint_cmd(
    samplesheet: Annotated[
        Path,
        typer.Argument(
            help="Path to YAML or JSON samplesheet file to validate against the JSON Schema.",
        ),
    ],
    loglevel: LogLevelOption = None,
) -> None:
    """Validate a YAML or JSON samplesheet against the PREFACE JSON Schema."""
    from preface.io import load_samplesheet

    if loglevel is not None:
        configure_logging(loglevel)

    try:
        entries = load_samplesheet(samplesheet)
        logger.info(
            "Samplesheet '%s' is valid (%d samples).", samplesheet, len(entries)
        )
    except (FileNotFoundError, OSError, ValueError) as exc:
        logger.error("%s", exc)
        raise typer.Exit(code=1) from exc


@app.command("train")
def train_cmd(
    target: Annotated[
        TrainTarget,
        typer.Argument(
            case_sensitive=False,
            help="Training target: 'all' (unified FFY + PREFACE), 'ffy' (FFY only), or 'preface' / 'ffx' (PREFACE + FFX only).",
        ),
    ] = TrainTarget.ALL,
    samplesheet: Annotated[
        Path,
        typer.Option(
            "-s",
            "--samplesheet",
            help="Path to training samplesheet (.yaml, .yml, or .json).",
        ),
    ] = ...,
    prefix: Annotated[
        Path,
        typer.Option(
            "-p",
            "--prefix",
            help="Output prefix for trained <prefix>.onnx model, samplesheet, statistics, and plots.",
        ),
    ] = ...,
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
            help="Include female fetuses in supervised learning phase (requires proportional ff labels).",
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
    binsize: Annotated[
        int,
        typer.Option(
            "--binsize",
            help="Genomic bin width in base pairs for autosomal/chrX feature extraction.",
        ),
    ] = 100_000,
    ffy_binsize: Annotated[
        int,
        typer.Option(
            "--ffy-binsize",
            help="Genomic bin width in base pairs on chrY when training FFY.",
        ),
    ] = 5_000,
    mapq: Annotated[
        int,
        typer.Option(
            "--mapq",
            help="Minimum mapping quality for aligned BAM/CRAM/SAM reads.",
        ),
    ] = 1,
    yfrac: Annotated[
        float | None,
        typer.Option(
            "--yfrac",
            help="Explicit normalized chrY read depth threshold for male/female sex classification.",
        ),
    ] = None,
    lasso_lambda: Annotated[
        float,
        typer.Option(
            "--lambda",
            help="Lasso regularization penalty lambda for chrY bin weighting (PMC6771918 Eq. 3).",
        ),
    ] = 1e-4,
    no_female_gate: Annotated[
        bool,
        typer.Option(
            "--no-female-gate",
            help="Disable gating FFY predictions to 0.0 for female samples.",
        ),
    ] = False,
    reference: Annotated[
        Path | None,
        typer.Option(
            "--reference",
            help="Reference FASTA file path for CRAM alignment decoding.",
        ),
    ] = None,
    male_ref: Annotated[
        list[Path] | None,
        typer.Option(
            "--male-ref",
            help="Paths to 100% male reference BAM/CRAM/SAM samples (LBs) when training FFY.",
        ),
    ] = None,
    cache_counts: Annotated[
        Path | None,
        typer.Option(
            "--cache-counts",
            help="Optional .npz file path to save/reuse single-pass binned read counts.",
        ),
    ] = None,
    device: Annotated[
        str,
        typer.Option("--device", help="Compute device to use (auto, cpu, cuda, mps)."),
    ] = "auto",
    seed: Annotated[
        int,
        typer.Option("--seed", help="Random seed for reproducibility."),
    ] = 1,
    loglevel: LogLevelOption = None,
) -> None:
    """Train FFY, PREFACE/FFX, or unified FFY+PREFACE models in a single pass and export to ONNX."""
    from preface.pipeline import run_train

    if loglevel is not None:
        configure_logging(loglevel)

    try:
        run_train(
            samplesheet_path=samplesheet,
            prefix=prefix,
            target=target,
            n_feat=nfeat,
            hidden=hidden,
            cpus=cpus,
            femprop=femprop,
            olm=olm,
            noskewcorrect=noskewcorrect,
            binsize=binsize,
            ffy_binsize=ffy_binsize,
            mapq=mapq,
            yfrac=yfrac,
            lasso_lambda=lasso_lambda,
            no_female_gate=no_female_gate,
            reference=reference,
            male_ref_paths=male_ref,
            cache_counts=cache_counts,
            device_name=device,
            seed=seed,
        )
    except (FileNotFoundError, OSError, ValueError, RuntimeError) as exc:
        logger.error("%s", exc)
        raise typer.Exit(code=1) from exc


@app.command("predict")
def predict_cmd(
    infile: Annotated[
        Path,
        typer.Option(
            "--infile", help="Path to sample aligned read file (.bam, .cram, .sam)."
        ),
    ],
    model: Annotated[
        Path,
        typer.Option("--model", help="Path to trained .onnx model file."),
    ],
    prefix: Annotated[
        Path | None,
        typer.Option(
            "-p",
            "--prefix",
            help="Output prefix for JSON predictions (0..1 floats). Omit to print plain-text percentages to stdout.",
        ),
    ] = None,
    mapq: Annotated[
        int,
        typer.Option(
            "--mapq",
            help="Minimum mapping quality for aligned BAM/CRAM/SAM reads.",
        ),
    ] = 1,
    reference: Annotated[
        Path | None,
        typer.Option(
            "--reference",
            help="Reference FASTA file path for CRAM alignment decoding.",
        ),
    ] = None,
    no_female_gate: Annotated[
        bool | None,
        typer.Option(
            "--no-female-gate/--female-gate",
            help="Disable or enable gating FFY predictions to 0.0 for female samples (defaults to model setting).",
        ),
    ] = None,
    device: Annotated[
        str,
        typer.Option(
            "--device", help="Execution device for ONNX Runtime (auto, cpu, cuda, mps)."
        ),
    ] = "auto",
    loglevel: LogLevelOption = None,
) -> None:
    """Predict fetal fraction (sex, FFY, FFX, and PREFACE) for an aligned read file using an ONNX model."""
    from preface.pipeline import run_predict

    if loglevel is not None:
        configure_logging(loglevel)

    try:
        run_predict(
            infile=infile,
            model_path=model,
            prefix=prefix,
            mapq=mapq,
            reference=reference,
            no_female_gate=no_female_gate,
            device_name=device,
        )
    except (FileNotFoundError, OSError, ValueError, RuntimeError) as exc:
        logger.error("%s", exc)
        raise typer.Exit(code=1) from exc


def main() -> None:
    """CLI entrypoint."""
    app()


if __name__ == "__main__":
    main()
