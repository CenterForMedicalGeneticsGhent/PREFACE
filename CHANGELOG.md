# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0] - 2026-11-23

### Breaking Changes

- Rewrote PREFACE in Python (requires Python >= 3.14), replacing the legacy `PREFACE.R` script with a PyTorch training backend, ONNX Runtime inference/export, and Typer CLI (`preface`) (#19).
- Replaced BED input files and tabular text configuration files with direct aligned sequencing read inputs (`.bam`, `.cram`, `.sam`) and a unified YAML/JSON samplesheet validated against `src/preface/samplesheet.schema.json` (#23).
- Unified `FFY`, `FFX`, and `PREFACE` training and prediction into single-pass workflows (`preface train [all|ffy|preface|ffx]` and `preface predict`) exporting and consuming a single unified `.onnx` model (#24).

### Added

- Automated baseline chromosome Y fetal fraction (`FFY`) training with two-component Gaussian mixture model sex classification, 10-fold cross-validated non-negative Lasso over informative `chrY` bins, optional female gating (`--no-female-gate`), and structured logging (`--loglevel`) (#22, #24).
- `preface lint` subcommand for JSON Schema samplesheet validation, file existence/readability checks, and optional alignment index paths (`.gzi`, `.bai`, `.csi`, `.crai`) (#23).
- Single-pass multi-binsize alignment binning with optional `.npz` count caching (`--cache-counts`), parallel sample loading (`--cpus`), and accelerator selection (`--device auto|cpu|cuda|mps`) (#19, #24).
- Example alignment files (`examples/sample.bam`, `examples/sample.bam.bai`) and YAML/JSON samplesheets (`examples/samplesheet.yaml`, `examples/samplesheet.json`) (#23).
- Pixi workspace and `pixi-build-python` packaging configuration for `linux-64` and `osx-arm64` (#18).
- Multi-stage `Dockerfile`, GitHub Actions CI/CD workflow (lint, test, wheel build, and Quay.io container publishing), Ruff formatting/linting, `prek` pre-commit configuration, and `pytest` test suite (#19, #20, #21, #26).

### Fixed

- Added `[tool.pixi.package]` metadata in `pyproject.toml` to support `pixi global install`.
- Removed the `dev` suffix from `preface.__version__` so CLI, package metadata, and exported ONNX metadata consistently report `0.2.0`.
- Cleaned up CI workflow artifact retention (#26).

## [0.1.2] - 2021-01-07

### Changed

- Decreased PCA runtime by adding the `irlba` dependency (#3).
- Increased PCA components for plotting purposes (#3).

### Fixed

- Improved robustness against small training sets for testing purposes (#3).

## [0.1.1] - 2019-06-27

### Changed

- Enabled compatibility with `neuralnet` v1.44.2.
- Updated plot symbols, README documentation, output formatting, and added publication reference.

### Fixed

- Fixed `--nfeat` parameter handling in plotting.

## [0.1.0] - 2019-02-22

### Added

- Initial non-dev release of PREFACE (R implementation) for fetal fraction prediction from shallow whole-genome NIPT sequencing data.
- Added `--hidden` neural network option, improved non-convergence handling, and support for fetal fraction measures other than `FFY`.

[0.2.0]: https://github.com/CenterForMedicalGeneticsGhent/PREFACE/compare/v0.1.2...v0.2.0
[0.1.2]: https://github.com/CenterForMedicalGeneticsGhent/PREFACE/compare/v0.1.1...v0.1.2
[0.1.1]: https://github.com/CenterForMedicalGeneticsGhent/PREFACE/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/CenterForMedicalGeneticsGhent/PREFACE/releases/tag/v0.1.0
