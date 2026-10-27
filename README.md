# Background

A share of all cell-free DNA fragments isolated from maternal plasma during pregnancy is fetal-derived. This amount is referred to as the 'fetal fraction' and represents an important estimate during routine noninvasive prenatal testing (NIPT). Its most essential role is informing geneticists whether an assay is conclusive: if the fetal fraction is insufficient (this limit has often been debated to be 4%) claims on fetal aneuploidies cannot be made accurately. Several techniques exist to deduce this figure, but the far most require additional experimental procedures, which impede routine execution. Therefore, we set out to develop PREFACE, a software to accurately predict fetal fraction based on solely shallow-depth whole-genome sequencing data, which is the fundamental base of a default NIPT assay. In contrast to previous efforts, PREFACE enables user-friendly model training with a limited amount of retrospective data, which eliminates between-laboratory bias. For sets of roughly 1100 male NIPT samples, a cross-validated correlation of 0.9 between predictions and fetal fractions according to Y chromosomal read counts was noted (FFY). Our approach enables training with both male and unlabeled female fetuses: using our complete cohort (nfemale=2468, nmale=2723), the correlation metric reached 0.94. In addition, PREFACE provides the fetal fraction based on the copy number state of chromosome X (FFX). The presented statistics indirectly predict mixed multiple pregnancies, the source of observed events and sex chromosomal aneuploidies. All details can be found in our [corresponding paper](https://www.ncbi.nlm.nih.gov/pubmed/31219182).

# Installation

We implemented PREFACE in Python with PyTorch for model training, ONNX Runtime for model export and inference, and Typer for the command-line interface. [Pixi](https://pixi.sh) manages dependencies and environments on `linux-64` and `osx-arm64`.

```bash
pixi install
```

Verify the CLI installation:

```bash
pixi run preface --version
pixi run preface --help
```

# Manual

## Required files

### Aligned sequencing read files (`.bam`, `.cram`, `.sam`)

PREFACE operates directly on aligned sequencing read files for baseline `FFY` calculation, model training, and inference:
- Primary, non-duplicate, non-QC-failed alignments with mapping quality at or above `--mapq` (default: `1`) are binned into fixed genomic windows of size `--binsize` (`5000` bp default for `preface ffy`; `100000` bp default for `preface train`).
- `.bam` files are read via `pysam` when installed or Python's built-in BGZF/BAM parser. `.cram` files require `--reference` (FASTA) and either `pysam` or `samtools` in `PATH`.
- Chromosome names may either omit or include the `chr` prefix (`1` until `22`, `X`, and `Y`, or `chr1` until `chr22`, `chrX`, and `chrY`). Chromosomes `13`, `18`, and `21` are excluded during diploid autosomal read-depth normalization and feature extraction.

### Samplesheet (`.yaml`, `.yml`, or `.json`)

Cohort commands (`preface ffy` and `preface train`) consume a unified YAML or JSON samplesheet containing a root-level array of sample objects:

- Example: `./examples/samplesheet.yaml`
- Every samplesheet is automatically validated against the bundled JSON Schema (`src/preface/samplesheet.schema.json`).
- Each sample object supports the following properties (`additionalProperties: false`):
    - `id` (required `string`): Unique identifier for the sample.
    - `filepath` (required `string`): Path to the sample's aligned read file (`.bam`, `.cram`, or `.sam`). Must exist on disk and be readable.
    - `index` (optional `string` or `null`): Path to the alignment index file (`.gzi`, `.bai`, `.csi`, or `.crai`). Must exist on disk and be readable when provided.
    - `sex` (optional `"M"`, `"F"`, or `null`): Fetal sex (`M` for male, `F` for female). Automatically inferred by `preface ffy` (or `preface train`) when omitted or `null`.
    - `ff` (optional `number` in `[0.0, 1.0]` or `null`): Fetal fraction expressed as a proportion between `0.0` and `1.0`. Automatically estimated by `preface ffy` (or `preface train`) when omitted or `null`.

Example YAML samplesheet:

```yaml
- id: ID_0
  filepath: /path/to/ID_0.bam
  index: /path/to/ID_0.bam.bai
  sex: M
  ff: 0.10470662
- id: ID_1
  filepath: /path/to/ID_1.bam
  sex: F
  ff: 0.0
```

### Linting samplesheets (`preface lint`)

Validate any YAML or JSON samplesheet against the PREFACE JSON Schema:

```bash
pixi run preface lint examples/samplesheet.yaml
```

## Calculating baseline FFY (`preface ffy`)

`preface ffy` implements the two-stage baseline `FFY` calculation described in Section 3.4 of the [PREFACE paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC6771918/):

1. Each chromosome `Y` bin count is normalized by the sample's mappable diploid autosomal read density (excluding chromosomes `13`, `18`, and `21`). Pseudoautosomal regions (`PAR1` $\le 2.78\text{ Mb}$, `PAR2`) and `Yq12` heterochromatin ($\ge 28.8\text{ Mb}$) are masked when computing the summary chromosome `Y` depth.
2. In cohort mode, fetal `sex` (`M` or `F`) is taken from the samplesheet when present or inferred via a two-component Gaussian mixture model (or `--yfrac` threshold). Prior `FFY` is scaled between the female cohort median (`0.0`) and the 100% male baseline (either empirical adult male liquid biopsies passed via `--male-ref` or theoretical diploid autosomal scaling $2 \times Y_{\text{norm}}$). A 10-fold cross-validated non-negative Lasso model ($\lambda = 10^{-4}$) is then fit across the chromosome `Y` bins to weight informative male-specific loci and produce the final `ff` value (`0.0` for females, clamped `[0.0, 1.0]` for males).

### Cohort mode (populate `sex` and `ff` in a samplesheet)

```bash
pixi run preface ffy -s path/to/samplesheet.yaml --outfile path/to/samplesheet_with_ffy.yaml [optional arguments]
```

### Single-sample mode

```bash
pixi run preface ffy --infile path/to/sample.bam [--model path/to/ffy_model.onnx] [--json]
```

In single-sample mode, `preface ffy` prints `sex` and `FFY` (as a `0.0`–`1.0` fraction) to stdout (or writes JSON when `--json` is passed).

| Optional argument | Function |
| :--- | :--- |
| `-s`, `--samplesheet x` | Cohort samplesheet file (`.yaml`, `.yml`, or `.json`). |
| `--infile x` | Single-sample `.bam`, `.cram`, or `.sam` alignment file. |
| `--outfile x` | Destination path (`.yaml`, `.yml`, or `.json`) for the populated samplesheet in cohort mode. (default: stdout) |
| `--male-ref x` | One or more 100% male reference alignment files (liquid biopsies) for empirical denominator scaling. Defaults to theoretical diploid autosomal scaling ($2 \times Y_{\text{norm}}$) when omitted. |
| `--save-model x` | Export the calibrated cohort Lasso FFY model and baselines to an ONNX file. |
| `--model x` | Load a calibrated FFY ONNX model in single-sample mode. |
| `--binsize x` | Genomic bin size in bp when binning aligned reads. (default: `x=5000`) |
| `--mapq x` | Minimum mapping quality for aligned reads. (default: `x=1`) |
| `--yfrac x` | Explicit normalized chrY depth threshold for `M`/`F` sex classification (overrides 2-component GMM). |
| `--lambda x` | Lasso regularization penalty $\lambda$ for chrY bin weighting. (default: `x=1e-4`) |
| `--reference x` | Reference genome FASTA file for `.cram` decoding. |
| `--json [x]` | In single-sample mode, write predictions as JSON to stdout or to file `x`. |
| `--cpus x` | Number of CPU threads for parallel sample loading. (default: all available cores) |
| `--seed x` | Random seed for reproducible cross-validation. (default: `x=1`) |
| `--loglevel x` | Logging verbosity level (`DEBUG`, `INFO`, `WARNING`, `ERROR`, `CRITICAL`). (default: `x=INFO`) |

## Model training (`preface train`)

```bash
pixi run preface train -s path/to/samplesheet.yaml --outdir path/to/dir/ [optional arguments]
```

If `sex` or `ff` is omitted (`null`) in the input samplesheet, `preface train` automatically runs cohort `FFY` estimation prior to training.

| Optional argument | Function |
| :--- | :--- |
| `-s`, `--samplesheet x` | Training samplesheet file (`.yaml`, `.yml`, or `.json`). |
| `--outdir x` | Output directory for trained `model.onnx` and diagnostic plots. |
| `--nfeat x` | Number of principal components to use during modeling. (default: `x=50`) |
| `--hidden x` | Number of hidden units used in neural network. Use with caution. (default: `x=2`) |
| `--cpus x` | Number of CPU threads for parallel file loading and linear algebra. (default: all available cores) |
| `--femprop` | Include female fetuses in supervised learning phase (requires proportional `ff` labels for females). |
| `--olm` | Use an ordinary linear model instead of a neural network. |
| `--noskewcorrect` | Disable linear skew correction on cross-validated predictions. |
| `--binsize x` | Genomic bin width in bp for autosomal and chrX feature extraction. (default: `x=100000`) |
| `--ffy-binsize x` | Genomic bin width in bp on chrY when auto-estimating FFY. (default: `x=5000`) |
| `--mapq x` | Minimum mapping quality for aligned reads. (default: `x=1`) |
| `--reference x` | Reference genome FASTA file for `.cram` decoding. |
| `--male-ref x` | One or more 100% male reference alignment files when auto-estimating FFY. |
| `--device x` | Compute accelerator to use (`auto`, `cpu`, `cuda`, `mps`). (default: `auto`) |
| `--seed x` | Random seed for reproducible shuffling, cross-validation, and initialization. (default: `x=1`) |
| `--loglevel x` | Logging verbosity level (`DEBUG`, `INFO`, `WARNING`, `ERROR`, `CRITICAL`). (default: `x=INFO`) |

Training writes `model.onnx` into `--outdir` as a self-contained ONNX graph that embeds bin size metadata, NA mean-imputation, PCA centering and projection, the neural network or linear model weights, skew correction, and FFX regression parameters. It also writes `overall_performance.png`, per-fold plots under `training_repeats/repeat_1.png` through `repeat_10.png`, `FFX.png`, and `training_statistics.txt`.

## Predicting (`preface predict`)

```bash
pixi run preface predict --infile path/to/sample.bam --model path/to/model.onnx [optional arguments]
```

| Optional argument | Function |
| :--- | :--- |
| `--infile x` | Sample `.bam`, `.cram`, or `.sam` alignment file. |
| `--model x` | Trained `model.onnx` file. |
| `--mapq x` | Minimum mapping quality for aligned reads. (default: `x=1`) |
| `--reference x` | Reference genome FASTA file for `.cram` decoding. |
| `--json [x]` | Predictions are written to stdout. Use `--json` for JSON format on stdout, or provide `--json path/to/output.json` to write to a JSON file. |
| `--device x` | Execution device for ONNX Runtime (`auto`, `cpu`, `cuda`, `mps`). (default: `auto`) |
| `--loglevel x` | Logging verbosity level (`DEBUG`, `INFO`, `WARNING`, `ERROR`, `CRITICAL`). (default: `x=INFO`) |

## Model optimization

- The most important parameter is `--nfeat`:
    - It represents the number of principal components (PCs) that will be used as features during model training. Depending on the bin size and the number of training samples, it might have different optimal values. In general, it is recommended to train a model using the default parameters. The output will contain a plot that enables you to review the selected `--nfeat`. Two parts should be seen in the proportion of variance across the principal components (indexed in order of importance):
        - A 'random' phase (representing PCs that explain variance caused by, inter alia, fetal fraction).
        - A 'non-random' phase (representing PCs that explain variance caused by natural Gaussian noise).
    - An optimal `--nfeat` captures the 'random' phase (as shown in the example at `./examples/overall_performance.png`). Capturing too much of the 'non-random' phase could lead to convergence problems during modeling.
    - If you are not satisfied with the performance of your model or with the position of `--nfeat`, re-run with a different number of features.
- Note that the final model will probably be a bit more accurate than what is claimed by the performance statistics. This is because PREFACE uses a cross-validation strategy where 10% of the (male) samples are excluded from training, after which these 10% serve as validation cases. This process is repeated 10 times. Therefore, the final performance measurements are based on models trained with only 90% of the (male) fetuses, yet the resulting `model.onnx` is trained with all provided cases.
