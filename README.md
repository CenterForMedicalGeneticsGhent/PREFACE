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

PREFACE operates directly on aligned sequencing read files for baseline `FFY` calculation, model training, and unified inference:
- Each alignment file is read in a single pass across all configured bin sizes (`--ffy-binsize`, default `5000` bp on `chrY`/autosomes; `--binsize`, default `100000` bp on autosomes/`chrX`).
- Primary, non-duplicate, non-QC-failed alignments with mapping quality at or above `--mapq` (default: `1`) are binned into fixed genomic windows.
- `.bam` files are read via `pysam` when installed or Python's built-in BGZF/BAM parser. `.cram` files require `--reference` (FASTA) and either `pysam` or `samtools` in `PATH`.
- Chromosome names may either omit or include the `chr` prefix (`1` until `22`, `X`, and `Y`, or `chr1` until `chr22`, `chrX`, and `chrY`). Chromosomes `13`, `18`, and `21` are excluded during diploid autosomal read-depth normalization and feature extraction.

### Samplesheet (`.yaml`, `.yml`, or `.json`)

`preface train` consumes a unified YAML or JSON samplesheet containing a root-level array of sample objects:

- Example: `./examples/samplesheet.yaml`
- Every samplesheet is automatically validated against the bundled JSON Schema (`src/preface/samplesheet.schema.json`).
- Each sample object supports the following properties (`additionalProperties: false`):
    - `id` (required `string`): Unique identifier for the sample.
    - `filepath` (required `string`): Path to the sample's aligned read file (`.bam`, `.cram`, or `.sam`). Must exist on disk and be readable.
    - `index` (optional `string` or `null`): Path to the alignment index file (`.gzi`, `.bai`, `.csi`, or `.crai`). Must exist on disk and be readable when provided.
    - `sex` (optional `"M"`, `"F"`, or `null`): Fetal sex (`M` for male, `F` for female). Automatically inferred by the `FFY` model when omitted or `null`.
    - `ff` (optional `number` in `[0.0, 1.0]` or `null`): Fetal fraction expressed as a proportion between `0.0` and `1.0`. Automatically estimated by the `FFY` model when omitted or `null`.

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

## Model training (`preface train [all|ffy|preface|ffx]`)

```bash
pixi run preface train [all|ffy|preface|ffx] -s path/to/samplesheet.yaml --prefix path/to/prefix [optional arguments]
```

`preface train` supports three training targets (defaulting to `all` when omitted):

1. `preface train` / `preface train all`:
    - Reads all cohort alignment files in a **single pass** (or loads/saves pre-binned counts via `--cache-counts`).
    - Trains the baseline `FFY` model (two-component Gaussian mixture model sex classification + 10-fold cross-validated non-negative Lasso across informative `chrY` bins per Section 3.4 of the [PREFACE paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC6771918/)) and writes the populated samplesheet config to `<prefix>.samplesheet.yaml`.
    - Trains the `PREFACE` (PCA + neural network or linear model + OLS skew correction) and `FFX` (Huber robust linear model) models using the `FFY` outputs.
    - Integrates the `FFY`, `FFX`, and `PREFACE` models into a **single unified ONNX file** at `<prefix>.onnx`, alongside `<prefix>.training_statistics.txt`, `<prefix>.overall_performance.png`, `<prefix>.FFX.png`, and `<prefix>.training_repeats/repeat_1..10.png`.
2. `preface train ffy`:
    - Trains only the `FFY` model and writes `<prefix>.ffy_statistics.txt`, `<prefix>.samplesheet.yaml`, and `<prefix>.onnx`.
3. `preface train preface` (alias: `preface train ffx`):
    - Trains the `PREFACE` + `FFX` model from a labeled samplesheet (such as `<prefix>.samplesheet.yaml` produced by `preface train ffy`), writing `<prefix>.training_statistics.txt`, plots, and `<prefix>.onnx`.

If `--prefix` ends with `/` or points to an existing directory, outputs are placed inside that directory (`<prefix>/model.onnx`, `<prefix>/samplesheet.yaml`, `<prefix>/training_statistics.txt`, etc.).

| Argument / Option | Function |
| :--- | :--- |
| `[all\|ffy\|preface\|ffx]` | Optional positional target (`all` by default; `ffx` is an alias for `preface`). |
| `-s`, `--samplesheet x` | Training samplesheet file (`.yaml`, `.yml`, or `.json`). |
| `-p`, `--prefix x` | Output path prefix for trained `<prefix>.onnx`, samplesheet, statistics, and plots. |
| `--nfeat x` | Number of principal components to use during modeling. (default: `x=50`) |
| `--hidden x` | Number of hidden units used in neural network. Use with caution. (default: `x=2`) |
| `--cpus x` | Number of CPU threads for parallel file loading and linear algebra. (default: all available cores) |
| `--femprop` | Include female fetuses in supervised learning phase (requires proportional `ff` labels for females). |
| `--olm` | Use an ordinary linear model instead of a neural network. |
| `--noskewcorrect` | Disable linear skew correction on cross-validated predictions. |
| `--binsize x` | Genomic bin width in bp for autosomal and chrX feature extraction. (default: `x=100000`) |
| `--ffy-binsize x` | Genomic bin width in bp on chrY when training FFY. (default: `x=5000`) |
| `--mapq x` | Minimum mapping quality for aligned reads. (default: `x=1`) |
| `--yfrac x` | Explicit normalized chrY depth threshold for `M`/`F` sex classification (overrides 2-component GMM). |
| `--lambda x` | Lasso regularization penalty $\lambda$ for chrY bin weighting. (default: `x=1e-4`) |
| `--no-female-gate` | Disable gating FFY predictions to `0.0` for female samples. |
| `--reference x` | Reference genome FASTA file for `.cram` decoding. |
| `--male-ref x` | One or more 100% male reference alignment files (liquid biopsies) when calibrating FFY. |
| `--cache-counts x` | Optional `.npz` path to save or reuse single-pass binned read counts across runs. |
| `--device x` | Compute accelerator to use (`auto`, `cpu`, `cuda`, `mps`). (default: `auto`) |
| `--seed x` | Random seed for reproducible shuffling, cross-validation, and initialization. (default: `x=1`) |
| `--loglevel x` | Logging verbosity level (`DEBUG`, `INFO`, `WARNING`, `ERROR`, `CRITICAL`). (default: `x=INFO`) |

## Predicting (`preface predict`)

```bash
pixi run preface predict --infile path/to/sample.bam --model path/to/model.onnx [--prefix path/to/output]
```

`preface predict` reads the sample alignment file in a single pass and executes the unified `FFY` + `FFX` + `PREFACE` ONNX graph in a single inference run:
- When `--prefix` is **omitted**, predictions are written in plain text to `stdout` with fetal fractions formatted as percentages (`sex = M`, `FFY = 12.34%`, `FFX = 11.89%`, `PREFACE = 12.05%`).
- When `--prefix` is **provided**, predictions are written as JSON to `<prefix>.json` with fetal fractions formatted as float proportions between `0.0` and `1.0` (`{"sex": "M", "FFY": 0.1234, "FFX": 0.1189, "PREFACE": 0.1205}`).

| Optional argument | Function |
| :--- | :--- |
| `--infile x` | Sample `.bam`, `.cram`, or `.sam` alignment file. |
| `--model x` | Trained `.onnx` model file. |
| `-p`, `--prefix x` | Output path prefix for JSON predictions (`0..1` floats). Omit to print plain-text percentages to `stdout`. |
| `--mapq x` | Minimum mapping quality for aligned reads. (default: `x=1`) |
| `--reference x` | Reference genome FASTA file for `.cram` decoding. |
| `--no-female-gate` / `--female-gate` | Override model setting to disable or enable gating `FFY` to `0.0` for female samples. |
| `--device x` | Execution device for ONNX Runtime (`auto`, `cpu`, `cuda`, `mps`). (default: `auto`) |
| `--loglevel x` | Logging verbosity level (`DEBUG`, `INFO`, `WARNING`, `ERROR`, `CRITICAL`). (default: `x=INFO`) |

## Model optimization

- The most important parameter is `--nfeat`:
    - It represents the number of principal components (PCs) that will be used as features during model training. Depending on the bin size and the number of training samples, it might have different optimal values. In general, it is recommended to train a model using the default parameters. The output will contain a plot that enables you to review the selected `--nfeat`. Two parts should be seen in the proportion of variance across the principal components (indexed in order of importance):
        - A 'random' phase (representing PCs that explain variance caused by, inter alia, fetal fraction).
        - A 'non-random' phase (representing PCs that explain variance caused by natural Gaussian noise).
    - An optimal `--nfeat` captures the 'random' phase (as shown in the example at `./examples/overall_performance.png`). Capturing too much of the 'non-random' phase could lead to convergence problems during modeling.
    - If you are not satisfied with the performance of your model or with the position of `--nfeat`, re-run with a different number of features (using `--cache-counts` to skip re-reading alignment files).
- Note that the final model will probably be a bit more accurate than what is claimed by the performance statistics. This is because PREFACE uses a cross-validation strategy where 10% of the (male) samples are excluded from training, after which these 10% serve as validation cases. This process is repeated 10 times. Therefore, the final performance measurements are based on models trained with only 90% of the (male) fetuses, yet the resulting `.onnx` model is trained with all provided cases.
