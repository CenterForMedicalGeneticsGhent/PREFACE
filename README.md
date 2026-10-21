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

The original R script (`PREFACE.R`) remains in the repository for reference.

# Manual

## Installation

Install the dependencies using `pixi`:

```bash
pixi install
```

and run the script using `pixi run`:

```bash
pixi run PREFACE.R <train/predict> ...
```

## Required files

### Copy number alteration .bed files

Each sample (whether it is used for training or for predicting) should be passed to PREFACE in the format shown below. During benchmarking, using a bin size of 100 kb (others might work equally well), copy number normalization was performed by [WisecondorX](https://github.com/CenterForMedicalGeneticsGhent/WisecondorX/), yet PREFACE is not limited to any copy number alteration software, however, the default output of WisecondorX is directly interpretable by PREFACE.

- Example: `./examples/infile.bed`
- Tab-separated file with at least four columns.
- The name of these columns (passed as a header) must be `chr`, `start`, `end` and `ratio`.
    - The possible values of `chr` are `1` until `22`, and `X` and `Y` (uppercase).
    - The `ratio` column contains the log2-transformed ratio between the observed and expected copy number.
    - The ratio can be unknown at certain loci (e.g. often seen at centromeres). Here, values should be expressed as `NaN` or `NA`.
- The order of rows does not matter. Yet, it is paramount that, for a certain line, file x deals with the same locus as file y. This implies, of course, that all copy number alteration files have the same number of lines.

### Binned read count or alignment files (for FFY calculation)

To calculate baseline chromosome Y fetal fraction (`FFY`) and infer fetal sex (`gender`) with `preface ffy`, each sample can be provided either as a binned read count file or as an aligned sequencing file:

1. Binned read count `.bed` / `.tsv` files
    - Example: `./examples/counts.bed`
    - Tab-separated file with at least four columns representing chromosome, start coordinate, end coordinate, and read count.
    - Both headered files and headerless 4-column BED files (such as outputs from `samtools bedcov`, `mosdepth`, or `bedtools multicov`) are accepted.
        - When a header is present, `preface ffy` matches column names `chr` / `chrom` / `#chrom`, `start` / `chromStart`, `end` / `chromEnd`, and `count` / `counts` / `reads` / `read_count` / `depth` / `coverage` / `ratio` (falling back to the first four columns by position).
        - Chromosome names may either omit or include the `chr` prefix (`1` until `22` and `Y`, or `chr1` until `chr22` and `chrY`).
        - The count column must contain non-negative raw or depth-normalized read counts (not log2 copy number ratios, which contain negative values and typically mask chromosome `Y`).
        - Each file must contain both autosomal bins (chromosomes `1` until `22`, where `13`, `18`, and `21` are automatically excluded during diploid read-depth normalization) and chromosome `Y` bins. A bin size of 5 kb on chromosome `Y` is recommended to match Section 3.4 of the PREFACE paper.
2. Aligned read files (`.bam`, `.cram`, `.sam`)
    - `preface ffy` bins primary, non-duplicate, non-QC-failed alignments with mapping quality at or above `--mapq` (default: `1`) into fixed genomic windows of size `--binsize` (default: `5000` bp) across autosomes and chromosome `Y`.
    - `.bam` files are read directly using Python's built-in BGZF/BAM parser (or `pysam` when installed). `.cram` files require `--reference` (FASTA) and either `pysam` or `samtools` in `PATH`.

### PREFACE's config.txt

For training, PREFACE requires a config file.

- Example: `./examples/config.txt`
- Tab-separated file with at least four columns.
- The name of these columns (passed as a header) must be `ID`, `filepath`, `gender` and `FF`.
    - `ID` is used to specify a mandatory unique identifier to each of the samples.
    - The `filepath` column holds the full absolute path of the training copy number alteration files.
    - The possible values for `gender` are either `M` (male) or `F` (female), representing fetal gender. Twins/triplets/... can be included if they are all male or all female.
    - The `FF` column contains the response variable (the 'true' fetal fraction). One can use any method he/she believes performs best at quantifying the actual fetal fraction. PREFACE was benchmarked using the number of mapped Y-reads, referred to as FFY. As FFY is not informative for female fetuses, this measure is ignored for cases labeled with `F`, unless the `--femprop` flag is given (see below).

## Calculating baseline FFY (`preface ffy`)

`preface ffy` implements the two-stage baseline `FFY` calculation described in Section 3.4 of the [PREFACE paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC6771918/):

1. Each chromosome `Y` bin count is normalized by the sample's mappable diploid autosomal read density (excluding chromosomes `13`, `18`, and `21`). Pseudoautosomal regions (`PAR1` $\le 2.78\text{ Mb}$, `PAR2`) and `Yq12` heterochromatin ($\ge 28.8\text{ Mb}$) are masked when computing the summary chromosome `Y` depth.
2. In cohort mode, fetal sex (`M` or `F`) is taken from the input table when present or inferred via a two-component Gaussian mixture model (or `--yfrac` threshold). Prior `FFY` is scaled between the female cohort median (`0%`) and the 100% male baseline (either empirical adult male liquid biopsies passed via `--male-ref` or theoretical diploid autosomal scaling $2 \times Y_{\text{norm}}$). A 10-fold cross-validated non-negative Lasso model ($\lambda = 10^{-4}$) is then fit across the chromosome `Y` bins to weight informative male-specific loci and produce the final `FF` column (`0.0` for females, clamped `[0, 100]` `%` for males).

### Cohort mode (generate `config.txt` for `preface train`)

```bash
pixi run preface ffy --config path/to/manifest.txt --outfile path/to/config.txt [optional arguments]
```

The input `--config` manifest is a tab-separated file with required columns `ID` and `filepath`, plus optional columns `reads_filepath` and `gender`:
- If `reads_filepath` is present, `preface ffy` reads the `.bam`/`.cram` or binned read count `.bed` file (such as `./examples/counts.bed`) from `reads_filepath` while keeping `filepath` (the copy number alteration `.bed` file for `preface train`, such as `./examples/infile.bed`) unchanged in the output `config.txt`.
- If `reads_filepath` is omitted, `preface ffy` reads the counts or alignments directly from `filepath`.
- If `gender` (`M` or `F`) is omitted or `NA`, `preface ffy` fills it in automatically.

### Single-sample mode

```bash
pixi run preface ffy --infile ./examples/counts.bed [--model path/to/ffy_model.onnx] [--json]
```

In single-sample mode, `preface ffy` prints `gender` and `FFY` to stdout (or writes JSON when `--json` is passed). When `--model` is provided, it applies the cohort-calibrated ONNX Lasso weights and sex threshold saved via `--save-model`; otherwise, it computes `FFY` directly from the PAR-masked and heterochromatin-masked diploid autosomal ratio.

| Optional argument | Function |
| :--- | :--- |
| `--config x` | Cohort TSV file with `ID`, `filepath`, and optional `reads_filepath` and `gender` columns. |
| `--infile x` | Single-sample `.bam`/`.cram` or 4-column binned read count `.bed`/`.tsv` file (e.g. `./examples/counts.bed`). |
| `--outfile x` | Destination path for the generated `config.txt` in cohort mode. (default: stdout) |
| `--male-ref x` | One or more 100% male reference samples (liquid biopsies) for empirical denominator scaling. Defaults to theoretical diploid autosomal scaling ($2 \times Y_{\text{norm}}$) when omitted. |
| `--save-model x` | Export the calibrated cohort Lasso FFY model and baselines to an ONNX file. |
| `--model x` | Load a calibrated FFY ONNX model in single-sample mode. |
| `--binsize x` | Genomic bin size in bp when binning aligned `.bam`/`.cram` files. (default: `x=5000`) |
| `--mapq x` | Minimum mapping quality for aligned reads. (default: `x=1`) |
| `--yfrac x` | Explicit normalized chrY depth threshold for `M`/`F` classification (overrides 2-component GMM). |
| `--lambda x` | Lasso regularization penalty $\lambda$ for 5 kb chrY bin weighting. (default: `x=1e-4`) |
| `--reference x` | Reference genome FASTA file for `.cram` decoding. |
| `--json [x]` | In single-sample mode, write predictions as JSON to stdout or to file `x`. |
| `--cpus x` | Number of CPU threads for parallel sample loading. (default: all available cores) |
| `--seed x` | Random seed for reproducible cross-validation. (default: `x=1`) |
| `--loglevel x` | Logging verbosity level (`DEBUG`, `INFO`, `WARNING`, `ERROR`, `CRITICAL`). (default: `x=INFO`) |

## Model training

```bash
pixi run preface train --config path/to/config.txt --outdir path/to/dir/ [optional arguments]
```

| Optional argument | Function |
| :--- | :--- |
| `--nfeat x` | Number of principal components to use during modeling. (default: `x=50`) |
| `--hidden x` | Number of hidden units used in neural network. Use with caution. (default: `x=2`) |
| `--cpus x` | Number of CPU threads for parallel file loading and linear algebra. (default: all available cores) |
| `--femprop` | When using FFY as FF (recommended), FF labels for female fetuses are irrelevant, and should be ignored in the supervised learning phase (default). If this behavior is not desired, use this flag, which demands that the given FFs for female fetuses are proportional to their actual FF. |
| `--olm` | It might be possible the neural network does not converge; or for your kind of data/sample size, an ordinary linear model might be a better option. In these cases, use this flag. |
| `--noskewcorrect` | This flag ascertains the best fit for most (instead of all) of the data is generated. Mostly not recommended. |
| `--device x` | Compute accelerator to use (`auto`, `cpu`, `cuda`, `mps`). Selects CUDA or Apple Silicon MPS when available. (default: `auto`) |
| `--seed x` | Random seed for reproducible shuffling, cross-validation, and initialization. (default: `x=1`) |
| `--loglevel x` | Logging verbosity level (`DEBUG`, `INFO`, `WARNING`, `ERROR`, `CRITICAL`). (default: `x=INFO`) |

Training writes `model.onnx` into `--outdir` as a self-contained ONNX graph that embeds NA mean-imputation, PCA centering and projection, the neural network or linear model weights, skew correction, and FFX regression parameters. It also writes `overall_performance.png`, per-fold plots under `training_repeats/repeat_1.png` through `repeat_10.png`, `FFX.png`, and `training_statistics.txt`.

## Predicting

```bash
pixi run preface predict --infile path/to/infile.bed --model path/to/model.onnx [optional arguments]
```

| Optional argument | Function |
| :--- | :--- |
| `--json [x]` | Predictions are written to stdout. Use `--json` for JSON format on stdout, or provide `--json path/to/output.json` to write to a JSON file. |
| `--device x` | Execution device for ONNX Runtime (`auto`, `cpu`, `cuda`, `mps`). (default: `auto`) |
| `--loglevel x` | Logging verbosity level (`DEBUG`, `INFO`, `WARNING`, `ERROR`, `CRITICAL`). (default: `x=INFO`) |

## Model optimization

- The most important parameter is `--nfeat`:
    - It represents the number of principal components (PCs) that will be used as features during model training. Depending on the used copy number alteration software, bin size and the number of training samples, it might have different optimal values. In general, it is recommended to train a model using the default parameters. The output will contain a plot that enables you to review the selected `--nfeat`. Two parts should be seen in the proportion of variance across the principal components (indexed in order of importance):
        - A 'random' phase (representing PCs that explain variance caused by, inter alia, fetal fraction).
        - A 'non-random' phase (representing PCs that explain variance caused by natural Gaussian noise).
    - An optimal `--nfeat` captures the 'random' phase (as shown in the example at `./examples/overall_performance.png`). Capturing too much of the 'non-random' phase could lead to convergence problems during modeling.
    - If you are not satisfied with the performance of your model or with the position of `--nfeat`, re-run with a different number of features.
- Note that the final model will probably be a bit more accurate than what is claimed by the performance statistics. This is because PREFACE uses a cross-validation strategy where 10% of the (male) samples are excluded from training, after which these 10% serve as validation cases. This process is repeated 10 times. Therefore, the final performance measurements are based on models trained with only 90% of the (male) fetuses, yet the resulting `model.onnx` is trained with all provided cases.
