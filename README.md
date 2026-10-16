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

### PREFACE's config.txt

For training, PREFACE requires a config file.

- Example: `./examples/config.txt`
- Tab-separated file with at least four columns.
- The name of these columns (passed as a header) must be `ID`, `filepath`, `gender` and `FF`.
    - `ID` is used to specify a mandatory unique identifier to each of the samples.
    - The `filepath` column holds the full absolute path of the training copy number alteration files.
    - The possible values for `gender` are either `M` (male) or `F` (female), representing fetal gender. Twins/triplets/... can be included if they are all male or all female.
    - The `FF` column contains the response variable (the 'true' fetal fraction). One can use any method he/she believes performs best at quantifying the actual fetal fraction. PREFACE was benchmarked using the number of mapped Y-reads, referred to as FFY. As FFY is not informative for female fetuses, this measure is ignored for cases labeled with `F`, unless the `--femprop` flag is given (see below).

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

Training writes `model.onnx` into `--outdir` as a self-contained ONNX graph that embeds NA mean-imputation, PCA centering and projection, the neural network or linear model weights, skew correction, and FFX regression parameters. It also writes `overall_performance.png`, per-fold plots under `training_repeats/repeat_1.png` through `repeat_10.png`, `FFX.png`, and `training_statistics.txt`.

## Predicting

```bash
pixi run preface predict --infile path/to/infile.bed --model path/to/model.onnx [optional arguments]
```

| Optional argument | Function |
| :--- | :--- |
| `--json [x]` | Predictions are written to stdout. Use `--json` for JSON format on stdout, or provide `--json path/to/output.json` to write to a JSON file. |
| `--device x` | Execution device for ONNX Runtime (`auto`, `cpu`, `cuda`, `mps`). (default: `auto`) |

## Model optimization

- The most important parameter is `--nfeat`:
    - It represents the number of principal components (PCs) that will be used as features during model training. Depending on the used copy number alteration software, bin size and the number of training samples, it might have different optimal values. In general, it is recommended to train a model using the default parameters. The output will contain a plot that enables you to review the selected `--nfeat`. Two parts should be seen in the proportion of variance across the principal components (indexed in order of importance):
        - A 'random' phase (representing PCs that explain variance caused by, inter alia, fetal fraction).
        - A 'non-random' phase (representing PCs that explain variance caused by natural Gaussian noise).
    - An optimal `--nfeat` captures the 'random' phase (as shown in the example at `./examples/overall_performance.png`). Capturing too much of the 'non-random' phase could lead to convergence problems during modeling.
    - If you are not satisfied with the performance of your model or with the position of `--nfeat`, re-run with a different number of features.
- Note that the final model will probably be a bit more accurate than what is claimed by the performance statistics. This is because PREFACE uses a cross-validation strategy where 10% of the (male) samples are excluded from training, after which these 10% serve as validation cases. This process is repeated 10 times. Therefore, the final performance measurements are based on models trained with only 90% of the (male) fetuses, yet the resulting `model.onnx` is trained with all provided cases.
