# Survival NPDR

**Nearest-Neighbor Projected Distance Regression for Biomarker Discovery in Gene-Expression Survival Data**

> Patrick Kampmeyer, Bryan A. Dawkins, Ryan J. Urbanowicz, Kia Kazemi-Nia, Brett A. McKinney

Survival NPDR extends [NPDR](https://github.com/insilico/npdr) to right-censored survival outcomes. It identifies prognostic biomarkers in high-dimensional gene-expression data by regressing pairwise survival-time differences on pairwise attribute differences within a censor-aware nearest-neighbor neighborhood.

---

## Repository Structure

```
survivalNPDR/
├── sNPDR/                        # Core R package — install this first
│   ├── R/
│   │   ├── npdr_surv.R           # npdr_surv_lm(), npdr_surv_lm_km(), npdr_surv_binomial()
│   │   ├── npdr_surv_symmetric.R # Symmetric pair variant
│   │   ├── nearestNeighbors.R    # ReliefF / MultiSURF neighborhood construction
│   │   ├── adaptive.R            # Adaptive k selection (MultiSURF)
│   │   ├── KM_estimate.R         # Kaplan-Meier conditional mean imputation
│   │   ├── regain.R              # ReGAIN pairwise interaction network
│   │   └── utils.R
│   └── DESCRIPTION
│
├── scripts/
│   ├── simulation/               # Figs 3 & 6, Tables 1 & 2
│   │   ├── main_effect_size_sweep_rsurv.r        # Main-effect benchmark + effect-size sweep
│   │   └── interaction_effect_size_sweep_rsurv.r # Interaction benchmark + effect-size sweep
│   │
│   ├── test_k/                   # Fig 7 — neighborhood size study
│   │   ├── test_k_lognormal.R   # runs both main + interaction k sweeps
│   │   └── plot_k_combined_lognormal.R
│   │
│   ├── real_data/                # Tables 3 & 4, Fig 8
│   │   ├── preprocess_breast_cancer.R  # Download + preprocess GEO breast datasets
│   │   ├── ovarian_risk_stratification.R                   # Ovarian cancer sNPDR + KM risk stratification
│   │   ├── breast_model_performance.R         # Performance metrics (C-index, logrank, HR) for GSE2034
│   │   └── breast_survival_curve.R         # KM descriptive figure for GSE9893
│   │
│   ├── network/                  # Figs 9 & 10 — interaction network analysis
│   │   ├── breast_cancer_network.r
│   │   └── ovarian_cancer_network.r
│   │
│   ├── figures/                  # Assemble publication figures
│   │   ├── plot_lognormal_paper_figures.R
│   │   ├── plot_ef_vs_effect_size.R
│   │   └── plot_mrr_vs_effect_size.R
│   │
│   └── analysis_statistical.R   # Friedman + Wilcoxon significance tests
│
├── data/                         # Processed real datasets
│   └── GSE{9893,2034,2990,9891,32062,13876,17260}/
│
├── results/                      # Pre-computed paper result tables (CSV)
└── figures/                      # Pre-computed paper figures (PNG)
```

---

## Installation

### 1. Install R dependencies

```r
install.packages(c(
  "survival", "dplyr", "rlang", "PRROC", "glmnet",
  "caret", "survcomp", "survivalsvm", "ranger",
  "R.utils", "ggplot2", "tibble", "rsurv", "parallel"
))

# Bioconductor (for real-data preprocessing)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("Biobase", "GEOquery", "curatedOvarianData"))
```

### 2. Install the sNPDR package

```r
# From local source (recommended)
devtools::install("sNPDR")

# Or from GitHub
devtools::install_github("pak611/survivalNPDR", subdir = "sNPDR")
```

---

## Quick Start

Your data should be a `data.frame` with a `time` column, a `status` column (1 = event, 0 = censored), and all other columns as features.

```r
library(sNPDR)

# ── Main-effect detection — signed (asymmetric), adaptive neighborhood ─────────
res <- npdr_surv(
  outcome   = c(time_var = "time", status_var = "status"),
  dataset   = my_data,
  diff.type = "signed"       # signed ΔT ~ signed ΔX
)
head(res[order(res$pval.att), ])
```

```r
# ── Interaction detection — absolute (symmetric), fixed k = 20 ────────────────
res <- npdr_surv(
  outcome    = c(time_var = "time", status_var = "status"),
  dataset    = my_data,
  diff.type  = "absolute",   # |ΔT| ~ |ΔX|
  nbd.method = "relieff",
  knn        = 20
)
```

```r
# ── With KM imputation and covariate adjustment ────────────────────────────────
res <- npdr_surv(
  outcome    = c(time_var = "time", status_var = "status"),
  dataset    = my_data,
  diff.type  = "signed",
  km.impute  = TRUE,
  covariates = c("age", "sex", "BMI")
)
```

```r
# ── Logistic regression with LASSO regularization ─────────────────────────────
res <- npdr_surv(
  outcome    = c(time_var = "time", status_var = "status"),
  dataset    = my_data,
  regression = "binomial",
  regularize = TRUE,
  alpha      = 1             # 1 = LASSO, 0 = ridge
)
```

### Visualize top features

```r
library(ggplot2)

top <- res |> dplyr::arrange(pval.att) |> head(20)

ggplot(top, aes(x = reorder(att, -abs(beta.raw.att)), y = beta.raw.att)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(x = "Feature", y = "sNPDR coefficient",
       title = "Top 20 survival-associated features")
```

---

## Reproducing Paper Results

All scripts use `root <- getwd()`. Run from the repository root:

```bash
cd /path/to/survivalNPDR
```

### Figs 3 & 6 / Tables 1 & 2 — Simulation benchmarks

Each sweep script runs all comparison methods (Cox, LASSO-Cox, survival SVM, RSF, sNPDR-signed, sNPDR-absolute) across a range of effect sizes in parallel.

```bash
# Main-effect sweep (default beta_main = 0.00, 0.07, 0.13, 0.20, 0.27, 0.33, 0.40)
Rscript scripts/simulation/main_effect_size_sweep_rsurv.r

# Run at a single effect size only (e.g. beta_main = 0.5 for Table 1)
BETA_MAIN_SWEEP="0.5" Rscript scripts/simulation/main_effect_size_sweep_rsurv.r

# Interaction-effect sweep (default beta_int = 0.00, 0.03, 0.07, 0.10, 0.13, 0.17, 0.20)
Rscript scripts/simulation/interaction_effect_size_sweep_rsurv.r

# Run at a single effect size only (e.g. beta_int = 1.0 for Table 2)
BETA_INT_SWEEP="1.0" Rscript scripts/simulation/interaction_effect_size_sweep_rsurv.r

# Assemble figures
Rscript scripts/figures/plot_lognormal_paper_figures.R
Rscript scripts/figures/plot_ef_vs_effect_size.R
```

### Fig 7 — Choice of neighborhood size k

```bash
# Both sweeps run in parallel (uses all available cores − 2)
Rscript scripts/test_k/test_k_lognormal.R   # runs both main + interaction k sweeps
Rscript scripts/test_k/plot_k_combined_lognormal.R
# → figures/auPRC_K_combined_lognormal_no_multisurf.png
```

### Tables 3 & 4 — Real cancer dataset benchmarks

```bash
# Breast cancer: download and preprocess from GEO (requires GEOquery)
Rscript scripts/real_data/preprocess_breast_cancer.R

# Ovarian cancer sNPDR + KM risk stratification
Rscript scripts/real_data/ovarian_risk_stratification.R

# Performance metrics (C-index, logrank, HR) for GSE2034 → Tables 3 & 4
Rscript scripts/real_data/breast_model_performance.R
```

### Fig 8 — Kaplan-Meier risk stratification

```bash
Rscript scripts/real_data/breast_survival_curve.R
```

### Figs 9 & 10 — Interaction networks (ReGAIN)

```bash
Rscript scripts/network/breast_cancer_network.r
Rscript scripts/network/ovarian_cancer_network.r
```

### Statistical significance tests

```bash
Rscript scripts/analysis_statistical.R
```

---

## Data

Processed expression matrices and clinical annotations are in `data/`:

| Dataset | Cancer type | Samples | Notes |
|---------|-------------|---------|-------|
| GSE9893 | Breast | 155 | ERα+ |
| GSE2034 | Breast | 286 | Lymph-node negative |
| GSE2990 | Breast | 189 | Mixed |
| GSE9891 | Ovarian | 260 | Stage III/IV |
| GSE32062 | Ovarian | 260 | High-grade serous |
| GSE13876 | Ovarian | 157 | High-grade serous |
| GSE17260 | Ovarian | 110 | High-grade serous |

Raw data can be re-downloaded from [GEO](https://www.ncbi.nlm.nih.gov/geo/) using `scripts/real_data/preprocess_breast_cancer.R` (breast) or `curatedOvarianData` (ovarian).

---

## Pre-computed Results

`results/` contains pre-computed CSV tables from the paper. `figures/` contains the corresponding plots. Re-running any script above will overwrite them.

---

## Key sNPDR Options

All functionality is accessed through a single function:

```r
npdr_surv(outcome, dataset, diff.type, regression, ...)
```

| Option | Parameter | Values |
|--------|-----------|--------|
| **Difference type** | `diff.type` | `"absolute"` — symmetric \|ΔT\| ~ \|ΔX\| (best for interactions); `"signed"` — asymmetric ΔT ~ ΔX (best for main effects) |
| **Regression model** | `regression` | `"lm"` (default) — linear regression; `"binomial"` — logistic regression |
| **KM imputation** | `km.impute` | `TRUE` — replace censored times with KM conditional mean before ΔT; `FALSE` (default) |
| **Covariate adjustment** | `covariates` | e.g. `c("age", "sex", "BMI")` — pairwise differences included as regressors |
| **Regularization** | `regularize`, `alpha` | `regularize = TRUE` enables elastic-net; `alpha = 1` → LASSO, `alpha = 0` → ridge |
| **Neighborhood** | `nbd.method`, `knn` | `"multisurf"` — adaptive k (default); `"relieff"` + `knn = 20` — fixed k |

---

## Citation

```bibtex
@article{kampmeyer2025survnpdr,
  title   = {Survival {NPDR}: Nearest-Neighbor Projected Distance Regression
             for Biomarker Discovery in Gene-Expression Survival Data},
  author  = {Kampmeyer, Patrick and Dawkins, Bryan A. and Urbanowicz, Ryan J.
             and Kazemi-Nia, Kia and McKinney, Brett A.},
  year    = {2026}
}
```

---

## Related

- [NPDR R package](https://github.com/insilico/npdr) — original NPDR for non-survival outcomes
- [rsurv](https://CRAN.R-project.org/package=rsurv) — lognormal AFT simulation used in benchmarks
