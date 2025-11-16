# Survival NPDR

Survival extension of Nearest-neighbor Projected Distance Regression (NPDR) for ranking features by association with time-to-event outcomes.

## Key Functions
`npdr_surv_binomial()` (univariate) and `npdr_surv_binomial_glm()` (penalized) score features using neighbor-based difference regression and survival time + status.

## Script Summary (`tests/`)
- `main/main_effects.r`: Benchmark on simulated main-effect data.
- `main/interaction_effects.r`: Benchmark with interaction-heavy data.
- `test_k/test_k_main.R` & `test_k/test_k_interaction.R`: Effect of neighborhood size k.
- `sensitivity_analysis/*`: Sweeps over effect size, interaction strength, censoring.
- `real_data_analysis/ovarian_main.r`: Ovarian cancer datasets; KM plots and feature ranking.
- `real_data_analysis/real_model_performance.r`: Aggregate real-data metrics.
- `regain_network/*_network.r`: Network construction using ranked features.


## Simulation
Simulated datasets for main effects and interaction effects are included under `data/`:

- `data/main_effects1.csv`
- `data/interaction_effects1.csv`
- `data/simulated_surv.csv`

You can run the benchmark scripts in `tests/main/` directly against these files.

## Dependencies
R packages: survival, glmnet, PRROC, dplyr, ggplot2, survcomp, survivalsvm, ranger, Biobase.

## Usage
```r
root <- getwd()            # project root
library(devtools)
devtools::load_all(file.path(root, "survival_NPDR/sNPDR"))
# run a model
df <- read.csv(file.path(root, "survival_NPDR/data/simulated_surv.csv"))

res <- npdr_surv_binomial(outcome = c(time_var="time", status_var="status"),
                          dataset = df,
                          nbd.method="multisurf", knn=20)
head(res)
```

## Quick example with included simulated data

```r
library(readr)
library(devtools)
devtools::load_all("sNPDR")

# Load a tiny simulated dataset (30 samples, 6 features)
sim_path <- file.path("data", "simulated_surv.csv")
dat <- readr::read_csv(sim_path, show_col_types = FALSE)

# Univariate per-feature GLMs
res_glm <- npdr_surv_binomial(
    outcome = c(time_var = "time", status_var = "status"),
    dataset = dat,
    nbd.method = "relieff",
    knn = 10,
    KM.weight = FALSE
)
head(res_glm)

# Elastic-net (joint) model; returns Feature and beta
res_net <- npdr_surv_binomial_glm(
    outcome = c(time_var = "time", status_var = "status"),
    dataset = dat,
    nbd.method = "relieff",
    knn = 10,
    use.glmnet = TRUE,
    glmnet.lam = "lambda.min"
)
head(res_net)
```

## Contact


