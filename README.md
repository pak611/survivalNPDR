# sNPDR

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
Simulation utilities used in our experiments come from the third‑party `coxed` package (functions like `sim.survdata`). We do not include that code here. Please refer to the `coxed` package and its vignettes for data generation, or install it locally and point scripts to your simulated datasets.

## Dependencies
R packages: survival, glmnet, PRROC, dplyr, ggplot2, survcomp, survivalsvm, ranger, Biobase.

## Usage
```r
root <- getwd()            # project root
library(devtools)
devtools::load_all(file.path(root, "sNPDR"))
# run a model
res <- npdr_surv_binomial(outcome = c(time_var="time", status_var="status"),
                          dataset = your_data_frame,
                          nbd.method="multisurf", knn=20)
head(res)
```

## Contact
