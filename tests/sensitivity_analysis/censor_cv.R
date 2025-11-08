# Load necessary libraries
root <- getwd()
devtools::load_all(file.path(root, "survival_NPDR/cox_epistasis"))
devtools::load_all(file.path(root, "survival_NPDR/sNPDR"))
library(dplyr)
library(PRROC)
library(survival)
library(survivalsvm)

# Initialize variables
num_iterations <- 30
inner_iterations <- 3
censor_rates <- seq(0.05, 0.5, length.out = num_iterations) # Censoring rates from 0% to 50%
auc_results <- data.frame(
    iteration = 1:num_iterations,
    censor_rate = censor_rates,
    npdr = numeric(num_iterations),
    standard_cox = numeric(num_iterations),
    survival_svm = numeric(num_iterations)
)

# Set seeds for reproducibility
seeds <- c(1234, 5678, 2468, 1357, 9876, 4321, 8765, 13579, 24680, 10101, 11111, 22222, 33333, 44444, 55555, 66666, 77777, 88888, 99999, 12121, 23232, 34343, 45454, 56565, 67676, 78787, 89898, 90909, 101010, 111111)

# Main loop for simulations
for (i in 1:num_iterations) {

  # Initialize variables to store AUC values for inner loop
  npdr_auc <- numeric(inner_iterations)
  standard_cox_auc <- numeric(inner_iterations)
  survival_svm_auc <- numeric(inner_iterations)

  for (j in 1:inner_iterations) {

    # Set simulation parameters
    num_features <- 500
    num_samples <- 200
    num_functional <- 50
    censor <- censor_rates[i]
    functional_vars <- paste0("simvar", 1:num_functional)

    # Simulate data
    simdata <- simul.int(seeds[i] + j, n = num_samples, p = num_features,
                         n.main = num_functional, n.int = 2,
                         beta.main = 0.70, beta.int = 0.01,
                         censparam = censor, lambda = 1/5)
    dat <- simdata$data

    # ------------------ survNPDR ------------------
    survNPDR.model <- sNPDR::npdr_surv_binomial(
      outcome = c("time_var" = "time", "status_var" = "status"),
      dataset = dat,
      attr.diff.type = "standard",
      nbd.method = "multisurf",
      nbd.metric = "manhattan",
      knn = 20,
      msurf.sd.frac = 0.5,
      glmnet.alpha = 1,
      model.type = "binomial",
      KM.weight = FALSE,
      KM.kernel.type = "gaussian",
      KM.kernel.sigma = 1.5
    )
    func_betas <- survNPDR.model$beta[survNPDR.model$Feature %in% functional_vars]
    neg_betas <- survNPDR.model$beta[!survNPDR.model$Feature %in% functional_vars]
    pr_curve_npdr <- PRROC::pr.curve(scores.class0 = abs(func_betas), scores.class1 = abs(neg_betas), curve = TRUE)
    npdr_auc[j] <- pr_curve_npdr$auc.integral

    # ------------------ Standard Cox Regression ------------------
    attr_mat <- select(dat, -c(time, status))
    cox.model <- lapply(colnames(attr_mat), function(x) {
      formula <- as.formula(paste0("Surv(time, status) ~ ", x))
      mod_tmp <- survival::coxph(formula, data = dat)
      as.data.frame(summary(mod_tmp)$coefficients)
    })
    cox.model <- do.call(rbind, cox.model) |>
      tibble::rownames_to_column(var = "Feature") |>
      rename(beta = "coef")
    func_betas <- cox.model$beta[cox.model$Feature %in% functional_vars]
    neg_betas <- cox.model$beta[!cox.model$Feature %in% functional_vars]
    pr_curve_cox <- PRROC::pr.curve(scores.class0 = abs(func_betas), scores.class1 = abs(neg_betas), curve = TRUE)
    standard_cox_auc[j] <- pr_curve_cox$auc.integral

    # ------------------ Survival SVM ------------------
    svm.fit <- survivalsvm(formula = Surv(time, status) ~ ., data = dat, type = "regression", gamma.mu = 0.2, opt.meth = "quadprog", kernel = "add_kernel")
    support_vectors <- svm.fit$model.fit$SV  
    coefficients <- svm.fit$model.fit$Beta
    weights <- t(support_vectors) %*% coefficients
    svm.positives <- data.frame(Feature = svm.fit$var.names, beta = weights) |> arrange(desc(abs(beta)))

    # Calculate AUC for SVM
    func_betas <- svm.positives$beta[svm.positives$Feature %in% functional_vars]
    neg_betas <- svm.positives$beta[!svm.positives$Feature %in% functional_vars]
    pr_curve_svm <- PRROC::pr.curve(scores.class0 = abs(func_betas), scores.class1 = abs(neg_betas), curve = TRUE)
    survival_svm_auc[j] <- pr_curve_svm$auc.integral
  }

  # Calculate average AUC for each algorithm
  auc_results$npdr[i] <- mean(npdr_auc)
  auc_results$standard_cox[i] <- mean(standard_cox_auc)
  auc_results$survival_svm[i] <- mean(survival_svm_auc)
}

# Print and save results
print(auc_results)
saveRDS(auc_results, "auc_results_censor.rds")
