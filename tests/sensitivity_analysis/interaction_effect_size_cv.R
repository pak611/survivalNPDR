# Load necessary libraries and your custom functions
devtools::load_all("C:/Users/patri/OneDrive/Desktop/snpdr_update2/snpdr_update/simSurvData/cox_epistasis")
devtools::load_all("C:/Users/patri/OneDrive/Desktop/snpdr_update2/snpdr_update/sNPDR")
library(dplyr)
library(PRROC)
library(survival)
library(survivalsvm)

# Initialize an empty data frame to store AUC values for each algorithm and iteration
num_iterations <- 20
inner_iterations <- 3
effect_sizes <- seq(1.0, 5.0, length.out = num_iterations) # Effect sizes from 0.2 to 1.0
auc_results <- data.frame(
    iteration = 1:num_iterations,
    effect_size = effect_sizes,
    npdr = numeric(num_iterations),
    standard_cox = numeric(num_iterations),
    survival_svm = numeric(num_iterations)
)

# Define explicit seeds for reproducibility
seeds <- c(1234, 5678, 91011, 121314, 151617, 181920, 212223, 242526, 272829, 303132, 333435, 363738, 394041, 424344, 454647, 484950, 515253, 545556, 575859)

# Main loop for simulations
for (iter in 1:num_iterations) {

    # Initialize variables to store AUC values for inner loop
    npdr_auc <- numeric(inner_iterations)
    standard_cox_auc <- numeric(inner_iterations)
    survival_svm_auc <- numeric(inner_iterations)

    for (j in 1:inner_iterations) {

        # Print the iteration number
        print(paste0("Iteration: ", iter, ", Inner Iteration: ", j))

        # Set parameters for the simulation
        beta_int <- effect_sizes[iter]

        # Simulate data
        simdata <- simul.int(seeds[iter] + j, n = 200, p = 500,
                             n.main = 2,
                             n.int = 50,
                             beta.main = 0.1, 
                             beta.int = beta_int, 
                             censparam = 1/4, 
                             lambda = 1/10)
        dat <- simdata$data

        # ------------------ survNPDR ------------------
        survNPDR.model <- sNPDR::npdr_surv_binomial(outcome = c("time_var" = "time", "status_var" = "status"),
                          dataset = dat, 
                          attr.diff.type = "standard",
                          nbd.method = "relieff",
                          nbd.metric = "manhattan",
                          knn = 5, 
                          msurf.sd.frac = 0.5, 
                          glmnet.alpha = 1, 
                          model.type = "binomial", 
                          KM.weight = FALSE,
                          KM.kernel.type = "gaussian",
                          KM.kernel.sigma = 1.5)
        # Arrange results by beta coefficient
        survNPDR.model <- survNPDR.model |> arrange(desc(beta))

        # Plot Precision-Recall curve for survNPDR
        functional_vars <- NULL
        interaction_vars <- grep("intervar", colnames(dat), value = TRUE)
        idx_func <- which(survNPDR.model$Feature %in% interaction_vars)
        func_betas <- survNPDR.model$beta[idx_func]
        non_func_betas <- survNPDR.model$beta[-idx_func]
        pr_curve_survNPDR <- PRROC::pr.curve(scores.class0 = func_betas, 
                            scores.class1 = non_func_betas, 
                            curve = TRUE)
        npdr_auc[j] <- pr_curve_survNPDR$auc.integral

        # ------------------ Standard Cox Regression ------------------
        attr_mat <- dat[, -c(1, 2)]
        cox.model <- lapply(colnames(attr_mat), function(x) {
                formula <- as.formula(paste0("Surv(time, status) ~ ", x))
                mod_tmp <- survival::coxph(formula, data = dat)
                mod_summary <- summary(mod_tmp)
                as.data.frame(mod_summary$coefficients)
        })
        cox.model <- do.call(rbind, cox.model) |> 
        tibble::rownames_to_column(var = "Feature") |> 
        rename(beta = "coef", HR = "exp(coef)", std.err = "se(coef)", 
            beta.z = "z", p.value = "Pr(>|z|)") |> 
        mutate(p.adj = p.adjust(p.value, method = "bonferroni", n = ncol(attr_mat))) |> 
        arrange(p.value)

        # Plot Precision-Recall curve for Cox regression
        idx_func <- which(cox.model$Feature %in% interaction_vars)
        func_betas <- cox.model$beta[idx_func]
        non_func_betas <- cox.model$beta[-idx_func]

        pr_curve_standard_cox <- PRROC::pr.curve(scores.class0 = abs(func_betas),
                                        scores.class1 = abs(non_func_betas), 
                                        curve = TRUE)
        standard_cox_auc[j] <- pr_curve_standard_cox$auc.integral

        # ------------------ Survival SVM ------------------
        svm.fit <- survivalsvm(formula = Surv(time, status) ~ ., data = dat, type = "regression", gamma.mu = 0.2, opt.meth = "quadprog", kernel = "add_kernel")
        support_vectors <- svm.fit$model.fit$SV  
        coefficients <- svm.fit$model.fit$Beta
        weights <- t(support_vectors) %*% coefficients
        svm.positives <- data.frame(Feature = svm.fit$var.names, beta = weights) |> arrange(desc(abs(beta)))

        # Calculate AUC for SVM
        idx_func <- which(svm.positives$Feature %in% interaction_vars)
        func_betas <- svm.positives$beta[idx_func]
        non_func_betas <- svm.positives$beta[-idx_func]

        pr_curve_svm <- PRROC::pr.curve(scores.class0 = abs(func_betas), 
                                        scores.class1 = abs(non_func_betas), 
                                        curve = TRUE)
        survival_svm_auc[j] <- pr_curve_svm$auc.integral
    }

    # Calculate average AUC for each algorithm
    auc_results$npdr[iter] <- mean(npdr_auc)
    auc_results$standard_cox[iter] <- mean(standard_cox_auc)
    auc_results$survival_svm[iter] <- mean(survival_svm_auc)
}

# Print the data frame with AUC results for each algorithm and iteration
print(auc_results)

# Save the auc_results
saveRDS(auc_results, file = "auc_results__interactions_effect_size.rds")