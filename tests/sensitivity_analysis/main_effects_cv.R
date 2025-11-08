# Load necessary libraries and your custom functions
root <- getwd()
devtools::load_all(file.path(root, "survival_NPDR/cox_epistasis"))
devtools::load_all(file.path(root, "survival_NPDR/sNPDR"))
library(dplyr)
library(PRROC)
library(survival)
library(survivalsvm)

# Initialize an empty data frame to store AUC values for each algorithm and iteration
num_iterations <- 20
auc_results <- data.frame(
    iteration = 1:num_iterations,
    npdr = numeric(num_iterations),
    npdr_glmnet = numeric(num_iterations),
    standard_cox = numeric(num_iterations),
    survival_svm = numeric(num_iterations)
)

seeds <- c(1234, 5678, 2468, 1357, 9876, 4321, 8765, 2468, 1357, 9876, 1234, 5678, 2468, 1357, 9876, 4321, 8765, 2468, 1357, 9876)

# Run the analysis 20 times for each algorithm
for (i in 1:num_iterations) {

    #browser()

    # Parameters for the simulation
    num_features <- 500
    num_samples <- 200
    num_functional <- 50
    censor <- 1/5

    functional.vars <- paste0("simvar", 1:num_functional)

    # Simulate data
    simdata <- simul.int(seeds[i], n = num_samples, p = num_features,
                         n.main = num_functional,
                         n.int = 2,
                         beta.main = 0.70, 
                         beta.int = 0.01, 
                         censparam = censor, 
                         lambda = 1/5)

    dat <- simdata$data

    # ------------------ survNPDR ------------------
    survNPDR.model <- sNPDR::npdr_surv_binomial(outcome = c("time_var" = "time", "status_var" = "status"),
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
                                                KM.kernel.sigma = 1.5)
    # Arrange results by beta coefficient
    survNPDR.model <- survNPDR.model |> arrange(desc(beta))

    # Calculate AUC for survNPDR
    func_betas <- survNPDR.model$beta[survNPDR.model$Feature %in% functional.vars]
    neg_betas <- survNPDR.model$beta[!survNPDR.model$Feature %in% functional.vars]
    pr_curve_npdr <- PRROC::pr.curve(scores.class0 = abs(func_betas), scores.class1 = abs(neg_betas), curve = TRUE)
    auc_results$npdr[i] <- pr_curve_npdr$auc.integral

    # ------------------ survNPDR with glmnet ------------------
    survNPDR_glmnet.model <- sNPDR::npdr_surv_binomial_glm(outcome = c("time_var" = "time", "status_var" = "status"),
                                                           dataset = dat, 
                                                           attr.diff.type = "standard",
                                                           nbd.method = "relieff", 
                                                           nbd.metric = "manhattan",
                                                           knn = num_samples - 1, 
                                                           msurf.sd.frac = 0.5, 
                                                           glmnet.alpha = 0.1, 
                                                           glmnet.lower = -Inf,
                                                           glmnet.lam = 0.01,
                                                           use.glmnet = TRUE,
                                                           model.type = "binomial", 
                                                           KM.weight = FALSE,
                                                           KM.kernel.type = "gaussian",
                                                           KM.kernel.sigma = 1.0)
    # Randomize the order of entries with zero beta coefficient
    zero_beta_indices <- which(survNPDR_glmnet.model$beta == 0)
    non_zero_beta_indices <- which(survNPDR_glmnet.model$beta != 0)
    randomized_zero_beta_indices <- sample(zero_beta_indices)
    randomized_indices <- c(non_zero_beta_indices, randomized_zero_beta_indices)
    survNPDR_glmnet.model <- survNPDR_glmnet.model[randomized_indices, ]

    # Calculate AUC for survNPDR with glmnet
    func_betas_glmnet <- survNPDR_glmnet.model$beta[survNPDR_glmnet.model$Feature %in% functional.vars]
    neg_betas_glmnet <- survNPDR_glmnet.model$beta[!survNPDR_glmnet.model$Feature %in% functional.vars]
    pr_curve_npdr_glmnet <- PRROC::pr.curve(scores.class0 = abs(func_betas_glmnet), scores.class1 = abs(neg_betas_glmnet), curve = TRUE)
    auc_results$npdr_glmnet[i] <- pr_curve_npdr_glmnet$auc.integral

    # ------------------ Standard Cox Regression ------------------
    attr_mat <- select(dat, -c(time, status))
    cox.model <- lapply(colnames(attr_mat), function(x) {
        formula <- as.formula(paste0("Surv(time, status) ~ ", x))
        mod_tmp <- survival::coxph(formula, data = dat)
        as.data.frame(summary(mod_tmp)$coefficients)
    })
    cox.model <- do.call(rbind, cox.model) |> 
        tibble::rownames_to_column(var = "Feature") |> 
        rename(beta = "coef", HR = "exp(coef)", std.err = "se(coef)", 
               beta.z = "z", p.value = "Pr(>|z|)") |> 
        mutate(p.adj = p.adjust(p.value, method = "bonferroni", n = ncol(attr_mat))) |> 
        arrange(p.value)

    # Calculate AUC for Standard Cox
    func_betas <- cox.model$beta[cox.model$Feature %in% functional.vars]
    neg_betas <- cox.model$beta[!cox.model$Feature %in% functional.vars]
    pr_curve_standard_cox <- PRROC::pr.curve(scores.class0 = abs(func_betas), scores.class1 = abs(neg_betas), curve = TRUE)
    auc_results$standard_cox[i] <- pr_curve_standard_cox$auc.integral

    # ------------------ Survival SVM ------------------
    svm.fit <- survivalsvm(formula = Surv(time, status) ~ ., data = dat, type = "regression", gamma.mu = 0.2, opt.meth = "quadprog", kernel = "add_kernel")
    support_vectors <- svm.fit$model.fit$SV  
    coefficients <- svm.fit$model.fit$Beta
    weights <- t(support_vectors) %*% coefficients
    svm.positives <- data.frame(Feature = svm.fit$var.names, beta = weights) |> arrange(desc(abs(beta)))

    # Calculate AUC for SVM
    func_betas <- svm.positives$beta[svm.positives$Feature %in% functional.vars]
    neg_betas <- svm.positives$beta[!svm.positives$Feature %in% functional.vars]
    pr_curve_svm <- PRROC::pr.curve(scores.class0 = abs(func_betas), scores.class1 = abs(neg_betas), curve = TRUE)
    auc_results$survival_svm[i] <- pr_curve_svm$auc.integral
}

# Print the data frame with AUC results for each algorithm and iteration
print(auc_results)

# Save the auc_results
saveRDS(auc_results, file = "auc_results.rds")