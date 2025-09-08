# Load necessary libraries
library(survival)
library(dplyr)
library(tibble)
library(glmnet)
library(caret)
library(survcomp)
library(survivalsvm)
library(PRROC)
library(R.utils)
library(ggplot2)
library(gridExtra)

# Load custom simulation function
devtools::load_all("C:/Users/patri/OneDrive/Desktop/snpdr_update2/snpdr_update/simSurvData/cox_epistasis")
devtools::load_all("C:/Users/patri/OneDrive/Desktop/snpdr_update2/snpdr_update/sNPDR")

# Simulation parameters
sim_seed <- 2467
n <- 200
p <- 500
n_main <- 50
n_int <- 2
beta_main <- 0.7
beta_int <- 0.1
censparam <- 1/4
lambda <- 1/100

# Simulate data
simdata <- simul.int(sim_seed, n = n, p = p,
                     n.main = n_main,
                     n.int = n_int,
                     beta.main = beta_main, 
                     beta.int = beta_int, 
                     censparam = censparam, 
                     lambda = lambda)

dat <- simdata$data

# Create 10-fold cross-validation splits
set.seed(123)
k_folds <- createFolds(dat$status, k = 10, list = TRUE, returnTrain = TRUE)

# Initialize error trackers
errors <- data.frame(
    method = character(),
    fold = integer(),
    c_index = numeric(),
    auc = numeric(),
    runtime = numeric(),
    n = integer(),
    p = integer(),
    n_main = integer(),
    n_int = integer(),
    beta_main = numeric(),
    beta_int = numeric(),
    censparam = numeric(),
    lambda = numeric(),
    stringsAsFactors = FALSE
)

# Cross-validation loop
for (fold_idx in seq_along(k_folds)) {
  cat("Processing Fold", fold_idx, "\n")
  #browser()

  # Split data into training and testing
  train_idx <- k_folds[[fold_idx]]
  train_data <- dat[train_idx, ]
  test_data <- dat[-train_idx, ]

  # Normalize training data
  attr_train <- train_data %>% select(-time, -status) %>% scale()
  attr_test <- test_data %>% select(-time, -status) %>% scale(center = attr_train %>% attr("scaled:center"), scale = attr_train %>% attr("scaled:scale"))

  # ----------------- survNPDR (KM.weight = FALSE) -----------------
  survNPDR_runtime <- system.time({
    survNPDR.model <- withTimeout({
      sNPDR::npdr_surv_binomial(outcome = c("time_var" = "time", "status_var" = "status"),
                                dataset = train_data, 
                                attr.diff.type = "standard",
                                #nbd.method = "relieff", 
                                nbd.method = "multisurf",
                                nbd.metric = "manhattan",
                                knn = 20, 
                                msurf.sd.frac = 0.5, 
                                glmnet.alpha = 1, 
                                model.type = "binomial", 
                                KM.weight = FALSE,
                                KM.kernel.type = "gaussian",
                                KM.kernel.sigma = 0.25)
    }, timeout = 20, onTimeout = "silent")
  })[3]  

  if (!is.null(survNPDR.model)) {
    threshold <- quantile(survNPDR.model$p.adj, 0.90)
    survNPDR.model <- survNPDR.model %>% filter(p.adj < threshold) %>% arrange(desc(abs(beta)))

    top.features <- survNPDR.model$Feature[1:50]
    top.coefficients <- survNPDR.model$beta[1:50]

    X_test <- test_data %>% select(one_of(top.features)) %>% as.matrix()
    relative.risk <- exp(X_test %*% as.numeric(top.coefficients))

    # Calculate C-index
    c.ind.survNPDR <- concordance.index(x = relative.risk, surv.time = test_data$time, surv.event = test_data$status)

    # Generate AUC
    functional.vars <- grep("simvar", colnames(dat), value = TRUE)
    idx_func <- which(survNPDR.model$Feature %in% functional.vars)
    func_betas <- survNPDR.model$beta[idx_func]
    neg_betas <- survNPDR.model$beta[-idx_func]
    pr_curve_survNPDR <- PRROC::pr.curve(scores.class0 = abs(func_betas), scores.class1 = abs(neg_betas), curve = TRUE)
    auc_survNPDR <- pr_curve_survNPDR$auc.integral

    # Save results
    errors <- rbind(errors, data.frame(method = "sNPDR", fold = fold_idx, c_index = c.ind.survNPDR$c.index,
                                       auc = auc_survNPDR, runtime = survNPDR_runtime,
                                       n = n, p = p, n_main = n_main, n_int = n_int, beta_main = beta_main, beta_int = beta_int, censparam = censparam, lambda = lambda))
  }

  npdr_surv_binomial <- function(outcome = c("time_var" = "time", 
                                           "status_var" = "status"),
                                dataset, 
                                attr.diff.type = "numeric-abs",
                                nbd.method     = "relieff", 
                                nbd.metric     = "manhattan",
                                knn            = 0, 
                                msurf.sd.frac  = 0.5, 
                                glmnet.alpha   = 1, 
                                glmnet.lower   = -Inf, 
                                lambda.grid    = seq(0.01, 1, by = 0.01),
                                model.type     = "multinomial", 
                                KM.weight      = FALSE, 
                                KM.kernel.type = "gaussian", 
                                KM.kernel.sigma = 1.0,
                                feature.subsample.frac = 1.0) {  # Fraction of features to subsample

    # Record the start time
    start_time <- Sys.time()
    
    # Extract time and status variable names
    time_var <- as.character(outcome["time_var"])
    status_var <- as.character(outcome["status_var"])

    # Ensure dataset is a data frame
    if (!is.data.frame(dataset)) {
        dataset <- as.data.frame(dataset)
    }

    # Extract the attribute matrix by removing outcome columns
    attr_mat <- dataset |> select(-all_of(as.character(outcome)))
    num_samp <- nrow(dataset)

    # Check and adjust knn value if it's too large
    if (knn > num_samp - 1) {
        warning("The number of nearest neighbors `knn` is too large! Setting `knn = nrow(dataset) - 1`.")
        knn <- num_samp - 1
    }

    # Find nearest neighbors using sNPDR package
    knn_edge_list <- sNPDR::nearestNeighbors(attr.mat   = attr_mat, 
                                             sd.frac    = msurf.sd.frac, 
                                             k          = knn,
                                             nbd.method = nbd.method,
                                             nbd.metric = nbd.metric)

    Ri.pheno.vals <- dataset[knn_edge_list$Ri_idx, c(time_var, status_var)]
    NN.pheno.vals <- dataset[knn_edge_list$NN_idx, c(time_var, status_var)]

    Ri.KM.est <- computeKM(Ri.pheno.vals) # Compute the Kaplan-Meier Estimate for Ri
    NN.KM.est <- computeKM(NN.pheno.vals) # Compute the Kaplan-Meier Estimate for NN

    # Initialize the npdr design matrix
    num_attr  <- ncol(attr_mat)
    num_pairs <- nrow(knn_edge_list)
    npdr_diff_mat <- matrix(0, nrow = num_pairs, ncol = num_attr)

    # Compute attribute differences and apply KM weights if specified
    for (i in seq(1, num_attr)) {
        attr_vec <- attr_mat[, i]
        attr_Ris <- attr_vec[knn_edge_list$Ri_idx]
        attr_NNs <- attr_vec[knn_edge_list$NN_idx]
        attr_diff_vec <- sNPDR::npdrDiff(attr_Ris, attr_NNs, diff.type = attr.diff.type)
        
        if (KM.weight) {
            if (KM.kernel.type == 'linear') {
                KM.kernel <- Ri.KM.est$surv_prob * NN.KM.est$surv_prob
            } else if (KM.kernel.type == 'gaussian') {
                KM.kernel <- exp((Ri.KM.est$surv_prob - NN.KM.est$surv_prob)^2 / (2 * KM.kernel.sigma^2))
            }
            attr_diff_vec <- attr_diff_vec * KM.kernel
        }
        npdr_diff_mat[, i] <- attr_diff_vec
    }

    colnames(npdr_diff_mat) <- colnames(attr_mat)

    # Add time and status information to the knn_edge_list
    knn_edge_list <- knn_edge_list |> 
        mutate(Ri_time = dataset[Ri_idx, time_var],
               NN_time = dataset[NN_idx, time_var],
               Ri_status = dataset[Ri_idx, status_var],
               NN_status = dataset[NN_idx, status_var])

    # Filter out pairs where NN_time < Ri_time and NN_status == 0
    filter_indices <- which(knn_edge_list$NN_time <= knn_edge_list$Ri_time & knn_edge_list$NN_status == 0)
    npdr_diff_mat_full <- knn_edge_list[-filter_indices, ]
    npdr_diff_mat <- npdr_diff_mat[-filter_indices, ]

    # Create a class column for the multinomial model
    npdr_diff_mat_full <- npdr_diff_mat_full |> 
        mutate(class = case_when(
            NN_time > Ri_time ~ 1,
            NN_time == Ri_time ~ 1,
            NN_time < Ri_time ~ 0
        ))

    # Merge the npdr_diff_mat_full with the npdr_diff_mat
    npdr_diff_mat_full <- cbind(npdr_diff_mat_full, npdr_diff_mat) |> 
        select(colnames(attr_mat), class)

    # Initialize an empty list to store the results
    results_list <- list()

    # Get the feature names
    features <- colnames(npdr_diff_mat_full)[colnames(npdr_diff_mat_full) != "class"]

    # Apply feature subsampling
    num_features_to_sample <- ceiling(length(features) * feature.subsample.frac)
    sampled_features <- sample(features, num_features_to_sample, replace = FALSE)

    # Iterate over each sampled feature
    for (feature in sampled_features) {
        # Subset the data to remove rows with NA in the current feature
        subset_data <- npdr_diff_mat_full |>
            select(class, all_of(feature)) |>
            na.omit()  # Remove rows with NA in the current feature

        # Fit the binomial logistic regression model for the current feature
        binom_fit <- tryCatch(
            {
                glm(class ~ ., data = subset_data, family = binomial(link = "logit"))
            },
            error = function(e) {
                return(NULL)  # Handle cases where the model cannot be fit
            }
        )

        # Check if the model was successfully fitted
        if (!is.null(binom_fit)) {
            # Extract and format the coefficients
            coef_summary <- as.data.frame(summary(binom_fit)$coefficients)
            coef_summary <- coef_summary[2, ]  # Extract the coefficient for the feature (excluding intercept)
            coef_summary$Feature <- feature

            # Store the result
            results_list[[feature]] <- coef_summary
        }
    }

    # Combine all results into a single data frame
    results_df <- do.call(rbind, results_list)

    # Rename 'Estimate' to 'beta'
    colnames(results_df)[colnames(results_df) == "Estimate"] <- "beta"

    # Add adjusted p-values
    results_df$p.adj <- p.adjust(results_df$`Pr(>|z|)`, method = "bonferroni", n = ncol(attr_mat))

    # Sort the results by the absolute value of beta (importance)
    results_df <- results_df[order(abs(results_df$beta), decreasing = TRUE), ]

    # Return the results
    return(results_df)
}

  # ----------------- survNPDR with KM weighting (KM.weight = TRUE) -----------------
    # survNPDR + GLM
  survNPDR_glm_runtime <- system.time({
    survNPDR.model.glm <- withTimeout({
      sNPDR::npdr_surv_binomial_glm(
        outcome = c("time_var" = "time", "status_var" = "status"),
        dataset = train_data, 
        attr.diff.type = "standard",
        #nbd.method = "relieff", 
        nbd.method = "multisurf",
        nbd.metric = "manhattan",
        knn = 20, 
        msurf.sd.frac = 0.5, 
        glmnet.alpha = 1.0, 
        glmnet.lower = -Inf,
        glmnet.lam = 0.005,
        use.glmnet = TRUE,
        model.type = "binomial", 
        KM.weight = FALSE,
        KM.kernel.type = "gaussian",
        KM.kernel.sigma = 1.0)
    }, timeout = 20, onTimeout = "silent")
  })[3]  # Extract elapsed time

  if (is.null(survNPDR.model.glm)) {
    cat("survNPDR + GLM model timed out\n")
    next
  }

  #browser()

  # Randomize the order of entries with zero beta coefficient
  zero_beta_indices <- which(survNPDR.model.glm$beta == 0)
  non_zero_beta_indices <- which(survNPDR.model.glm$beta != 0)
  randomized_zero_beta_indices <- sample(zero_beta_indices)
  randomized_indices <- c(non_zero_beta_indices, randomized_zero_beta_indices)
  survNPDR.model.glm <- survNPDR.model.glm[randomized_indices, ]

  # Calculate C-index for survNPDR.model.glm
  top.features.glm <- survNPDR.model.glm$Feature[1:50]
  top.coefficients.glm <- survNPDR.model.glm$beta[1:50]

  X_train_glm <- train_data %>% select(one_of(top.features.glm)) %>% as.matrix()
  X_test_glm <- test_data %>% select(one_of(top.features.glm)) %>% as.matrix()

  linear.predictor.glm <- X_test_glm %*% as.numeric(top.coefficients.glm)
  relative.risk.glm <- exp(linear.predictor.glm)

  c.ind.survNPDR.glm <- concordance.index(x = relative.risk.glm, surv.time = test_data$time, surv.event = test_data$status)

  # Generate Precision-Recall curve for survNPDR + GLM
  idx_func <- which(survNPDR.model.glm$Feature %in% functional.vars)
  func_betas <- survNPDR.model.glm$beta[idx_func]
  neg_betas <- survNPDR.model.glm$beta[-idx_func]
  pr_curve_glm <- PRROC::pr.curve(scores.class0 = abs(func_betas), scores.class1 = abs(neg_betas), curve = TRUE)
  glm_auc <- pr_curve_glm$auc.integral

  # Add survNPDR + GLM results to errors
  errors <- rbind(errors, data.frame(method = "sNPDRL", fold = fold_idx, c_index = c.ind.survNPDR.glm$c.index,
                                     auc = glm_auc, runtime = survNPDR_glm_runtime,
                                     n = n, p = p, n_main = n_main, n_int = n_int, beta_main = beta_main, beta_int = beta_int, censparam = censparam, lambda = lambda))

    # Cox Regression
  cox_runtime <- system.time({
    cox.model <- withTimeout({
      lapply(colnames(attr_train), function(x) {
        formula <- as.formula(paste0("Surv(time, status) ~ ", x))
        mod_tmp <- survival::coxph(formula, data = train_data)
        mod_summary <- summary(mod_tmp)
        as.data.frame(mod_summary$coefficients)
      })
    }, timeout = 20, onTimeout = "silent")
  })[3]  # Extract elapsed time

  if (is.null(cox.model)) {
    cat("Cox model timed out\n")
    next
  }

  cox.model <- do.call(rbind, cox.model) |> 
    tibble::rownames_to_column(var = "Feature") |> 
    rename(beta = "coef", HR = "exp(coef)", std.err = "se(coef)", 
           beta.z = "z", p.value = "Pr(>|z|)") |> 
    mutate(p.adj = p.adjust(p.value, method = "bonferroni", n = ncol(attr_train))) |> 
    arrange(p.value)

  top.features.cox <- cox.model$Feature[1:50]
  top.coefficients.cox <- cox.model$beta[1:50]

  X_train_cox <- train_data %>% select(one_of(top.features.cox)) %>% as.matrix()
  X_test_cox <- test_data %>% select(one_of(top.features.cox)) %>% as.matrix()

  linear.predictor.cox <- X_test_cox %*% as.numeric(top.coefficients.cox)
  relative.risk.cox <- exp(linear.predictor.cox)

  # Calculate C-index for Cox Regression
  c.ind.cox <- concordance.index(x = relative.risk.cox, surv.time = test_data$time, surv.event = test_data$status)

  # Generate PRC curve for Cox AUC
  idx_func <- which(cox.model$Feature %in% functional.vars)
  func_betas <- cox.model$beta[idx_func]
  neg_betas <- cox.model$beta[-idx_func]
  pr_curve_cox <- PRROC::pr.curve(scores.class0 = abs(func_betas), scores.class1 = abs(neg_betas), curve = TRUE)
  auc_cox <- pr_curve_cox$auc.integral

  # Add Cox results to errors
  errors <- rbind(errors, data.frame(method = "Cox", fold = fold_idx, c_index = c.ind.cox$c.index,
                                     auc = auc_cox, runtime = cox_runtime,
                                     n = n, p = p, n_main = n_main, n_int = n_int, beta_main = beta_main, beta_int = beta_int, censparam = censparam, lambda = lambda))

  # Survival SVM
  svm_runtime <- system.time({
    svm.fit <- withTimeout({
      survivalsvm(formula = Surv(train_data$time, train_data$status) ~ ., 
                  data = as.data.frame(attr_train), 
                  type = "regression", 
                  gamma.mu = 0.2, 
                  opt.meth = "quadprog", 
                  kernel = "add_kernel")
    }, timeout = 20, onTimeout = "silent")
  })[3]  # Extract elapsed time

  if (is.null(svm.fit)) {
    cat("SVM model timed out\n")
    next
  }

  #browser()

  support_vectors <- svm.fit$model.fit$SV  
  coefficients <- svm.fit$model.fit$Beta
  weights <- t(support_vectors) %*% coefficients

  # need to reassign the feature names
  names(weights) <- colnames(attr_train)

  # Construct a dataframe with feature names and their corresponding weights
  svm.model <- data.frame(Feature = names(weights), Weight = weights) %>% 
    arrange(desc(abs(Weight))) %>% 
    mutate(p.adj = p.adjust(Weight, method = "bonferroni", n = ncol(attr_train))) %>% 
    arrange(p.adj)


  top.features.svm <- svm.model$Feature[1:50]
  top.coefficients.svm <- svm.model$Weight[1:50]

  X_train_svm <- train_data %>% select(one_of(top.features.svm)) %>% as.matrix()
  X_test_svm <- test_data %>% select(one_of(top.features.svm)) %>% as.matrix()

  linear.predictor.svm <- X_test_svm %*% as.numeric(top.coefficients.svm)
  relative.risk.svm <- exp(linear.predictor.svm) %>% as.vector()



  # Calculate C-index for SVM
  c.ind.svm <- concordance.index(x = relative.risk.svm, surv.time = test_data$time, surv.event = test_data$status)

  # Generate PRC curve for SVM AUC
  #browser()
  idx_func <- which(svm.model$Feature %in% functional.vars)
  func_betas <- svm.model$Weight[idx_func]
  neg_betas <- svm.model$Weight[-idx_func]
  pr_curve_svm <- PRROC::pr.curve(scores.class0 = abs(func_betas), scores.class1 = abs(neg_betas), curve = TRUE)
  auc_svm <- pr_curve_svm$auc.integral

  # Add SVM results to errors
  errors <- rbind(errors, data.frame(method = "SVM", fold = fold_idx, c_index = c.ind.svm$c.index,
                                     auc = auc_svm, runtime = svm_runtime,
                                     n = n, p = p, n_main = n_main, n_int = n_int, beta_main = beta_main, beta_int = beta_int, censparam = censparam, lambda = lambda))
}

# Summarize cross-validation results
results_summary <- errors %>% 
  group_by(method, n, p, n_main, n_int, beta_main, beta_int, censparam, lambda) %>% 
  summarise(mean_c_index = mean(c_index, na.rm = TRUE), 
            median_c_index = median(c_index, na.rm = TRUE),
            sd_c_index = sd(c_index, na.rm = TRUE),
            mean_auc = mean(auc, na.rm = TRUE),
            median_auc = median(auc, na.rm = TRUE),
            sd_auc = sd(auc, na.rm = TRUE),
            mean_runtime = mean(runtime, na.rm = TRUE),
            median_runtime = median(runtime, na.rm = TRUE),
            sd_runtime = sd(runtime, na.rm = TRUE))

print(results_summary)

# Save the results summary to a CSV file
write.csv(results_summary, "paper_tables/results_summary_main.csv", row.names = FALSE)

# ----------------- PLOTTING C-INDEX AND AUC -----------------
# Load necessary library for color palette
library(RColorBrewer)

# Define a unified y-axis range
y_axis_range <- range(c(errors$c_index, errors$auc), na.rm = TRUE)

# Create the C-index boxplot
c_index_plot <- ggplot(errors, aes(x = method, y = c_index)) +
  geom_boxplot(aes(fill = method), outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 3, alpha = 0.6) +
  stat_summary(fun = mean, geom = "text", aes(label = round(..y.., 2)), 
               vjust = -3.5, size = 10, color = "black") +
  labs(title = NULL, x = NULL, y = "C-index") +
  scale_fill_brewer(palette = "Set2") +
  coord_cartesian(ylim = y_axis_range) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 30, face = "bold"),
    axis.text = element_text(size = 28),
    axis.text.x = element_text(size = 36, face = "bold"),
    axis.text.y = element_text(size = 36),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, size = 1.5)
  )

# Create the AUC boxplot
auc_plot <- ggplot(errors, aes(x = method, y = auc)) +
  geom_boxplot(aes(fill = method), outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 3, alpha = 0.6) +
  stat_summary(fun = mean, geom = "text", aes(label = round(..y.., 2)), 
               vjust = -2.5, size = 10, color = "black") +
  labs(title = NULL, x = NULL, y = "AUC") +
  scale_fill_brewer(palette = "Set2") +
  coord_cartesian(ylim = y_axis_range) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 36, face = "bold"),
    axis.text = element_text(size = 34),
    axis.text.x = element_text(size = 36, face = "bold"),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, size = 1.5)
  )

# Arrange the plots side by side
combined_plot <- grid.arrange(c_index_plot, auc_plot, ncol = 2)

# Save the combined plot
ggsave("paper_graphics/c_index_auc_boxplot_updated.png", plot = combined_plot, width = 20, height = 20)
