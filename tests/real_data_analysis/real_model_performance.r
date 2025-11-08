# Load required libraries
library(survcomp)
library(dplyr)
library(devtools)
library(survival)
library(survivalsvm)

# Load sNPDR package and data
library(rprojroot)
root <- getwd()
devtools::load_all(file.path(root, "survival_NPDR/sNPDR"))
dat <- readRDS(file.path(root, "survival_NPDR/data", "GSE2034", "dat.rds"))
dat <- as.data.frame(dat)
dat <- dat %>% filter(!is.na(time) & !is.na(status))

# Preprocess the dataset
outcome <- dat[, c("time", "status")]
attr_mat <- dat[, setdiff(names(dat), c("time", "status"))]
variance_values <- apply(attr_mat, 2, var)
low_variance_genes <- order(variance_values)[1:500]
attr_mat <- attr_mat[, low_variance_genes]

# Initialize results table
results_table <- data.frame(
  Model = character(),
  C_Index = numeric(),
  Logrank_Chi2 = numeric(),
  Hazard_Ratio = numeric(),
  Num_Weights = numeric()
)

# Function to calculate metrics
calculate_metrics <- function(dat, risk_scores, beta_coefficients, model_name) {
  # Group by median risk score
  dat$group <- ifelse(risk_scores > median(risk_scores), "High", "Low")
  
  # Logrank χ²
  logrank_test <- survdiff(Surv(time, status) ~ group, data = dat)
  logrank_chi2 <- logrank_test$chisq
  
  # Hazard ratio
  cox_model <- coxph(Surv(time, status) ~ risk_scores, data = dat)
  hazard_ratio <- exp(coef(cox_model))
  
  # Number of weights
  num_weights <- sum(abs(beta_coefficients) > 1e-8)
  
  # C-index
  c_index <- concordance.index(
    x = risk_scores,
    surv.time = dat$time,
    surv.event = dat$status
  )$c.index
  
  # Return results as a list
  return(data.frame(Model = model_name, C_Index = c_index, Logrank_Chi2 = logrank_chi2, 
                    Hazard_Ratio = hazard_ratio, Num_Weights = num_weights))
}

#------------------------------------ survNPDR ------------------------------------
survNPDR <- sNPDR::npdr_surv_binomial(outcome = c("time_var" = "time", "status_var" = "status"),
                                       dataset = dat, attr.diff.type = "standard",
                                       nbd.method = "multisurf", nbd.metric = "manhattan",
                                       knn = 20, msurf.sd.frac = 0.5, glmnet.alpha = 1,
                                       glmnet.lower = -Inf, model.type = "binomial")

# save the model
write.csv(survNPDR, "survNPDR_top_genes2_GSE2034.csv")
# read the model
survNPDR <- read.csv("survNPDR_top_genes2_GSE2034.csv")
# Filter top features
threshold <- quantile(survNPDR$p.adj, 0.90)
survNPDR <- survNPDR %>% filter(p.adj < threshold) %>% arrange(desc(abs(beta)))
top_features <- survNPDR$Feature[1:30]
top_coefficients <- survNPDR$beta[1:30]

# Calculate risk scores
X <- dat %>% select(all_of(top_features)) %>% as.matrix()
risk_scores <- X %*% as.numeric(top_coefficients)

# Get survNPDR metrics
results_table <- rbind(results_table, calculate_metrics(dat, risk_scores, top_coefficients, "survNPDR"))

# ------------------------------------ survNPDR.glm ------------------------------------
survNPDR_glm <- sNPDR::npdr_surv_binomial_glm(outcome = c("time_var" = "time", "status_var" = "status"),
                                              dataset = dat, attr.diff.type = "standard",
                                              nbd.method = "relieff", nbd.metric = "manhattan",
                                              knn = 20, msurf.sd.frac = 0.5, glmnet.alpha = 0.1,
                                              glmnet.lower = -Inf, glmnet.lam = "lambda.1se",
                                              use.glmnet = TRUE, model.type = "binomial")

# save the model
write.csv(survNPDR_glm, "survNPDR_glm_top_genes2_GSE2034.csv")
# read the model
survNPDR_glm <- read.csv("survNPDR_glm_top_genes2_GSE2034.csv")
# Filter top features
# Build predictive model using the top 30 features
top.features <- survNPDR.glm$Feature[1:30]
top.coefficients <- survNPDR.glm$beta[1:30]


# Calculate risk scores
X <- dat %>% select(one_of(top_features)) %>% as.matrix()
risk_scores <- X %*% as.numeric(top_coefficients)

# Get survNPDR.glm metrics
results_table <- rbind(results_table, calculate_metrics(dat, risk_scores, top_coefficients, "survNPDR.glm"))

# ------------------------------------ Cox Regression ------------------------------------
# Run univariate Cox regression for each gene
cox.model <- lapply(colnames(dat |> select(-c(time, status))), function(x) {
  tryCatch({
    formula <- as.formula(paste("Surv(time, status) ~ `", x, "`", sep = ""))
    coxFit <- coxph(formula, data = dat)
    summary(coxFit)
  }, error = function(e) {
    message(paste("Error in column:", x, "-", e$message))
    return(NULL)
  })
})


# Remove NULL elements and extract results
cox.model <- cox.model[!sapply(cox.model, is.null)]
results <- lapply(cox.model, function(summary_obj) {
  coef_table <- summary_obj$coefficients
  data.frame(Feature = rownames(coef_table)[1], 
             Beta = coef_table[1, "coef"], 
             P.Value = coef_table[1, "Pr(>|z|)"])
})

# Combine results into a data frame
results_df <- do.call(rbind, results)
results_df <- as.data.frame(results_df)  # Ensure results_df is a data frame

# Arrange by absolute value of Beta
results_df <- results_df %>% arrange(desc(abs(Beta)))

# # Save the model
# write.csv(results_df, "cox_regression_top_genes2_GSE2034.csv")

# # Read the model
# results_df <- read.csv("cox_regression_top_genes2_GSE2034.csv")

# Top features and coefficients
top_features <- results_df$Feature[1:30]
top_coefficients <- results_df$Beta[1:30]

# Select only the features that exist in the dataset
existing_features <- intersect(top_features, colnames(dat))
missing_features <- setdiff(top_features, existing_features)

# Note the number of genes that could not be selected
num_missing_features <- length(missing_features)
message(paste("Number of missing features:", num_missing_features))

# Calculate risk scores
X <- dat %>% select(all_of(existing_features)) %>% as.matrix()
risk_scores <- X %*% as.numeric(top_coefficients[match(existing_features, top_features)])

# Get Cox regression metrics
results_table <- rbind(results_table, calculate_metrics(dat, risk_scores, top_coefficients, "Cox Regression"))

# ------------------------------------ Survival SVM ------------------------------------
svm.model <- survivalsvm(
  formula = Surv(time, status) ~ ., 
  data = dat, 
  type = "regression", 
  gamma.mu = 0.2, 
  opt.meth = "quadprog", 
  kernel = "add_kernel"
)

# Extract and save SVM results
support_vectors <- svm.model$model.fit$SV  
coefficients <- svm.model$model.fit$Beta
weights <- t(support_vectors) %*% coefficients
svm.results <- data.frame("Feature" = svm.model$var.names, beta = weights) %>% 
               arrange(desc(abs(beta)))


# Calculate and print concordance index for SVM model
top.features <- svm.results$Feature[1:30]
top.coefficients <- svm.results$beta[1:30]

# save the model
write.csv(as.data.frame(top.features, top.coefficients), "survival_svm_top_genes2_GSE2034.csv")

# Select only the features that exist in the dataset
existing.features <- intersect(top.features, colnames(dat))
missing.features <- setdiff(top.features, existing.features)

# Note the number of genes that could not be selected
num.missing.features <- length(missing.features)
message(paste("Number of missing features:", num.missing.features))

# Print missing features for debugging
if (num.missing.features > 0) {
    message("Missing features: ", paste(missing.features, collapse = ", "))
}

# Calculate risk scores
X <- dat %>% select(all_of(existing_features)) %>% as.matrix()
risk_scores <- X %*% as.numeric(top.coefficients[match(existing.features, top.features)])

# Get Survival SVM metrics
results_table <- rbind(results_table, calculate_metrics(dat, risk_scores, top_coefficients, "Survival SVM"))


# Save results table
write.csv(results_table, "model_performance_metrics_2034.csv")
print(results_table)
