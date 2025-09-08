library(Biobase)
library(MetaGxOvarian)
library(survcomp)
library(dplyr)
library(devtools)
library(survival)
library(survivalsvm)
library(survminer)

dataset <- "GSE9891"
# # Load ovarian datasets
# esetsAndDups <- loadOvarianEsets()
# dataSets <- esetsAndDups$esets

# Select datasets
# GSE9891 <- dataSets[["GSE9891"]]
# GSE32062 <- dataSets[["GSE32062"]]
# GSE13876 <- dataSets[["GSE13876"]]
# GSE17260 <- dataSets[["GSE17260"]]

# # save the datasets
# save(GSE9891, file = "C:/Users/patri/OneDrive/Desktop/snpdr_update2/snpdr_update/data/GSE_9891/GSE9891.rda")
# save(GSE32062, file = "C:/Users/patri/OneDrive/Desktop/snpdr_update2/snpdr_update/data/GSE_32062/GSE32062.rda")
# save(GSE13876, file = "C:/Users/patri/OneDrive/Desktop/snpdr_update2/snpdr_update/data/GSE_13876/GSE13876.rda")
# save(GSE17260, file = "C:/Users/patri/OneDrive/Desktop/snpdr_update2/snpdr_update/data/GSE_17260/GSE17260.rda")

dataset <- "GSE9891"

# Load the datasets
load(paste0("C:/Users/patri/OneDrive/Desktop/snpdr_update2/snpdr_update/data/", dataset, "/", dataset, ".rda")) # both survival time and recurrence status
#load("C:/Users/patri/OneDrive/Desktop/snpdr_update2/snpdr_update/data/GSE_32062/GSE32062.rda") # just survival time, no recurrence status
#load("C:/Users/patri/OneDrive/Desktop/snpdr_update2/snpdr_update/data/GSE_13876/GSE13876.rda") # just survival time, no recurrence status
#load("C:/Users/patri/OneDrive/Desktop/snpdr_update2/snpdr_update/data/GSE_17260/GSE17260.rda") # both survival time and recurrence status

# read the data


# Choose dataset for analysis
gset <- get(dataset)
expr_data <- exprs(gset)
feature_data <- fData(gset)

# Map probe IDs to gene symbols
probe_to_gene <- feature_data[, c("probeset", "gene")]
colnames(probe_to_gene) <- c("ID", "gene")  # Rename column for consistency

# Convert expression data to a data frame and add row names as a column
expr_data_df <- as.data.frame(expr_data)
expr_data_df$ID <- rownames(expr_data_df)

# Merge expression data with gene names
expr_data_with_genes <- merge(expr_data_df, probe_to_gene, by = "ID", all.x = TRUE)

# Ensure unique gene names
expr_data_with_genes$gene <- make.unique(as.character(expr_data_with_genes$gene))

# Set row names to gene symbols
rownames(expr_data_with_genes) <- expr_data_with_genes$gene

# Remove ID and gene columns
expr_data_with_genes <- expr_data_with_genes[, !(names(expr_data_with_genes) %in% c("ID", "gene"))]

# Print a preview of the updated expression data
print(head(expr_data_with_genes))

# Extract phenotype data
pheno_data <- pData(gset)

#pheno_data <- pheno_data[, c("days_to_tumor_recurrence", "recurrence_status")]
pheno_data <- pheno_data[, c("days_to_death", "vital_status")]

# Convert recurrence_status to numeric (1 = recurrence, 0 = no recurrence)
#pheno_data$recurrence_status <- ifelse(pheno_data$recurrence_status == "recurrence", 1, 0)

# Convert vital_status to numeric (1 = deceased, 0 = alive)
pheno_data$vital_status <- ifelse(pheno_data$vital_status == "deceased", 1, 0)

# Ensure sample names match between phenotype and expression data
sample_names <- colnames(expr_data_with_genes)
pheno_data <- pheno_data[match(sample_names, rownames(pheno_data)), ]


# Combine phenotype data with expression data
combined_data <- cbind(pheno_data, t(expr_data_with_genes))



# Load sNPDR package and data
devtools::load_all("C:/Users/patri/OneDrive/Desktop/snpdr_update2/snpdr_update/sNPDR")
# dat <- readRDS("C:/Users/patri/OneDrive/Desktop/snpdr_update2/snpdr_update/data/GSE_2034/dat.rds")
#dat <- readRDS("C:/Users/patri/OneDrive/Desktop/snpdr_update2/snpdr_update/data/GSE_2990/dat.rds")
#dat <- readRDS("C:/Users/patri/OneDrive/Desktop/snpdr_update2/snpdr_update/data/GSE_9893/dat.rds")
dat <- as.data.frame(combined_data)


# Preprocess the dataset
# Preprocess the dataset
# dat$time <- dat$days_to_tumor_recurrence
# dat$status <- dat$recurrence_status
dat$time <- dat$days_to_death
dat$status <- dat$vital_status

# Remove the original columns
#dat <- dat[, setdiff(names(dat), c("days_to_tumor_recurrence", "recurrence_status"))]
dat <- dat[, setdiff(names(dat), c("days_to_death", "vital_status"))]
# remove rows with NA for the time or status
dat <- dat[complete.cases(dat[, c("time", "status")]), ]
# Ensure all columns in attr_mat are numeric
attr_mat <- dat[, setdiff(names(dat), c("time", "status"))]
attr_mat <- attr_mat[, sapply(attr_mat, is.numeric)]

# Calculate variance and select top 50th percentile of highest variance genes
variance_values <- apply(attr_mat, 2, var, na.rm = TRUE)
#top_variance_genes <- names(sort(variance_values, decreasing = TRUE)[1:1000])
top_variance_genes <- names(sort(variance_values, decreasing = TRUE)[1:20000])
#top_variance_genes <- names(sort(variance_values, decreasing = TRUE)[1:42400])
attr_mat <- attr_mat[, top_variance_genes]


# Combine the high variance genes with the time and status columns
dat <- cbind(attr_mat, time = dat$time, status = dat$status)


# outcome <- dat[, c("time", "status")]
# attr_mat <- dat[, setdiff(names(dat), c("time", "status"))]
# variance_values <- apply(attr_mat, 2, var)
# low_variance_genes <- order(variance_values)[1:500]
# attr_mat <- attr_mat[, low_variance_genes]

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

# # save the model
write.csv(survNPDR, paste0("C:/Users/patri/OneDrive/Desktop/snpdr_update2/snpdr_update/data/", dataset, "/", dataset, "NPDR_RANKING", ".csv"))
# # read the model
# survNPDR <- read.csv("survNPDR_top_genes2_GSE2034.csv")
# # Filter top features
threshold <- quantile(survNPDR$p.adj, 0.90)
survNPDR <- survNPDR %>% filter(p.adj < threshold) %>% arrange(desc(abs(beta)))
top_features <- survNPDR$Feature[1:30]
top_coefficients <- survNPDR$beta[1:30]

# Calculate risk scores
X <- dat %>% select(all_of(top_features)) %>% as.matrix()
risk_scores_npdr <- X %*% as.numeric(top_coefficients)

# Get survNPDR metrics
results_table <- rbind(results_table, calculate_metrics(dat, risk_scores_npdr, top_coefficients, "survNPDR"))

# ------------------------------------ survNPDR.glm ------------------------------------
# survNPDR.glm <- sNPDR::npdr_surv_binomial_glm(outcome = c("time_var" = "time", "status_var" = "status"),
#                                               dataset = dat, attr.diff.type = "standard",
#                                               nbd.method = "relieff", nbd.metric = "manhattan",
#                                               knn = 20, msurf.sd.frac = 0.5, glmnet.alpha = 0.1,
#                                               glmnet.lower = -Inf, glmnet.lam = "lambda.1se",
#                                               use.glmnet = TRUE, model.type = "binomial")

# # # save the model
# # write.csv(survNPDR_glm, "survNPDR_glm_top_genes2_GSE2034.csv")
# # # read the model
# # survNPDR_glm <- read.csv("survNPDR_glm_top_genes2_GSE2034.csv")
# # Filter top features
# # Build predictive model using the top 30 features
# top.features <- survNPDR.glm$Feature[1:30]
# top.coefficients <- survNPDR.glm$beta[1:30]


# # Calculate risk scores
# X <- dat %>% select(one_of(top_features)) %>% as.matrix()
# risk_scores_glm <- X %*% as.numeric(top_coefficients)

# # Get survNPDR.glm metrics
# results_table <- rbind(results_table, calculate_metrics(dat, risk_scores, top_coefficients, "survNPDR.glm"))

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

# Save the model
write.csv(results_df, paste0("C:/Users/patri/OneDrive/Desktop/snpdr_update2/snpdr_update/data/", dataset, "/", dataset, "COX_RANKING", ".csv"))
#write.csv(results_df, "C:/Users/patri/OneDrive/Desktop/snpdr_update2/snpdr_update/data/GSE_32062/cox_regression_top_genes_GSE32062.csv")

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
risk_scores_cox <- X %*% as.numeric(top_coefficients[match(existing_features, top_features)])

# Get Cox regression metrics
results_table <- rbind(results_table, calculate_metrics(dat, risk_scores_cox, top_coefficients, "Cox Regression"))

# ------------------------------------ Survival SVM ------------------------------------
svm.model <- survivalsvm(
  formula = Surv(time, status) ~ ., 
  data = dat, 
  type = "regression", 
  gamma.mu = 1, 
  opt.meth = "quadprog", 
  kernel = "add_kernel")

# Extract and save SVM results
support_vectors <- svm.model$model.fit$SV  
coefficients <- svm.model$model.fit$Beta
weights <- t(support_vectors) %*% coefficients
svm.results <- data.frame("Feature" = svm.model$var.names, beta = weights) %>% 
               arrange(desc(abs(beta)))


# Calculate and print concordance index for SVM model
top.features <- svm.results$Feature[1:30]
top.coefficients <- svm.results$beta[1:30]

# Save the model with dynamic file naming
write.csv(as.data.frame(top.features, top.coefficients), paste0("C:/Users/patri/OneDrive/Desktop/snpdr_update2/snpdr_update/data/", dataset, "/", dataset, "SVM_RANKING", ".csv"))
#write.csv(as.data.frame(top.features, top.coefficients), "C:/Users/patri/OneDrive/Desktop/snpdr_update2/snpdr_update/data/GSE_32062/survival_svm_top_genes_GSE32062.csv")

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
risk_scores_svm <- X %*% as.numeric(top.coefficients[match(existing.features, top.features)])

# Get Survival SVM metrics
results_table <- rbind(results_table, calculate_metrics(dat, risk_scores_svm, top_coefficients, "Survival SVM"))


# # Save results table
#write.csv(results_table, "C:/Users/patri/OneDrive/Desktop/snpdr_update2/snpdr_update/data/GSE_17260/model_performance_metrics_17260.csv")
#write.csv(results_table, "C:/Users/patri/OneDrive/Desktop/snpdr_update2/snpdr_update/data/GSE_32062/model_performance_metrics_32062.csv")
print(results_table)


# ------------------------------------ Kaplan-Meier Plot ------------------------------------
# Categorize patients into high and low risk groups based on risk scores

#risk_scores <- risk_scores_cox  # Choose risk scores from the desired model
risk_scores <- risk_scores_npdr
# risk_scores <- risk_scores_glm
#risk_scores <- risk_scores_svm

dat$group <- cut(risk_scores, breaks = quantile(risk_scores, probs = c(0, 0.4, 0.6, 1)), labels = c("Low", "Intermediate", "High"), include.lowest = TRUE)

# Filter out the intermediate group
dat_km <- dat %>% filter(group != "Intermediate")

# Create a Kaplan-Meier plot
km_fit <- survfit(Surv(time, status) ~ group, data = dat_km)
ggsurvplot(km_fit, data = dat_km, pval = TRUE, 
           title = "Kaplan-Meier Survival Curve survNPDR", 
           xlab = "Time (days)", 
           ylab = "Survival Probability", 
           legend.title = "Risk Group", 
           legend.labs = c("Low Risk", "High Risk"),
           font.main = c(20, "bold"), 
           font.x = c(20, "bold"), 
           font.y = c(20, "bold"), 
           font.tickslab = c(18, "plain"), 
           font.legend = c(18, "plain"), 
           font.risk.table = c(18, "plain"),
           linetype = "solid", 
           size = 1.5,
           risk.table = FALSE)


# Now read in the saved models
survNPDR <- read.csv("C:/Users/patri/OneDrive/Desktop/snpdr_update2/snpdr_update/data/GSE9891/NPDR_RANKING.csv")
#survNPDR_glm <- read.csv("C:/Users/patri/OneDrive/Desktop/snpdr_update2/snpdr_update/data/GSE_9891/survNPDR_glm_top_genes_GSE9891.csv")
cox_regression <- read.csv("C:/Users/patri/OneDrive/Desktop/snpdr_update2/snpdr_update/data/GSE9891/GSE9891COX_RANKING.csv")
#survival_svm <- read.csv("C:/Users/patri/OneDrive/Desktop/snpdr_update2/snpdr_update/data/GSE_9891/survival_svm_top_genes_GSE9891.csv")


top_features_survNPDR <- survNPDR$Feature[1:30]
top_features_cox <- cox_regression$Feature[1:30]

# Calculate the overlap between the top features
overlap_features <- intersect(top_features_survNPDR, top_features_cox)

print(overlap_features)
