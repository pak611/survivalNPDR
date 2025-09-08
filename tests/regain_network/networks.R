# Load sNPDR package and prepare data
devtools::load_all("C:/Users/patri/OneDrive/Desktop/snpdr_update2/snpdr_update/sNPDR")

# Read the data file with specified encoding
library(data.table)
dat <- fread("C:/Users/patri/OneDrive/Desktop/snpdr_update2/snpdr_update/tests/processed_data.txt", 
             sep = "\t", 
             header = TRUE, 
             encoding = "UTF-8")

outcome <- dat[, .(time, status)]
attr_mat <- dat[, !c("time", "status"), with = FALSE]

#-------------------------- PREPROCESSING ------------------------------------------------

# Filter to keep the 500 genes with the lowest variance
variance_values <- apply(attr_mat, 2, var) # Calculate variance for each column
low_variance_genes <- order(variance_values)[1:500]
filtered_attr_mat <- attr_mat[, low_variance_genes, with = FALSE]


# Initialize a matrix for interaction data
num_cols <- ncol(filtered_attr_mat)
inter.data <- matrix(0, nrow = nrow(filtered_attr_mat), ncol = num_cols * num_cols)

# Initialize a vector to store column names
col_names <- vector("character", num_cols * num_cols)

# Generate pairwise interactions
for (i in 1:num_cols) {
    for (j in 1:num_cols) {
        col_index <- (i - 1) * num_cols + j
        inter.data[, col_index] <- filtered_attr_mat[[i]] * filtered_attr_mat[[j]]
        col_names[col_index] <- paste0("V", i, "V", j)
    }
}

# Assign the generated column names to the interaction data
colnames(inter.data) <- col_names

# Convert to data frame for easier handling
inter.data <- as.data.frame(inter.data)

# Add back on time and status
inter.data$time <- outcome$time
inter.data$status <- outcome$status

# Print the first few rows of the interaction data
head(inter.data)

# ------------------------------------ Univariate Cox Regression ----------------------------------

library(survival)
# Run univariate Cox regression for each gene
cox.model <- lapply(colnames(inter.data %>% select(-time, -status)), function(x) {
  tryCatch({
    formula <- as.formula(paste("Surv(time, status) ~ `", x, "`", sep = ""))
    coxFit <- coxph(formula, data = inter.data)
    print(paste0("Feature", x))
    summary(coxFit)
  }, error = function(e) {
    message(paste("Error in column:", x, "-", e$message))
    return(NULL)
  })
})

cox.model.df <- cox.model

# Remove NULL elements and extract results
cox.model <- cox.model[!sapply(cox.model, is.null)]
results <- lapply(cox.model, function(summary_obj) {
  coef_table <- summary_obj$coefficients
  data.frame(Feature = rownames(coef_table)[1], 
             Beta = coef_table[1, "coef"], 
             P.Value = coef_table[1, "Pr(>|z|)"])
})

results_df <- do.call(rbind, results)
results_df <- results_df %>% arrange(Beta)

gene_names <- colnames(filtered_attr_mat)

# Function to extract gene names by index and create the interaction string
extract_gene_interactions <- function(features, gene_names) {
  sapply(features, function(feature) {
    # Extract the indices from each interaction (e.g., "V156V341" to 156 and 341)
    indices <- as.numeric(unlist(regmatches(feature, gregexpr("\\d+", feature))))
    # Retrieve gene names by indices and form the interaction string
    paste(gene_names[indices[1]], gene_names[indices[2]], sep = "_")
  })
}

# Apply the function to the features in results_df and create the new column
results_df$Gene_Interactions <- extract_gene_interactions(results_df$Feature, gene_names)

# View the updated results_df
head(results_df, 20)

# Save the results dataframe as a CSV file
write.csv(results_df, "univariate_cox_results.csv", row.names = FALSE)
