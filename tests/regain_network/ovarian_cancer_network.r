
library(rprojroot)
root <- rprojroot::find_root(rprojroot::has_file("README.md"))
pkg_path <- file.path(root, "sNPDR")
## Load sNPDR package relative to project root
devtools::load_all(pkg_path)

library(gridExtra)
library(ggplot2)
library(patchwork)
library(survminer)
library(survivalsvm)
library(Biobase)

datasets <- c("GSE9891", "GSE32062", "GSE13876")
plot_list <- list()

devtools::load_all(pkg_path)

dataset <- datasets[2]  # Change index to select different dataset
load(file.path(root, "data", dataset, paste0(dataset, ".rda")))
gset <- get(dataset)
expr_data <- exprs(gset)
feature_data <- fData(gset)
probe_to_gene <- feature_data[, c("probeset", "gene")]
colnames(probe_to_gene) <- c("ID", "gene")
expr_data_df <- as.data.frame(expr_data)
expr_data_df$ID <- rownames(expr_data_df)
expr_data_with_genes <- merge(expr_data_df, probe_to_gene, by = "ID", all.x = TRUE)
expr_data_with_genes$gene <- make.unique(as.character(expr_data_with_genes$gene))
rownames(expr_data_with_genes) <- expr_data_with_genes$gene
expr_data_with_genes <- expr_data_with_genes[, !(names(expr_data_with_genes) %in% c("ID", "gene"))]
pheno_data <- pData(gset)
pheno_data <- pheno_data[, c("days_to_death", "vital_status")]
pheno_data$vital_status <- ifelse(pheno_data$vital_status == "deceased", 1, 0)
sample_names <- colnames(expr_data_with_genes)
pheno_data <- pheno_data[match(sample_names, rownames(pheno_data)), ]
combined_data <- cbind(pheno_data, t(expr_data_with_genes))
dat <- as.data.frame(combined_data)
dat$time <- dat$days_to_death
dat$status <- dat$vital_status
dat <- dat[, setdiff(names(dat), c("days_to_death", "vital_status"))]
dat <- dat[complete.cases(dat[, c("time", "status")]), ]
attr_mat <- dat[, setdiff(names(dat), c("time", "status"))]
attr_mat <- attr_mat[, sapply(attr_mat, is.numeric)]
variance_values <- apply(attr_mat, 2, var, na.rm = TRUE)
top_variance_genes <- names(sort(variance_values, decreasing = TRUE)[1:15000])
#top_variance_genes <- names(sort(variance_values, decreasing = TRUE)[1:100])
attr_mat <- attr_mat[, top_variance_genes]
dat <- cbind(attr_mat, time = dat$time, status = dat$status)

outcome <- dat[, c("time", "status")]


#-------------------------- PREPROCESSING ------------------------------------------------

# Filter to keep the 500 genes with the lowest variance
variance_values <- apply(attr_mat, 2, var) # Calculate variance for each column
low_variance_genes <- order(variance_values)[1:500]
filtered_attr_mat <- attr_mat[, low_variance_genes]


# Initialize a matrix for interaction data
num_cols <- ncol(filtered_attr_mat)
inter.data <- matrix(0, nrow = nrow(filtered_attr_mat), ncol = num_cols * num_cols)

# Initialize a vector to store column names
col_names <- vector("character", num_cols * num_cols)

# # Generate pairwise interactions
# for (i in 1:num_cols) {
#     for (j in 1:num_cols) {
#         col_index <- (i - 1) * num_cols + j
#         inter.data[, col_index] <- filtered_attr_mat[[i]] * filtered_attr_mat[[j]]
#         col_names[col_index] <- paste0("V", i, "V", j)
#     }
# }

# # Assign the generated column names to the interaction data
# colnames(inter.data) <- col_names

# # Convert to data frame for easier handling
# inter.data <- as.data.frame(inter.data)

# # Add back on time and status
# inter.data$time <- outcome$time
# inter.data$status <- outcome$status

# # Print the first few rows of the interaction data
# head(inter.data)

# # ------------------------------------ Univariate Cox Regression ----------------------------------

# library(survival)
# library(parallel)

# # Run univariate Cox regression for each gene in parallel
# library(parallel)
# num_cores <- detectCores() - 1  # Leave one core free

# cox.model <- mclapply(
#   colnames(inter.data %>% dplyr::select(-time, -status)),
#   function(x) {
#     tryCatch({
#       formula <- as.formula(paste("Surv(time, status) ~ `", x, "`", sep = ""))
#       coxFit <- coxph(formula, data = inter.data)
#       print(paste0("Feature", x))
#       summary(coxFit)
#     }, error = function(e) {
#       message(paste("Error in column:", x, "-", e$message))
#       return(NULL)
#     })
#   },
#   mc.cores = num_cores
# )

# cox.model.df <- cox.model

# # Remove NULL elements and extract results
# cox.model <- cox.model[!sapply(cox.model, is.null)]
# results <- lapply(cox.model, function(summary_obj) {
#   coef_table <- summary_obj$coefficients
#   data.frame(Feature = rownames(coef_table)[1], 
#              Beta = coef_table[1, "coef"], 
#              P.Value = coef_table[1, "Pr(>|z|)"])
# })

# results_df <- do.call(rbind, results)
# results_df <- results_df %>% arrange(Beta)

# gene_names <- colnames(filtered_attr_mat)

# # Function to extract gene names by index and create the interaction string
# extract_gene_interactions <- function(features, gene_names) {
#   sapply(features, function(feature) {
#     # Extract the indices from each interaction (e.g., "V156V341" to 156 and 341)
#     indices <- as.numeric(unlist(regmatches(feature, gregexpr("\\d+", feature))))
#     # Retrieve gene names by indices and form the interaction string
#     paste(gene_names[indices[1]], gene_names[indices[2]], sep = "_")
#   })
# }

# # Apply the function to the features in results_df and create the new column
# results_df$Gene_Interactions <- extract_gene_interactions(results_df$Feature, gene_names)

# # View the updated results_df
# head(results_df, 20)

# # Save the results dataframe as a CSV file
# write.csv(results_df, "univariate_cox_results.csv", row.names = FALSE)


#------------------------------------ survival NPDR on Ovarian Data ------------------------------------------------

devtools::load_all(pkg_path)

library(gridExtra)
library(ggplot2)
library(patchwork)
library(survminer)
library(survivalsvm)
library(Biobase)

datasets <- c("GSE9891", "GSE32062", "GSE13876")
plot_list <- list()

devtools::load_all(pkg_path)


dataset <- datasets[1]  # Change index to select different dataset
load(file.path(root, "data", dataset, paste0(dataset, ".rda")))
gset <- get(dataset)
expr_data <- exprs(gset)
feature_data <- fData(gset)
probe_to_gene <- feature_data[, c("probeset", "gene")]
colnames(probe_to_gene) <- c("ID", "gene")
expr_data_df <- as.data.frame(expr_data)
expr_data_df$ID <- rownames(expr_data_df)
expr_data_with_genes <- merge(expr_data_df, probe_to_gene, by = "ID", all.x = TRUE)
expr_data_with_genes$gene <- make.unique(as.character(expr_data_with_genes$gene))
rownames(expr_data_with_genes) <- expr_data_with_genes$gene
expr_data_with_genes <- expr_data_with_genes[, !(names(expr_data_with_genes) %in% c("ID", "gene"))]
pheno_data <- pData(gset)
pheno_data <- pheno_data[, c("days_to_death", "vital_status")]
pheno_data$vital_status <- ifelse(pheno_data$vital_status == "deceased", 1, 0)
sample_names <- colnames(expr_data_with_genes)
pheno_data <- pheno_data[match(sample_names, rownames(pheno_data)), ]
combined_data <- cbind(pheno_data, t(expr_data_with_genes))
dat <- as.data.frame(combined_data)
dat$time <- dat$days_to_death
dat$status <- dat$vital_status
dat <- dat[, setdiff(names(dat), c("days_to_death", "vital_status"))]
dat <- dat[complete.cases(dat[, c("time", "status")]), ]
attr_mat <- dat[, setdiff(names(dat), c("time", "status"))]
attr_mat <- attr_mat[, sapply(attr_mat, is.numeric)]
variance_values <- apply(attr_mat, 2, var, na.rm = TRUE)
top_variance_genes <- names(sort(variance_values, decreasing = TRUE)[1:15000])
#top_variance_genes <- names(sort(variance_values, decreasing = TRUE)[1:100])
attr_mat <- attr_mat[, top_variance_genes]
dat <- cbind(attr_mat, time = dat$time, status = dat$status)

outcome <- dat[, c("time", "status")]

# Now run survival NPDR on the interaction data
snpdr.model <- sNPDR::npdr_surv_binomial(
  outcome = c("time_var" = "time", "status_var" = "status"),
  dataset = dat, 
  attr.diff.type = "standard",
  nbd.method = "relieff", 
  nbd.metric = "manhattan",
  knn = 20, 
  msurf.sd.frac = 0.5,
  glmnet.alpha = 1,
  model.type = "binomial",
  KM.weight = FALSE,
  KM.kernel.type = "gaussian",
  KM.kernel.sigma = 1.5
)


head(snpdr.model, 20)

# Save the results dataframe as a CSV file
write.csv(snpdr.model[1:5000, ], file.path(root, "paper_tables", paste0("snpdr_", dataset, "_results.csv")), row.names = FALSE)


#-------------------------------- Cox Univariate (not pairwise) ------------------------------------------------


# --- Cox Univariate Regression Feature Extraction ---
library(survival)

# Assume your data frame is called 'dat' with columns: time, status, feature1, feature2, ...
# If your data frame is named differently, change 'dat' below.
features <- setdiff(colnames(dat), c("time", "status"))
results <- data.frame(Feature = character(), Beta = numeric(), P.Value = numeric(), stringsAsFactors = FALSE)


for (feature in features) {
  # Skip features not present or all NA
  if (!(feature %in% colnames(dat))) next
  if (all(is.na(dat[[feature]]))) next
  # Use backticks to handle special characters in feature names
  formula <- as.formula(paste("Surv(time, status) ~ `", feature, "`", sep = ""))
  fit <- survival::coxph(formula, data = dat)
  summary_fit <- summary(fit)
  results <- rbind(results, data.frame(
    Feature = feature,
    Beta = summary_fit$coefficients[1, "coef"],
    P.Value = summary_fit$coefficients[1, "Pr(>|z|)"]
  ))
}

# View results
print(head(results))

write.csv(results, file.path(root, "paper_tables", paste0("cox_univariate_results_", dataset, ".csv")), row.names = FALSE)


#########################
# Load Data and Libraries
#########################

# Load required libraries
library(survcomp)
library(dplyr)
library(stringr)
library(tidyr)
library(data.table)
library(igraph)

# Load data
#dat <- suppressWarnings(fread("./snpdr_update/tests/processed_data.txt", sep = "\t", header = TRUE, encoding = "UTF-8"))

#########################
# Preprocessing
#########################

# dat <- read.csv(file.path(root, "paper_tables", "top_genes_GSE2034.csv"))

# # Extract outcome and attributes
# outcome <- dat[, .(time, status)]
# attr_mat <- dat[, !c("time", "status"), with = FALSE]

# # Filter to keep the 500 genes with the lowest variance
# variance_values <- apply(attr_mat, 2, var)
# low_variance_genes <- order(variance_values)[1:500]
# attr_mat <- attr_mat[, low_variance_genes, with = FALSE]

# # Combine data
# dats <- cbind(attr_mat, status = outcome$status)

#########################
# Process Cox Results
#########################

# Read in univariate Cox interaction results
cox_interactions <- read.csv(file.path(root, "paper_tables", "interaction_cox_results_GSE9891.csv"))
cox_univariate <- read.csv(file.path(root, "paper_tables", "cox_univariate_results_GSE9891.csv"))
# cox_interactions <- read.csv(file.path(root, "paper_tables", "interaction_cox_results_GSE32062.csv"))
# cox_univariate <- read.csv(file.path(root, "paper_tables", "cox_univariate_results_GSE32062.csv"))

# Extract unique genes from interactions
gene_interactions <- cox_interactions$Gene_Interactions
gene_list <- str_split(gene_interactions, "_")
unique_genes <- unique(unlist(gene_list))

# Subset Cox univariate results to include only the unique genes
cox_univariate_subset <- cox_univariate[cox_univariate$Feature %in% unique_genes, ]

# Split interactions into individual genes and recompute unique genes
cox_interactions_split <- cox_interactions %>%
  separate(Gene_Interactions, into = c("Gene1", "Gene2"), sep = "_") %>%
  select(-Feature)
unique_genes <- unique(c(cox_interactions_split$Gene1, cox_interactions_split$Gene2))

# Create a beta matrix
num_genes <- length(unique_genes)
betas <- matrix(0, nrow = num_genes, ncol = num_genes, dimnames = list(unique_genes, unique_genes))

# Fill diagonal with univariate beta values
for (i in 1:nrow(cox_univariate_subset)) {
  gene <- cox_univariate_subset$Feature[i]
  if (gene %in% unique_genes) {
    betas[gene, gene] <- cox_univariate_subset$beta[i]
  }
}

# Fill off-diagonal with interaction beta values
for (i in 1:nrow(cox_interactions_split)) {
  gene1 <- cox_interactions_split$Gene1[i]
  gene2 <- cox_interactions_split$Gene2[i]
  if (gene1 %in% unique_genes && gene2 %in% unique_genes) {
    betas[gene1, gene2] <- cox_interactions_split$Beta[i]
    betas[gene2, gene1] <- cox_interactions_split$Beta[i]  # Assuming symmetry
  }
}

#########################
# Filter Based on Interactions
#########################

main_betas <- diag(betas)
# Set the diagonal to zero to exclude main effects
diag(betas) <- 0


# Threshold for interactions
interaction_threshold <- 0.075  # for dataset 1 (GSE 9891)
#interaction_threshold <- 0.175  # for dataset 2 (GSE 32062)
#interaction_threshold <- 0.025  # for dataset 3 (GSE 13876)
regain.nodiag.adj <- as.matrix(abs(betas) > interaction_threshold) + 0  # Threshold for edges
regain.nodiag.weight <- regain.nodiag.adj * betas

# Remove genes with no significant off-diagonal interactions
degree_mask <- rowSums(regain.nodiag.adj) > 0
regain.nodiag.adj <- regain.nodiag.adj[degree_mask, degree_mask]
regain.nodiag.weight <- regain.nodiag.weight[degree_mask, degree_mask]

# Keep only main effects for retained genes
main_effects <- main_betas[degree_mask]

# Create igraph object
A.adj <- graph_from_adjacency_matrix(regain.nodiag.adj, mode = "undirected")
A.weight <- graph_from_adjacency_matrix(regain.nodiag.weight, mode = "undirected", weight = TRUE)

#########################
# Define Node and Edge Properties
#########################
# Use Kamada-Kawai layout for better hub separation
my.layout <- layout_with_kk(A.adj)

# Node properties: highlight hubs
deg <- degree(A.adj)
hub_threshold <- quantile(deg, 0.9)  # Top 10% degree as hubs
vertex_sizes <- ifelse(deg > hub_threshold, abs(main_effects) * 800, abs(main_effects) * 400)
vertex_colors <- ifelse(deg > hub_threshold, "red", ifelse(main_effects > 0, "limegreen", "gray"))

# Edge properties
edge_weights <- E(A.weight)$weight
max_weight <- max(abs(edge_weights))  # Maximum weight for normalization

# Map interaction strengths to colors (gradient from red to gray to green)
color_palette <- colorRampPalette(c("red", "gray", "limegreen"))
edge_colors <- color_palette(100)[as.numeric(cut(edge_weights, breaks = 100, labels = FALSE))]
E(A.adj)$color <- edge_colors
E(A.adj)$width <- abs(edge_weights) * 2  # Scale edge width by interaction weight

#########################
# Plot the Network
#########################

# Increase the font size of the gene names
vertex_label_cex <- 0.3  # Adjust the font size as needed

# Save the plot as a high-resolution PNG file
png(file.path(root, "network_plot_GSE9891.png"), width = 2000, height = 2000, res = 300)
# png(file.path(root, "network_plot_GSE32062.png"), width = 2000, height = 2000, res = 300)
# png(file.path(root, "network_plot_GSE13876.png"), width = 2000, height = 2000, res = 300)

# Plot the Network
plot.igraph(
  A.adj,
  layout = my.layout,
  vertex.color = vertex_colors,
  vertex.size = vertex_sizes,
  edge.color = E(A.adj)$color,
  edge.width = E(A.adj)$width,
  vertex.label.dist = 0.5,
  vertex.label.color = "black",
  vertex.label.cex = vertex_label_cex
)

# # Add the legend
# legend(
#   x = -1.0, y = -1.0,  # Adjust these values to move the legend further down and to the right
#   legend = c("Strong Negative", "Neutral", "Strong Positive"),
#   col = c("red", "gray", "limegreen"),
#   lwd = 3,
#   title = "Edge Strength"
# )

# Close the PNG device
dev.off()




write_graph(A.adj, "network.graphml", format = "graphml")

