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
dat <- suppressWarnings(fread("./snpdr_update/tests/processed_data.txt", sep = "\t", header = TRUE, encoding = "UTF-8"))

#########################
# Preprocessing
#########################

# Extract outcome and attributes
outcome <- dat[, .(time, status)]
attr_mat <- dat[, !c("time", "status"), with = FALSE]

# Filter to keep the 500 genes with the lowest variance
variance_values <- apply(attr_mat, 2, var)
low_variance_genes <- order(variance_values)[1:500]
attr_mat <- attr_mat[, low_variance_genes, with = FALSE]

# Combine data
dats <- cbind(attr_mat, status = outcome$status)

#########################
# Process Cox Results
#########################

# Read in univariate Cox interaction results
cox_interactions <- read.csv("univariate_cox_interaction_results.csv")
cox_univariate <- read.csv("univariate_cox_results.csv")

# Extract unique genes from interactions
gene_interactions <- cox_interactions$Gene_Interactions
gene_list <- str_split(gene_interactions, "_")
unique_genes <- unique(unlist(gene_list))

# Subset Cox univariate results to include only the unique genes
cox_univariate_subset <- cox_univariate[cox_univariate$Feature %in% unique_genes, ] %>% select(-X)

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
    betas[gene, gene] <- cox_univariate_subset$Beta[i]
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
interaction_threshold <- 0.3  # Adjust as needed
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

# Define layout
my.layout <- layout_with_fr(A.adj, niter = 100)

# Node properties
vertex_sizes <- abs(main_effects) * 20  # Scale node size by main effect magnitude
vertex_colors <- ifelse(main_effects > 0, "limegreen", "gray")  # Color by positive/negative main effect

# Edge properties
edge_weights <- E(A.weight)$weight
max_weight <- max(abs(edge_weights))  # Maximum weight for normalization

# Map interaction strengths to colors (gradient from red to gray to green)
color_palette <- colorRampPalette(c("red", "gray", "limegreen"))
edge_colors <- color_palette(100)[as.numeric(cut(edge_weights, breaks = 100, labels = FALSE))]
E(A.adj)$color <- edge_colors
E(A.adj)$width <- abs(edge_weights) * 4  # Scale edge width by interaction weight

#########################
# Plot the Network
#########################

# Increase the font size of the gene names
vertex_label_cex <- 0.7  # Adjust the font size as needed

# Save the plot as a high-resolution PNG file
png("network_plot.png", width = 2000, height = 2000, res = 300)

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
