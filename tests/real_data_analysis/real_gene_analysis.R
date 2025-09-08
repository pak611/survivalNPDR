# Load necessary library
library(readr)
library(dplyr)

# Read in the CSV file
survNPDR_top_genes2 <- read_csv("survNPDR_top_genes2.csv")

# Arrange the data by absolute value of beta in descending order
survNPDR_top_genes2 <- survNPDR_top_genes2 |> arrange(desc(abs(beta)))

# Put the top 50 genes into a list
top_genes_survNPDR <- survNPDR_top_genes2$Feature[1:50]

# Read the univariate Cox results
cox_univariate_top_genes <- read_csv("univariate_cox_results.csv")

# Filter genes with P.Value < 0.05 and arrange by absolute value of Beta in descending order
cox_univariate_top_genes <- cox_univariate_top_genes |>
  filter(P.Value < 0.05) |>
  arrange(desc(abs(Beta)))

# Put the top 50 genes into a list
top_genes_cox_univariate <- cox_univariate_top_genes$Feature[1:50]

# Read the interaction Cox results
cox_interactions_top_genes <- read_csv("univariate_cox_interaction_results.csv") |> as.data.frame()

# Remove genes with P.Value > 0.05, drop the Feature column, and rename Gene_Interactions to Feature
cox_interactions_top_genes <- cox_interactions_top_genes |>
  filter(P.Value < 0.05) |>
  select(-Feature) |>
  rename(Feature = Gene_Interactions)

# Arrange the data by absolute value of Beta in descending order
cox_interactions_top_genes <- cox_interactions_top_genes |> arrange(desc(abs(Beta)))

# Put the top 50 genes into a list
top_genes_cox_interactions <- cox_interactions_top_genes$Feature[1:50]

# Calculate and print the overlap between each list
overlap_survNPDR_cox_univariate <- intersect(top_genes_survNPDR, top_genes_cox_univariate)
overlap_survNPDR_cox_interactions <- intersect(top_genes_survNPDR, top_genes_cox_interactions)
overlap_cox_univariate_cox_interactions <- intersect(top_genes_cox_univariate, top_genes_cox_interactions)

print("Overlap between survNPDR and Cox Univariate:")
print(overlap_survNPDR_cox_univariate)

print("Overlap between survNPDR and Cox Interactions:")
print(overlap_survNPDR_cox_interactions)

print("Overlap between Cox Univariate and Cox Interactions:")
print(overlap_cox_univariate_cox_interactions)

# Create a summary dataframe with the top 50 genes for each method
summary_df <- data.frame(
  Rank = 1:50,
  survNPDR = top_genes_survNPDR,
  Cox_Univariate = top_genes_cox_univariate,
  Cox_Interactions = top_genes_cox_interactions
)

# Print the summary dataframe
print(summary_df)

# Save the summary dataframe as a CSV file
write.csv(summary_df, "top_50_genes_summary.csv", row.names = FALSE)










