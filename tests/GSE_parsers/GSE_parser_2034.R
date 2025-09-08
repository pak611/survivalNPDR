# Version info: R 4.2.2, Biobase 2.58.0, GEOquery 2.66.0, limma 3.54.0
################################################################
#   Differential expression analysis with limma
# Load necessary libraries
library(GEOquery)
library(limma)
library(umap)
library(Biobase)

# Load series and platform data from GEO
gset <- getGEO("GSE2034", GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]]

# Inspect the gset object
print(gset)

# Expression data
expr_data <- exprs(gset)
print(dim(expr_data))  # Print dimensions of the expression data
print(head(expr_data))  # Print the first few rows of the expression data

# Log2 normalization
expr_data <- log2(expr_data + 1)

# Transpose the expression data so that samples are rows and features are columns
expr_data <- t(expr_data)

# Read in the phenotype data from the text file
pheno_data_file <- "C:/Users/patri/OneDrive/Desktop/snpdr_update2/snpdr_update/data/GSE_2034/phenotype_data.txt"
pheno_data <- read.table(pheno_data_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Print the first few rows of the phenotype data
print(head(pheno_data))

# Ensure the sample names in the phenotype data match the sample names in the gset object
sample_names <- sampleNames(gset)
pheno_data <- pheno_data[match(sample_names, pheno_data$GEO.asscession.number), ]

print(head(pheno_data))

# Check if the sample names match and handle missing values
if (any(is.na(pheno_data$geo_accn)) || !all(sample_names == pheno_data$GEO.asscession.number, na.rm = FALSE)) {
  stop("Sample names in the phenotype data do not match the sample names in the gset object.")
}

# Assign the phenotype data to the gset object
pData(gset) <- pheno_data

# Extract the platform annotation
platform <- annotation(gset)
print(platform)

# Download the annotation file for the platform
annotation_data <- getGEO(platform, AnnotGPL = TRUE)

# Inspect the annotation data
print(head(annotation_data@dataTable@table))

# Map probe set IDs to gene names
probe_to_gene <- annotation_data@dataTable@table[, c("ID", "Gene symbol")]
print(head(probe_to_gene))

# change Gene symbol to Gene.symbol
colnames(probe_to_gene) <- c("ID", "GeneSymbol")

# Transpose the expression data so that probe set IDs are in the rows
expr_data_t <- t(expr_data)

# Convert the transposed expression data to a data frame and add row names as a column
expr_data_df <- as.data.frame(expr_data_t)
expr_data_df$ID <- rownames(expr_data_df)

# Merge the expression data with the gene names
expr_data_with_genes <- merge(expr_data_df, probe_to_gene, by = "ID", all.x = TRUE)

# Handle missing gene symbols
expr_data_with_genes$GeneSymbol[expr_data_with_genes$GeneSymbol == ""] <- paste0("Gene_", seq_len(sum(expr_data_with_genes$GeneSymbol == "")))

# Handle duplicate gene symbols by making them unique
expr_data_with_genes$GeneSymbol <- make.unique(expr_data_with_genes$GeneSymbol)

# Set the row names to the gene symbols
rownames(expr_data_with_genes) <- expr_data_with_genes$GeneSymbol

# Remove the ID and Gene.symbol columns
expr_data_with_genes <- expr_data_with_genes[, !(names(expr_data_with_genes) %in% c("ID", "Gene.symbol"))]

# Transpose the data back to the original orientation
expr_data_with_genes <- t(expr_data_with_genes)


# Print the updated expression data with gene names
print(head(expr_data_with_genes))

# Convert to numeric without getting rid of the row names
expr_data_with_genes_num <- apply(expr_data_with_genes,2, as.numeric)
rownames(expr_data_with_genes_num) <- rownames(expr_data_with_genes)

# set all NA values to 0
expr_data_with_genes_num[is.na(expr_data_with_genes_num)] <- 0

# Filter low variation genes
gene_variances <- apply(expr_data_with_genes_num, 2, function(x) var(x))



high_variance_genes <- gene_variances > quantile(gene_variances, 0.25)
expr_data_filtered <- expr_data_with_genes_num[, high_variance_genes]

# Print the filtered expression data
print(dim(expr_data_filtered))  # Print dimensions of the filtered expression data
print(head(expr_data_filtered))  # Print the first few rows of the filtered expression data

# drop the row with name GeneSymbol
expr_data_filtered <- expr_data_filtered[-1,]

# Extract RFS time and status
rfs_time <- pheno_data$time.to.relapse.or.last.follow.up..months.
rfs_status <- pheno_data$relapse..1.True.

# Ensure the lengths of rfs_time and rfs_status match the number of samples
if (length(rfs_time) != nrow(expr_data_filtered) || length(rfs_status) != nrow(expr_data_filtered)) {
  stop("The lengths of rfs_time and rfs_status do not match the number of samples in expr_data_filtered.")
}

# Append RFS time and status to the expression data
expr_data_final <- cbind(expr_data_filtered, rfs_time, rfs_status)

# Set column names to gene names
colnames(expr_data_final) <- c(colnames(expr_data_filtered), "time", "status")

# Print the updated expression data
print(dim(expr_data_final))  # Print dimensions of the updated expression data
print(head(expr_data_final))  # Print the first few rows of the updated expression data

# Save the final dataset as an RDS file
saveRDS(expr_data_final, file = "C:/Users/patri/OneDrive/Desktop/snpdr_update2/snpdr_update/data/GSE_2034/dat.rds")

# Save the final dataset as a text file
Sys.setlocale("LC_ALL", "C")
write.table(expr_data_final, "processed_data.txt", sep = "\t", row.names = FALSE)
