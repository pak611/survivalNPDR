# Install required packages
install.packages("BiocManager")
BiocManager::install("TCGAbiolinks")
BiocManager::install("SummarizedExperiment")


# Query and download the TCGA-OV RNA-Seq Dataset

library(TCGAbiolinks)

# Define the query
query <- GDCquery(
    project = "TCGA-OV",
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
)

# Download the data
GDCdownload(query)

# Prepare the data for analysis
data <- GDCprepare(query)

# View the structure of the downloaded data
str(data)


# Extract Expression Data

library(SummarizedExperiment)

# Extract the expression matrix
expr_matrix <- assay(data)

# Convert to a data frame
expr_df <- as.data.frame(expr_matrix)

# Add gene names (rownames → first column)
expr_df <- cbind(Gene = rownames(expr_df), expr_df)
rownames(expr_df) <- NULL

# Preview the expression data
head(expr_df)


# Download TCGA Clinical Data (Including Survival Information)

clinical_query <- GDCquery(
    project = "TCGA-OV",
    data.category = "Clinical",
    file.type = "xml"
)

GDCdownload(clinical_query)

clinical_data <- GDCprepare_clinic(clinical_query, clinical.info = "patient")
head(clinical_data)


# Merge Expression & Clinical Data

# Ensure patient barcodes match between datasets
rownames(expr_df) <- gsub("\\..*", "", rownames(expr_df))  # Remove extra IDs
clinical_data$patient <- gsub("-", ".", clinical_data$bcr_patient_barcode)  # Format patient IDs

# Merge datasets (matching on patient ID)
merged_data <- merge(clinical_data, t(expr_df), by.x = "patient", by.y = "Gene")
head(merged_data)


# Save the data as CSV for future use

write.csv(expr_df, "TCGA_OV_expression.csv", row.names = FALSE)
write.csv(clinical_data, "TCGA_OV_clinical.csv", row.names = FALSE)
write.csv(merged_data, "TCGA_OV_merged.csv", row.names = FALSE)
