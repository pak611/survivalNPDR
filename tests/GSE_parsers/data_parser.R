# Load necessary libraries
library(preprocessCore)
library(Biobase)
library(gcrma)
library(affy)
library(affyPLM)  # for normalization
library(dplyr)

# ------------------------------------ Set Paths & Load Files ------------------------------------

# Set the path to the GSM_data folder
gsmDataFolder <- "C:/Users/patri/OneDrive/Desktop/snpdr_update2/snpdr_update/data/GSM_data"

# Get the list of .gpr files in the GSM_data folder
fileList <- list.files(gsmDataFolder, pattern = "\\.gpr$", full.names = TRUE)

# Initialize lists to store data and sample IDs
dataList <- list()
sampleIDs <- list()

# ------------------------------------ Process GPR Files -----------------------------------------

# Loop over the files and process the data
for (i in seq_along(fileList)) {
  gprDataRaw <- read.table(fileList[i], sep = "\t", header = TRUE, skip = 31)
  
  # Aggregate by gene names (mean of duplicates)
  aggregatedData <- gprDataRaw %>%
    dplyr::group_by(Name) %>%
    dplyr::summarise(MeanValue = mean(F635.Mean, na.rm = TRUE))
  
  # Store the aggregated data
  dataList[[i]] <- aggregatedData$MeanValue
  geneIDs <- aggregatedData$Name
  sampleIDs[[i]] <- basename(fileList[i])
}

# ------------------------------------ Data Quality Control --------------------------------------

# Ensure all files have the same number of genes
numGenes <- length(dataList[[1]])
keepIndices <- which(sapply(dataList, length) == numGenes)
dataList <- dataList[keepIndices]
sampleIDs <- sampleIDs[keepIndices]

# Combine the data into a single matrix
dataMatrix <- do.call(cbind, dataList)
rownames(dataMatrix) <- geneIDs

# Rename columns with sample IDs and clean up the names
colnames(dataMatrix) <- sub("\\.gpr$", "", sampleIDs)

# ------------------------------------ Load Phenotype Data ---------------------------------------

# Load phenotype data
txtFilePath <- "C:/Users/patri/OneDrive/Desktop/snpdr_update/snpdr_update/data/GSM_data/GSE9893_clinicalData.txt"
phenotypeData <- read.table(txtFilePath, sep = "\t", header = TRUE, fill = TRUE)

# Transpose and clean up phenotype data
transposedPhenoData <- as.data.frame(t(phenotypeData))
colnames(transposedPhenoData) <- phenotypeData$"Tumor.sample"
transposedPhenoData <- transposedPhenoData[-which(rownames(transposedPhenoData) == "Tumor.sample"), ]

# Update the phenotypeData object
phenotypeData <- transposedPhenoData

# ------------------------------------ Map GSM to EB Identifiers ---------------------------------

# Load the GSM to EB identifier data
identifierPath <- "C:/Users/patri/OneDrive/Desktop/snpdr_update2/snpdr_update/data/GSM_data/identifier.csv"
identifierData <- read.csv(identifierPath, header = FALSE, stringsAsFactors = FALSE)
colnames(identifierData) <- c("GSM", "EB")

# Replace GSM with EB numbers in data matrix and match the order
colnames(dataMatrix) <- identifierData$EB[match(colnames(dataMatrix), identifierData$GSM)]
dataMatrix <- dataMatrix[, match(colnames(phenotypeData), colnames(dataMatrix))]
dataMatrix <- dataMatrix[, !is.na(colnames(dataMatrix))]

# Subset phenotype data to match data matrix columns
phenotypeData <- phenotypeData[, colnames(dataMatrix)]

# ------------------------------------ Create ExpressionSet --------------------------------------

# Convert the processed data into an ExpressionSet object
phenoData <- new("AnnotatedDataFrame", data = as.data.frame(t(phenotypeData)))
eset <- ExpressionSet(assayData = dataMatrix, phenoData = phenoData)

# Save the ExpressionSet object to a file
save(eset, file = "eset.RData")

# Load the ExpressionSet object if needed
load("eset.RData")

# ------------------------------------ Normalize Expression Data ---------------------------------

# Normalize expression data using quantile normalization and log2 transformation
ExprData <- exprs(eset)
ExprData_quantile <- normalize.quantiles(ExprData)
ExprData_quantileLog2 <- log2(ExprData_quantile)

# ------------------------------------ Extract and Process Survival Data -------------------------

# Extract phenotype and survival data
pheno_data <- pData(eset)
time <- as.numeric(trimws(pheno_data[["Decease.delay.after.surgery..months."]]))
max_time <- max(time, na.rm = TRUE)
time[is.na(time)] <- max_time
outcome <- ifelse(pheno_data[["State.of.health"]] == "deceased", 1, 0)

# Handle missing values and filter low-expressing genes
attr.mat <- t(ExprData_quantileLog2)
attr.mat <- t(apply(attr.mat, 1, function(x) {
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  return(x)
}))

# Filter low-expressing genes based on threshold
threshold <- 1
keep_genes <- rowMeans(attr.mat) > threshold

attr.mat_filtered <- attr.mat[keep_genes, ]
gene.names_filtered <- rownames(ExprData)[keep_genes]

# ------------------------------------ Finalize and Save Dataset ---------------------------------

# Finalize the dataset
dat <- cbind(time = time, status = outcome, attr.mat_filtered)
dat <- as.data.frame(dat)
colnames(dat) <- c("time", "status", gene.names_filtered)

# Save as rds file
saveRDS(dat, file = "dat.rds")


Sys.setlocale("LC_ALL", "C")
write.table(dat, "processed_data.txt", sep = "\t", row.names = FALSE)
