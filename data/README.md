# Data Directory

The original datasets (~1.7 GB) are omitted from the repository to keep the clone lightweight.

Contents that were present locally:
- GEO expression and clinical datasets: GSE13876, GSE1456, GSE17260, GSE2034, GSE2990, GSE32062, GSE9891, GSE9893
- Aggregated objects: `dat.rds`, `eset.RData`, `processed_data.txt`, `identifier.csv`
- Supplementary: `GSE9893_clinicalData.txt`, `GSE1456_suppl_info.txt`, `ID_REF_to_GSM.txt`

## Reconstitution
1. Download GEO series matrices from NCBI GEO using accession IDs above (e.g. via GEOquery):
```r
library(GEOquery)
ids <- c("GSE13876","GSE1456","GSE17260","GSE2034","GSE2990","GSE32062","GSE9891","GSE9893")
raw <- lapply(ids, getGEO)
```
2. Harmonize probe IDs and clinical metadata to match scripts in `tests/real_data_analysis/`.
3. Save processed ExpressionSets or data frames into this directory with the same filenames expected by the scripts.

## Optional Network Inputs
If running network scripts under `tests/regain_network/`, ensure any precomputed network or intermediate CSV files are regenerated or copied here before execution.

## Note
If you need the exact original processed files for reproducibility beyond what GEO reconstruction permits, host an external archive (e.g. Zenodo, Figshare) and document the DOI here.
