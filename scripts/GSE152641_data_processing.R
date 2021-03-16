# GSE152641_data_processing.R
#
# Purpose: Clean up and normalise the count matrix under GSE152641 on GEOdb
# Version: 1.0
# Date: 2021-03-13
# Author: Feifei Li <ff.li@mail.utoronto.ca>
#
# Output: A TMM normalised count matrix
# Dependencies:
#     Bioconductor 3.11
#     GEOquery
#     annotate
#     org.Hs.eg.db
#     edgeR
#
# ToDo: 
# Notes: This script only functions properly with Bioconductor version 3.11!
#
# ==============================================================================


# ====  PARAMETERS  ============================================================

SERIES <- "GSE152641"

# ====  PACKAGES  ==============================================================
# Check that required packages have been installed. Install if needed.

# Run it if Bioconductor is not installed or version is not 3.11
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install(version = "3.11",
#                      ask     = FALSE)
# For downloading data file ====================================================
if (!requireNamespace("GEOquery", quietly = TRUE))
    BiocManager::install("GEOquery", ask = FALSE)
# For HGNC symbol mapping ======================================================
# Use org.Hs.eg.db database to map identifiers
if (!requireNamespace("annotate", quietly = TRUE))
    BiocManager::install("annotate", ask = FALSE)
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
    BiocManager::install("org.Hs.eg.db", ask = FALSE)
library(org.Hs.eg.db)       # Attach database for entrez to HGNC symbol mapping
# For data processing ==========================================================
if (!requireNamespace("edgeR", quietly = TRUE))
    BiocManager::install("edgeR", ask = FALSE)


# ====  PROCESS  ===============================================================
# Where we put all the data files in ===========================================
if(!dir.exists("./data")) {
    dir.create("./data", showWarnings = FALSE)
}
# Get file name from GEOdb =====================================================
fname <- GEOquery::getGEOSuppFiles(GEO           = SERIES,
                                   fetch_files   = FALSE,
                                   makeDirectory = FALSE)$fname
# Download dataset if it doesn't exist in the working directory ================
if (!file.exists(file.path(getwd(), "data", fname))) {
    GEOquery::getGEOSuppFiles(GEO           = SERIES,
                              baseDir       = "./data",
                              makeDirectory = FALSE)
}
# Read in data file ============================================================
gene_counts <- read.csv(
    file        = file.path(getwd(), "data", fname),
    header      = TRUE,
    check.names = FALSE
)
colnames(gene_counts)[1] <- "entrezgene_id" # Genes were labeled w/ Entrez ID's

# define sample groups
samples_by_group           <- grepl(pattern = "^ESC",
                                    colnames(gene_counts[2:87]))
samples_by_group           <- data.frame(samples_by_group)
rownames(samples_by_group) <- colnames(gene_counts[2:87])
colnames(samples_by_group) <- "group"
samples_by_group$group     <- ifelse(samples_by_group$group, "COVID-19", "HC")
samples_by_group           <- data.frame(
    sample = rownames(samples_by_group),
    group  = samples_by_group$group
)
# remove the duplicate for gene with Entrez ID 463
gene_counts <- gene_counts[-which(gene_counts$entrezgene_id == "388289"), ]
# Normalising count ============================================================
gene_counts_matrix <- as.matrix(gene_counts[2:87])  # DGEList uses matrix input
rownames(gene_counts_matrix) <- gene_counts$entrezgene_id # Label w/ Entrez IDs
# Remove outliers by edgeR protocol:
keep <- (rowSums(edgeR::cpm(gene_counts_matrix) > 1) >= 24)
gene_counts_matrix <- gene_counts_matrix[keep, ]       # remove low-count genes
gene_counts_DGE    <- edgeR::DGEList(
    counts = gene_counts_matrix,
    group  = samples_by_group$group
)
d <- edgeR::calcNormFactors(gene_counts_DGE)     # compute normalization factor
gene_counts_norm <- edgeR::cpm(d)                   # extract normalised matrix
# Gene name mapping ============================================================
entrez_to_gname <- annotate::getSYMBOL(
    x    = rownames(gene_counts_norm),
    data = 'org.Hs.eg')
# Manually label the unmapped gene with updated its HGNC symbol
entrez_to_gname["285464"] <- "CRIPAK"
rownames(gene_counts_norm) <- entrez_to_gname[rownames(gene_counts_norm)]

# [END]
