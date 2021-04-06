# diff_expr_qlf.R
#
# Purpose: Perform differential expression analysis by QLF for GSE152641 RNA-seq
# Version: 1.0
# Date:2021-03-30
# Author: Feifei Li <ff.li@mail.utoronto.ca>
#
# Input: ./GSE152641_data_processing.R
# Output: tophits_qlf, upregulated_genes, downregulated_genes
# Dependencies: edgeR
#
# ToDo:
# Notes:
#
# ==============================================================================


# ====  PARAMETERS  ============================================================

source("./GSE152641_data_processing.R") # read in normalised RNA-seq matrix

# ====  PACKAGES  ==============================================================

# Required for performing differential expression analysis
if (! requireNamespace("edgeR", quietly = TRUE)) {
    BiocManager::install("edgeR")
}

# ====  PROCESS  ===============================================================

## differential expression analysis procedure
model_design    <- model.matrix(~group, data = samples_by_group)
d               <- edgeR::estimateDisp(d, model_design)
fit_qlm         <- edgeR::glmQLFit(d, model_design)
qlf_SARS2vsHC   <- edgeR::glmQLFTest(fit_qlm, coef = colnames(model_design)[2])
qlf_output_hits <- edgeR::topTags(object         = qlf_SARS2vsHC,
                                  sort.by        = "PValue",
                                  n              = nrow(gene_counts_norm))

## HGNC symbols of significant genes deferentially expressed
tophits_qlf <- entrez_to_gname[
    rownames(qlf_output_hits$table[
        which(qlf_output_hits$table$FDR < 0.05)
        , ])
]

## HGNC symbols of up-regulated significant genes
upregulated_genes <- entrez_to_gname[rownames(qlf_output_hits$table[which(
    qlf_output_hits$table$FDR < 0.05 & qlf_output_hits$table$logFC > 0
),])]

## HGNC symbols of down-regulated significant genes
downregulated_genes <- entrez_to_gname[rownames(qlf_output_hits$table[which(
    qlf_output_hits$table$FDR < 0.05 & qlf_output_hits$table$logFC < 0
),])]

# [END]