# plotters.R
#
# Purpose: Helper functions to plot conveniently
# Version: 1.0
# Date: 2021-03-13
# Author: Feifei Li <ff.li@mail.utoronto.ca>
#
# Dependencies: Bioconductor, ComplexHeatmap, circlize
#
# ToDo: Modify plotHeatMap s.t. it annotates multiple groups rather than only 2
# Notes:
#
# ==============================================================================

# ====  PACKAGES  ==============================================================
# Dependency for ComplexHeatmap
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
    BiocManager::install(version = "3.11", # release version is optional here
                         ask     = FALSE)
}
# For ploting heatmap ==========================================================
if (!requireNamespace("ComplexHeatmap", quietly = TRUE))
    BiocManager::install("ComplexHeatmap")
if (!requireNamespace("circlize", quietly = TRUE))
    install.packages("circlize")

# ====  FUNCTIONS  =============================================================


heatmap_colours <- function(mat) {
    # Purpose:
    #     A helper function to produce a colour gradient for a heatmap
    # Parameters:
    #     mat: A RNA-seq count matrix
    # Value:
    #     heatmap_cols: a vector of colour gradient indicating expression level
    
    if (min(mat) == 0) {
        heatmap_cols <- circlize::colorRamp2(
            c(min(mat), max(mat)),
            c("white", "red")  # hotter colours for more upregulated expression
        )
    } else {
        heatmap_cols <- circlize::colorRamp2(
            c(min(mat), 0, max(mat)),
            c("blue", "white", "red") # white for same level
        )       # colder colours for more downregulated
    }
    
    return(heatmap_cols)
}


plotHeatMap <- function(tophits = integer(0),
                        m,
                        samples,
                        title = "",
                        group_cols = c("#51af84", "#b7371d")) {
    # Purpose:
    #     Standardise values in the input count matrix to Z-scores, and
    #     generate a heatmap plot for gene expression with samples ordered
    #     by their groups and annotated.
    # Parameters:
    #     tophits:
    #         Top-hit genes by differential expression analysis, optional
    #     m      : RNA-seq count matrix
    #     samples:
    #              A dataframe specifying to which group a sample belongs.
    #              Note that it only supports samples divided into 2 groups:
    #              control and experimental
    #     title:   Title for the heatmap
    #     group_cols:
    #              Colours to annotate the heatmap by groups.
    #              Use preset colours if not supplied.
    # Value:
    #     gene_heatmap: A heatmap object
    
    require(ComplexHeatmap)
    
    # Row normalisation/standardisation to Z-scores
    if (length(tophits) == 0) {
        m <- t(
            scale(
                t(m)            # transpose matrix => row = samples, col = genes
            )                  # convert values of each column(gene) to Z scores
        )       # transpose the transposed matrix => row =  genes, col = samples
    } else {
        m <- t(scale(t(m[which(rownames(m) %in% tophits), ])))
    }
    
    # Get group names of samples
    groups <- unique(samples$group)
    
    # Sample order divided by group
    group_order <- c(which(samples$group == groups[2]),
                     which(samples$group == groups[1]))
    
    # Annotate heatmap by groups
    ha_cols <- group_cols
    names(ha_cols) <- groups
    ha <- ComplexHeatmap::HeatmapAnnotation(
        df  = data.frame(Group = samples$group[group_order]),
        col = list(Group = ha_cols)
    )
    # Generate a heatmap
    gene_heatmap <- ComplexHeatmap::Heatmap(
        matrix              = m[ , group_order],
        column_title        = title,
        name                = "Expr lvl",
        cluster_rows        = TRUE,
        cluster_columns     = FALSE,
        show_row_dend       = TRUE,
        show_column_dend    = TRUE,
        col                 = heatmap_colours(m),
        show_column_names   = FALSE,
        show_row_names      = FALSE,
        show_heatmap_legend = TRUE,
        top_annotation      = ha
    )
    
    return(gene_heatmap)
}


plotVolcano <- function(df, title, gene_of_interest = integer(0)) {
    # Purpose:
    #     Plot a vocanol plot conveniently with highlighted genes of interest.
    # Parameters:
    #     df: A dataframe dedicated to this dataset,
    #         must include the following ordered columns, regardless of names:
    #         column 1: HGNC gene symbol
    #         column 2: Fold change on log2 scale
    #         column 3: negated P-value on log10 scale
    #     title: title of the plot
    #     gene_of_interst:
    #         A vector of genes to highlight, optional
    # Value:
    #     result: A volcano plot with genes of interest highlighted.
    
    cols <- c("insig"     = "grey",
              "sig"       = "#2cbe88",
              "highlight" = "#d03b41")
    
    df$colour <- cols["insig"]
    
    sig <- which(df[ , 3] > -log10(0.05))
    df[sig, "colour"] <- cols["sig"]
    
    highlight <- integer(0)
    if (length(gene_of_interest) > 0) {
        highlight <- which(df[ , 1] %in% gene_of_interest)
    }
    x <- df[ ,2]
    y <- df[ ,3]
    if (length(highlight) > 0) {
        df[highlight, "colour"] <- cols["highlight"]
    }
    
    plot(x    = df[-highlight, 2],
         y    = df[-highlight, 3],
         col  = df[-highlight, "colour"],
         xlab = expression("Log"[2]*" fold change"),
         ylab = expression("-Log"[10]*" P"),
         main = title)
    if (length(highlight) > 0) {
        points(x      = df[highlight, 2],
               y      = df[highlight, 3],
               col    = df[highlight, "colour"],
               pch    = 8,
               cex    = 1.5,
               lwd    = 2)
        text(x        = df[highlight, 2],
             y        = df[highlight, 3],
             labels   = df[highlight, 1],
             cex      = 0.75)
        legend(x      = min(df[ , 2]),
               y      = max(df[ , 3]),
               legend = c("has evidence of DE",
                          "insignificant",
                          "gene of interest"),
               col    = c(cols["sig"], cols["insig"], cols["highlight"]),
               pch    = c(1, 1, 8),
               cex    = 0.75)
        
    } else {
        legend(x      = min(df[ , 2]),
               y      = max(df[ , 3]),
               legend = c("has evidence of DE", "insignificant"),
               col    = c(cols["sig"], cols["insig"]),
               cex    = 0.75)
    }
    abline(h    = -log10(0.05),
           col  = "#0c82b7",
           lty  = 2,
           lwd  = 1.5)
    text(x      = min(df[ , 2]) + 0.25,
         y      = -log10(0.05)  + 0.5,
         labels = "p = 0.05",
         col    = "#0c82b7")
    
}


# ====  TESTS  =================================================================
# function tests...


# [END]