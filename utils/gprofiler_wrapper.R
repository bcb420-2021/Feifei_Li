# gprofiler_wrapper.R
#
# Purpose:
#     Wrappers for the g:profile R interface functions, with predefined values
# Version: 1.0
# Date: 2021-03-13
# Author: Feifei Li <ff.li@mail.utoronto.ca>
#
# Dependencies: gprofiler2
#
# ToDo: Add unit tests for the wrappers
# Notes:
#
# ==============================================================================

# ====  PACKAGES  ==============================================================

# For enrichment analysis/ORA =================================================
if (!requireNamespace("gprofiler2", quietly = TRUE))
    install.packages("gprofiler2")

# ====  FUNCTIONS  =============================================================

gostQuery <- function(genelist) {
    # Purpose:
    #     g:Profiler Query Builder. Customise query parameters in one wrapper
    # Parameters:
    #     genelist: A list of genes to query
    #     sources: Data source of pathways, optional with default values
    # Value:
    #     gostres: A list containning the query result from g:Profiler
    
    gostres <- gprofiler2::gost(
        query = genelist,
        organism = "hsapiens",
        ordered_query = F, # this option is TRUE for the ranked list
        multi_query = F,
        significant = F, # want ALL the results, not just the significant
        exclude_iea = T, # human-curated result only for the first analysis
        evcodes = F, # too long to parse
        correction_method = "fdr", # Benjamini-Hochberg for false positive
        domain_scope = "annotated",
        measure_underrepresentation = F, # default
        user_threshold = 0.05,
        sources = c("GO:BP", "REAC", "WP", "KEGG")
    )
    
    return(gostres)
}


summariseGost <- function(gostres) {
    # Purpose:
    #     g:Profiler query result brief summary
    # Parameters:
    #     gostres: a list returned by gprofiler2::gost()
    # Value:
    #     result:
    #         a dataframe with the number of genesets returned by the query,
    #         the number of genesets whose pathways pass the p-value threshold,
    #         ambiguous genes and genes that couldn't be used for ORA
    
    result <- data.frame(
        "Genesets Returned"     = nrow(gostres$result),
        "Significant Genesets " = nrow(
            gostres$result[gostres$result$significant == T, ]
        ),
        "Ambiguous Genes"       = length(
            gostres$meta$genes_metadata$ambiguous
        ),
        "Failed Genes"          = length(
            gostres$meta$genes_metadata$failed),
        check.names                      = F
    )
    
    return(result)
}


termTable <- function(gostres, lower = 1, upper = Inf, sig = TRUE, source) {
    # Purpose:
    #     Filter the g:Profiler query result by term size, significance, source
    # Parameters:
    #     gostres: A list returned by gprofileR2 query
    #     lower: smallest size of geneset to filter, default 1, optional
    #     upper: greatest size of geneset to filter, default no limit, optonal
    #     sig: Filter genesets that are significant. FALSE return all
    #     source: Get genesets by annotation sources.
    # Value:
    #     result: ...
    
    require(gprofiler2)
    
    bySource <- (gostres$result$source %in% source)
    result_source <- gostres$result[bySource, ]
    
    limit <- (result_source$term_size >= lower &
              result_source$term_size <= upper)
    
    isSig <- (result_source$significant == sig |
              result_source$significant == TRUE)
    
    result <- result_source[limit & isSig, ]

    gosttable <- data.frame(
        "Term name" = result$term_name,
        "Term ID"   = result$term_id,
        "p-value"   = result$p_value,
        check.names = FALSE
    )
    
    gosttable <- gosttable[order(gosttable[ , "p-value"]), ]
    
    return(gosttable)
}

# ====  TESTS  =================================================================
# function tests...


# [END]