library(gprofiler2)
library(tidyverse)

#Reads in the count table
counts_df <- read.csv(snakemake@input[["preranked_genes"]])

#Convert gene name column to vector
gene_vector <- c(counts_df$symbol)

# converted_symbols_df <- gconvert(query = gene_vector,
#                                 organism = "foxysporum",
#                                 target = "GO",
#                                 mthreshold = Inf,
#                                 filter_na = TRUE)

#Function Enrichment Anaylsis of preranked gene set
enriched_functions <- gost(query = gene_vector,
                             organism = "foxysporum",
                             ordered_query = TRUE,
                             multi_query = FALSE,
                             significant = TRUE,
                             exclude_iea = FALSE,
                             measure_underrepresentation = FALSE,
                             evcodes = FALSE,
                             user_threshold = 0.05,
                             correction_method = "g_SCS",
                             domain_scope = "annotated",
                             custom_bg = NULL,
                             numeric_ns = "",
                             sources = c("GO:BP", "KEGG", "REAC"),
                             as_short_link = FALSE,
                             highlight = FALSE
                             )

res <- enriched_functions$result

mod_res <- res[, c("source",
                    "term_id",
                    "term_name",
                    "term_size",
                    "intersection_size",
                    "p_value"
                    )]

write.csv(mod_res, snakemake@output[["enriched_functions"]])