library(DESeq2)
library(tidyverse)

#Contrasting conditions info
contrast_info <- read.csv(snakemake@input[["contrast"]])

dds <- readRDS(snakemake@input[["deseq_dataset"]])

#Convert signficant (padj = 0.05) results to dataframe
contrast_conditions <- function(treatment, control){
    res <- data.frame(results(dds, contrast = c("condition",
                                        treatment,
                                        control)))
    res <- rownames_to_column(res, var = "symbol")
    res$row <- NULL

    #Order the results based on log2FC (high->low)
    res <- res[order(res$log2FoldChange,
                            decreasing = TRUE,
                            na.last = TRUE), ]
}

for (pair in seq_len(nrow(contrast_info))) {
    pair_vector <- as.character(contrast_info[pair,])
    conc_pairs <- paste(pair_vector, collapse = "_vs_")
    contrast_res <- contrast_conditions(pair_vector[1], pair_vector[2])
    write.csv(contrast_res, gsub(" ", "",
                            paste("results/DEG_analysis/",
                                    snakemake@params[["run_id"]],
                                    "/",
                                    snakemake@params[["run_id"]],
                                    "_",
                                    conc_pairs,
                                    ".csv")))
}
