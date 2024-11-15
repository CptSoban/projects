library(DESeq2)
library(tidyverse)

# Read in contrast information and DESeq2 dataset
contrast_info <- read.csv(snakemake@input[["contrast"]])
dds <- readRDS(snakemake@input[["deseq_dataset"]])

# Function to perform differential expression analysis for given treatment and control
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
    return(res)
}

# Loop through each pair in the contrast information and generate results
for (pair in seq_len(nrow(contrast_info))) {
    pair_vector <- as.character(contrast_info[pair, ])
    conc_pairs <- paste(pair_vector, collapse = "_vs_")
    contrast_res <- contrast_conditions(pair_vector[1], pair_vector[2])

    # Write results to CSV file
    write.csv(contrast_res, snakemake@output[["mod_deg_table"]])
}
