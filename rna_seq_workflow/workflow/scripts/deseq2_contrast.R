library(DESeq2)
library(tidyverse)


# Get the contrast pair
contrast_pair <- snakemake@params[["contrast_pair"]]

# Read in DESeq2 dataset
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

# Perform analysis
contrast_res <- contrast_conditions(contrast_pair[1], contrast_pair[2])

# Write results to the specified output file
write.csv(contrast_res, snakemake@output[["mod_deg_table"]], row.names = FALSE)
