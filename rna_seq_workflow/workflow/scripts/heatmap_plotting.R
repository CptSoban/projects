#Snakemake variables for column_title, row_title, number of top_res, sample name colum in sample info
library(ggplot2)
library(tidyverse)
library(ComplexHeatmap)

#Ranked DeSeq2 results
mod_res <- read.csv(snakemake@input[["mod_deg_table"]])

#Top ranked DeSeq2 results
top_res <- (c(head(mod_res, 50)))
top_res_df <- data.frame(top_res)

#IDs as rownames
top_res_df <- column_to_rownames(top_res_df, var = "symbol")

#DeSeq2 normalized counts 
normalized_counts <- read.csv(snakemake@input[["normalized_counts"]])
normalized_counts_df <- data.frame(normalized_counts)
#IDs as rownames
normalized_counts <- column_to_rownames(normalized_counts, var = "X")

#Calculate Z-scores
norm_counts_z <- t(apply(normalized_counts, 1, scale))
# Extract Z-scores for Top ranked results 
norm_counts_z <- norm_counts_z[match(rownames(top_res_df), rownames(norm_counts_z)), ]

#Metadata (Sample name, treatment)
sample_info <- read.csv(snakemake@input[["sample_info"]], sep = "\t")
#Sample name as rowname
sample_info <- column_to_rownames(sample_info, var = "sample")
#Reset column names to sample names 
colnames(norm_counts_z) <- rownames(sample_info)

#HEATMAP

#Build Heatmap
pdf(snakemake@output[["heatmap_plot"]])
Heatmap(norm_counts_z,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        column_labels = colnames(norm_counts_z),
        name = "Z-score",
        row_labels = rownames(norm_counts_z),
        )
dev.off()