library(DESeq2)
library(tidyverse)
library(ggplot2)
library(apeglm)

#Count table
counts_data <- read.csv(snakemake@input[["prep_counts_table"]])

#Sample info
sample_info <- read.csv(snakemake@input[["sample_info"]], sep = "\t")

#Converst first column to rownames
counts_data <- column_to_rownames(counts_data, var = "Geneid")


#Builds DESeq2 dataset from count table and sample info
dds <- DESeqDataSetFromMatrix(countData = counts_data,
                                colData = sample_info,
                                design = ~condition)

#Sets control samples as "baseline" to compare the treated samples against
# dds$condition <- relevel(dds$condition, ref = "control")

#Run DESeq2
dds <- DESeq(dds)

#Variance stabilization transformation
vsdata <- vst(dds, blind = FALSE)

#Extract normalized counts
normalized_counts <- counts(dds, normalize = TRUE)
#Create .csv from results
write.csv(normalized_counts, snakemake@output[["normalized_counts"]])


#DISPERSION
pdf(snakemake@output[["disp_plot"]])
plotDispEsts(dds)
dev.off()

#PCA
svg(snakemake@output[["pca_plot"]])
plotPCA(vsdata) + geom_text(aes(label=name),vjust=0,hjust=2)+ theme_bw()
dev.off()

saveRDS(dds, file = snakemake@output[["deseq_dataset"]])

#Replaces the log2FC standard error column with the shrinked logFC values (column name stays the same)
# res_lfc <- data.frame(lfcShrink(dds,
#                                 coef = "condition_treated_vs_control",
#                                 type = "apeglm"))

#Rename the log2FC column with the shrinked log2FC values
# res_lfc <- rename(res_lfc, log2FC_shrinked = log2FoldChange)