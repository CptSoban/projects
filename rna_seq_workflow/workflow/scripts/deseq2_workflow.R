library(DESeq2)
library(tidyverse)
library(ggplot2)
library(apeglm)

#Reads in the count table
counts_data <- read.csv(snakemake@input[["prep_counts_table"]])

#Converst first column to rownames
counts_data <- column_to_rownames(counts_data, var = "Geneid")

#Reads in sample info ("treated","control", etc.)
sample_info <- read.csv(snakemake@input[["sample_info"]], sep = "\t")

#Builds DESeq2 dataset from count table and sample info
dds <- DESeqDataSetFromMatrix(countData = counts_data,
                                colData = sample_info,
                                design = ~condition)

#Sets control samples as "baseline" to compare the treated samples against
dds$condition <- relevel(dds$condition, ref = "control")

#Run DESeq2
dds <- DESeq(dds)

#Variance stabilization transformation
vsdata <- vst(dds, blind=FALSE)

#Extract normalized counts
normalized_counts <- counts(dds, normalize = TRUE)
#Create .csv from results
write.csv(normalized_counts, snakemake@output[["normalized_counts"]])

#Convert results to dataframe
res <- data.frame(results(dds, tidy=TRUE))
#Create .csv from results
write.csv(res, snakemake@output[["deg_table"]])

#Replaces the log2FC standard error column with the shrinked logFC values (column name stays the same)
res_lfc <- data.frame(lfcShrink(dds,
                                coef = "condition_treated_vs_control",
                                type = "apeglm"))

#Read in file where IDs are mapped to symbols
#map_id <- read.csv(snakemake@input[["mapID"]])
#Add a column to results containing corresponding gene symbols
#res <- add_column(res, symbol = map_id$gene)

#Filter out non-significant expression changes
#res_sig <- results(dds, alpha = 0.05)
#Create .csv from filtered results
#write.csv(res_sig, snakemake@output[["filt_deg_table"]])

#Modifications to the results for more convinient plotting

#If no corresponding symbol exists fill in ID instead
#mod_res$symbol <- ifelse(mod_res$symbol == "", rownames(mod_res), mod_res$symbol)
#Add a column with ranking values (padj/log2FC)
#mod_res <- add_column(res, rank = log10(res$padj)/res$log2FoldChange)
#mod_res <- add_column(res_lfc, rank = log10(res$padj)/res$log2FoldChange)

#Rename the log2FC column with the shrinked log2FC values
res_lfc <- rename(res_lfc, log2FC_shrinked = log2FoldChange)
res_lfc <- rownames_to_column(res_lfc, var="symbol")
#Order the results based on shrinked log2FC (high->low)
res_lfc <- res_lfc[order(res_lfc$log2FC_shrinked,
                        decreasing = TRUE,
                        na.last = TRUE), ]

write.csv(res_lfc, snakemake@output[["mod_deg_table"]])

#DISPERSION
pdf(snakemake@output[["disp_plot"]])
plotDispEsts(dds)
dev.off()

#PCA
svg(snakemake@output[["pca_plot"]])
plotPCA(vsdata) + geom_text(aes(label=name),vjust=0,hjust=2)+ theme_bw()
dev.off()
