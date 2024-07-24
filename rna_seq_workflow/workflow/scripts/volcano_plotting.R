library(ggplot2)
library(EnhancedVolcano)
library(tidyverse)

#Ranked DeSeq2 results (log2FC/padj)
mod_res <- read.csv(snakemake@input[["mod_deg_table"]])

top_over <- head(mod_res$log2FoldChange, 5)

#top_under <- tail(mod_res[order(mod_res$log2FC_shrinked), ], 10)

#lab_italics <- paste0("italic('", rownames(mod_res), "')")
# selectLab_italics <- paste0(
#     "italic('",
#     c('FOXG_12350', 'FOXG_03770', 'FOXG_09804', 'FOXG_20229', 'FOXG_11545'),
#     "')")
#Top ranked DeSeq2 results
# top_res <- c('FOXG_12350', 'FOXG_03770', 'FOXG_09804', 'FOXG_20229', 'FOXG_11545')

#VOLCANO PLOT

#Build Volcano Plot
pdf(snakemake@output[["volcano_plot"]])
EnhancedVolcano(mod_res,
                x = "log2FoldChange",
                y = "padj",
                lab = mod_res$symbol,
                #selectLab = top_res,
                pCutoff = 0.05,
                FCcutoff = 2,
                pointSize = 3.0,
                labSize = 4.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                colAlpha = 4/5,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                legendPosition = 'bottom',
                legendLabSize = 12,
                legendIconSize = 3.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black') + coord_flip()
dev.off()