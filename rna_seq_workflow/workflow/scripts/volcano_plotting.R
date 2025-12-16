library(ggplot2)
library(EnhancedVolcano)
library(tidyverse)

#Ranked DeSeq2 results (log2FC/padj)
mod_res <- read.csv(snakemake@input[["mod_deg_table"]])

top_over <- head(mod_res$log2FoldChange, 5)


#VOLCANO PLOT

#Build Volcano Plot
pdf(snakemake@output[["volcano_plot"]])
EnhancedVolcano(mod_res,
                x = "log2FoldChange",
                y = "padj",
                lab = rep('', nrow(mod_res)),
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