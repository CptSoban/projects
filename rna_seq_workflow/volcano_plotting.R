library(ggplot2)
library(EnhancedVolcano)
library(tidyverse)


deseq2_res <- read.csv(snakemake@input[["mod_deg_table"]])

top_res <- c(head(deseq2_res$symbol, 10))

#VOLCANO PLOT

#Build Volcano Plot
pdf(snakemake@output[["volcano_plot"]])
EnhancedVolcano(deseq2_res,
                x = "log2FoldChange",
                y = "padj",
                lab = deseq2_res$symbol,
                selectLab = top_res,
                pCutoff = 0.05,
                FCcutoff = 1,
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