library(WGCNA)
library(tidyverse)
library(CorLevelPlot)
library(gridExtra)

norm_data <- read.csv(snakemake@input[["normalized_counts"]])

#Vector of soft-thresholding powers
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

#Network topology anaylsis function to determine adequate power value
sft <- pickSoftThreshold(norm_data,
                        powerVector = power,
                        networkType = "signed",
                        verbose = 5 )


sft_data <- sft$fitIndices

#Get the first power value that crosses the R^2 threshold
sft_data$threshold_crossed <- ifelse(sft_data$SFT.R.sq >= 0.8,
                                    sft_data$Power,
                                    NA)

sft_power <- min(sft_data$threshold_crossed, na.rm = TRUE)

bwnet <- blockwiseModules(norm_data,
                        maxBlockSize = 15000,
                        TOMType = "signed",
                        power = sft_power,
                        mergeCutHeight = 0.25,
                        randomSeed = 27,
                        verbose = 3)