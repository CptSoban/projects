library(ggplot2)
library(tidyverse)
library(ComplexHeatmap)

# Read normalized counts
normalized_counts <- read.csv(snakemake@input[["normalized_counts"]])
normalized_counts <- column_to_rownames(normalized_counts, var = "X")

# Z-score normalization
norm_counts_z <- t(apply(normalized_counts, 1, scale))

# Metadata
sample_info <- read.csv(snakemake@input[["sample_info"]])
sample_info <- column_to_rownames(sample_info, var = "sample")
colnames(norm_counts_z) <- rownames(sample_info)

# Initial clustering to assign row clusters
heatmap <- Heatmap(norm_counts_z,
                   show_row_names = FALSE,
                   show_row_dend = TRUE,
                   cluster_rows = TRUE,
                   cluster_columns = TRUE,
                   name = "Z-score")

ht <- draw(heatmap)
row_dend <- row_dend(ht)
row_clusters <- cutree(as.hclust(row_dend), k = 9)

# Extract gene lists per cluster
cluster1_genes <- rownames(norm_counts_z)[row_clusters == 1]
cluster2_genes <- rownames(norm_counts_z)[row_clusters == 2]
cluster3_genes <- rownames(norm_counts_z)[row_clusters == 3]
cluster4_genes <- rownames(norm_counts_z)[row_clusters == 4]
cluster5_genes <- rownames(norm_counts_z)[row_clusters == 5]
cluster6_genes <- rownames(norm_counts_z)[row_clusters == 6]
cluster7_genes <- rownames(norm_counts_z)[row_clusters == 7]
cluster8_genes <- rownames(norm_counts_z)[row_clusters == 8]
cluster9_genes <- rownames(norm_counts_z)[row_clusters == 9]



# Save gene lists if desired
write.csv(cluster1_genes, "cluster1_genes.csv", row.names = FALSE)
write.csv(cluster2_genes, "cluster2_genes.csv", row.names = FALSE)
write.csv(cluster3_genes, "cluster3_genes.csv", row.names = FALSE)
write.csv(cluster4_genes, "cluster4_genes.csv", row.names = FALSE)
write.csv(cluster5_genes, "cluster5_genes.csv", row.names = FALSE)
write.csv(cluster6_genes, "cluster6_genes.csv", row.names = FALSE)
write.csv(cluster7_genes, "cluster7_genes.csv", row.names = FALSE)
write.csv(cluster8_genes, "cluster8_genes.csv", row.names = FALSE)
write.csv(cluster9_genes, "cluster9_genes.csv", row.names = FALSE)


# Add annotation using the **correct matching order**
row_ha <- rowAnnotation(
  cluster = row_clusters,
  col = list(cluster = c(
    "1" = "red",
    "2" = "blue",
    "3" = "green",
    "4" = "purple",
    "5" = "orange",
    "6" = "brown",
    "7" = "pink",
    "8" = "yellow",
    "9" = "gray"
  ))
)


# Reuse the same row clustering to keep annotations aligned
ht_colored <- Heatmap(norm_counts_z,
                      name = "Z-score",
                      show_row_names = FALSE,
                      cluster_rows = as.hclust(row_dend),  # reuse the clustering
                      cluster_columns = TRUE,
                      left_annotation = row_ha)  # correct annotation placement

# Save to SVG
svg(snakemake@output[["heatmap_plot"]])
draw(ht_colored)
dev.off()



# #Build Heatmap
# pdf(snakemake@output[["heatmap_plot"]])
# Heatmap(norm_counts_z,
#         show_row_names = FALSE,
#         show_row_dend = TRUE,
#         cluster_rows = TRUE,
#         cluster_columns = TRUE,
#         column_labels = colnames(norm_counts_z),
#         name = "Z-score",
#         row_labels = rownames(norm_counts_z),
#         )
# dev.off()
