# Load necessary libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(snakemake)

# Input and output files from Snakemake
input_file <- snakemake@input[[1]]
output_file <- snakemake@output[[1]]

# Read the input gene list
gene_list <- read.csv(input_file, header = TRUE, stringsAsFactors = FALSE)
ko_ids <- gene_list$ko_id

# Perform KEGG pathway analysis using KO ids
kegg_result <- enrichKEGG(gene = ko_ids, organism = 'fox', keyType = 'kegg')

# Save the results to a file
write.csv(as.data.frame(kegg_result), file = output_file)
# Extract log2 fold change values
log2fc <- gene_list$log2FoldChange

# Create a named vector of log2 fold changes with KO ids as names
log2fc_named <- setNames(log2fc, ko_ids)

# Perform KEGG pathway analysis using KO ids and log2 fold changes
kegg_result <- gseKEGG(geneList = log2fc_named, organism = 'fox', keyType = 'kegg')

# Save the results to a file
write.csv(as.data.frame(kegg_result), file = output_file)