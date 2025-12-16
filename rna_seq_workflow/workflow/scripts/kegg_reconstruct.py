import pandas as pd
import numpy as np

# Read the combined annotations file
combined_annotations = pd.read_csv(snakemake.input["combined_annotations"])
trim_df = combined_annotations[["symbol","log2FoldChange","KEGG_ko"]].set_index("symbol")

### Reconstruct the DataFrame with KEGG IDs

# Remove rows where 'KEGG_ko' is NaN or '-'
reconstruct_df = trim_df[trim_df['KEGG_ko'] != '-'].dropna()

# Convert 'KEGG_ko' column to a list of KEGG IDs
reconstruct_df["KEGG_ko"] = reconstruct_df["KEGG_ko"].str.split(',')

# Explode the 'KEGG_ko' column to have one row per KEGG ID
reconstruct_df =reconstruct_df.explode("KEGG_ko")

# Write the reconstructed DataFrame to a CSV file
reconstruct_df["KEGG_ko"].to_csv(snakemake.output["annotate_id_with_ko"], sep="\t")


### Create DataFrame for color mapping

# Calculate the mean log2FoldChange for each KEGG ID
color_df = reconstruct_df.groupby('KEGG_ko', as_index=False).mean()

# Rename the columns for clarity
color_df = color_df.sort_values(by=["log2FoldChange"], ascending=False)

# Create a new DataFrame with KEGG IDs and their corresponding log2FoldChange
color_df = color_df.set_index("KEGG_ko")

# Create a new column for color mapping based on log2FoldChange
color_df.to_csv(snakemake.output["annotate_ko_with_log2fc"], sep="\t", header=None)

