import pandas as pd
import numpy as np
import gseapy as gp
import matplotlib.pyplot as plt
import seaborn as sns
import json

### Load and preprocess the data ###
df = pd.read_csv(snakemake.input["functional_annotations"])

# Prepare the KEGG KO column
df = df[['KEGG_ko', 'log2FoldChange']]
df['KEGG_ko'] = df['KEGG_ko'].str.split(',')
df_exploded = df.explode('KEGG_ko')
df_exploded['KEGG_ko'] = df_exploded['KEGG_ko'].str.strip()
df_exploded = df_exploded.replace('-', np.nan).dropna()

# Load the cached KEGG pathway data
with open(snakemake.input["kegg_cache"], "r") as cache_file:
    pathway_cache = json.load(cache_file)
    
# Map KO IDs to pathways using the cache
df_exploded['Pathway_Names'] = df_exploded['KEGG_ko'].map(pathway_cache)
df_exploded = df_exploded.explode("Pathway_Names").dropna()

# Count the total number of KOs and recovered pathway names
total_kos = df_exploded["KEGG_ko"].count()
total_pathways = df_exploded["Pathway_Names"].count()
with open(snakemake.log["pathway_recovery"], "w") as log_file:
    log_file.write(f"Total number of KOs:   {total_kos}\nRecovered pathways: {total_pathways}")


### Prepare data for GSEA ###

# Create a ranked list of genes (ranked by log2FoldChange)
ranked_list = (
    df_exploded.groupby("KEGG_ko")["log2FoldChange"]
    .mean()
    .sort_values(ascending=False)
)
ranked_list = ranked_list.reset_index()
ranked_list_dict = dict(zip(ranked_list["KEGG_ko"], ranked_list["log2FoldChange"]))
df_rnk = pd.DataFrame(ranked_list_dict.items())

# Prepare a dictionary of gene sets (pathway -> KOs)
gene_sets = (
    df_exploded.set_index('Pathway_Names')['KEGG_ko']
    .groupby("Pathway_Names")
    .apply(list)
    .to_dict()
)

### Perform GSEA ###
# Run GSEA using gseapy
enrichment_results = gp.prerank(
    rnk=ranked_list,
    gene_sets=gene_sets,
    outdir=None,  # Set to a folder if you want to save results
    permutation_num=10000,  # Increase for more robust results
    min_size=3,  # Minimum size of gene sets to include in analysis
    max_size=500,  # Maximum size of gene sets
)

### Process and save results ###
# Convert GSEA results to a DataFrame
gsea_results_df = pd.DataFrame(enrichment_results.res2d)
gsea_results_df['FDR q-val'] = pd.to_numeric(gsea_results_df['FDR q-val'], errors='coerce')
gsea_results_df = gsea_results_df.dropna(subset=['FDR q-val'])

# Save GSEA results to CSV
gsea_results_df.to_csv(snakemake.output["gsea_results"], index=False)

# Filter the results by FDR q-value
gsea_results_df = gsea_results_df[gsea_results_df["FDR q-val"] < 0.1]

# Select the significant pathways and sort by NES
top_pathways = gsea_results_df[["Term", "NES"]].sort_values("NES", ascending=False)

# Create a bar plot
plt.figure(figsize=(10, 6))
sns.barplot(
    data=top_pathways,
    x="NES",  # Normalized Enrichment Score
    y="Term",
    color="red",
)

# Add labels and title
plt.xlabel("Normalized Enrichment Score (NES)")
plt.ylabel("Pathway")
plt.title("Enriched KEGG pathways in GSEA (FDR < 0.05)")
plt.tight_layout()

# Save the plot
plt.savefig(snakemake.output["gsea_plot"], dpi=300)
