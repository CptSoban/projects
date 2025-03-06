import pandas as pd
import numpy as np
import gseapy as gp
import matplotlib.pyplot as plt
import seaborn as sns
import json

# Added function to filter pathways
def filter_pathways(df, exclude_terms):
    """Remove pathways that contain any of the exclude_terms."""
    pattern = '|'.join(exclude_terms)  # Create regex pattern from list
    return df[~df['Term'].str.contains(pattern, case=False, na=False)]

### Load and preprocess the data ###
df = pd.read_csv(snakemake.input["functional_annotations"])

# Filter out non-significant genes
df = df[['KEGG_ko', 'log2FoldChange']]
df['KEGG_ko'] = df['KEGG_ko'].str.split(',')
df_exploded = df.explode('KEGG_ko')
df_exploded['KEGG_ko'] = df_exploded['KEGG_ko'].str.strip()
df_exploded = df_exploded.replace('-', np.nan).dropna()

# Load KEGG pathway cache
with open(snakemake.input["kegg_cache"], "r") as cache_file:
    pathway_cache = json.load(cache_file)

# Map KEGG KOs to pathway names
df_exploded['Pathway_Names'] = df_exploded['KEGG_ko'].map(pathway_cache)
df_exploded = df_exploded.explode("Pathway_Names").dropna()

# Log the number of KOs and pathways recovered
total_kos = df_exploded["KEGG_ko"].count()
total_pathways = df_exploded["Pathway_Names"].count()
with open(snakemake.log["pathway_recovery"], "w") as log_file:
    log_file.write(f"Total number of KOs:   {total_kos}\nRecovered pathways: {total_pathways}")

### Prepare data for GSEA ###

# Create ranked list
ranked_list = (
    df_exploded.groupby("KEGG_ko")["log2FoldChange"]
    .mean()
    .sort_values(ascending=False)
)
ranked_list_dict = dict(zip(ranked_list.index, ranked_list.values))
df_rnk = pd.DataFrame(ranked_list_dict.items())

# Create gene sets
gene_sets = (
    df_exploded.set_index('Pathway_Names')['KEGG_ko']
    .groupby("Pathway_Names")
    .apply(list)
    .to_dict()
)

### Perform GSEA ###

# Run GSEA
enrichment_results = gp.prerank(
    rnk=ranked_list,
    gene_sets=gene_sets,
    outdir=None,
    permutation_num=10000,
    min_size=10,
    max_size=250,
)

# Save GSEA results
gsea_results_df = pd.DataFrame(enrichment_results.res2d)
gsea_results_df['FDR q-val'] = pd.to_numeric(gsea_results_df['FDR q-val'], errors='coerce')
gsea_results_df = gsea_results_df.dropna(subset=['FDR q-val'])

gsea_results_df.to_csv(snakemake.output["gsea_results"], index=False)

gsea_results_df = gsea_results_df[gsea_results_df["FDR q-val"] < 0.05]

# Filtering of unwanted pathways
exclude_terms = snakemake.params["filter_terms"]  #terms to exclude
gsea_results_df = filter_pathways(gsea_results_df, exclude_terms)

top_pathways = gsea_results_df[["Term", "NES"]].sort_values("NES", ascending=False)

# Plot the top pathways
plt.figure(figsize=(10, 6))
sns.barplot(
    data=top_pathways,
    x="NES",
    y="Term",
    color="red",
)

plt.xlabel("Normalized Enrichment Score (NES)")
plt.ylabel("Pathway")
plt.title("Enriched KEGG pathways in GSEA (FDR < 0.05)")
plt.tight_layout()

plt.savefig(snakemake.output["gsea_plot"], dpi=300)

