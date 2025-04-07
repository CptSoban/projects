import pandas as pd
import numpy as np
import gseapy as gp
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
import json

# Function to filter pathways/modules
def filter_terms(df, exclude_terms):
    """Remove pathways/modules that contain any of the exclude_terms."""
    pattern = '|'.join(exclude_terms)  # Create regex pattern
    return df[~df['Term'].str.contains(pattern, case=False, na=False)]

### Load and preprocess the data ###
df = pd.read_csv(snakemake.input["functional_annotations"])

# Extract KEGG KOs and log2FoldChange
df = df[['KEGG_ko', 'log2FoldChange']]
df['KEGG_ko'] = df['KEGG_ko'].str.split(',')
df_exploded = df.explode('KEGG_ko')
df_exploded['KEGG_ko'] = df_exploded['KEGG_ko'].str.strip()
df_exploded = df_exploded.replace('-', np.nan).dropna()

# Load KEGG pathway & module cache
with open(snakemake.input["kegg_cache"], "r") as cache_file:
    kegg_cache = json.load(cache_file)

# Extract pathways and modules separately
df_exploded["Pathways"] = df_exploded["KEGG_ko"].map(lambda ko: kegg_cache.get(ko, {}).get("pathways", []))
df_exploded["Modules"] = df_exploded["KEGG_ko"].map(lambda ko: kegg_cache.get(ko, {}).get("modules", []))

# Explode pathways and modules into separate rows
df_pathways = df_exploded.explode("Pathways").dropna(subset=["Pathways"])
df_modules = df_exploded.explode("Modules").dropna(subset=["Modules"])

# Log the number of KOs, pathways, and modules recovered
with open(snakemake.log[0], "w") as log_file:
    log_file.write(f"Total number of KOs: {df_exploded['KEGG_ko'].nunique()}\n")
    log_file.write(f"Recovered pathways: {df_pathways['Pathways'].nunique()}\n")
    log_file.write(f"Recovered modules: {df_modules['Modules'].nunique()}\n")


### Prepare data for GSEA ###
def prepare_gsea_data(df, term_column):
    """Prepare ranked list and gene sets for GSEA analysis."""
    ranked_list = (
        df.groupby("KEGG_ko")["log2FoldChange"]
        .mean()
        .sort_values(ascending=False)
    )
    gene_sets = (
        df.set_index(term_column)['KEGG_ko']
        .groupby(term_column)
        .apply(list)
        .to_dict()
    )
    return ranked_list, gene_sets

# Prepare GSEA input for pathways and modules
ranked_list_pathways, gene_sets_pathways = prepare_gsea_data(df_pathways, "Pathways")
ranked_list_modules, gene_sets_modules = prepare_gsea_data(df_modules, "Modules")

### Perform GSEA ###
def run_gsea(ranked_list, gene_sets, term_type):
    """Run GSEA analysis and return filtered results."""
    results = gp.prerank(
        rnk=ranked_list,
        gene_sets=gene_sets,
        outdir=None,
        permutation_num=10000,
        min_size=10,
        max_size=250,
    )
    
    results_df = pd.DataFrame(results.res2d)
    results_df['FDR q-val'] = pd.to_numeric(results_df['FDR q-val'], errors='coerce')
    results_df = results_df.dropna(subset=['FDR q-val'])
    
    # Filter results by FDR threshold
    results_df = results_df[results_df["FDR q-val"] < 0.1]
    
    # Apply filtering for unwanted terms
    exclude_terms = snakemake.params["filter_terms"]
    results_df = filter_terms(results_df, exclude_terms)
    
    # Save results
    results_df.to_csv(snakemake.output[f"gsea_results_{term_type}"], index=False)
    

    
    return results_df

# Run GSEA for pathways and modules
gsea_results_pathways = run_gsea(ranked_list_pathways, gene_sets_pathways, "pathways")
gsea_results_modules = run_gsea(ranked_list_modules, gene_sets_modules, "modules")

### Visualization ###
def plot_gsea_results(results_df, term_type):
    """Plot the top enriched pathways/modules with color based on FDR q-value."""
    top_terms = results_df[["Term", "NES", "FDR q-val"]].sort_values("NES", ascending=False)

    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Normalize FDR q-value for color mapping
    norm = plt.Normalize(top_terms["FDR q-val"].min(), top_terms["FDR q-val"].max())
    
    # Create a colormap from blue to purple to red
    blue_purple_red = LinearSegmentedColormap.from_list("red_purple_blue", ["red", "purple", "blue"])
    cmap = blue_purple_red
    #cmap = plt.cm.plasma 
    bar_colors = cmap(norm(top_terms["FDR q-val"].values))

    # Create barplot
    sns.barplot(data=top_terms, x="NES", y="Term", palette=bar_colors, ax=ax)

    ax.set_xlabel("Normalized Enrichment Score (NES)")
    ax.set_ylabel(term_type.capitalize())
    ax.set_title(f"Enriched KEGG {term_type} in GSEA (FDR < 0.1)")

    # Create colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    # Shrink colorbar and reverse direction (high FDR at bottom)
    cbar = fig.colorbar(sm, ax=ax, aspect=30)
    cbar.ax.invert_yaxis()  # invert the colorbar to show low FDR at the top
    cbar.set_label("FDR q-value")

    plt.tight_layout()
    plt.savefig(snakemake.output[f"gsea_plot_{term_type}"], dpi=300)
    plt.close()



# Plot results
plot_gsea_results(gsea_results_pathways, "pathways")
plot_gsea_results(gsea_results_modules, "modules")