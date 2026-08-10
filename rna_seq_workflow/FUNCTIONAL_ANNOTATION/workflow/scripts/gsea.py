import pandas as pd
import numpy as np
import gseapy as gp
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
import json

### Prepare data for GSEA ###
def prepare_gsea_data(df, term_column, baseMean_floor=1.0):
    """Prepare ranked list and gene sets for GSEA analysis."""
    # Create a copy to avoid modifying the original dataframe
    df = df.copy()
    
    # Handle baseMean values to avoid issues with zero or near-zero values
    df["baseMean_safe"] = df["baseMean"].clip(lower=baseMean_floor)

    df["weight"] = df["baseMean_safe"]

    # Calculate weighted log2FoldChange
    weighted_sum = (
        df.assign(weighted_lfc=df["log2FoldChange"] * df["weight"])
          .groupby("KEGG_ko")[["weighted_lfc", "weight"]]
          .sum()
    )
    
    # Create ranked list 
    ranked_list = (
        weighted_sum["weighted_lfc"] / weighted_sum["weight"]
    ).sort_values(ascending=False)
    
    # Create gene sets mapping terms to KOs
    gene_sets = (
        df.set_index(term_column)['KEGG_ko']
        .groupby(term_column)
        .apply(list)
        .to_dict()
    )
    
    return ranked_list, gene_sets

# Function to filter pathways/modules
def filter_terms(df, exclude_terms):
    """Filter pathways or modules based on exclude terms."""
    pattern = '|'.join(exclude_terms)  # Create regex pattern
    return df[~df['Term'].str.contains(pattern, case=False, na=False)]

### Perform GSEA ###
def run_gsea(ranked_list, gene_sets, term_type):
    """Run GSEA analysis and return filtered results."""
    results = gp.prerank(
        rnk=ranked_list,
        gene_sets=gene_sets,
        outdir=None,
        permutation_num=10000,
        min_size=10,
        max_size=200,
    )
    
    results_df = pd.DataFrame(results.res2d)
    results_df['FDR q-val'] = pd.to_numeric(results_df['FDR q-val'], errors='coerce')
    results_df = results_df.dropna(subset=['FDR q-val'])
    
    # Filter results by FDR threshold
    results_df = results_df[results_df["FDR q-val"] <= snakemake.params["fdr_threshold"]]
    # Apply filtering for unwanted terms
    exclude_terms = snakemake.params["filter_terms"]
    results_df = filter_terms(results_df, exclude_terms)
    # Save results
    results_df.to_csv(snakemake.output[f"gsea_results_{term_type}"], index=False)
    
    
    return results_df

### Visualization ###
def plot_gsea_results(results_df, term_type):
    """Plot the top enriched pathways/modules with color based on FDR q-value."""
    top_terms = results_df[["Term", "NES", "FDR q-val"]].sort_values("NES", ascending=False)
    num_bars = len(top_terms)

    # --- Fixed axes dimensions (in inches) ---
    BAR_HEIGHT   = 0.4   # height per bar
    AX_WIDTH     = 6.0   # fixed axes box width
    AX_HEIGHT    = BAR_HEIGHT * num_bars
    LEFT_PAD     = 2.5   # room for y-tick labels
    RIGHT_PAD    = 1.2   # room for colorbar
    TOP_PAD      = 0.3
    BOTTOM_PAD   = 0.7   # room for x-axis label

    fig_width  = LEFT_PAD + AX_WIDTH + RIGHT_PAD
    fig_height = TOP_PAD + AX_HEIGHT + BOTTOM_PAD

    fig = plt.figure(figsize=(fig_width, fig_height))

    # Place axes at exact position: [left, bottom, width, height] as fractions
    ax = fig.add_axes([
        LEFT_PAD   / fig_width,
        BOTTOM_PAD / fig_height,
        AX_WIDTH   / fig_width,
        AX_HEIGHT  / fig_height,
    ])

    # Normalize FDR q-value for color mapping
    fdr_threshold = snakemake.params["fdr_threshold"]
    norm = plt.Normalize(vmin=0, vmax=fdr_threshold)

    blue_purple_red = LinearSegmentedColormap.from_list(
        "red_purple_blue", ["red", "purple", "blue"]
    )
    cmap = blue_purple_red
    bar_colors = cmap(norm(top_terms["FDR q-val"].values))

    sns.barplot(data=top_terms, x="NES", y="Term", palette=bar_colors, ax=ax)
    
    # --- Fixed x-axis range and ticks ---
    X_MIN   = -3.0   # adjust to your expected NES range
    X_MAX   =  3.0
    X_TICKS = [-3, -2, -1, 0, 1, 2, 3]

    ax.set_xlim(X_MIN, X_MAX)
    ax.set_xticks(X_TICKS)

    ax.set_xlabel("Normalized Enrichment Score (NES)", fontsize=16)
    ax.set_ylabel(term_type.capitalize(), fontsize=16)
    ax.tick_params(axis="both", which="major", labelsize=12)

    for spine in ax.spines.values():
        spine.set_linewidth(1.2)

    ax.set_ylim(-0.5, num_bars - 0.5)

    # Colorbar: manually placed to the right of the axes box
    cbar_left   = (LEFT_PAD + AX_WIDTH + 0.15) / fig_width
    cbar_bottom = (BOTTOM_PAD + AX_HEIGHT * 0.3) / fig_height
    cbar_width  = 0.08 / fig_width
    cbar_height = (AX_HEIGHT * 0.4) / fig_height

    cbar_ax = fig.add_axes([cbar_left, cbar_bottom, cbar_width, cbar_height])
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cbar_ax)
    cbar.ax.invert_yaxis()
    cbar.set_label("FDR q-value", fontsize=14)
    cbar.ax.tick_params(labelsize=12)

    # No tight_layout — axes are already manually positioned
    plt.savefig(
        snakemake.output[f"gsea_plot_{term_type}"],
        dpi=300,
        transparent=True,
        bbox_inches="tight",   # catches any label overflow without resizing the axes box
    )
    plt.close()

### Load and preprocess the data ###
df = pd.read_csv(snakemake.input["functional_annotations"])

# Extract KEGG KOs and log2FoldChange
df = df[['KEGG_ko', 'log2FoldChange','baseMean','padj']]
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

# Prepare GSEA input for pathways and modules
ranked_list_pathways, gene_sets_pathways = prepare_gsea_data(df_pathways, "Pathways")
ranked_list_modules, gene_sets_modules = prepare_gsea_data(df_modules, "Modules")

# Run GSEA for pathways and modules
gsea_results_pathways = run_gsea(ranked_list_pathways, gene_sets_pathways, "pathways")
gsea_results_modules = run_gsea(ranked_list_modules, gene_sets_modules, "modules")

# Plot results
plot_gsea_results(gsea_results_pathways, "pathways")
plot_gsea_results(gsea_results_modules, "modules")