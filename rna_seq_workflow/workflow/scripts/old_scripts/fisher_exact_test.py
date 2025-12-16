import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import json


### KEGG Pathway Enrichment Analysis ###

### Load the data ###

# Read in the functional annotations file
df = pd.read_csv(snakemake.input["functional_annotations"])

# Filter out the rows that have no KEGG annotation
df = df[['KEGG_ko', 'log2FoldChange']]
df['KEGG_ko'] = df['KEGG_ko'].str.split(',')

# Explode the KEGG_ko column
df_exploded = df.explode('KEGG_ko')
df_exploded['KEGG_ko'] = df_exploded['KEGG_ko'].str.strip()
df_exploded = df_exploded.replace('-', np.nan)
df_exploded = df_exploded.dropna()

# Split the data into positive (upregulated) and background (non-upregulated) subsets
df_exploded_positive = df_exploded[df_exploded["log2FoldChange"] > 1.5]
df_exploded_background = df_exploded[df_exploded["log2FoldChange"] < 1.5]


### Fetch the pathways for the KOs ###

# Load the cached KEGG pathway data
with open(snakemake.input["kegg_cache"], "r") as cache_file:
    pathway_cache = json.load(cache_file)
    
# Map KO IDs to pathways using the cache
df_exploded_positive['Pathway_Names'] = df_exploded_positive['KEGG_ko'].map(pathway_cache)
df_exploded_background['Pathway_Names'] = df_exploded_background['KEGG_ko'].map(pathway_cache)


### Calculate the enrichment of pathways in the positive, negative, and background subsets ###

# Function to calculate the enrichment of pathways in a subset compared to a background
def calculate_enrichment(subset, background):
    # Calculate the total number of genes in the subset and background
    total_subset = len(subset)
    total_background = len(background)
    # Count the occurrences of each pathway in the subset and background
    subset_pathway_counts = subset["Pathway_Names"].explode().value_counts()
    background_pathway_counts = background["Pathway_Names"].explode().value_counts()
    pseudo_count = 0.5
    # Add pseudo counts to avoid division by zero
    subset_pathway_counts += pseudo_count
    background_pathway_counts += pseudo_count
    
    # Initialize a list to store the enrichment results
    enrichment = []
    # Iterate over the pathway IDs in the subset
    for pathway_name, count in subset_pathway_counts.items():
        # Get the count of the pathway in the background
        background_count = background_pathway_counts.get(pathway_name, 0)
        # Calculate the contingency table for the Fisher's exact test
        non_pathway_count = total_subset - count
        non_background_count = total_background - background_count
        contingency_table = [[count, background_count], [non_pathway_count, non_background_count]]
        # Perform the Fisher's exact test
        oddsratio, pvalue = fisher_exact(contingency_table)
        # Append the results to the enrichment list
        enrichment.append((pathway_name, count, background_count, oddsratio, pvalue))
    
    # Create a DataFrame from the enrichment list
    enrichment_df = pd.DataFrame(enrichment, columns=["Pathway_Name", "Count", "Background_Count", "Odds_Ratio", "P_Value"])
    
    # Perform multiple testing correction on the p-values
    enrichment_df["Adjusted_P_Value"] = multipletests(enrichment_df["P_Value"], method="fdr_bh")[1]
    
    # Replace infinite odds ratios with a high value for visualization
    enrichment_df["Odds_Ratio"] = enrichment_df["Odds_Ratio"].replace(np.inf, 1000)

        
    return enrichment_df

# Calculate the enrichment of pathways in the high, low, and background subsets
positive_enrichment = calculate_enrichment(df_exploded_positive, df_exploded_background)

# Save the enriched pathways to CSV files
positive_enrichment.to_csv(snakemake.output["positive_enrichment"])


# Create a forest plot for the enriched pathways
def create_forest_plot(enrichment_df, output_file):
    # Sort the DataFrame by Adjusted P-Value (smallest first)
    enrichment_df = enrichment_df.sort_values(by="Adjusted_P_Value")
    
    # Keep only the significant pathways (Adjusted P-Value < 0.05)
    significant_enrichment = enrichment_df[enrichment_df["Adjusted_P_Value"] < 0.05]
    
    # If there are no significant pathways, exit
    if significant_enrichment.empty:
        print("No significant pathways to plot.")
        return
    
    # Extract relevant data
    pathway_names = significant_enrichment["Pathway_Name"]
    odds_ratios = significant_enrichment["Odds_Ratio"]
    
    # Calculate confidence intervals for Odds Ratios
    ci_lower = []
    ci_upper = []
    pseudo_count = 0.5
    for _, row in significant_enrichment.iterrows():
        a = row["Count"]
        b = row["Background_Count"] + pseudo_count
        c = df_exploded_positive.shape[0] - a # Not in pathway (positive subset)
        d = df_exploded_background.shape[0] - b # Not in pathway (background)
        se_ln_or = np.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
        z = 1.96  # 95% CI
        ln_or = np.log(row["Odds_Ratio"])
        ci_lower.append(np.exp(ln_or - z * se_ln_or))
        ci_upper.append(np.exp(ln_or + z * se_ln_or))
    
    # Add confidence intervals to the DataFrame
    significant_enrichment["CI_Lower"] = ci_lower
    significant_enrichment["CI_Upper"] = ci_upper
    
    # Plot the forest plot
    plt.figure(figsize=(8, len(pathway_names) * 0.5))
    plt.errorbar(
        odds_ratios, 
        range(len(pathway_names)), 
        xerr=[np.array(odds_ratios) - np.array(ci_lower), np.array(ci_upper) - np.array(odds_ratios)],
        fmt='o', 
        color='blue', 
        capsize=5, 
        label="95% CI"
    )
    plt.axvline(1, color='red', linestyle='--', label='No Enrichment (OR=1)')
    plt.yticks(range(len(pathway_names)), pathway_names)
    plt.xlabel("Odds Ratio (log scale)")
    plt.xscale("log")  # Log scale for Odds Ratios
    plt.title("Forest Plot of Enriched Pathways")
    plt.tight_layout()
    plt.legend()
    
    # Save the plot to the output file
    plt.savefig(output_file)
    plt.close()

# Create the forest plot
create_forest_plot(positive_enrichment, snakemake.output["forest_plot"])



# ### Plot the enriched pathways ###

# def plot_stacked_bar_counts(positive_enrichment, negative_enrichment, top_n, title):
#     """
#     Create a stacked horizontal bar plot to visualize the top enriched pathways 
#     from positive and negative enrichment datasets, showing both pathway counts and background counts.
#     The pathways are sorted by the ratio (Count / Background_Count) in descending order 
#     within each group, with a visual separation between groups.

#     """
    
#     # Filter pathways with significant enrichment (FDR < 0.05)
#     significant_high = positive_enrichment[positive_enrichment["Adjusted_P_Value"] < 0.01]
#     significant_low = negative_enrichment[negative_enrichment["Adjusted_P_Value"] < 0.01]
    
#     # Compute ratio and select top pathways based on the ratio (Count / Background_Count)
#     top_high = (
#         significant_high
#         .assign(Ratio=significant_high["Count"] / significant_high["Background_Count"])
#         .nlargest(top_n, "Ratio")
#         .assign(Color="skyblue", Group="Positive Enrichment")
#     )
#     top_low = (
#         significant_low
#         .assign(Ratio=significant_low["Count"] / significant_low["Background_Count"])
#         .nlargest(top_n, "Ratio")
#         .assign(Color="salmon", Group="Negative Enrichment")
#     )
    
#     # Combine data into a single DataFrame but keep groups visually separated
#     top_high = top_high.sort_values("Ratio", ascending=False)
#     top_low = top_low.sort_values("Ratio", ascending=False)
    
#     # Add a gap row for visual separation
#     gap_row = pd.DataFrame({
#         "Pathway_Names": [" "], 
#         "Count": [0], 
#         "Background_Count": [0], 
#         "Color": "white", 
#         "Group": "Gap", 
#         "Ratio": [0]
#     })
    
#     combined_data = pd.concat([top_high, gap_row, top_low], ignore_index=True)
    
#     # Prepare data for plotting
#     pathway_names = combined_data["Pathway_Names"]
#     counts = combined_data["Count"]
#     background_counts = combined_data["Background_Count"]
#     colors = combined_data["Color"]

#     # Plot settings
#     plt.figure(figsize=(10, 8))
#     bar_positions = range(len(pathway_names))

#     # Plot the background counts
#     plt.barh(bar_positions, background_counts, color="lightgray", label="Background Count")

#     # Overlay the counts for high and low enrichment
#     for i, (count, bg_count, color) in enumerate(zip(counts, background_counts, colors)):
#         plt.barh(i, count, color=color, left=bg_count)

#     # Custom legend
#     legend_patches = [
#         mpatches.Patch(color="lightgray", label="Background"),
#         mpatches.Patch(color="skyblue", label="Positive Enrichment"),
#         mpatches.Patch(color="salmon", label="Negative Enrichment"),
#     ]
#     plt.legend(handles=legend_patches, loc="best")

#     # Add labels and formatting
#     plt.yticks(bar_positions, pathway_names)
#     plt.xlabel("Counts")
#     plt.ylabel("Pathways")
#     plt.title(title)
#     plt.tight_layout()
#     plt.gca().invert_yaxis()

#     # Save the plot
#     plt.savefig(snakemake.output["stacked_bar_plot"], dpi=300, bbox_inches="tight")


# # Call the function to plot
# plot_stacked_bar_counts(positive_enrichment, negative_enrichment, 10, "Enriched KEGG Pathways (FDR < 0.01)")
