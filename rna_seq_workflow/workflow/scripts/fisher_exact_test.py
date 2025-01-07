import pandas as pd
import numpy as np
import requests
from itertools import chain
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


### KEGG Pathway Enrichment Analysis ###


### Load the data ###

# Read in the functional annotations file
df = pd.read_csv(snakemake.input["functional_annotations"])

# Filter out the rows that have no KEGG annotation
df = df[['KEGG_ko', 'log2FoldChange']]
df['KEGG_ko'] = df['KEGG_ko'].str.split(',')

# Explode the KEGG_ko column
df_exploded = df.explode('KEGG_ko')
df_exploded['KEGG_ko'] = df_exploded['KEGG_ko'].str.strip().str.replace('ko:', '')
df_exploded = df_exploded.replace('-', np.nan)
df_exploded = df_exploded.dropna()

# Filter the data based on the log2FoldChange values
filter_mask_high = (df_exploded["log2FoldChange"] > 2)
filter_mask_background = ((df_exploded["log2FoldChange"] < 1) & (df_exploded["log2FoldChange"] > -1))
filter_mask_low = (df_exploded["log2FoldChange"] < -2)

# Apply the filters
df_exploded_positive = df_exploded[filter_mask_high]
df_exploded_background = df_exploded[filter_mask_background]
df_exploded_negative = df_exploded[filter_mask_low]


### Fetch the pathways for the KOs ###

# Function to fetch the pathways for a given KO
def fetch_pathways(ko):
    try:
        # Make a request to the KEGG API
        url = f"https://rest.kegg.jp/link/pathway/{ko}"
        response = requests.get(url)
        # Check if the response is successful
        if response.status_code == 200:
            if response.text.strip():
                pathways = []
                # Parse the response and extract the pathway IDs
                for line in response.text.splitlines():
                    pathway_list = line.split('\t')[1].strip().split(',')
                    ko_paths = [pathway for pathway in pathway_list if pathway.startswith('path:ko')]
                    # Append the pathway IDs to the list
                    if ko_paths:
                        pathways.append(ko_paths)
                # Flatten the list of pathways
                return list(chain.from_iterable(pathways))
            # If no pathways are found, return an empty list
            else:
                print("No pathways found for KO: ", ko)
                return []
        # If the response is not successful, print the error code    
        else:
            print("Error: ", response.status_code)
    # Handle exceptions
    except Exception as e:
        print(e)
    # Return an empty list if no pathways are found
    return []

# Apply the fetch_pathways function to the KEGG_ko column
df_exploded_positive.loc[:, "Pathway_IDs"] = df_exploded_positive["KEGG_ko"].apply(fetch_pathways)
df_exploded_background.loc[:, "Pathway_IDs"] = df_exploded_background["KEGG_ko"].apply(fetch_pathways)
df_exploded_negative.loc[:, "Pathway_IDs"] = df_exploded_negative["KEGG_ko"].apply(fetch_pathways)


### Calculate the enrichment of pathways in the positive, negative, and background subsets ###

# Function to calculate the enrichment of pathways in a subset compared to a background
def calculate_enrichment(subset, background):
    # Calculate the total number of genes in the subset and background
    total_subset = len(subset)
    total_background = len(background)
    # Count the occurrences of each pathway in the subset and background
    subset_pathway_counts = subset.explode("Pathway_IDs")["Pathway_IDs"].value_counts()
    background_pathway_counts = background.explode("Pathway_IDs")["Pathway_IDs"].value_counts()
    enrichment = []
    # Iterate over the pathway IDs in the subset
    for pathway_id, count in subset_pathway_counts.items():
        # Get the count of the pathway in the background
        background_count = background_pathway_counts.get(pathway_id, 0)
        # Calculate the contingency table for the Fisher's exact test
        non_pathway_count = total_subset - count
        non_background_count = total_background - background_count
        contingency_table = [[count, background_count], [non_pathway_count, non_background_count]]
        # Perform the Fisher's exact test
        oddsratio, pvalue = fisher_exact(contingency_table)
        # Append the results to the enrichment list
        enrichment.append((pathway_id, count, background_count, oddsratio, pvalue))
    
    # Create a DataFrame from the enrichment list
    enrichment_df = pd.DataFrame(enrichment, columns=["Pathway_ID", "Count", "Background_Count", "Odds_Ratio", "P_Value"])
    
    # Perform multiple testing correction on the p-values
    enrichment_df["Adjusted_P_Value"] = multipletests(enrichment_df["P_Value"], method="fdr_bh")[1]
        
    return enrichment_df

# Calculate the enrichment of pathways in the high, low, and background subsets
positive_enrichment = calculate_enrichment(df_exploded_positive, df_exploded_background)
negative_enrichment = calculate_enrichment(df_exploded_negative, df_exploded_background)


### Fetch the pathway names for the enriched pathways ###

pathway_cache = {}

# Function to fetch the pathway names for a given pathway ID
def fetch_pathway_names(pathway_id):
    pathway_names = []
    # Check if the pathway ID is already in the cache
    if pathway_id in pathway_cache:
        pathway_names.append(pathway_cache[pathway_id])
    else:
        try:
            # Make a request to the KEGG API
            url = f"https://rest.kegg.jp/get/{pathway_id}"
            response = requests.get(url)
            # Check if the response is successful
            if response.status_code == 200:
                # Parse the response and extract the pathway name
                for line in response.text.splitlines():
                    # Extract the name of the pathway
                    if line.startswith("NAME"):
                        name = line.replace("NAME", "").strip()
                        pathway_cache[pathway_id] = name
                        pathway_names.append(name)
                        break
            # If the response is not successful, print the error code        
            else:
                print("Error: ", response.status_code)
        # Handle exceptions
        except Exception as e:
            print(e)
    # Return a semicolon-separated list of pathway names            
    return ";".join(set(pathway_names))

# Apply the fetch_pathway_names function to the Pathway_IDs column
positive_enrichment.loc[:, "Pathway_Names"] = positive_enrichment["Pathway_ID"].apply(fetch_pathway_names)
negative_enrichment.loc[:, "Pathway_Names"] = negative_enrichment["Pathway_ID"].apply(fetch_pathway_names)

# Save the enriched pathways to CSV files
positive_enrichment.to_csv(snakemake.output["positive_enrichment"], positive_enrichment)
negative_enrichment.to_csv(snakemake.output["negative_enrichment"], negative_enrichment)


### Plot the enriched pathways ###

def plot_stacked_bar_counts(positive_enrichment, negative_enrichment, top_n, title):
    """
    Create a stacked horizontal bar plot to visualize the top enriched pathways 
    from high and low enrichment datasets, showing both pathway counts and background counts.
    The pathways are sorted by the ratio (Count / Background_Count) in descending order 
    within each group, with a visual separation between groups.

    """
    
    # Filter pathways with significant enrichment (FDR < 0.05)
    significant_high = positive_enrichment[positive_enrichment["Adjusted_P_Value"] < 0.01]
    significant_low = negative_enrichment[negative_enrichment["Adjusted_P_Value"] < 0.01]
    
    # Compute ratio and select top pathways based on the ratio (Count / Background_Count)
    top_high = (
        significant_high
        .assign(Ratio=significant_high["Count"] / significant_high["Background_Count"])
        .nlargest(top_n, "Ratio")
        .assign(Color="skyblue", Group="MOT vs OT")
    )
    top_low = (
        significant_low
        .assign(Ratio=significant_low["Count"] / significant_low["Background_Count"])
        .nlargest(top_n, "Ratio")
        .assign(Color="salmon", Group="OT vs MOT")
    )
    
    # Combine data into a single DataFrame but keep groups visually separated
    top_high = top_high.sort_values("Ratio", ascending=False)
    top_low = top_low.sort_values("Ratio", ascending=False)
    
    # Add a gap row for visual separation
    gap_row = pd.DataFrame({
        "Pathway_Names": [" "], 
        "Count": [0], 
        "Background_Count": [0], 
        "Color": "white", 
        "Group": "Gap", 
        "Ratio": [0]
    })
    
    combined_data = pd.concat([top_high, gap_row, top_low], ignore_index=True)
    
    # Prepare data for plotting
    pathway_names = combined_data["Pathway_Names"]
    counts = combined_data["Count"]
    background_counts = combined_data["Background_Count"]
    colors = combined_data["Color"]

    # Plot settings
    plt.figure(figsize=(10, 8))
    bar_positions = range(len(pathway_names))

    # Plot the background counts
    plt.barh(bar_positions, background_counts, color="lightgray", label="Background Count")

    # Overlay the counts for high and low enrichment
    for i, (count, bg_count, color) in enumerate(zip(counts, background_counts, colors)):
        plt.barh(i, count, color=color, left=bg_count)

    # Custom legend
    legend_patches = [
        mpatches.Patch(color="lightgray", label="Background"),
        mpatches.Patch(color="skyblue", label="MOT vs OT"),
        mpatches.Patch(color="salmon", label="OT vs MOT")
    ]
    plt.legend(handles=legend_patches, loc="best")

    # Add labels and formatting
    plt.yticks(bar_positions, pathway_names)
    plt.xlabel("Counts")
    plt.ylabel("Pathways")
    plt.title(title)
    plt.tight_layout()
    plt.gca().invert_yaxis()

    # Save the plot
    plt.savefig(snakemake.output["stacked_bar_plot"], dpi=300, bbox_inches="tight")


# Call the function to plot
plot_stacked_bar_counts(positive_enrichment, negative_enrichment, 10, "Enriched KEGG Pathways (FDR < 0.01)")
