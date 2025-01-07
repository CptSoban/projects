import pandas as pd
import numpy as np
import requests
from itertools import chain
import gseapy as gp
import matplotlib.pyplot as plt
import seaborn as sns

### Load and preprocess the data ###
df = pd.read_csv(snakemake.input["functional_annotations"])

# Prepare the KEGG KO column
df = df[['KEGG_ko', 'log2FoldChange']]
df['KEGG_ko'] = df['KEGG_ko'].str.split(',')
df_exploded = df.explode('KEGG_ko')
df_exploded['KEGG_ko'] = df_exploded['KEGG_ko'].str.strip().str.replace('ko:', '')
df_exploded = df_exploded.replace('-', np.nan).dropna()

# Fetch pathway mappings for each KEGG KO
def fetch_pathways(ko):
    try:
        url = f"https://rest.kegg.jp/link/pathway/{ko}"
        response = requests.get(url)
        if response.status_code == 200:
            return [
                line.split('\t')[1].replace("path:", "")
                for line in response.text.splitlines()
                if line.startswith("path:ko")
            ]
    except Exception as e:
        print(f"Error fetching pathways for KO {ko}: {e}")
    return []

df_exploded["Pathway_IDs"] = df_exploded["KEGG_ko"].apply(fetch_pathways)
df_exploded = df_exploded.explode("Pathway_IDs").dropna()

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
df_exploded["Pathway_Name"] = df_exploded["Pathway_IDs"].apply(fetch_pathway_names)
df_exploded = df_exploded.drop("Pathway_IDs", axis=1)

### Prepare data for GSEA ###
# Create a ranked list of genes (ranked by log2FoldChange)
ranked_list = (
    df_exploded.groupby("KEGG_ko")["log2FoldChange"]
    .mean()
    .sort_values(ascending=False)
)
ranked_list = ranked_list.reset_index()
ranked_list.columns = ["KEGG_ko", "Rank"]
ranked_list_dict = dict(zip(ranked_list["KEGG_ko"], ranked_list["Rank"]))

# Prepare a dictionary of gene sets (pathway -> KOs)
gene_sets = (
    df_exploded.groupby("Pathway_Name")["KEGG_ko"]
    .apply(list)
    .to_dict()
)

### Perform GSEA ###
# Run GSEA using gseapy
enrichment_results = gp.prerank(
    rnk=ranked_list,
    gene_sets=gene_sets,
    outdir=None,  # Set to a folder if you want to save results
    permutation_num=1000,  # Increase for more robust results
    min_size=5,  # Minimum size of gene sets to include in analysis
    max_size=500,  # Maximum size of gene sets
)

### Process and save results ###
# Convert GSEA results to a DataFrame
gsea_results_df = pd.DataFrame(enrichment_results.res2d)
gsea_results_df.reset_index(inplace=True)
gsea_results_df.rename(columns={"index": "Pathway_Name"}, inplace=True)

# Save GSEA results to CSV
gsea_results_df.to_csv(snakemake.output["gsea_results"], index=False)


# Select the top enriched pathways
top_pathways = gsea_results_df.nsmallest(10, "FDR.q-val")  # Top 10 pathways by FDR

# Create a bar plot
plt.figure(figsize=(10, 6))
sns.barplot(
    data=top_pathways,
    x="NES",  # Normalized Enrichment Score
    y="Pathway_Name",
    hue="FDR.q-val",  # Use FDR as a color gradient
    palette="viridis",
)

# Add labels and title
plt.xlabel("Normalized Enrichment Score (NES)")
plt.ylabel("Pathway")
plt.title("Top 10 Enriched Pathways (GSEA)")
plt.legend(title="FDR.q-val")
plt.tight_layout()

# Save the plot
plt.savefig(snakemake.output["gsea_plot"], dpi=300)
