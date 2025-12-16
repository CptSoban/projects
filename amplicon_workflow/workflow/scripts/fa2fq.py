import pandas as pd
import numpy as np
from Bio import SeqIO
import re
import gzip

# Load the no-chimera records
nochim_records = SeqIO.to_dict(SeqIO.parse(snakemake.input["chimera_filtered_clusters"], "fasta"))
# Load the clustering data
clust_df = pd.read_csv(snakemake.input["clustering_table"], sep="\t", header=None)
# Load the dereplication data
derep_df = pd.read_csv(snakemake.input["dereplication_table"], sep="\t", header=None)

# Remove the cluster records
clust_df[[0]] = np.where(clust_df[[0]] == "C", np.nan, clust_df[[0]])
clust_df = clust_df.dropna(axis=0)

# Select the last two columns and rename them
clust_df = clust_df.iloc[:, -2:]
clust_df = clust_df.rename(columns={clust_df.columns[0]: "derep_id", clust_df.columns[1]: "cluster_id"})
# Replace '*' in cluster_id with derep_id
clust_df["cluster_id"] = np.where(clust_df["cluster_id"] == "*", clust_df["derep_id"], clust_df["cluster_id"])

# Extract the part before the semicolon in derep_id
def extract_before_semicolon(cell):
    if isinstance(cell, str):
        match = re.match(r'^(.*?);', cell)
        if match:
            return match.group(1)
    return cell  # Return original if no match or not a string

# Apply the function to the derep_id column
clust_df = clust_df.map(extract_before_semicolon)
# Group by cluster_id and aggregate derep_id into a set
clust_df = clust_df.groupby("cluster_id")["derep_id"].apply(set)
# Convert the DataFrame to a dictionary
clust_dict = clust_df.to_dict()
# Filter the clust_dict to only include keys that are in nochim_records
filtered_records = {}
for key, value in clust_dict.items():
    if key in nochim_records:
        filtered_records[key] = value

# Select the first two columns and rename them
derep_df = derep_df.iloc[:, :2]
derep_df = derep_df.rename(columns={derep_df.columns[0]: "fq_id", derep_df.columns[1]: "fa_id"})

# Group by fa_id and aggregate fq_id into a list
derep_df = derep_df.groupby("fa_id")["fq_id"].apply(list)
# Convert the DataFrame to a dictionary
derep_dict = derep_df.to_dict()

# Create a set to hold the original IDs from derep_dict
filt_og_ids = set()
# Iterate over the filtered records to collect original IDs
for cluster, derep_ids in filtered_records.items():
    # For each cluster, collect the corresponding original_ids from derep_dict
    original_ids = []
    for id_ in derep_ids:
        if id_ in derep_dict:
            original_ids.append(",".join(derep_dict[id_]))
    # Add the collected original_ids to the result dictionary
    filt_og_ids.update(original_ids)

# Write the filtered original IDs to a new FASTQ file
with gzip.open(snakemake.input["chopper_filt_reads"], "rt") as out_handle:
    subset_records = (record for record in SeqIO.parse(out_handle, "fastq") if record.id in filt_og_ids)
    with gzip.open(snakemake.output["clean_reads"], "wt") as out_handle:
        SeqIO.write(subset_records, out_handle, "fastq")
