import pandas as pd
import requests
import json
import time

# Collect all unique KO IDs from input files
input_files = snakemake.input["functional_annotations"]  # List of input files
all_kos = set()

for file in input_files:
    df = pd.read_csv(file)
    df['KEGG_ko'] = df['KEGG_ko'].str.split(',')
    kos = df['KEGG_ko'].explode().str.strip().dropna()
    all_kos.update(kos)

# Query KEGG API for all unique KO IDs
def fetch_pathway_ids(all_kos):
    # Create a cache to store the pathway IDs for each KO
    pathway_ids_cache = {}
    for ko in all_kos:
        try:
            # Make a request to the KEGG API
            url = f"https://rest.kegg.jp/link/pathway/{ko}"
            response = requests.get(url)
            # Check if the response is successful
            if response.status_code == 200:
                pathway_ids = []
                # Parse the response and extract the pathway IDs
                for line in response.text.splitlines():
                    pathway_ids_list = line.split('\t')[1].strip().split(',')
                    ko_path_ids = [pathway_id for pathway_id in pathway_ids_list if pathway_id.startswith('path:ko')]
                    if ko_path_ids:
                        print(f"Pathway ID found for KO {ko}: {ko_path_ids}")
                        pathway_ids.extend(ko_path_ids)
                pathway_ids_cache[ko] = pathway_ids
                time.sleep(0.5)  # Add a delay between requests
            else:
                print(f"Error {response.status_code} for KO {ko}")
        except Exception as e:
            print(f"Exception: {e}")
    return pathway_ids_cache

# Fetch the pathway IDs for all unique KO IDs
ko_pathway_ids = fetch_pathway_ids(all_kos)

# Extract all unique pathway IDs from the cache
unique_pathway_ids = set(pathway_id for pathway_ids in ko_pathway_ids.values() for pathway_id in pathway_ids)

# Query KEGG API for pathway names
def fetch_pathway_names(unique_pathway_ids):
    pathway_names_cache = {}
    for pathway_id in unique_pathway_ids:
        try:
            # Make a request to the KEGG API
            url = f"https://rest.kegg.jp/get/{pathway_id}"
            response = requests.get(url)
            # Check if the response is successful
            if response.status_code == 200:
                # Parse the response and extract the pathway name
                for line in response.text.splitlines():
                    if line.startswith("NAME"):
                        pathway_name = line.replace("NAME", "").strip()
                        pathway_names_cache[pathway_id] = pathway_name
                        print(f"Pathway Name found for {pathway_id}: {pathway_name}")       
                time.sleep(0.5)  # Add a delay between requests
            else:
                print(f"Error {response.status_code} for pathway {pathway_id}")
        except Exception as e:
            print(f"Exception: {e}")
    return pathway_names_cache

# Fetch the pathway names for all unique pathway IDs
pathway_names_cache = fetch_pathway_names(unique_pathway_ids)

# Create a dictionary to map KO IDs to pathway names
ko_pathway_names = {
    key: [pathway_names_cache[pathway_id] for pathway_id in pathway_id_list if pathway_id in pathway_names_cache]
    for key, pathway_id_list in ko_pathway_ids.items()}

# Save the pathway cache to a JSON file
with open(snakemake.output["kegg_cache"], "w") as cache_file:
    json.dump(ko_pathway_names, cache_file)