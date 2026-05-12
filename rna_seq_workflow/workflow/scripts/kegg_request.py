import pandas as pd
import requests
import json
import time

# Query KEGG API for all unique KO IDs
def fetch_kegg_data(all_kos):
    
    session = requests.Session()  # Optimize API requests
    
    # Create a cache to store the pathway IDs for each KO
    ids_cache = {}
    for ko in all_kos:
        ids_cache[ko] = {"pathway": [], "module": []}
        try:
            # Make a request to the KEGG API for pathway IDs
            url_pathway = f"https://rest.kegg.jp/link/pathway/{ko}"
            response = session.get(url_pathway)
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
                ids_cache[ko]["pathway"] = pathway_ids
                #time.sleep(0.5)  # Add a delay between requests
                
            # Make a request to the KEGG API for module IDs
            url_module = f"https://rest.kegg.jp/link/module/{ko}"
            response = requests.get(url_module)
            # Check if the response is successful
            if response.status_code == 200:
                module_ids = []
                # Parse the response and extract the pathway IDs
                for line in response.text.splitlines():
                    module_ids_list = line.split('\t')[1].strip().split(',')
                    ko_module_ids = [module_id for module_id in module_ids_list if module_id.startswith('md:')]
                    if ko_module_ids:
                        print(f"Module ID found for KO {ko}: {ko_module_ids}")
                        module_ids.extend(ko_module_ids)
                ids_cache[ko]["module"] = module_ids
                #time.sleep(0.5)  # Add a delay between requests
            else:
                print(f"Error {response.status_code} for KO {ko}")
        except Exception as e:
            print(f"Exception: {e}")
    return ids_cache

# Function to fetch names for pathways & modules
def fetch_kegg_names(kegg_ids, category):
    """
    Fetch names for KEGG pathways or modules.

    :param kegg_ids: Set of KEGG pathway/module IDs.
    :param category: "pathway" or "module" (used for logging).
    :return: Dictionary mapping KEGG IDs to their names.
    """
    names_cache = {}
    session = requests.Session()  # Use a session for efficiency

    for kegg_id in kegg_ids:
        try:
            url = f"https://rest.kegg.jp/get/{kegg_id}"
            response = session.get(url)
            if response.status_code == 200:
                for line in response.text.splitlines():
                    if line.startswith("NAME"):
                        name = line.replace("NAME", "").strip()
                        names_cache[kegg_id] = name
                        print(f"{category.capitalize()} Name found for {kegg_id}: {name}")
                        break  # Stop after finding the name
            else:
                print(f"Error {response.status_code} for {category} {kegg_id}")

            time.sleep(0.5)  # Rate limiting

        except requests.RequestException as e:
            print(f"Exception fetching {category} {kegg_id}: {e}")

    return names_cache

# Collect all unique KO IDs from input files
input_files = snakemake.input["functional_annotations"]  # List of input files
all_kos = set()

for file in input_files:
    df = pd.read_csv(file)
    df['KEGG_ko'] = df['KEGG_ko'].str.split(',')
    kos = df['KEGG_ko'].explode().str.strip().dropna()
    all_kos.update(kos)

# Fetch the pathway and module IDs for all unique KO IDs
ko_data = fetch_kegg_data(all_kos)

# Extract unique pathway IDs
unique_pathway_ids = {pathway for pathways in ko_data.values() for pathway in pathways["pathway"]}

# Extract unique module IDs
unique_module_ids = {module for modules in ko_data.values() for module in modules["module"]}


# Fetch pathway and module names
pathway_names_cache = fetch_kegg_names(unique_pathway_ids, "pathway")
module_names_cache = fetch_kegg_names(unique_module_ids, "module")

# Map KO IDs to pathway and module names
ko_kegg_names = {
    ko_id: {
        "pathways": [pathway_names_cache[pathway_id] for pathway_id in ko_data[ko_id]["pathway"] if pathway_id in pathway_names_cache],
        "modules": [module_names_cache[module_id] for module_id in ko_data[ko_id]["module"] if module_id in module_names_cache]
    }
    for ko_id in ko_data
}

# Save the combined KEGG cache to a JSON file
with open(snakemake.output["kegg_cache"], "w") as cache_file:
    json.dump(ko_kegg_names, cache_file, indent=4)