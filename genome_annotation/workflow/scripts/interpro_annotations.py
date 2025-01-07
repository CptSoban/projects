import pandas as pd
import gffutils
from collections import defaultdict

# Load GFF file
interpro_gff = gffutils.FeatureDB(snakemake.input["interpro_gff"])

interpro_results = defaultdict(lambda: {'annotations': defaultdict(list), 'GO_terms': set()})

# Iterate through the features in the database
for feature in interpro_gff.all_features():
    # Check if the 'Target' (protein id) and 'signature_desc' (protein description) attributes exist
    if 'Target' in feature.attributes and 'signature_desc' in feature.attributes:
        source = feature.source
        descriptions = feature.attributes['signature_desc']
        targets = feature.attributes['Target']
        
        # Extract the protein ID without location coordinates
        target_id = targets[0].split(" ")[0]
        
        # Add descriptions to the appropriate source
        interpro_results[target_id]['annotations'][source].extend(descriptions)
        
        # Check for GO terms
        if 'Ontology_term' in feature.attributes:
            go_terms = feature.attributes['Ontology_term']
            interpro_results[target_id]['GO_terms'].update(go_terms)

# Convert defaultdict to regular dictionary and process GO terms
processed_results = {}
for protein_id, data in interpro_results.items():
    processed_results[protein_id] = {'GO_terms': ';'.join(data['GO_terms'])}
    for source, desc in data['annotations'].items():
        processed_results[protein_id][source] = list(set(desc))  # Remove duplicates

# Convert dictionary to DataFrame
interpro_results_df = pd.DataFrame.from_dict(processed_results, orient="index")

# Drop the 'MobiDBLite' column if it exists
if 'MobiDBLite' in interpro_results_df.columns:
    interpro_results_df = interpro_results_df.drop('MobiDBLite', axis=1)

# Save annotated results to CSV
interpro_results_df.to_csv(snakemake.output["annotated_results"])

