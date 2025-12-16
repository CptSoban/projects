import pandas as pd
import gffutils
from collections import defaultdict

# Load initial results data
first_results = pd.read_csv(snakemake.input["first_results"])
first_results = first_results.set_index("protein_ids")

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

# Join initial results with the annotations
annotated_results = first_results.join(interpro_results_df)
# Drop the 'MobiDBLite' column if it exists
if 'MobiDBLite' in annotated_results.columns:
    annotated_results = annotated_results.drop('MobiDBLite', axis=1)
# Sort results by 'log2FoldChange' and 'padj'
annotated_results = annotated_results.sort_values(by=["log2FoldChange", "padj"], ascending=False)

# Save annotated results to CSV
annotated_results.to_csv(snakemake.output["annotated_results"])


# import pandas as pd
# import numpy as np
# import gffutils
# from collections import defaultdict

# # Load initial results data
# first_results = pd.read_csv(snakemake.input["first_results"])
# first_results = first_results.set_index("protein_ids")

# #Load GFF file
# interpro_gff = gffutils.FeatureDB(snakemake.input["interpro_gff"])

# interpro_results = defaultdict(dict)

# # Iterate through the features in the database
# for feature in interpro_gff.all_features():
#   # Check if the 'Target'(=protein id) and 'signature_desc'(=protein description) attributes exist
#   if 'Target' in feature.attributes and 'signature_desc' in feature.attributes:
#     #Source = Database
#     source = feature.source
#     descriptions = feature.attributes['signature_desc']
#     targets = feature.attributes['Target']

#     #Split the protein id from the specified location coordinates
#     targets = targets[0].split(" ")[0]

#     #Check if source already exists in the nested dict as subkey
#     if source in interpro_results[targets].keys():
#       #If it exists add the new descritption as value, turn it into a set to remove duplicates and back into a list 
#       interpro_results[targets][source] = list(set(interpro_results[targets][source]+descriptions))
#     else:
#       #If it doesn't exist add it with its corresponding descritption as value
#       interpro_results[targets][source] = descriptions

# # Convert defaultdict to regular dictionary
# interpro_results = dict(interpro_results)
# # Convert dictionary to DataFrame
# interpro_results_df = pd.DataFrame.from_dict(interpro_results, orient="index")

# # Join initial results with the annotations
# annotated_results = first_results.join(interpro_results_df)
# # Drop the 'MobiDBLite' column
# annotated_results = annotated_results.drop('MobiDBLite', axis=1)
# # Sort results by 'log2FoldChange' and 'padj'
# annotated_results = annotated_results.sort_values(by=["log2FoldChange","padj"], ascending=False)

# # Save annotated results to CSV
# annotated_results.to_csv(snakemake.output["annotated_results"])
