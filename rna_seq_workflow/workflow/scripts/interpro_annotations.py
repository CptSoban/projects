import pandas as pd
import numpy as np
import gffutils
from collections import defaultdict

first_results = pd.read_csv(snakemake.input["first_results"])
# eggnog_results = pd.read_csv(snakemake.input["eggnog_results"], 
#                             sep="\t",
#                             header=4,
#                             index_col=0,
#                             skipfooter=3)

first_results = first_results.set_index("protein_ids")

# eggnog_results = eggnog_results[["Description","KEGG_ko"]]
# eggnog_results = eggnog_results.rename(columns={"Description": "eggNOG_description",})

# annotated_results = first_results.join(eggnog_results)

interpro_gff = gffutils.FeatureDB(snakemake.input["interpro_gff"])

source_desc = {}
interpro_results = defaultdict(dict)
# Iterate through the features in the database
for feature in interpro_gff.all_features():
  # Check if the 'Target'(=protein id) and 'signature_desc'(=protein description) attributes exist
  if 'Target' in feature.attributes and 'signature_desc' in feature.attributes:
    #Source = Database
    source = feature.source
    descriptions = feature.attributes['signature_desc']
    targets = feature.attributes['Target']
    #Split the protein id from the specified location coordinates
    targets = targets[0].split(" ")[0]
    #Check if source already exists in the nested dict as subkey
    if source in interpro_results[targets].keys():
      #If it exists add the new descritption as value, turn it into a set to remove duplicates and back into a list 
      interpro_results[targets][source] = list(set(interpro_results[targets][source]+descriptions))
    else:
      #If it doesn't exist add it with its corresponding descritption as value
      interpro_results[targets][source] = descriptions

interpro_results = dict(interpro_results)
interpro_results_df = pd.DataFrame.from_dict(interpro_results, orient="index")

annotated_results = first_results.join(interpro_results_df)
annotated_results = annotated_results.drop('MobiDBLite', axis=1)
annotated_results = annotated_results.sort_values(by=["log2FC_shrinked","padj"], ascending=False)

annotated_results.to_csv(snakemake.output["annotated_results"])
