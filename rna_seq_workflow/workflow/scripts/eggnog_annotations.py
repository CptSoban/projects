#!/usr/bin/python
import os
import pandas as pd

annotated_results = pd.read_csv(snakemake.input["annotated_results"],)

eggnog_results = pd.read_csv(snakemake.input["eggnog_results"], 
                            sep="\t",
                            header=4,
                            index_col=0,
                            skipfooter=3)

annotated_results = annotated_results.set_index("protein_ids")

eggnog_results = eggnog_results[["Description","KEGG_ko","PFAMs"]]
eggnog_results = eggnog_results.rename(columns={"Description": "eggNOG_description", "PFAMs": "eggnog_Pfam"})

annotated_results = annotated_results.join(eggnog_results)

# annotated_results = annotated_results[["symbol",
#                                        "log2FC_shrinked",
#                                        "padj",
#                                        "product",
#                                        "eggNOG_description",
#                                        "KEGG_ko",
#                                        "ProSitePatterns",
#                                        "ProSiteProfiles",
#                                        "CDD",
#                                        "Pfam",
#                                        "Hamap",
#                                        "FunFam",
#                                        "PRINTS",
#                                        "NCBIfam",
#                                        "SFLD"]]

annotated_results = annotated_results.sort_values(by=["log2FoldChange","padj"], ascending=False)

annotated_results.to_csv(snakemake.output["combined_annotations"])

# Try to delete the interpro file.
try:
    os.remove(snakemake.input["annotated_results"])
except OSError as e:
    # If it fails, inform the user.
    print("Error: %s - %s." % (e.filename, e.strerror))