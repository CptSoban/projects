import pandas as pd
import numpy as np


mod_res = pd.read_csv(snakemake.input["mod_deg_table"])

#Only keep those columns
gene_ranking = mod_res[["symbol", "log2FC_shrinked"]]

#gene_ranking = gene_ranking.drop_duplicates(subset = "symbol", keep="first")

#Only keep rows (=genes) where the log2FC is positive (=Increased expression)
gene_ranking = gene_ranking[gene_ranking["log2FC_shrinked"] > 0]

#Sort the rows (=genes) by their shrinked log2FC value (=Highest expression change first)
#gene_ranking = gene_ranking.sort_values(by=["log2FC_shrinked"], ascending=False)

#Remove any empthy rows
gene_ranking.dropna(inplace = True)

#Output new table without header row and index value as csv
gene_ranking.to_csv(snakemake.output["preranked_genes"], index=False)
