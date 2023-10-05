import pandas as pd
import numpy as np
import gseapy as gp
from gseapy.plot import gseaplot

gene_ranking = pd.read_csv(snakemake.input["preranked_genes"], header=None)
#gene_ranking = gene_ranking.drop([0])
print(gene_ranking)
pre_res = gp.prerank(rnk = gene_ranking, gene_sets = snakemake.params["gene_set"], seed = snakemake.params["seed"], min_size=2, max_size=1500, verbose=True)

out = []
for term in list(pre_res.results):
    out.append([term, pre_res.results[term]["fdr"], 
                pre_res.results[term]["es"], 
                pre_res.results[term]["nes"]])

out_df = pd.DataFrame(out, columns = ["Term", "fdr", "es", "nes"]).sort_values("fdr").reset_index(drop=True)

out_df.to_csv(snakemake.output["gsea_results"])