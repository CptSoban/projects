import pandas as pd
import numpy as np
import gseapy as gp
from gseapy.plot import gseaplot

mod_res = pd.read_csv(snakemake.input["mod_deg_table"])

gene_ranking = mod_res[["symbol", "rank"]]
gene_ranking.dropna(inplace = True)
gene_ranking = gene_ranking.drop_duplicates(subset = "symbol", keep="first")
print(gene_ranking)
pre_res = gp.prerank(rnk = gene_ranking, gene_sets = 'GO_Biological_Process_2023', seed = snakemake.params["seed"], max_size=10000)

out = []
for term in list(pre_res.results):
    out.append([term, pre_res.results[term]["fdr"], 
                pre_res.results[term]["es"], 
                pre_res.results[term]["nes"]])

out_df = pd.DataFrame(out, columns = ["Term", "fdr", "es", "nes"]).sort_values("fdr").reset_index(drop=True)

out_df.to_csv(snakemake.output["gsea_results"])