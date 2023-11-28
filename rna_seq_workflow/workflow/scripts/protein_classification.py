import pandas as pd

positive_DEG_df = pd.read_csv(snakemake.input[["positive_DEGs"]])

product_annotation_df = pd.read_csv(snakemake.input[["product_annotation"]])

#Annotate ranked gene list with products
positive_DEG_df = positive_DEG_df.set_index("symbol").join(product_annotation_df.set_index("Geneid"))

#Filter out genes with a shrinked log2FC of less than 2
positive_DEG_df_filt = positive_DEG_df[positive_DEG_df.iloc[:, 0] >= 2]



