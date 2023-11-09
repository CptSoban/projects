import pandas as pd

positive_DEG_df = pd.read_csv(snakemake.input[["positive_DEGs"]])

positive_DEG_df = positive_DEG_df.set_index("symbol").join(proteinID_annotation_df.set_index("Geneid"))

positive_DEG_df_filt = positive_DEG_df[positive_DEG_df.iloc[:, 0] >= 2]

