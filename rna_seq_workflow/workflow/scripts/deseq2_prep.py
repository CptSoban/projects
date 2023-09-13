import pandas as pd 
import numpy as np
from pathlib import Path

#Read in counts table
counts_df = pd.read_csv(snakemake.input["counts_table"], sep= "\t", header=1)

#Filters out genes with less than 10 counts across all samples
unfiltered_counts = len(counts_df.index)
counts_df = counts_df[counts_df.iloc[:, 7:].sum(axis=1) >= 10] # Removes genes with a total count of less than 10
filtered_counts = len(counts_df.index)
with open(snakemake.log["filtered"], "w") as log_file:
    log_file.write(f"From {unfiltered_counts} genes, {filtered_counts} occured at least 10 times across samples and where used for subsequent analysis")
        
#Makes a copy of the "Geneid" and "gene" column, which can be used to convert Gene IDs to Gene symbols for GO-Analysis and GSEA)
#symbols_df = counts_df.filter(items=["Geneid","gene"])
#symbols_df.to_csv(snakemake.output["mapID"], index=False)

#Drop unnecessary columns that would stop DESeq2 from correctly identifying the count table
counts_df = counts_df.drop(["Chr","Start","End","Strand","Length"], axis=1)
       
#Replace "realpath" for sample id as column name
for col_name in counts_df.columns:
    new_col_name = Path(col_name).stem.removesuffix(".sortedByCoord.out")
    counts_df.rename(columns={col_name: new_col_name}, inplace=True)

#Read in sample information (e.g. "treated" or "control")
sample_info_df = pd.read_csv(snakemake.input["sample_info"], sep="\t", header=0)
#Put sample ids from sample information in a list
sorted_sample_list = sample_info_df["sample"].tolist()
#Geneid is a column name in the countstable, an has to be included in the sorted list
sorted_sample_list.insert(0, "Geneid") 
#Orders the columns based on the order of sample ids in the sample information (Requiered by DESeq2)
counts_df = counts_df[sorted_sample_list]

counts_df.to_csv(snakemake.output["prep_counts_table"], index=False)