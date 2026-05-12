import pandas as pd
import numpy as np
from Bio import SeqIO
import gffutils

#Load gtf database
gtf_db = gffutils.FeatureDB(snakemake.input["gtf_database"])

gffutils.constants.always_return_list = True

#DEG analysis results table 
all_genes_df = pd.read_csv(snakemake.input["mod_deg_table"])

#Only keep those columns
all_genes_df = all_genes_df[["symbol", "log2FoldChange", "baseMean", "padj"]]

#Remove any empthy rows
all_genes_df.dropna(inplace = True)

#Filter lowly expressed genes
most_DEG_df = all_genes_df.loc[all_genes_df["baseMean"] >= snakemake.params["base_mean_threshold"]]

#Filter for most significant DEGs
# most_DEG_df = most_DEG_df.loc[(most_DEG_df["padj"] <= 0.05) &
#                               ((most_DEG_df["log2FoldChange"] >= 0.5) | (most_DEG_df["log2FoldChange"] <= -0.5))]


#Annotates GTF product information to the created dataframe
most_DEG_df = most_DEG_df.set_index("symbol")

#Returns string instead of list
gffutils.constants.always_return_list = False
gene_ids = most_DEG_df.index.tolist()
protein_list = []
all_protein_ids = []
protein_dict = {}

#Builds a dict with gene_ids as key and protein_ids as value
for gene_id in gene_ids:
    #Gene entry in GTF
    gene = gtf_db[gene_id]
    #The CDS entry contains the protein_id so here for every CDS feature ("children") of a gene of interest 
    #all corresponding protein_ids are appended to a list
    for i in gtf_db.children(gene, featuretype="CDS"):
        try:
            protein_list.append(i["protein_id"])
        #If genome was annotated using BRAKER pipeline
        except:
            protein_list.append(i["transcript_id"])
    
    #To remove duplicates the list is turned into a set
    protein_dict[gene_id] = set(protein_list)
    #The list for the current gene_id is emptied, so that it can be filled with the protein_ids for the next gene
    protein_list = []

#Convert dict to dataframe. Keeps protein_ids in a list in a single column
protein_df = pd.DataFrame([protein_dict])

#Gene_ids as index and renaming new protein_ids column
protein_df = protein_df.transpose().rename(columns={0:"protein_ids"})

protein_df = protein_df.explode('protein_ids')

#Flattens nested list of all protein_ids of interest
ids_to_extract = [item for sublist in list(protein_dict.values()) for item in sublist]

#Parse the fasta file
records = list(SeqIO.parse(snakemake.input["ref_protein_fasta"], "fasta"))


#Joining protein_ids with main DEGs dataframe
most_DEG_df = most_DEG_df.join(protein_df).sort_values(by=['log2FoldChange','padj'], ascending=False)

most_DEG_df.to_csv(snakemake.output["main_results"])

# #Table containing positive DEGs
# positive_DEG_df = most_DEG_df[most_DEG_df["log2FoldChange"] >= 2.0]

# #Mask to filter non-significant positive DEGs
# filter_mask_padj = (positive_DEG_df["padj"] <= 0.05)

# #Table containing significant positive DEGs
# positive_DEG_df = positive_DEG_df[filter_mask_padj]

# positive_DEG_df.to_csv(snakemake.output["positive_DEGs"])

# #Table containing negative DEGs
# negative_DEG_df = most_DEG_df[most_DEG_df["log2FoldChange"] <= -2.0]

# #Mask to filter non-significant negative DEGs
# filter_mask_padj = (negative_DEG_df["padj"] <= 0.05)

# #Table containing significant negative DEGs
# negative_DEG_df = negative_DEG_df[filter_mask_padj]

# negative_DEG_df.to_csv(snakemake.output["negative_DEGs"])
