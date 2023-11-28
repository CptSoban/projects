import pandas as pd
import numpy as np
from Bio import SeqIO
import gffutils

mod_res = pd.read_csv(snakemake.input["mod_deg_table"])

#Only keep those columns
positive_DEG_df = mod_res[["symbol", "log2FC_shrinked", "padj"]]

#Only keep rows (=genes) where the log2FC equal or greater than 2 (=Increased expression)
positive_DEG_df_filt = positive_DEG_df[positive_DEG_df["log2FC_shrinked"] >= 2]

#Sort the rows (=genes) by their shrinked log2FC value (=Highest expression change first)
positive_DEG_df = positive_DEG_df.sort_values(by=["log2FC_shrinked"], ascending=False)

#Remove any empthy rows
positive_DEG_df.dropna(inplace = True)

#Contains gene product information contained in the corresponding GTF
product_annotation_df = pd.read_csv(snakemake.input[["product_annotation"]])

#Creates a database of the GTF/GFF3 file
gffutils.constants.always_return_list = True
gtf_db = gffutils.create_db(snakemake.input[["gtf"]], snakemake.output[["gtf_database"]], merge_strategy="merge", keep_order=True)

#Join both dfs based on their index(=Gene ID/Symbol)
positive_DEG_df = positive_DEG_df.set_index("symbol").join(product_annotation_df.set_index("Geneid"))

#Returns string instead of list
gffutils.constants.always_return_list = False
gene_ids = positive_DEG_df_filt.index.tolist()
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
        protein_list.append(i["protein_id"])
    
    #To remove duplicates the list is turned into a set
    protein_dict[gene_id] = set(protein_list)
    #The list for the current gene_id is emptied, so that it can be filled with the protein_ids for the next gene
    protein_list = []

#Convert dict to dataframe. Keeps protein_ids in a list in a single column
protein_df = pd.DataFrame([protein_dict])
#Gene_ids as index and renaming new protein_ids column
protein_df = protein_df.transpose().rename(columns={0:"protein_ids"})
#Joining protein_ids with main DEGs dataframe
positive_DEG_df_filt = positive_DEG_df_filt.join(protein_df)

positive_DEG_df_filt.to_csv(snakemake.output[["main_results"]])

#Flattens nested list of all protein_ids of interest
ids_to_extract = [item for sublist in list(protein_dict.values()) for item in sublist]

# Parse the fasta file
records = SeqIO.parse(snakemake.input[["ref_protein_fasta"]], "fasta")

# Write selected sequences to a new fasta file
with open(snakemake.output[["extracted_protein_sequences"]], "w") as output_handle:
    SeqIO.write((record for record in records if record.id in ids_to_extract), output_handle, "fasta")
