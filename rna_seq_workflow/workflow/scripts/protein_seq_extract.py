import pandas as pd
import numpy as np
from Bio import SeqIO
import gffutils

#Creates a database of the GTF/GFF3 file

try:
    gtf_db = gffutils.FeatureDB(snakemake.output["gtf_database"])
except:
    gtf_db = gffutils.create_db(snakemake.input["gtf"], snakemake.output["gtf_database"])


gffutils.constants.always_return_list = True

#DEG analysis results table 
all_genes_df = pd.read_csv(snakemake.input["mod_deg_table"])

#Only keep those columns
all_genes_df = all_genes_df[["symbol", "log2FC_shrinked", "padj"]]

#Remove any empthy rows
all_genes_df.dropna(inplace = True)

#Sort the rows (=genes) by their shrinked log2FC value (=Highest expression change first)
all_genes_df = all_genes_df.sort_values(by=["log2FC_shrinked"], ascending=False)

#Mask to extract highly increaded and decreased DEGs into ONE table
filter_mask = (all_genes_df["log2FC_shrinked"] >= 2) | (all_genes_df["log2FC_shrinked"] <= -2)

#Table containing highly increaded and decreased DEGs
most_DEG_df = all_genes_df[filter_mask]

#Contains gene product information orignally contained in the corresponding GTF
product_annotation_df = pd.read_csv(snakemake.input["product_annotation"])

#Annotates GTF product information to the created dataframe
most_DEG_df = most_DEG_df.set_index("symbol").join(product_annotation_df.set_index("Geneid"))


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
most_DEG_df = most_DEG_df.join(protein_df)

most_DEG_df.to_csv(snakemake.output["main_results"])

#Table containing highly increased DEGs (Treatment->Control)
positive_DEG_df = most_DEG_df[most_DEG_df["log2FC_shrinked"] >= 2]

positive_DEG_df.to_csv(snakemake.output["positive_DEGs"])

#Table containing highly decreased DEGs (Treatment->Control)
negative_DEG_df = most_DEG_df[most_DEG_df["log2FC_shrinked"] <= 2]

negative_DEG_df.to_csv(snakemake.output["negative_DEGs"])

#Flattens nested list of all protein_ids of interest
ids_to_extract = [item for sublist in list(protein_dict.values()) for item in sublist]

#Parse the fasta file
records = SeqIO.parse(snakemake.input["ref_protein_fasta"], "fasta")

#Write selected sequences to a new fasta file
with open(snakemake.output["extracted_protein_sequences"], "w") as output_handle:
    SeqIO.write((record for record in records if record.id in ids_to_extract), output_handle, "fasta")
