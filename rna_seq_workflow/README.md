RNA-seq workflow information

This workflow is designed to process mRNA sequencing data from Illumina sequencing runs.

It includes steps for:  Quality control (fastqc),
                        Differential gene expression analysis (hisat2->featurecounts->DESeq2),
                        Visualization of expression changes (Volcano plot, heatmap, PCA),
                        Function annotation (InterproScan, eggNOG),
                        Pathway enrichment analysis (GSEA)

To get eggNOG annotations take your AA.fa (from reference or from BRAKER) and upload it to http://eggnog-mapper.embl.de/ with default settings.
When finished download the resulst in .csv format to a folder inside the resources directory, ideally named after your run id.
Include the path in your config.yaml and your are ready to start the workflow.

1.  Activate the snakemake conda activate:

    conda activate rna_seq

2.  Fill in the fields in the config.yaml with your information and save the file as config.yaml in the directory projects/rna_seq_workflow/config/.
    You can keep a copy of your config file in config/old_configs for reference.
    
3.  To test the workflow, perform a dry run use the following command within projects/rna_seq_workflow/:

    snakemake -np

4.  To run the workflow, use the following command within projects/rna_seq_workflow/:

    snakemake --use-conda --conda-frontend conda --cores all
    
    Only use --cores all if only one workflow is running. If multiple workflows or other computational intensive systems (sequencing) are running, specify the number of cores to use (e.g., --cores 20 (max. 32)).