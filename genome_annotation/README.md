Genome Annotation Workflow Information 

This workflow is designed to annotate eukaryotic genomes with coordinate information (exons, introns, etc.).

It includes steps for:  Modeling and masking of repeat regions (Repeatmodeler, Repeatmasker),  
                        Alignment of RNA reads (HISAT2), 
                        Genome annotation (BRAKER)

For improved predictions the workflow relies on either your own RNA-seq reads or RNA-seq data stored in SRA.

1.  Activate the snakemake conda activate:

    conda activate snakemake

2.  Fill in the fields in the config.yaml with your information and save the file as config.yaml in the directory projects/genome_annotation/config/.
    You can keep a copy of your config file in config/old_configs for reference.
    
3.  To test the workflow, perform a dry run with the command in projects/genome_annotation/:

    snakemake -np

4.  To run the workflow, use the command in projects/genome_annotation/:

    snakemake --use-conda --conda-frontend conda --cores all
    
    This workflow does not benefit from >8 cores.