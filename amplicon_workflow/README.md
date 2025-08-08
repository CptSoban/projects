Amplicon Analyis Workflow Information 

This workflow is designed to taxonomically classify amplicons (e.g. full-length 16S rRNA gene).
You have to run the duplex_basecalling workflow before starting this one!

It includes steps for:  - Quality and length filtering (chopper v0.10.0)
                        - Dereplication (vsearch v2.30.0)
                        - Clustering (vsearch)
                        - Chimera filtering (vsearch)
                        - Rereplication 
                        - Taxonomic classification (Emu v3.5.1)
                        
Note: Trimming of barcodes and sequencing adapter is usually already performed by MinKnow/dorado. If not visit https://github.com/nanoporetech/dorado for the necessary commands (dorado trim).
Trimming of PCR primers can be run manually after the duplex_basecalling pipeline. However, it is not a requisite as the aligner utilized here (minimap2) would automatically soft-clip these regions anyway.

1.  Fill in the fields in the config/config.yaml with your information and save it as {NAME}.yaml.

2.  Move to the directory amplicon_workflow/

3.  Activate the snakemake environment:

    conda activate snakemake
    
4.  To test the worfklow, perform a dry run with the command:

    snakemake -np

5.  To run the workflow, use the command:

    snakemake --use-conda --conda-frontend conda --cores all

6. If the workflow is stopped you can restart it:

    snakemake --use-conda --conda-frontend conda --cores all --rerun-incomplete
    