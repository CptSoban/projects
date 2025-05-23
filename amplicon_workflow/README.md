Amplicon Workflow Information 

This workflow is designed to process amplicon sequencing data from MinION runs.

It includes steps for:  Concatenation, length and quality filtering (fastcat v0.22.0),  
                        Chimera removal (vsearch v2.30.0), 
                        Taxonomic classification using emu (emu v3.5.1)

Regarding Primer trimming: When setting up the sequencing run through MinKNOW it is highly recommended to activate the option to trim adapter sequences during basecalling. Trimming of PCR-Primer sequences is possible but will have virtually no effect on classification results, since modern aligners like minimap2 will "soft-clip" these sequences during alignment.

1.  Activate the snakemake conda activate:

    conda activate snakemake

2.  Fill in the fields in the config.yaml with your information and save the file as config.yaml in the directory projects/amplicon_workflow/config/.
    Keep a copy of your finished config file in config/old_configs for reference (recommended).

3.  Make sure to have dorado installed and available in your PATH. Check installation with:

    dorado --version
    
4.  To test the workflow, perform a dry run with the command in projects/amplicon_workflow/:

    snakemake -np

5.  To run the workflow, use the command in projects/amplicon_workflow/:

    snakemake --use-conda --conda-frontend conda --cores all
    
    Only use --cores all if only one workflow is running. If multiple workflows are running, specify the number of cores to use (e.g., --cores 20 (max=28)).