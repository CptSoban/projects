Duplex Basecalling Workflow Information 

This workflow is designed to perform duplex basecalling on raw sequencing signals in the .pod5 file format.

It includes steps for:  - Duplex Basecalling (dorado v1.0.2)
                        - Filter duplexed simplex reads (samtools v1.13)
                        - Demultiplexing (dorado)

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

