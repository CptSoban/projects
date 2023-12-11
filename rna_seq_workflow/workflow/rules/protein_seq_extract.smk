rule protein_seq_extract:
    input:
        mod_deg_table = "/home/marc/projects/rna_seq_workflow/results/DEG_analysis/{run_id}/{run_id}_mod_deseq2_results.csv",
        product_annotation = "/home/marc/projects/rna_seq_workflow/results/feature_counts/{run_id}/{run_id}_product_annotation.csv",
        gtf = config["gtf"],
        ref_protein_fasta = config["protein_fasta"],

    output:
        gtf_database = "/home/marc/projects/rna_seq_workflow/resources/{run_id}_gtf_db",
        main_results = "/home/marc/projects/rna_seq_workflow/results/DEG_analysis/{run_id}/{run_id}_main_results.csv",
        positive_DEGs = "/home/marc/projects/rna_seq_workflow/results/DEG_analysis/{run_id}/{run_id}_positive_DEGs.csv",
        negative_DEGs = "/home/marc/projects/rna_seq_workflow/results/DEG_analysis/{run_id}/{run_id}_negative_DEGs.csv",
        extracted_protein_sequences = "/home/marc/projects/rna_seq_workflow/results/DEG_analysis/{run_id}/{run_id}_protein.fasta",

    conda:
        "../envs/protein_seq_extract.yaml"

    script:
        "../scripts/protein_seq_extract.py"


