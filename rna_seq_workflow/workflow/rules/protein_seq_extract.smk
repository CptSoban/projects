rule protein_seq_extract:
    input:
        mod_deg_table = config["workflow_dir"]+"/results/DEG_analysis/{run_id}/{run_id}_mod_deseq2_results.csv",
        product_annotation = config["workflow_dir"]+"/results/feature_counts/{run_id}/{run_id}_product_annotation.csv",
        gtf_database = config["workflow_dir"]+"/resources/{run_id}_gtf_db",
        ref_protein_fasta = config["protein_fasta"],

    output:
        main_results = config["workflow_dir"]+"/results/DEG_analysis/{run_id}/{run_id}_main_results.csv",
        positive_DEGs = config["workflow_dir"]+"/results/DEG_analysis/{run_id}/{run_id}_positive_DEGs.csv",
        negative_DEGs = config["workflow_dir"]+"/results/DEG_analysis/{run_id}/{run_id}_negative_DEGs.csv",
        all_protein_sequences = config["workflow_dir"]+"/results/DEG_analysis/{run_id}/{run_id}_protein.fasta",
        pos_protein_sequences = config["workflow_dir"]+"/results/DEG_analysis/{run_id}/{run_id}_pos_protein.fasta",
        neg_protein_sequences = config["workflow_dir"]+"/results/DEG_analysis/{run_id}/{run_id}_neg_protein.fasta",

    conda:
        "../envs/protein_seq_extract.yaml"

    script:
        "../scripts/protein_seq_extract.py"


