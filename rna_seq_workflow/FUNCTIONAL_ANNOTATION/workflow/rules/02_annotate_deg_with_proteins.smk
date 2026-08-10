rule protein_seq_extract:
    input:
        deg_table = config["deg_results"]+"/DEG_analysis/deg_{contrast}.csv",
        gtf_database = "resources/{run_id}_gtf_db",
        ref_protein_fasta = config["protein_fasta"],

    output:
        main_results = "results/{run_id}/trimmed/trim_{contrast}.csv"

    params:
        base_mean_threshold = config["min_base_mean"]

    conda:
        "../envs/annotate_deg_with_proteins.yaml"

    script:
        "../scripts/annotate_deg_with_proteins.py"


