rule protein_seq_extract:
    input:
        mod_deg_table = "results/{run_id}/DEG_analysis/deg_{contrast}.csv",
        gtf_database = "resources/{run_id}_gtf_db",
        ref_protein_fasta = config["protein_fasta"],

    output:
        main_results = "results/{run_id}/DEG_analysis/trimmed/trim_{contrast}.csv",
        # positive_DEGs = "results/{run_id}/DEG_analysis/{contrast}/pos_{contrast}.csv",
        # negative_DEGs = "results/{run_id}/DEG_analysis/{contrast}/neg_{contrast}.csv",

    conda:
        "../envs/protein_seq_extract.yaml"

    script:
        "../scripts/protein_seq_extract.py"


