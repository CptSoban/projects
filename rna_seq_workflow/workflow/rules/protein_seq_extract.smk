rule protein_seq_extract:
    input:
        mod_deg_table = "results/DEG_analysis/{run_id}/deg_{contrast}.csv",
        gtf_database = "resources/{run_id}_gtf_db",
        ref_protein_fasta = config["protein_fasta"],

    output:
        main_results = "results/DEG_analysis/{run_id}/{contrast}/trim_{contrast}.csv",
        positive_DEGs = "results/DEG_analysis/{run_id}/{contrast}/pos_{contrast}.csv",
        negative_DEGs = "results/DEG_analysis/{run_id}/{contrast}/neg_{contrast}.csv",
        all_protein_sequences = "results/DEG_analysis/{run_id}/{contrast}/{contrast}.aa",
        pos_protein_sequences = "results/DEG_analysis/{run_id}/{contrast}/pos_{contrast}.aa",
        neg_protein_sequences = "results/DEG_analysis/{run_id}/{contrast}/neg_{contrast}.aa",

    conda:
        "../envs/protein_seq_extract.yaml"

    script:
        "../scripts/protein_seq_extract.py"


