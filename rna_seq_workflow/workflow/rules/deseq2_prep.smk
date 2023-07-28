rule deseq2_prep:
    input:
        counts_table = f"{results_dir}/feature_counts/{run_id}/{run_id}_counts.txt",
        sample_info = f"{seq_dir}/{run_id}_info.txt"
    output:
        prep_counts_table = f"{results_dir}/feature_counts/{run_id}/mod_{run_id}_counts.csv",
        mapID = f"{results_dir}/feature_counts/{run_id}/{run_id}_mapid.csv"
    log:
        filtered = f"{results_dir}/feature_counts/{run_id}/abundance_filter.log"
    conda:
        "deseq2_prep.yaml"
    script:
        "deseq2_prep.py"