rule deseq2_prep:
    input:
        counts_table = "results/feature_counts/{run_id}/{run_id}_counts.txt",
        sample_info = config["sample_info"]
    output:
        prep_counts_table = "results/feature_counts/{run_id}/mod_{run_id}_counts.csv",
    log:
        filtered = "reports/{run_id}/{run_id}_abundance_filter.log"
    conda:
        "../envs/deseq2_prep.yaml"
    script:
        "../scripts/deseq2_prep.py"