rule deseq2_prep:
    input:
        counts_table = config["workflow_dir"]+"/results/feature_counts/{run_id}/{run_id}_counts.txt",
        sample_info = config["sample_info"]
    output:
        prep_counts_table = config["workflow_dir"]+"/results/feature_counts/{run_id}/mod_{run_id}_counts.csv",
        product_annotation = config["workflow_dir"]+"/results/feature_counts/{run_id}/{run_id}_product_annotation.csv"
    log:
        filtered = config["workflow_dir"]+"/reports/{run_id}/{run_id}_abundance_filter.log"
    conda:
        "../envs/deseq2_prep.yaml"
    script:
        "../scripts/deseq2_prep.py"