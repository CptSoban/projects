rule deseq2_prep:
    input:
        counts_table = "/home/marc/projects/rna_seq_workflow/results/feature_counts/{run_id}/{run_id}_counts.txt",
        sample_info = config["sample_info"]
    output:
        prep_counts_table = "/home/marc/projects/rna_seq_workflow/results/feature_counts/{run_id}/mod_{run_id}_counts.csv",
        product_annotation = "/home/marc/projects/rna_seq_workflow/results/feature_counts/{run_id}/{run_id}_product_annotation.csv"
    log:
        filtered = "/home/marc/projects/rna_seq_workflow/reports/{run_id}/{run_id}_abundance_filter.log"
    conda:
        "../envs/deseq2_prep.yaml"
    script:
        "../scripts/deseq2_prep.py"