rule deseq2:
    input:
        prep_counts_table = "results/feature_counts/{run_id}/mod_{run_id}_counts.csv",
        sample_info = config["sample_info"]
        

    output:
        deseq_dataset = "results/DEG_analysis/{run_id}/{run_id}_dds.rds",
        normalized_counts = "results/DEG_analysis/{run_id}/{run_id}_normalized_counts.csv",
        disp_plot = "results/DEG_analysis/{run_id}/plots/{run_id}_dispersion_plot.pdf",
        pca_plot = "results/DEG_analysis/{run_id}/plots/{run_id}_pca_plot.svg",

    conda:
        "../envs/deg_analysis.yaml"
        
    script:
        "../scripts/deseq2.R"

rule deseq2_contrast:
    input:
        deseq_dataset = "results/DEG_analysis/{run_id}/{run_id}_dds.rds",

    output:
        mod_deg_table = "results/DEG_analysis/{run_id}/deg_{contrast}.csv",

    params:
        contrast_pair = lambda wildcards: wildcards.contrast.split("_vs_"),

    conda:
        "../envs/deg_analysis.yaml"

    script:
        "../scripts/deseq2_contrast.R"