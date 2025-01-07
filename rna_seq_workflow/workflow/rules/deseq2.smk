rule deseq2:
    input:
        prep_counts_table = "results/{run_id}/feature_counts/mod_{run_id}_counts.csv",
        sample_info = config["sample_info"]
        

    output:
        deseq_dataset = "results/{run_id}/DEG_analysis/{run_id}_dds.rds",
        normalized_counts = "results/{run_id}/DEG_analysis/{run_id}_normalized_counts.csv",
        disp_plot = "results/{run_id}/DEG_analysis/plots/{run_id}_dispersion_plot.pdf",
        pca_plot = "results/{run_id}/DEG_analysis/plots/{run_id}_pca_plot.svg",

    conda:
        "../envs/deg_analysis.yaml"
        
    script:
        "../scripts/deseq2.R"

rule deseq2_contrast:
    input:
        deseq_dataset = "results/{run_id}/DEG_analysis/{run_id}_dds.rds",

    output:
        mod_deg_table = "results/{run_id}/DEG_analysis/deg_{contrast}.csv",

    params:
        contrast_pair = lambda wildcards: wildcards.contrast.split("_vs_"),

    conda:
        "../envs/deg_analysis.yaml"

    script:
        "../scripts/deseq2_contrast.R"