rule deseq2:
    input:
        prep_counts_table = config["workflow_dir"]+"/results/feature_counts/{run_id}/mod_{run_id}_counts.csv",
        sample_info = config["sample_info"],
       # mapID = "results/feature_counts/{run_id}_mapid.csv"
    output:
        deg_table = config["workflow_dir"]+"/results/DEG_analysis/{run_id}/{run_id}_deseq2_results.csv",
        mod_deg_table = config["workflow_dir"]+"/results/DEG_analysis/{run_id}/{run_id}_mod_deseq2_results.csv",
        normalized_counts = config["workflow_dir"]+"/results/DEG_analysis/{run_id}/{run_id}_normalized_counts.csv",
        disp_plot = config["workflow_dir"]+"/results/DEG_analysis/{run_id}/plots/{run_id}_dispersion_plot.pdf",
        pca_plot = config["workflow_dir"]+"/results/DEG_analysis/{run_id}/plots/{run_id}_pca_plot.svg"
    conda:
        "../envs/deg_analysis.yaml"
    script:
        "../scripts/deseq2_workflow.R"