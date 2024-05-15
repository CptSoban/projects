rule heatmap:
    input:
        sample_info = config["sample_info"],
        mod_deg_table = config["workflow_dir"]+"/results/DEG_analysis/{run_id}/{run_id}_mod_deseq2_results.csv",
        normalized_counts = config["workflow_dir"]+"/results/DEG_analysis/{run_id}/{run_id}_normalized_counts.csv",
    output:
        heatmap_plot = config["workflow_dir"]+"/results/DEG_analysis/{run_id}/plots/{run_id}_heatmap.pdf",
    conda:
        "../envs/complex_heatmap.yaml"
    script:
        "../scripts/heatmap_plotting.R"