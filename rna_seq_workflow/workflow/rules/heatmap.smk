rule heatmap:
    input:
        sample_info = config["sample_info"],
        normalized_counts = "results/DEG_analysis/{run_id}/{run_id}_normalized_counts.csv",

    output:
        heatmap_plot = "results/DEG_analysis/{run_id}/plots/heatmap.pdf",

    conda:
        "../envs/complex_heatmap.yaml"
        
    script:
        "../scripts/heatmap_plotting.R"