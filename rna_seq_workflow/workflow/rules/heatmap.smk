rule heatmap:
    input:
        sample_info = config["sample_info"],
        normalized_counts = "results/{run_id}/DEG_analysis/{run_id}_normalized_counts.csv",

    output:
        heatmap_plot = "results/{run_id}/DEG_analysis/plots/heatmap.pdf",

    conda:
        "../envs/complex_heatmap.yaml"
        
    script:
        "../scripts/heatmap_plotting.R"