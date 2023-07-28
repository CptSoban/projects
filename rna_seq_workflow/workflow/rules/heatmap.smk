rule heatmap:
    input:
        sample_info = f"{seq_dir}/{run_id}_info.txt",
        mod_deg_table = f"{results_dir}/DEG_analysis/{run_id}/{run_id}_mod_deseq2_results.csv",
        normalized_counts = f"{results_dir}/DEG_analysis/{run_id}/{run_id}_normalized_counts.csv",
    output:
        heatmap_plot = f"{results_dir}/DEG_analysis/{run_id}/plots/{run_id}_heatmap.pdf",
    conda:
        "complex_heatmap.yaml"
    script:
        "heatmap_plotting.R"