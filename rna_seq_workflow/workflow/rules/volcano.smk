rule volcano:
    input:
        mod_deg_table = f"{results_dir}/DEG_analysis/{run_id}/{run_id}_mod_deseq2_results.csv",
    output:
        volcano_plot = f"{results_dir}/DEG_analysis/{run_id}/plots/{run_id}_volcano_plot.pdf",
    conda:
        "enhanced_volcano.yaml"
    script:
        "volcano_plotting.R"