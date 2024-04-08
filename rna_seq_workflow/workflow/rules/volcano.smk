rule volcano:
    input:
        mod_deg_table = config["workflow_dir"]+"/results/DEG_analysis/{run_id}/{run_id}_mod_deseq2_results.csv",
    output:
        volcano_plot = config["workflow_dir"]+"/results/DEG_analysis/{run_id}/plots/{run_id}_volcano_plot.pdf",
    conda:
        "../envs/enhanced_volcano.yaml"
    script:
        "../scripts/volcano_plotting.R"