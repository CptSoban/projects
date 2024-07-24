rule volcano:
    input:
        mod_deg_table = "results/DEG_analysis/{run_id}/{run_id}_{contrast}.csv",
    output:
        volcano_plot = "results/DEG_analysis/{run_id}/plots/{run_id}_{contrast}_volcano.pdf",
    conda:
        "../envs/enhanced_volcano.yaml"
    script:
        "../scripts/volcano_plotting.R"