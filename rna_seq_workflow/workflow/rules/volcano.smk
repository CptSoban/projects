rule volcano:
    input:
        mod_deg_table = "results/DEG_analysis/{run_id}/deg_{contrast}.csv",

    output:
        volcano_plot = "results/DEG_analysis/{run_id}/plots/{contrast}_volcano.pdf",

    conda:
        "../envs/enhanced_volcano.yaml"
        
    script:
        "../scripts/volcano_plotting.R"