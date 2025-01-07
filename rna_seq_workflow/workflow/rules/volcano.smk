rule volcano:
    input:
        mod_deg_table = "results/{run_id}/DEG_analysis/deg_{contrast}.csv",

    output:
        volcano_plot = "results/{run_id}/DEG_analysis/plots/{contrast}_volcano.pdf",

    conda:
        "../envs/enhanced_volcano.yaml"
        
    script:
        "../scripts/volcano_plotting.R"