rule volcano:
    input:
        mod_deg_table = "/home/marc/projects/rna_seq_workflow/results/DEG_analysis/{run_id}_mod_deseq2_results.csv",
    output:
        volcano_plot = "/home/marc/projects/rna_seq_workflow/results/DEG_analysis/plots/{run_id}_volcano_plot.pdf",
    conda:
        "../envs/enhanced_volcano.yaml"
    script:
        "../scripts/volcano_plotting.R"