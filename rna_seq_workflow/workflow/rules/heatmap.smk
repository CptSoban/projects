rule heatmap:
    input:
        sample_info = config["sample_info"],
        mod_deg_table = "/home/marc/projects/rna_seq_workflow/results/DEG_analysis/{run_id}/{run_id}_mod_deseq2_results.csv",
        normalized_counts = "/home/marc/projects/rna_seq_workflow/results/DEG_analysis/{run_id}/{run_id}_normalized_counts.csv",
    output:
        heatmap_plot = "/home/marc/projects/rna_seq_workflow/results/DEG_analysis/plots/{run_id}/{run_id}_heatmap.pdf",
    conda:
        "../envs/complex_heatmap.yaml"
    script:
        "../scripts/heatmap_plotting.R"