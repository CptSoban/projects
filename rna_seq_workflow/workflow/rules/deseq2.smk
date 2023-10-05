rule deseq2:
    input:
        prep_counts_table = "/home/marc/projects/rna_seq_workflow/results/feature_counts/mod_{run_id}_counts.csv",
        sample_info = config["sample_info"],
       # mapID = "/home/marc/projects/rna_seq_workflow/results/feature_counts/{run_id}_mapid.csv"
    output:
        deg_table = "/home/marc/projects/rna_seq_workflow/results/DEG_analysis/{run_id}_deseq2_results.csv",
        mod_deg_table = "/home/marc/projects/rna_seq_workflow/results/DEG_analysis/{run_id}_mod_deseq2_results.csv",
        normalized_counts = "/home/marc/projects/rna_seq_workflow/results/DEG_analysis/{run_id}_normalized_counts.csv",
        disp_plot = "/home/marc/projects/rna_seq_workflow/results/DEG_analysis/plots/{run_id}_dispersion_plot.pdf",
        pca_plot = "/home/marc/projects/rna_seq_workflow/results/DEG_analysis/plots/{run_id}_pca_plot.pdf"
    conda:
        "../envs/deg_analysis.yaml"
    script:
        "../scripts/deseq2_workflow.R"