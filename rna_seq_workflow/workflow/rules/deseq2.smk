rule deseq2:
    input:
        prep_counts_table = f"{results_dir}/feature_counts/{run_id}/mod_{run_id}_counts.csv",
        sample_info = f"{seq_dir}/{run_id}_info.txt",
        mapID = f"{results_dir}/feature_counts/{run_id}/{run_id}_mapid.csv"
    output:
        deg_table = f"{results_dir}/DEG_analysis/{run_id}/{run_id}_deseq2_results.csv",
        filt_deg_table = f"{results_dir}/DEG_analysis/{run_id}/{run_id}_sig_deseq2_results.csv",
        mod_deg_table = f"{results_dir}/DEG_analysis/{run_id}/{run_id}_mod_deseq2_results.csv",
        normalized_counts = f"{results_dir}/DEG_analysis/{run_id}/{run_id}_normalized_counts.csv",
        disp_plot = f"{results_dir}/DEG_analysis/{run_id}/plots/{run_id}_dispersion_plot.pdf",
        pca_plot = f"{results_dir}/DEG_analysis/{run_id}/plots/{run_id}_pca_plot.pdf"
    conda:
        "deg_analysis.yaml"
    script:
        "deseq2_workflow.R"