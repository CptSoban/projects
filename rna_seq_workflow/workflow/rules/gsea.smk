rule gsea:
    input:
        functional_annotations = "results/{run_id}/functional_annotations/{contrast}_combined_annotations.csv",
        kegg_cache = "resources/{run_id}/enrichment_analysis/kegg_request/{run_id}_kegg_cache.json",
    output:
        gsea_results_pathways = "results/{run_id}/enrichment_analysis/gsea/{contrast}_gsea_results_pathways.csv",
        gsea_plot_pathways = "results/{run_id}/enrichment_analysis/gsea/{contrast}_gsea_plot_pathways.png",
        gsea_results_modules = "results/{run_id}/enrichment_analysis/gsea/{contrast}_gsea_results_modules.csv",
        gsea_plot_modules = "results/{run_id}/enrichment_analysis/gsea/{contrast}_gsea_plot_modules.png"
    
    params:
        filter_terms = config["filter_terms"],
        fdr_threshold = config["fdr_threshold"]
    
    log:
        recovery_log = "reports/{run_id}/enrichment_analysis/gsea/{contrast}_recovery.log"

    conda:
        "../envs/gsea.yaml"

    script:
        "../scripts/gsea.py"
