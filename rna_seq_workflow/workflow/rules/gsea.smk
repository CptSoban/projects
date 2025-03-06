rule gsea:
    input:
        functional_annotations = "results/{run_id}/functional_annotations/{contrast}_combined_annotations.csv",
        kegg_cache = "resources/{run_id}/enrichment_analysis/kegg_request/{run_id}_kegg_cache.json",
    output:
        gsea_results = "results/{run_id}/enrichment_analysis/gsea/{contrast}_gsea_results.csv",
        gsea_plot = "results/{run_id}/enrichment_analysis/gsea/{contrast}_gsea_plot.png"
    
    params:
        filter_terms = config["filter_terms"]
    
    log:
        pathway_recovery = "reports/{run_id}/enrichment_analysis/gsea/{contrast}_pathway_recovery.log"

    conda:
        "../envs/gsea.yaml"

    script:
        "../scripts/gsea.py"