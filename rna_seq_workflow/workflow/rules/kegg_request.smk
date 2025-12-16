rule kegg_request:
    input:
        functional_annotations = expand("results/{run_id}/functional_annotations/{contrast}_combined_annotations.csv", run_id=RUN_ID, contrast=CONTRAST)
    output:
        kegg_cache = "resources/{run_id}/enrichment_analysis/kegg_request/{run_id}_kegg_cache.json",
    
    # log:
    #     pathway_recovery = "reports/{run_id}/enrichment_analysis/kegg_request/{contrast}_pathway_recovery.log"

    conda:
        "../envs/kegg_request.yaml"

    script:
        "../scripts/kegg_request.py"