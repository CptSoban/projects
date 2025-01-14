rule fisher_exact_test:
    input:
        functional_annotations = "results/{run_id}/functional_annotations/{contrast}_combined_annotations.csv",
        kegg_cache = "resources/{run_id}/enrichment_analysis/kegg_request/{run_id}_kegg_cache.json"
    output:
        positive_enrichment = "results/{run_id}/enrichment_analysis/fisher_exact/{contrast}_positive_enrichment.csv",
        forest_plot = "results/{run_id}/enrichment_analysis/fisher_exact/{contrast}_forest_plot.png"
    conda:
        "../envs/fisher_exact.yaml"
    script:
        "../scripts/fisher_exact_test.py"