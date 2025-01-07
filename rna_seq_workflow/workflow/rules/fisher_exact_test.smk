rule fisher_exact_test:
    input:
        functional_annotations = "results/{run_id}/functional_annotations/{contrast}_combined_annotations.csv"
    output:
        positive_enrichment = "results/{run_id}/enrichment_analysis/fisher_exact/{contrast}_positive_enrichment.csv",
        negative_enrichment = "results/{run_id}/enrichment_analysis/fisher_exact/{contrast}_negative_enrichment.csv",
        stacked_bar_plot = "results/{run_id}/enrichment_analysis/fisher_exact/{contrast}_stacked_bar_plot.png"
    conda:
        "../envs/enrichment_analysis.yaml"
    script:
        "../scripts/fisher_exact_test.py"