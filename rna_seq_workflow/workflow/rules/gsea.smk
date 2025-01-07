rule gsea:
    input:
        functional_annotations = "results/{run_id}/functional_annotations/{contrast}_combined_annotations.csv"
    output:
        gsea_results = "results/{run_id}/enrichment_analysis/gsea/{contrast}_gsea_results.csv",
        gsea_plot = "results/{run_id}/enrichment_analysis/gsea/{contrast}_gsea_plot.png"

    conda:
        "../envs/gsea.yaml"

    script:
        "../scripts/gsea.py"