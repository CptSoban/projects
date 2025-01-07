rule interpro_annotation:
    input:
        first_results = "results/{run_id}/DEG_analysis/{contrast}/trim_{contrast}.csv",
        interpro_gff = "resources/{run_id}/interpro_results_db"

    output:
        annotated_results = "results/{run_id}/functional_annotations/interproscan/{contrast}_interpro_annotations.csv",

    conda:
        "../envs/annotation_integration.yaml"

    script:
        "../scripts/interpro_annotations.py"  