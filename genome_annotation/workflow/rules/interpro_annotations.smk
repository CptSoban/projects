rule interpro_annotation:
    input:
        interpro_gff = "resources/{run_id}/interpro_results_db"

    output:
        annotated_results = "results/functional_annotations/interproscan/{run_id}/interpro_annotations.csv",

    conda:
        "../envs/annotation_integration.yaml"

    script:
        "../scripts/interpro_annotations.py"  