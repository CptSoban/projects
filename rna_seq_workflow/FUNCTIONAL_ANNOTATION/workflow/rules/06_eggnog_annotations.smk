rule eggnog_annotations:
    input:
        annotated_results = "results/{run_id}/interproscan/{contrast}_interpro_annotations.csv",
        eggnog_results = "results/{run_id}/eggnog/{run_id}.emapper.annotations",

    output:
        combined_annotations = "results/{run_id}/functional_annotations/{contrast}_combined_annotations.csv"

    conda:
        "../envs/annotation_integration.yaml"

    script:
        "../scripts/eggnog_annotations.py" 