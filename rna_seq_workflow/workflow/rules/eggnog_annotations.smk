rule eggnog_annotations:
    input:
        annotated_results = "results/{run_id}/functional_annotations/interproscan/{contrast}_interpro_annotations.csv",
        eggnog_results = config["eggnog_results"],

    output:
        combined_annotations = "results/{run_id}/functional_annotations/{contrast}_combined_annotations.csv"

    conda:
        "../envs/annotation_integration.yaml"

    script:
        "../scripts/eggnog_annotations.py" 