rule eggnog_annotations:
    input:
        annotated_results = "results/functional_annotations/interproscan/{run_id}/{contrast}_interpro_annotations.csv",
        eggnog_results = config["eggnog_results"],

    output:
        combined_annotations = "results/functional_annotations/{run_id}/{contrast}_combined_annotations.csv"

    conda:
        "../envs/annotation_integration.yaml"

    script:
        "../scripts/eggnog_annotations.py" 