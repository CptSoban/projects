rule eggnog_annotations:
    input:
        annotated_results = config["workflow_dir"]+"/results/functional_annotations/{run_id}/interproscan/{run_id}_interpro_annotations.csv",
        eggnog_results = config["eggnog_results"],

    output:
        combined_annotations = config["workflow_dir"]+"/results/functional_annotations/{run_id}/{run_id}_combined_annotations.csv"

    threads:
        config["threads"]

    conda:
        "../envs/annotation_integration.yaml"

    script:
        "../scripts/eggnog_annotations.py" 