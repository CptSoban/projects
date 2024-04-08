rule interpro_annotation:
    input:
        first_results = config["workflow_dir"]+"/results/DEG_analysis/{run_id}/{run_id}_main_results.csv",
        interpro_gff = config["workflow_dir"]+"/resources/{run_id}_interpro_results_db"

    output:
        annotated_results = config["workflow_dir"]+"/results/functional_annotations/{run_id}/interproscan/{run_id}_interpro_annotations.csv",

    threads:
        config["threads"]

    conda:
        "../envs/annotation_integration.yaml"

    script:
        "../scripts/interpro_annotations.py"  