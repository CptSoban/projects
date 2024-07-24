rule interpro_annotation:
    input:
        first_results = "results/DEG_analysis/{run_id}/{contrast}/main_{run_id}_{contrast}.csv",
        interpro_gff = "resources/{run_id}/{run_id}_{contrast}_interpro_results_db"

    output:
        annotated_results = "results/functional_annotations/interproscan/{run_id}/{contrast}/{run_id}_{contrast}_interpro_annotations.csv",

    threads:
        config["threads"]

    conda:
        "../envs/annotation_integration.yaml"

    script:
        "../scripts/interpro_annotations.py"  