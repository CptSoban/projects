rule interproscan_db:
    input:
        interpro_gff = config["workflow_dir"]+"/results/functional_annotations/{run_id}/interproscan/{run_id}.gff3",

    output:
        interpro_results_db = config["workflow_dir"]+"/resources/{run_id}_interpro_results_db",

    threads: config["threads"]
    
    conda:
        "../envs/gffutils_db.yaml"

    script:
       "../scripts/interproscan_db.py"