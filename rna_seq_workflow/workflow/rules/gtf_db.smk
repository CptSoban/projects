rule gtf_db:
    input:
        gtf = config["gtf"],

    output:
        gtf_database = "resources/{run_id}_gtf_db",

    threads: config["threads"]
    
    conda:
        "../envs/gffutils_db.yaml"

    script:
       "../scripts/gtf_gffutils_db.py"
    
