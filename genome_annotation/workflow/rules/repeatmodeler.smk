rule repeatmodeler:
    input:
        genome_assembly = config["genome_assembly"],

    output:
        repeat_db = "resources/{run_id}/repeatmodeler/{run_id}-families.fa"
    
    params:
        run_id = RUN_ID,

    threads: config["threads"]

    conda:
        "../envs/repeatmodeler.yml"

    shell: """BuildDatabase -name {run_id}_DB {input.genome_assembly} && \
            RepeatModeler -database {run_id}_DB -threads {threads} -LTRStruct""" 