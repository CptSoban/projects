rule repeatmodeler:
    input:
        genome_assembly = config["genome_assembly"],

    output:
        repeat_db = "results/{run_id}/repeatmodeler/{run_id}_DB-families.fa"
    
    params:
        run_id = RUN_ID,

    threads: config["threads"]

    conda:
        "../envs/repeatmodeler.yml"

    shell:  """
            mkdir -p results/{wildcards.run_id}/repeatmodeler && \
            BuildDatabase -name results/{wildcards.run_id}/repeatmodeler/{wildcards.run_id}_DB {input.genome_assembly} && \
            RepeatModeler -database results/{wildcards.run_id}/repeatmodeler/{wildcards.run_id}_DB \
              -threads {threads} \
              -LTRStruct \
"""

#)"""BuildDatabase -name {wildcards.run_id}_DB {input.genome_assembly} && \
 #           RepeatModeler -database {wildcards.run_id}_DB -threads {threads} -LTRStruct""" 

 #              -dir results/{wildcards.run_id}/repeatmodeler
