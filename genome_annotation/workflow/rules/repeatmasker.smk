
rule repeatmasker:
    input:
        genome_assembly = config["genome_assembly"],
        repeat_db = "results/{run_id}/repeatmodeler/{run_id}_DB-families.fa"

    output:
        masked_genome = "results/{run_id}/repeatmasker/"+ASSEMBLY_FILE+".masked",

    params:
        out_dir = directory("results/{run_id}/repeatmasker"),
        #reference_taxon = config["DFAM taxon"],

    conda:
        "../envs/repeatmasker.yaml"

    shell: """RepeatMasker -dir {params.out_dir} -gff -xsmall -e rmblast -lib {input.repeat_db} -s {input.genome_assembly}""" 