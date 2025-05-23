rule repeatmodeler:
    input:
        genome_assembly = config["genome_assembly"],

    output:
        masked_genome = "results/{run_id}/repeatmasker/"+ASSEMBLY_BASE+".fasta.masked",

    params:
        out_dir = directory("results/{run_id}/repeatmasker"),
        reference_taxon = config["DFAM taxon"],

    conda:
        "../envs/repeatmodeler.yml"

    shell: """RepeatMasker -species {params.reference_taxon} -dir {params.out_dir} -gff -xsmall -e ncbi -s {input.genome_assembly}""" 