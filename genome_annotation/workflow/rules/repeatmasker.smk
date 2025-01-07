

rule repeatmasker:
    input:
        genome_assembly = config["genome_assembly"],

    output:
        masked_genome = "results/repeatmasker/{run_id}/{run_id}.fasta.masked",

    params:
        out_dir = directory("results/repeatmasker/{run_id}"),
        reference_taxon = config["DFAM taxon"],

    conda:
        "../envs/repeatmasker.yaml"

    shell: """RepeatMasker -species {params.reference_taxon} -dir {params.out_dir} -gff -xsmall -e ncbi {input.genome_assembly}""" 