rule repeatmasker:
    input:
        genome_assembly = config["genome_assembly"],

    output:
        masked_genome = config["work_dir"]+"/workflow/results/repeatmasker/{run_id}/{run_id}.fasta.masked",

    params:
        out_dir = directory(config["work_dir"]+"/workflow/results/repeatmasker/{run_id}"),
        reference_taxon = config["DFAM taxon"]

    threads:
        config["threads"]

    conda:
        "../envs/repeatmasker.yaml"

    shell: """RepeatMasker -species {params.reference_taxon} -dir {params.out_dir} -gff -e ncbi -s {input.genome_assembly}""" 