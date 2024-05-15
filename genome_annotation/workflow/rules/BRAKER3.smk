rule BRAKER3:
    input:
        masked_genome = config["workflow_dir"]+"/results/repeatmasker/{run_id}/{run_id}.fasta.masked",
        rna_reads_ids = config["RNA reads directory"],
        rna_reads_dir = config["RNA reads directory"]

    output:
        config["workflow_dir"]+"/results/repeatmasker/{run_id}/braker.gtf"
    
    params:
        fungi = "--fungus" if config["Fungi"] == "yes" else "",
        species_name = config["Run ID"],
        out_dir = directory(config["workflow_dir"]+"/results/repeatmasker/{run_id}")

    container: 
        "docker://teambraker/braker3:latest",

    shell: """braker.pl {params.fungi} --species={params.species_name} --genome={input.masked_genome} --rnaseq_sets_ids=BAM_ID1,BAM_ID2 --rnaseq_sets_dirs={input.rna_reads_dir} ----output_dir={params.out_dir}"""
