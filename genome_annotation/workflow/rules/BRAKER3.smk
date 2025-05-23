rule AUGUSTUS_install:
    output:
        augustus_config = directory("resources/Augustus/config"),

    shell:  """git clone https://github.com/Gaius-Augustus/Augustus.git resources/Augustus"""  

rule build_braker3:
    output:
        "resources/braker3.sif"

    shell:  """apptainer build --force resources/braker3.sif docker://teambraker/braker3:latest"""

rule BRAKER3:
    input:
        rna_reads = expand(config["rna_reads_directory"]+"/{sample}_{group}.fastq.gz", sample=set(SAMPLES.sample), group=set(SAMPLES.group)),
        augustus_config = "resources/Augustus/config",
        masked_genome = "results/{run_id}/repeatmasker/"+ASSEMBLY_BASE+".fasta.masked",
        rna_reads_dir = config["rna_reads_directory"],
        container_file = "resources/braker3.sif",
    output:
        braker_gtf = "results/{run_id}/braker/braker.gtf",
        braker_aa = "results/{run_id}/braker/braker.aa",
    params:
        fungi = "--fungus" if config["Fungi"] == "yes" else "",
        species_name = config["Run ID"], ###?
        rna_ids = "cloroAOK1", #",".join(set(SAMPLES.sample)),
        out_dir = directory("results/{run_id}/braker/")
    

    shell:  """
            apptainer exec -B $(realpath {params.out_dir}) \
            {input.container_file} braker.pl \
            {params.fungi} \
            --species={params.species_name} \
            --genome={input.masked_genome} \
            --rnaseq_sets_ids={params.rna_ids} \
            --rnaseq_sets_dirs={input.rna_reads_dir} \
            --workingdir={params.out_dir} \
            --AUGUSTUS_CONFIG_PATH=$(realpath {input.augustus_config}) \
            --threads=8
            """

    # in case Augustus is updated and is not the same version as in the braker3.sif, you can instead copy the config folder from within the braker3.sif
    # apptainer exec --cleanenv {input.container_file} cp -r /opt/Augustus/config . |\
    # chmod -R 777 ./config |\