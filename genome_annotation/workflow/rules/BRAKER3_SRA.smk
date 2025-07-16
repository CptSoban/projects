rule AUGUSTUS_install:
    output:
        augustus_config = directory("resources/Augustus/config"),

    shell:  """git clone https://github.com/Gaius-Augustus/Augustus.git resources/Augustus"""  

# rule decompress_gz:
#     input:
#         config["rna_reads_directory"]+"/{sample}_{group}.fastq.gz"
#     output:
#         config["rna_reads_directory"]+"/{sample}_{group}.fastq"

#     shell:  """gzip -d {input}"""

rule build_braker3:
    output:
        "resources/braker3.sif"

    shell:  """apptainer build --force resources/braker3.sif docker://teambraker/braker3:v3.0.7.6"""

rule BRAKER3_SRA:
    input:
        augustus_config = "/home/nanopore/projects/genome_annotation/resources/Augustus/config", # ???
        masked_genome = "results/{run_id}/repeatmasker/{run_id}.fasta.masked",
        container_file = "resources/braker3.sif",

    output:
        braker_gtf = "results/{run_id}/braker/braker.gtf",
        braker_aa = "results/{run_id}/braker/braker.aa",
    params:
        fungi = "--fungus" if config["Fungi"] == "+" else "",
        species_name = config["Run ID"],
        out_dir = directory("results/{run_id}/braker"),
        rna_reads = config["SRA_accession"],

    threads: config["threads"],
    
    conda:
        "../envs/braker3.yaml"

    shell: """apptainer exec {input.container_file} braker.pl \
            --threads={threads} 
            {params.fungi} \
            --species={params.species_name} \
            --genome={input.masked_genome} \
            --rnaseq_sets_ids={params.rna_reads} \
            --workingdir={params.out_dir} \
            --AUGUSTUS_CONFIG_PATH={input.augustus_config}"""

# rule compress_gz:
#     input:
#         rna_reads = config["rna_reads_directory"]+"/{sample}_{group}.fastq",
#         expand("results/braker/{run_id}/braker.gtf", run_id=config["Run ID"])
#     output:
#         config["rna_reads_directory"]+"/{sample}_{group}.fastq.gz"

#     shell:  """gzip {input.rna_reads}"""