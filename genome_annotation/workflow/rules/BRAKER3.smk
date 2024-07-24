rule AUGUSTUS_install:
    output:
        augustus_config = directory(config["work_dir"]+"/resources/Augustus/config"),

    shell:  """git clone https://github.com/Gaius-Augustus/Augustus.git resources/"""  

rule decompress_gz:
    input:
        config["rna_reads_directory"]+"/{sample}_{group}.fastq.gz"
    output:
        config["rna_reads_directory"]+"/{sample}_{group}.fastq"

    shell:  """gzip -d {input}"""


rule BRAKER3:
    input:
        expand(config["rna_reads_directory"]+"/{sample}_{group}.fastq", sample=set(samples.sample), group=["1","2"]),
        augustus_config = config["work_dir"]+"/resources/Augustus/config",
        masked_genome = config["work_dir"]+"/workflow/results/repeatmasker/{run_id}/{run_id}.fasta.masked",
        rna_reads_dir = config["rna_reads_directory"],
    output:
        braker_gtf = config["work_dir"]+"/workflow/results/braker/{run_id}/braker.gtf"
    params:
        # rna_reads_ids = get_basenames(config["rna_reads_directory"]),
        fungi = "--fungus" if config["Fungi"] == "+" else "",
        species_name = config["Run ID"],
        out_dir = directory(config["work_dir"]+"/workflow/results/braker/{run_id}")

    container:
        "/home/marc/projects/genome_annotation/resources/braker3.sif"
    threads:
        config["Threads"]

    shell: """braker.pl --threads={threads} {params.fungi} --species={params.species_name} --genome={input.masked_genome} --rnaseq_sets_ids={wildcards.sample} --rnaseq_sets_dirs={input.rna_reads_dir} --workingdir={params.out_dir} --AUGUSTUS_CONFIG_PATH={input.augustus_config}"""

rule compress_gz:
    input:
        rna_reads = config["rna_reads_directory"]+"/{sample}_{group}.fastq",
        expand(config["work_dir"]+"/workflow/results/braker/{run_id}/braker.gtf", run_id=config["Run ID"])
    output:
        config["rna_reads_directory"]+"/{sample}_{group}.fastq.gz"

    shell:  """gzip {input.rna_reads}"""