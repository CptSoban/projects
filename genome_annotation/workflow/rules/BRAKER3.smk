
rule AUGUSTUS_install:
    output:
        augustus_config = directory("resources/Augustus/config"),

    shell:  """git clone https://github.com/Gaius-Augustus/Augustus.git resources/Augustus"""  

rule build_braker3:
    output:
        "resources/braker3.sif"

    shell:  """apptainer build --force resources/braker3.sif docker://teambraker/braker3:v3.0.7.6"""

rule merge_bams:
    input:
        expand("results/{run_id}/hisat2_align/{sample}.sortedByCoord.out.bam", run_id=RUN_ID, sample=SAMPLES.sample)
    output:
        "results/{run_id}/hisat2_align/merged.bam"
    
    shell:  """samtools merge {output} {input}"""
    
rule BRAKER3:
    input:
        merged_bam = "results/{run_id}/hisat2_align/merged.bam",
        augustus_config = "resources/Augustus/config",
        masked_genome = "results/{run_id}/repeatmasker/"+ASSEMBLY_FILE+".masked",
        container_file = "resources/braker3.sif",
    output:
        braker_gtf = "results/{run_id}/braker/braker.gtf",
        braker_aa = "results/{run_id}/braker/braker.aa",
    params:
        fungi = "--fungus" if config["Fungi"] == "yes" else "",
        species_name = RUN_ID,
        out_dir = directory("results/{run_id}/braker/")

    threads: config["threads"]

    shell:  """
            apptainer exec -B $(realpath {params.out_dir}) \
            {input.container_file} braker.pl \
            {params.fungi} \
            --species={params.species_name} \
            --genome={input.masked_genome} \
            --bam={input.merged_bam} \
            --workingdir={params.out_dir} \
            --AUGUSTUS_CONFIG_PATH=$(realpath {input.augustus_config}) \
            --threads={threads} 
            """

    # in case Augustus is updated and is not the same version as in the braker3.sif, you can instead copy the config folder from within the braker3.sif
    # apptainer exec --cleanenv {input.container_file} cp -r /opt/Augustus/config . |\
    # chmod -R 777 ./config |\