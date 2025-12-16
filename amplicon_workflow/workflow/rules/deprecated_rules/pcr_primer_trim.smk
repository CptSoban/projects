rule fastcat2:
    input:
        reads = lambda wildcards: expand(config['reads_dir']+"/{barcode}/{sample}.fastq.gz", barcode=wildcards.barcode, sample=samples_per_barcode[wildcards.barcode]),
    output:
        concat_trimm_reads = temp("results/{run_id}/trimming/{barcode}_trim.fastq")
    params:
        reads_dir = config["reads_dir"]+"/{barcode}",
        min_length = config["min_length"],
        max_length = config["max_length"],
        min_quality = config["min_quality"]
    
    log:    temp("reports/{run_id}/{barcode}_fastcat.log")
    conda:
        "../envs/trim.yaml"
    
    shell:  """
            fastcat \
                -a {params.min_length} \
                -b {params.max_length} \
                --min_qscore={params.min_quality} \
                --histograms {log} \
                {params.reads_dir} \
                > {output.concat_trimm_reads} \
            """

rule dorado_exists:
    output: touch("command_available.txt")
    params: cmd="dorado"
    run:
        import shutil
        if shutil.which(params.cmd) is None:
            raise ValueError(f"Command '{params.cmd}' not found in PATH. Please install dorado and make sure it is in your PATH.")
        shell("touch {output}")

rule primer_trim:
    input:
        concat_trimm_reads = "results/{run_id}/trimming/{barcode}_trim.fastq",
        primer = config["primer_fasta"]
    output:
        trim_reads = "results/{run_id}/trimming/{barcode}_finaltrim.fastq.gz"
    
    shell: """
            dorado trim --primer-sequences {input.primer} --sequencing-kit SQK-NBD114-96 --emit-fastq {input.concat_trimm_reads} | gzip > {output.trim_reads}
            """