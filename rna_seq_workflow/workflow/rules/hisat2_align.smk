rule hisat2_align:
    input:
        r1 = working_dir+"/cutadapt/trim_{sample}_R1.fastq" if trimming == "y" else seq_dir+"/{sample}_R1.fastq",
        r2 = working_dir+"/cutadapt/trim_{sample}_R2.fastq" if trimming == "y" else seq_dir+"/{sample}_R2.fastq"
    output:
        bam_file = working_dir+"/results/{aligner}_alignment/{sample}.sortedByCoord.out.bam"
    log:
        summary = working_dir+"/report/{aligner}_alignment/{sample}_summary.log",
        metrics = working_dir+"/report/{aligner}_alignment/{sample}_metrics.log"
    params:
        indexed_genome_prefix = working_dir+"/{aligner}/index/"+run_id+"/"+run_id
    threads:
        config["threads"]
    conda:
        f"{working_dir}/workflow/envs/hisat2.yaml"

    shell: """hisat2 -p {threads} -x {params.indexed_genome_prefix} -1 {input.r1} -2 {input.r2} --summary-file {log.summary} --met-file {log.metrics} | samtools sort -o {output.bam_file}"""
