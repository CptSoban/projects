rule star_align:
    input:
        r1 = results_dir+"/cutadapt/trim_{sample}_R1.fastq" if trimming == "y" else seq_dir+"/{sample}_R1.fastq",
        r2 = results_dir+"/cutadapt/trim_{sample}_R2.fastq" if trimming == "y" else seq_dir+"/{sample}_R2.fastq",
        indexed_genome_dir = results_dir+"/{aligner}/index/"+run_id
    output:
        bam_file = results_dir+"/{aligner}/alignment/"+run_id+"/{sample}.sortedByCoord.out.bam"
    params:
        output_prefix = results_dir+"/{aligner}/alignment/{sample}/{sample}"
    threads: config["threads"]
    conda:
        "star.yaml"

    shell: """STAR runMode alignReads --runThreadN {threads} --genomeDir {input.indexed_genome} --readFilesIn {input.r1}, {input.r2} --readFilesCommand zcat --outFileNamePrefix {params.output_prefix} --outSAMtype BAM SortedByCoordinate"""
