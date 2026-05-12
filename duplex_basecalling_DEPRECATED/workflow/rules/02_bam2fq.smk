rule bam2fq:
    input:
        filtered_bam = config["pod5_dir"] + "duplex_basecalled/duplex_filtered.bam",
    output:
        fastq = config["pod5_dir"] + "duplex_basecalled/duplex_filtered.fastq.gz",
    threads: config["threads"],

    conda:
        "../envs/bam2fq.yaml"

    shell:"""
        echo "Converting filtered BAM to FASTQ..." &&\
        samtools fastq -@ {threads} {input.filtered_bam} | \
        gzip > {output.fastq}
        """
