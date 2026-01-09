rule bam2fq:
    input:
        filtered_bam = config["pod5_dir"] + "duplex_basecalled/duplex_filtered.bam",
    output:
        fastq = config["pod5_dir"] + "duplex_basecalled/duplex_filt.fastq",
    threads: config["threads"],

    conda:
        "../envs/bam2fq.yaml"

    shell:"""
        echo "Converting filtered BAM to FASTQ..." &&\
        samtools fastq -@ {threads} {input.filtered_bam} > {output.fastq}
        """
