rule hisat2_align:
    input:
        r1 = config["rna_reads_directory"]+"/{sample}_1.fastq.gz",
        r2 = config["rna_reads_directory"]+"/{sample}_2.fastq.gz",
        index_files = "resources/{run_id}/hisat2_index/{run_id}.1.ht2"
    output:
        bam_file = "results/{run_id}/hisat2_align/{sample}.sortedByCoord.out.bam"

    log:
        summary = "reports/{run_id}/hisat2_align/{sample}_summary.log",
        metrics = "reports/{run_id}/hisat2_align/{sample}_metrics.log"

    params:
        basename = "resources/{run_id}/hisat2_index/{run_id}",
        strandness = config["Strandness"][0]

    conda:
        "../envs/hisat2.yaml"

    shell:  """hisat2 -x {params.basename} \
            -1 {input.r1} \
            -2 {input.r2} \
            --dta \
            --rna-strandness {params.strandness} \
            --summary-file {log.summary} \
            --met-file {log.metrics} | \
            samtools sort -o {output.bam_file}"""
