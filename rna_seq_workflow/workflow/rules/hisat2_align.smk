rule hisat2_align:
    input:
        r1 = config["seq_dir"]+"/{sample}_R1.fastq.gz",
        r2 = config["seq_dir"]+"/{sample}_R2.fastq.gz",
        index_files = "results/hisat2_index/{run_id}/{run_id}.1.ht2"
    output:
        bam_file = "results/hisat2_align/{run_id}/{sample}.sortedByCoord.out.bam"

    log:
        summary = "reports/hisat2_align/{run_id}/{sample}_summary.log",
        metrics = "reports/hisat2_align/{run_id}/{sample}_metrics.log"

    params:
        basename = "results/hisat2_index/{run_id}/{run_id}",
        strandness = config["Strandness"][0]
    threads:
        config["threads"]

    conda:
        "../envs/hisat2.yaml"

    shell:  """hisat2 -p {threads} \
            -x {params.basename} \
            -1 {input.r1} \
            -2 {input.r2} \
            --rna-strandness {params.strandness} \
            --summary-file {log.summary} \
            --met-file {log.metrics} | \
            samtools sort -o {output.bam_file}"""
