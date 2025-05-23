rule hisat2_align:
    input:
        r1 = config["seq_dir"]+"/{sample}_R1.fastq.gz",
        r2 = config["seq_dir"]+"/{sample}_R2.fastq.gz",
        index_files = "results/{run_id}/hisat2_index/{run_id}.1.ht2",
        fastqc_output = expand("results/{run_id}/fastqc/{sample}_{group}_fastqc.html", run_id=RUN_ID, sample=set(SAMPLES.sample), group=set(SAMPLES.group))
    output:
        bam_file = "results/{run_id}/hisat2_align/{sample}.sortedByCoord.out.bam"

    log:
        summary = "reports/{run_id}/hisat2_align/{sample}_summary.log",
        metrics = "reports/{run_id}/hisat2_align/{sample}_metrics.log"

    params:
        basename = "results/{run_id}/hisat2_index/{run_id}",
        strandness = config["Strandness"][0]

    conda:
        "../envs/hisat2.yaml"

    shell:  """hisat2 -x {params.basename} \
            -1 {input.r1} \
            -2 {input.r2} \
            --rna-strandness {params.strandness} \
            --summary-file {log.summary} \
            --met-file {log.metrics} | \
            samtools sort -o {output.bam_file}"""
