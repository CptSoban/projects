rule hisat2_align:
    input:
        r1 = config["workflow_dir"]+"/results/cutadapt/{run_id}/trim_{sample}_R1.fastq" if config["Trimming"] == "yes" else config["seq_dir"]+"/{sample}_R1.fastq",
        r2 = config["workflow_dir"]+"/results/cutadapt/{run_id}/trim_{sample}_R2.fastq" if config["Trimming"] == "yes" else config["seq_dir"]+"/{sample}_R2.fastq",
        index_files = config["workflow_dir"]+"/results/hisat2_index/{run_id}/{run_id}.1.ht2"
    output:
        bam_file = config["workflow_dir"]+"/results/hisat2_align/{run_id}/{sample}.sortedByCoord.out.bam"

    log:
        summary = config["workflow_dir"]+"/reports/hisat2_align/{run_id}/{sample}_summary.log",
        metrics = config["workflow_dir"]+"/reports/hisat2_align/{run_id}/{sample}_metrics.log"

    params:
        basename = config["workflow_dir"]+"/results/hisat2_index/{run_id}/{run_id}",
        strandness = config["Strandness"][0]
    threads:
        config["threads"]

    conda:
        "../envs/hisat2.yaml"

    shell: """hisat2 -p {threads} -x {params.basename} -1 {input.r1} -2 {input.r2} --rna-strandness {params.strandness} --summary-file {log.summary} --met-file {log.metrics} | samtools sort -o {output.bam_file}"""
