#Rule for aligning reads to the reference genome using HISAT2
if config["Strandness"] == "Forward":
    strandness = "FR"
elif config["Strandness"] == "Reverse":
    strandness = "RF"
elif config["Strandness"] == "Unstranded":
    strandness = ""

rule hisat2_align:
    input:
        r1 = "/home/marc/projects/rna_seq_workflow/results/cutadapt/{run_id}/trim_{sample}_R1.fastq" if config["Trimming"] == "yes" else config["seq_dir"]+"/{sample}_R1.fastq",
        r2 = "/home/marc/projects/rna_seq_workflow/results/cutadapt/{run_id}/trim_{sample}_R2.fastq" if config["Trimming"] == "yes" else config["seq_dir"]+"/{sample}_R2.fastq",
        index_files = "/home/marc/projects/rna_seq_workflow/results/hisat2_index/{run_id}/{run_id}.1.ht2"
    output:
        bam_file = "/home/marc/projects/rna_seq_workflow/results/hisat2_align/{run_id}/{sample}.sortedByCoord.out.bam"

    log:
        summary = "/home/marc/projects/rna_seq_workflow/reports/hisat2_align/{run_id}/{sample}_summary.log",
        metrics = "/home/marc/projects/rna_seq_workflow/reports/hisat2_align/{run_id}/{sample}_metrics.log"

    params:
        basename = "/home/marc/projects/rna_seq_workflow/results/hisat2_index/{run_id}/{run_id}",
        strandness = "RF"
    threads:
        config["threads"]

    conda:
        "../envs/hisat2.yaml"

    shell: """hisat2 -p {threads} -x {params.basename} -1 {input.r1} -2 {input.r2} --rna-strandness {params.strandness} --summary-file {log.summary} --met-file {log.metrics} | samtools sort -o {output.bam_file}"""
