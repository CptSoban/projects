#Rule to run cutadapt on all paired-end reads
rule cutadapt:
    input:
        seq_r1 = seq_dir+"/{sample}_R1.fastq" if trimming == "y" else "",
        seq_r2 = seq_dir+"/{sample}_R2.fastq" if trimming == "y" else ""
    output:
        trimmed_r1 = working_dir+"/results/cutadapt/trim_{sample}_R1.fastq",
        trimmed_r2 = working_dir+"/results/cutadapt/trim_{sample}_R2.fastq"
    params:
        adapter_r1 = config["adapter_r1"],
        adapter_r2 = config["adapter_r2"],
        quality_cutoff = "-q "+config["quality_cutoff"] if config["quality_cutoff"] else ""
    threads:
        config["threads"]
    conda:
        f"{working_dir}/workflow/envs/cutadapt.yaml"

    shell:  """cutadapt -j 2 {params.quality_cutoff} -a {params.adapter_r1} -A {params.adapter_r2} -o {output.trimmed_r1} -p {output.trimmed_r2} {input.seq_r1} {input.seq_r2}"""