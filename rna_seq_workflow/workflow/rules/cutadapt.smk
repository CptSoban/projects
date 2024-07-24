#Rule to run cutadapt on all paired-end reads

rule cutadapt:
    input:
        r1 = config["seq_dir"]+"/{sample}_R1.fastq",
        r2 = config["seq_dir"]+"/{sample}_R2.fastq"

    output:
        trimmed_r1 = "results/cutadapt/{run_id}/trim_{sample}_R1.fastq",
        trimmed_r2 = "results/cutadapt/{run_id}/trim_{sample}_R2.fastq"

    params:
        adapter_r1 = config["adapter_r1"],
        adapter_r2 = config["adapter_r2"],
        quality_cutoff = "-q "+config["quality_cutoff"] if config["quality_cutoff"] else ""

    threads:
        config["threads"]

    conda:
        "../envs/cutadapt.yaml"

    shell: """ cutadapt \
            -j 2 {params.quality_cutoff} \
            -a {params.adapter_r1} \
            -A {params.adapter_r2} \
            -o {output.trimmed_r1} \
            -p {output.trimmed_r2} \
            {input.r1} {input.r2}"""
       