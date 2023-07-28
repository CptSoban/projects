# Rule to run FastQC on all samples
rule fastqc:
    input: 
        seq = seq_dir+"/{sample}_{pos}.fastq"
    output:
        html = working_dir+"/workflow/report/{sample}_{pos}_fastqc.html",
        zip_file = working_dir+"/workflow/report/{sample}_{pos}_fastqc.zip"
    params:
        out_dir = f"{working_dir}/workflow/report"
    threads:
        config["threads"]
    conda:
        f"{working_dir}/workflow/envs/fastqc.yaml"

    shell:"""fastqc -t {threads} {input.seq} --outdir {params.out_dir}"""