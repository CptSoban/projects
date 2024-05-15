# Rule to run FastQC on all samples

# Needs out_dir to be created before exectuion. Somehow does not work with mkdir?
rule fastqc:
    input: 
        seq = config["seq_dir"]+"/{sample}_{group}.fastq",
    
    output:
        html = config["workflow_dir"]+"/results/fastqc/{run_id}/{sample}_{group}_fastqc.html",
        zip_file = config["workflow_dir"]+"/results/fastqc/{run_id}/{sample}_{group}_fastqc.zip"
    
    params:
        out_dir = config["workflow_dir"]+"/results/fastqc/{run_id}",
    
    threads:
        config["threads"]
    
    conda:  "../envs/fastqc.yaml"

    shell:"""mkdir -p {params.out_dir} && fastqc -t {threads} {input.seq} --outdir {params.out_dir}"""
