rule fastqc:
    input: 
        seq = config["seq_dir"]+"/{sample}_{group}.fastq.gz",
    
    output:
        html = "results/{run_id}/fastqc/{sample}_{group}_fastqc.html",
        zip_file = "results/{run_id}/fastqc/{sample}_{group}_fastqc.zip"
    
    params:
        out_dir = "results/{run_id}/fastqc",
    
    conda:  "../envs/fastqc.yaml"

    shell:  """mkdir -p \
            {params.out_dir} && \
            fastqc -t {threads} {input.seq} \
            --outdir {params.out_dir}"""
