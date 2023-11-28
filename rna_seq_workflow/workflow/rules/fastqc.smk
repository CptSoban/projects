# Rule to run FastQC on all samples
rule fastqc:
    input: 
        seq = config["seq_dir"]+"/{sample}_{group}.fastq",
    
    output:
        html = "/home/marc/projects/rna_seq_workflow/results/fastqc/{run_id}/{sample}_{group}fastqc.html",
        zip_file = "/home/marc/projects/rna_seq_workflow/results/fastqc/{run_id}/{sample}_{group}fastqc.zip"
    
    params:
        out_dir = "/home/marc/projects/rna_seq_workflow/results/{run_id}/fastqc",
    
    threads:
        config["threads"]
    
    conda:  "../envs/fastqc.yaml"

    shell:"""fastqc -t {threads} {input.seq} --outdir {params.out_dir}"""
