#Rule for indexing the reference genome using HISAT2
rule hisat2_index:
    input:
        ref_genome
    output:
        index_file = working_dir+"/results/{aligner}/index/"+run_id+"/"+run_id+".1.ht2"
    params:
        index_dir = directory(working_dir+"/HISAT2/index/"+run_id),
        basename = working_dir+"/HISAT2/index/"+run_id+"/"+run_id
    threads:
        config["threads"]
    conda:
        f"{working_dir}/workflow/envs/hisat2.yaml"

    shell: """mkdir {params.index_dir} && hisat2-build -p {threads} {input} {params.basename}""" 