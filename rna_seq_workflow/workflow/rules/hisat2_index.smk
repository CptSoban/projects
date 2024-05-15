#Rule for indexing the reference genome using HISAT2

rule hisat2_index:
    input:
        config["ref_genome"]

    output:
        index_files = config["workflow_dir"]+"/results/hisat2_index/{run_id}/{run_id}.1.ht2"

    params:
        index_dir = directory(config["workflow_dir"]+"/results/hisat2_index/{run_id}"),
        basename = config["workflow_dir"]+"/results/hisat2_index/{run_id}/{run_id}",

    threads:
        config["threads"]

    conda:
        "../envs/hisat2.yaml"

    shell: """mkdir -p {params.index_dir} && hisat2-build -p {threads} {input} {params.basename}""" 