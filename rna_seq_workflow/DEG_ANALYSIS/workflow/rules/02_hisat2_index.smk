#Rule for indexing the reference genome using HISAT2
rule hisat2_index:
    input:
        config["ref_genome"]

    output:
        index_files = "results/{run_id}/hisat2_index/{run_id}.1.ht2"

    params:
        index_dir = directory("results/{run_id}/hisat2_index"),
        basename = "results/{run_id}/hisat2_index/{run_id}",

    conda:
        "../envs/hisat2.yaml"

    shell: """mkdir -p {params.index_dir} && hisat2-build -p {threads} {input} {params.basename}""" 