#Rule for indexing the reference genome using HISAT2
rule hisat2_index:
    input:
        masked_genome = "results/{run_id}/repeatmasker/"+ASSEMBLY_FILE+".masked"

    output:
        index_files = "resources/{run_id}/hisat2_index/{run_id}.1.ht2"

    params:
        index_dir = directory("resources/{run_id}/hisat2_index"),
        basename = "resources/{run_id}/hisat2_index/{run_id}",
    
    threads: config["threads"]

    conda:
        "../envs/hisat2.yaml"

    shell: """mkdir -p {params.index_dir} && hisat2-build -p {threads} {input.masked_genome} {params.basename}""" 