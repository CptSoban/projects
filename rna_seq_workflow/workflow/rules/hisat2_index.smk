#Rule for indexing the reference genome using HISAT2

rule hisat2_index:
    input:
        config["ref_genome"]

    output:
        index_files = "/home/marc/projects/rna_seq_workflow/results/hisat2_index/"+config["Run ID"]+".1.ht2"

    params:
        index_dir = directory("/home/marc/projects/rna_seq_workflow/results/hisat2_index/"+config["Run ID"]),
        basename = "/home/marc/projects/rna_seq_workflow/results/hisat2_index/"+config["Run ID"]

    threads:
        config["threads"]

    conda:
        "../envs/hisat2.yaml"

    shell: """mkdir {params.index_dir} && hisat2-build -p {threads} {input} {params.basename}""" 