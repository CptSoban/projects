rule clustering:
    input:
        derep_fasta = "results/{run_id}/dereplication/{barcode}_derep.fa"
    output:
        centroids = "results/{run_id}/clustering/{barcode}_centroids.fasta",
    log:
        cluster_report = "reports/{run_id}/clustering/{barcode}_cluster_report.tsv"
    params:
        cluster_threshold = config["cluster_threshold"]

    threads: config["threads"]
    conda:
        "../envs/vsearch.yaml"
        
    shell: """
        vsearch --cluster_size {input.derep_fasta} \
            --id {params.cluster_threshold} \
            --centroids {output.centroids} \
            --uc {log.cluster_report} \
            --sizein \
            --sizeout \
            --threads {threads} \
    """