rule minimap2_align:
    input:
        fasta = "results/{run_id}/clustering/{barcode}_centroids.fasta",
        reads = "results/{run_id}/chimera_filtering/{barcode}_nochim.fastq"
    output:
        sam = temp("results/{run_id}/alignment/{barcode}_align.sam")
    
    params:
        threads = config["threads"]
    
    conda:
        "../envs/minimap2.yaml"
    
    shell: """
        minimap2 -ax map-ont -t {params.threads} {input.fasta} {input.reads} > {output.sam}
    """