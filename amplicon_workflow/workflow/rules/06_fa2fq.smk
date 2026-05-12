rule fa2fq:
    input:
        chimera_filtered_clusters = "results/{run_id}/chimera_filtering/{barcode}_nochim.fasta",
        clustering_table = "reports/{run_id}/clustering/{barcode}_cluster_report.tsv",
        dereplication_table = "reports/{run_id}/dereplication/{barcode}_derep.tsv",
        chopper_filt_reads = "results/{run_id}/filtering/{barcode}_filt.fastq.gz"
    output:
        clean_reads = "results/{run_id}/filtering/{barcode}_clean.fastq.gz",
    
    conda:
        "../envs/fa2fq.yaml"
    script:
        "../scripts/fa2fq.py"