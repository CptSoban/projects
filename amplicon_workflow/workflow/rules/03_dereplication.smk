rule dereplicate:
    input:
        filtered_reads = "results/{run_id}/filtering/{barcode}_filt.fastq.gz"
    output:
        derep_fasta = temp("results/{run_id}/dereplication/{barcode}_derep.fa"),

    log:
        derep_info = "reports/{run_id}/dereplication/{barcode}_derep.tsv"
    conda:
        "../envs/vsearch.yaml"
    shell: """
        vsearch --fastx_uniques {input.filtered_reads} \
            --fastaout {output.derep_fasta} \
            --sizeout \
            --relabel Derep \
            --tabbedout {log.derep_info}
    """