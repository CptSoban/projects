rule chimera_filtering:
    input:
        centroids = "results/{run_id}/clustering/{barcode}_centroids.fasta",
    output:
        chimera_filtered = "results/{run_id}/chimera_filtering/{barcode}_nochim.fasta",

    conda:
        "../envs/vsearch.yaml"
    shell: """
        vsearch --uchime_denovo {input.centroids} \
            --nonchimeras {output.chimera_filtered} \
            --abskew 8 \
            --minh 0.5 \
            --sizein \
            --fasta_width 0 \
            --xsize
    """