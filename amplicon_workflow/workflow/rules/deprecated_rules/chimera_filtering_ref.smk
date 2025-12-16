rule fq2fa:
    input:
        concat_filt_reads = "results/{run_id}/filtering/{barcode}_filt.fastq",
    output:
        concat_filt_reads_fa = temp("results/{run_id}/filtering/{barcode}_filt.fa")
    conda:
        "../envs/chimera_filtering.yaml"
    
    shell:"""
        vsearch --fastq_filter {input.concat_filt_reads} \
            --fastaout {output.concat_filt_reads_fa} \
            --fastq_qmax 93 \
            --fastq_maxee 1.0
        """

rule chimera_filtering:
    input:
        concat_filt_reads_fa = "results/{run_id}/filtering/{barcode}_filt.fa",
        chimera_db = config["chimera_db"]

    output:
        chim_filt_reads = temp("results/{run_id}/chimera_filtering/{barcode}_nochim.fa")

    log: "../reports/{run_id}/{barcode}_chimera_filtering.log"

    params: threads = config["threads"],

    conda:
        "../envs/chimera_filtering.yaml"

    shell:"""
        vsearch --uchime_ref {input.concat_filt_reads_fa} \
            --db {input.chimera_db} \
            --nonchimeras {output.chim_filt_reads} \
            --threads {params.threads} \
            --log {log}
        """

rule fa2fq:
    input:
        chim_filt_reads = "results/{run_id}/chimera_filtering/{barcode}_nochim.fa",
        concat_filt_reads = "results/{run_id}/filtering/{barcode}_filt.fastq"
    output:
        chim_filt_reads_fq = "results/{run_id}/chimera_filtering/{barcode}_nochim.fastq"
    conda:
        "../envs/chimera_filtering.yaml"
    
    shell:"""
        seqkit fa2fq -f {input.chim_filt_reads} \
        {input.concat_filt_reads}
        """
