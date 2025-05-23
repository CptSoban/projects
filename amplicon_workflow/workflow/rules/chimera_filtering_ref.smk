rule fq2fa:
    input:
        trimmed_reads = "results/{run_id}/trimming/{barcode}_finaltrim.fastq.gz",
    output:
        trimmed_reads_fa = temp("results/{run_id}/dorado_trim/{barcode}_trim.fa")
    conda:
        "../envs/chimera_filtering.yaml"
    
    shell:"""
        vsearch --fastq_filter {input.trim_reads} \
            --fastaout {output.trimmed_reads_fa} \
            --fastq_qmax 93 \
            --fastq_maxee 1.0
        """

rule chimera_filtering:
    input:
        trimmed_reads_fa = "results/{run_id}/dorado_trim/{barcode}_trim.fa",
        chimera_db = config["chimera_db"]

    output:
        chim_filt_reads = "results/{run_id}/chimera_filtering/{barcode}_nochim.fa"
    
    log: "../reports/{run_id}/{barcode}_chimera_filtering.log"

    params: threads = config["threads"],

    conda:
        "../envs/chimera_filtering.yaml"

    shell:"""
        vsearch --uchime_ref {input.trimmed_reads_fa} \
            --db {input.chimera_db} \
            --nonchimeras {output.chim_filt_reads} \
            --threads {params.threads} \
            --log {log}
        """

rule fa2fq:
    input:
        chim_filt_reads = "results/{run_id}/chimera_filtering/{barcode}_nochim.fa",
        trimmed_reads = "results/{run_id}/dorado_trim/{barcode}_trim.fastq.gz"
    output:
        chim_filt_reads_fq = "results/{run_id}/chimera_filtering/{barcode}_nochim.fastq.gz"
    conda:
        "../envs/chimera_filtering.yaml"
    
    shell:"""
        seqkit fa2fq -f {input.chim_filt_reads} \
        {input.trimmed_reads} |\
        gzip > {output.chim_filt_reads_fq}
        """