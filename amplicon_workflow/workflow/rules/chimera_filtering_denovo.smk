# rule fq2fa:
#     input:
#         trimmed_reads = "results/{run_id}/trimming/{barcode}_trim.fastq"
#     output:
#         trimmed_reads_fa = temp("results/{run_id}/trimming/{barcode}_trim.fa")
#     conda:
#         "../envs/chimera_filtering.yaml"
    
#     shell:"""
#         vsearch --fastq_filter {input.trimmed_reads} \
#             --fastaout {output.trimmed_reads_fa} \
                    
            # """


rule derep:
    input:
        trimmed_reads = "results/{run_id}/trimming/{barcode}_trim.fastq"
    output:
        derep_reads_fa = "results/{run_id}/trimming/{barcode}_trim_derep.fa"
    conda:
        "../envs/chimera_filtering.yaml"
    
    shell:"""
        vsearch --fastx_uniques {input.trimmed_reads} \
            --fastaout {output.derep_reads_fa} \
            --fastq_qmax 93 \
            --sizeout \
            --relabel Uniq
      """

rule chimera_filtering:
    input:
        derep_reads_fa = "results/{run_id}/trimming/{barcode}_trim_derep.fa",
    output:
        chim_filt_reads = "results/{run_id}/chimera_filtering/{barcode}_nochim.fa"
    
    log: "reports/{run_id}/{barcode}_chimera_filtering.log"

    params: threads = config["threads"],

    conda:
        "../envs/chimera_filtering.yaml"

    shell:"""
        vsearch --uchime_denovo {input.derep_reads_fa} \
            --nonchimeras {output.chim_filt_reads} \
            --fasta_width 0 \
            --threads {params.threads} \
            --log {log}
        """

rule fa2fq:
    input:
        chim_filt_reads = "results/{run_id}/chimera_filtering/{barcode}_nochim.fa",
        trimmed_reads = "results/{run_id}/trimming/{barcode}_trim.fastq"
    output:
        chim_filt_reads_fq = "results/{run_id}/chimera_filtering/{barcode}_nochim.fastq"

    conda:
        "../envs/chimera_filtering.yaml"
    
    shell:"""
        vsearch --search_exact {input.trimmed_reads} \
            --db {input.chim_filt_reads} \
            --strand plus \
            --matched {output.chim_filt_reads_fq}
        """

# rule gzip:
#     input:
#         chim_filt_reads_fq = "results/{run_id}/chimera_filtering/{barcode}_nochim.fastq"
#     output:
#         chim_filt_reads_fq_gz = "results/{run_id}/chimera_filtering/{barcode}_nochim.fastq.gz"

#     conda:
#         "../envs/chimera_filtering.yaml"
    
#     shell:"""
#         gzip {input.chim_filt_reads_fq} > {output.chim_filt_reads_fq_gz}
#         """
# rule fa2fq:
#     input:
#         chim_filt_reads = "results/{run_id}/chimera_filtering/{barcode}_nochim.fa",
#         trimmed_reads = "results/{run_id}/trimming/{barcode}_trim.fastq"
#     output:
#         chim_filt_reads_fq = "results/{run_id}/chimera_filtering/{barcode}_nochim.fastq.gz"
#     conda:
#         "../envs/chimera_filtering.yaml"
    
#     shell:"""
#         seqkit fa2fq -f {input.chim_filt_reads} \
#         {input.trimmed_reads} |\
#         gzip > {output.chim_filt_reads_fq}
#         """