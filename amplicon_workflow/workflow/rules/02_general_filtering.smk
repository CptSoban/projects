# rule cutadapt:
#     input:
#         reads = demux_filename,
#     output:
#         cutadapt_filt_reads = temp("results/{run_id}/filtering/{barcode}_trim.fastq"),

#     params:
#         forward_primer = config["Forward_primer"],
#         reverse_primer = config["Reverse_primer"],

#     threads: config["threads"]
    
#     conda:
#         "../envs/general_filtering.yaml"

#     shell: """
#         cutadapt \
#             -g {params.forward_primer} \
#             -a {params.reverse_primer} \
#             -o {output.cutadapt_filt_reads} \
#             {input.reads} \
#             -e 0.05 \
#             -rc \
#             -m 10 \
#             -O 10 \
#             --cores={threads}
#     """

rule chopper:
    input:
        cutadapt_filt_reads = demux_filename, #"results/{run_id}/filtering/{barcode}_trim.fastq",
    output:
        chopper_filt_reads = temp("results/{run_id}/filtering/{barcode}_filt.fastq.gz"),

    params:
        min_length = config["min_length"],
        max_length = config["max_length"],
        min_quality = config["min_quality"],
        contam = "resources/DCS.fasta",

    threads: config["threads"]
    
    conda:
        "../envs/general_filtering.yaml"

    shell: """
        chopper \
            -i {input.cutadapt_filt_reads} \
            -q {params.min_quality} \
            --minlength {params.min_length} \
            --maxlength {params.max_length} \
            -t {threads} \
            | gzip > {output.chopper_filt_reads}
    """

     #           --contam {params.contam} \