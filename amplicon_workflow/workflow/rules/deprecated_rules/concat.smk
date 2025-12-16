rule fastcat:
    output:
        concat_filt_reads = "results/{run_id}/filtering/{barcode}_filt.fastq",
        
    params:
        reports = directory("reports/{run_id}/{barcode}"),
        reads_dir = config["demultiplexed_reads"]+"/{barcode}",
        min_length = config["min_length"],
        max_length = config["max_length"],
        min_quality = config["min_quality"],
    
    conda:
        "../envs/concat_filt.yaml"

    shell: """
        fastcat \
            -a {params.min_length} \
            -b {params.max_length} \
            --min_qscore {params.min_quality} \
            --histograms {params.reports} \
            {params.reads_dir} \
            > {output.concat_filt_reads}
    """



# rule dorado_exists:
#     output: touch("command_available.txt")
#     params: cmd="dorado"
#     run:
#         import shutil
#         if shutil.which(params.cmd) is None:
#             raise ValueError(f"Command '{params.cmd}' not found in PATH. Please install dorado and make sure it is in your PATH.")
#         shell("touch {output}")

# rule seq_trim:
#     input:
#         concat_trimm_reads = "results/{run_id}/trimming/{barcode}_trim.fastq",
#         primer = config["primer_fasta"]
#     output:
#         trim_reads = "results/{run_id}/trimming/{barcode}_finaltrim.fastq.gz"
    
#     shell: """
#             dorado trim \
#             --sequencing-kit SQK-NBD114-96 \
#             --emit-fastq \
#             {input.concat_trimm_reads} |\
#             gzip > {output.trim_reads}
#             """


