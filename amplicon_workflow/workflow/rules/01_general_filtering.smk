rule chopper:
    input:
        reads = config['demultiplexed_reads']+"/{barcode}.fastq.gz",
    output:
        chopper_filt_reads = temp("results/{run_id}/filtering/{barcode}_filt.fastq.gz"),

    params:
        min_length = config["min_length"],
        max_length = config["max_length"],
        min_quality = config["min_quality"],

    threads: config["threads"]
    
    conda:
        "../envs/general_filtering.yaml"

    shell: """
        chopper \
            -i {input.reads} \
            -q {params.min_quality} \
            --minlength {params.min_length} \
            --maxlength {params.max_length} \
            -t {threads} \
            | gzip > {output.chopper_filt_reads}
    """
