checkpoint demultiplex:
    input:
        filtered_bam = config["pod5_dir"] + "duplex_basecalled/duplex_filtered.bam",
    output:
        output_dir = directory(config["pod5_dir"] + "duplex_basecalled/fastq/"),
    params:
        kit_name = config["kit_name"],

    threads: config["threads"],

    script:
        "../scripts/demultiplex.sh"