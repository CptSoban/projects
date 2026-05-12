rule dorado_exists:
    output: touch("command_available.txt")
    params: cmd="dorado"
    run:
        import shutil
        if shutil.which(params.cmd) is None:
            raise ValueError(f"Command '{params.cmd}' not found in PATH. Please install dorado and make sure it is in your PATH.")
        shell("touch {output}")


rule simplex_basecalling:
    output:
        basecalled_reads = config["pod5_dir"] + "/simplex_basecalled/simplex.bam",

    params:
        pod5_dir = config["pod5_dir"],
        min_quality = config["min_quality"],

    threads: config["threads"],

    shell:  """
            dorado basecaller sup {params.pod5_dir} \
                --min-qscore {params.min_quality} \
                --no-trim \
                > {output.basecalled_reads}
            """

checkpoint demultiplex:
    input:
        basecalled_reads = config["pod5_dir"] + "/simplex_basecalled/simplex.bam",
    output:
        output_dir = directory(config["pod5_dir"] + "/simplex_basecalled/fastq/"),
    params:
        kit_name = config["kit_name"],

    threads: config["threads"],

    script:
        "../scripts/demultiplex.sh"