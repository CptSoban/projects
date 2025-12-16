rule nanoplot:
    input:
        demux_input,
    output:
        output_dir = directory(config["pod5_dir"] + "duplex_basecalled/fastq/{barcode}_nanoplot"),
        #stats = config["pod5_dir"] + "duplex_basecalled/fastq/{barcode}_nanoplot/NanoStats.txt",

    threads: config["threads"],
    conda:
        "../envs/nanoplot.yaml"
    shell: """
        NanoPlot \
            --fastq {input[0]} \
            --outdir {output.output_dir} \
            -t {threads} \
            --plots dot \
            --legacy hex
        """