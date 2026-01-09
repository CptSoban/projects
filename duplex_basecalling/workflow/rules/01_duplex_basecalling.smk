rule dorado_exists:
    output: touch("command_available.txt")
    params: cmd="dorado"
    run:
        import shutil
        if shutil.which(params.cmd) is None:
            raise ValueError(f"Command '{params.cmd}' not found in PATH. Please install dorado and make sure it is in your PATH.")
        shell("touch {output}")


rule duplex_basecalling:
    output:
        basecalled_reads = config["pod5_dir"] + "duplex_basecalled/duplex.bam",

    params:
        pod5_dir = config["pod5_dir"],

    threads: config["threads"],

    shell:"""
        dorado duplex sup -t {threads} {params.pod5_dir} > {output.basecalled_reads}
        """

rule filter_bam:
# Removes simplex reads with an duplex equivalent
    input:
        basecalled_reads = config["pod5_dir"] + "duplex_basecalled/duplex.bam",

    output:
        filtered_bam = config["pod5_dir"] + "duplex_basecalled/duplex_filtered.bam",
        tag_file = temp("tag.txt"),

    shell:"""
        echo "Cleaning up duplex basecalled reads..." &&\
        printf "1\n0\n" > {output.tag_file} &&\
        samtools view -b -h -D dx:{output.tag_file} {input.basecalled_reads} > {output.filtered_bam}
        """
