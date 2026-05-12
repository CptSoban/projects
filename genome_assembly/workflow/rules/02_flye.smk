rule flye:
    input:
        reads = expand(config["nanopore_reads_directory"]+"/{sample}.fastq", sample=config["Nanopore Samples"]),
"results/{run_id}/filtering/{barcode}_filt.fastq.gz"
    output:
        assembly = "results/{run_id}/flye/assembly.fasta",

    params:
        read_quality = "--nano-hq" if config["q20_plus_reads"] == "yes" else "--nano-raw",
        out_dir = directory("results/{run_id}/flye"),

    threads: config["threads"],

    conda:
        "../envs/flye.yaml"

    shell: """flye --{params.read_quality} {input.reads} \
            --out-dir {params.out_dir} \
            --threads {threads} """