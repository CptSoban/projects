rule emu_classify:
    input:
        clean_reads = "results/{run_id}/filtering/{barcode}_clean.fastq.gz",
    output:
        rel_abundance = "results/{run_id}/emu/{barcode}_rel-abundance.tsv",

    params:
        output_dir = directory("results/{run_id}/emu"),
        ref_db = config["reference_db"] if config["reference_db"] else "resources/emu_db",
        min_abundance = config["min_abundance"],
        seq_type = "lr:hq" if int(config["min_quality"]) >= 20 else "map-ont",

    threads: config["threads"]

    conda:
        "emu"

    shell:"""
        emu abundance {input.clean_reads} \
            --db {params.ref_db} \
            --min-abundance {params.min_abundance} \
            --N 75 \
            --output-basename {wildcards.barcode} \
            --output-dir {params.output_dir} \
            --keep-counts \
            --type {params.seq_type} \
            --threads {threads}
        """
    
rule emu_combine:
    input:
        rel_abundance = expand("results/{run_id}/emu/{barcode}_rel-abundance.tsv", run_id=RUN_ID, barcode=SAMPLES.barcode),
        
    output:
        combined_rel_abundance = "results/{run_id}/emu/emu-combined-"+config["taxonomic_rank"]+".tsv",
    
    params:
        output_dir = directory("results/{run_id}/emu"),
        taxonomic_rank = config["taxonomic_rank"],

    conda:
        "emu"

    shell:"""
        rm {params.output_dir}/*threshold-*.tsv ;\
        emu combine-outputs {params.output_dir} {params.taxonomic_rank}
        """