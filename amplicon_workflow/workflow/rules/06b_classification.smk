rule k2_build_silva:
    output:
        rel_abundance = "results/{run_id}/emu/{barcode}_rel-abundance.tsv",

    params:
        silva_dir = "resources/kraken2_db/silva",

    threads: config["threads"]

    conda:
        "../envs/kraken2.yaml"

    shell:"""
        k2 build \
            --db {params.silva_dir} \
            --special silva \
            --threads {threads}
        """
    
rule k2_classify:
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