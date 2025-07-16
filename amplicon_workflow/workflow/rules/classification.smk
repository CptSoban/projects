rule emu_classify:
    input:
        chim_filt_reads_fq = "results/{run_id}/chimera_filtering/{barcode}_nochim.fastq"
    
    output:
        rel_abundance = "results/{run_id}/emu/{barcode}_rel-abundance.tsv",
    
    params:
        output_dir = directory("results/{run_id}/emu"),
        ref_db = directory(config["reference_db_emu"])

    threads: config["threads"],

    conda:
        "../envs/classification.yaml"

    shell:"""
        emu abundance {input.chim_filt_reads_fq} \
            --db {params.ref_db} \
            --output-basename {wildcards.barcode} \
            --output-dir {params.output_dir} \
            --keep-counts \
            --type map-ont \
            --threads {threads} \
        """

rule emu_combine:
    input:
        rel_abundance = lambda wildcards: expand("results/{run_id}/emu/{barcode}_rel-abundance.tsv", barcode=wildcards.barcode),
    
    output:
        combined_rel_abundance = "results/{run_id}/emu/combined_rel-abundance.tsv",
    
    params:
        output_dir = directory("results/{run_id}/emu"),
        taxonomic_rank = config["taxonomic_rank"]

    conda:
        "../envs/classification.yaml"

    shell:"""
        emu combine-outputs {params.output_dir} {params.taxonomic_rank}
        """