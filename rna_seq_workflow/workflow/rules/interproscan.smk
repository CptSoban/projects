rule interproscan:
    input:
        extracted_protein_sequences = config["workflow_dir"]+"/results/DEG_analysis/{run_id}/{run_id}_protein.fasta",
        interproscan_dir = config["interproscan_dir"]

    output:
        interpro_gff = config["workflow_dir"]+"/results/functional_annotations/{run_id}/interproscan/{run_id}.gff3",

    params:
        interpro_output = config["workflow_dir"]+"/results/functional_annotations/{run_id}/interproscan/{run_id}",

    shell:  """{input.interproscan_dir}/interproscan.sh -cpu 4 -i {input.extracted_protein_sequences} -b {params.interpro_output}"""

