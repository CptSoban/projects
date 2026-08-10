rule deeploc:
    input:
        mod_deg_table = "results/{run_id}/DEG_analysis/trimmed/trim_{contrast}.csv",
        ref_protein_fasta = config["protein_fasta"],
        deeploc_model = config["deeploc_model"]

    output:
        main_results = "results/{run_id}/functional_annotation/deeploc/{contrast}_deeploc_results.csv"

    conda:
        "deeploc"
    
    shell:
        """
        deeploc2 -f {input.ref_protein_fasta} -m Accurate -o {output.main_results}
        """
