rule eggnog:
    input:
        cleaned_aa = "/home/nanopore/projects/rna_seq_workflow/FUNCTIONAL_ANNOTATION/results/20260603_CLORO/cleaned.aa",
        interpro_results_db = "resources/{run_id}/interpro_results_db"
    output:
        eggnog_results ="results/{run_id}/eggnog/{run_id}.emapper.annotations",
    params:
        eggnog_dir = directory("results/{run_id}/eggnog"),
        eggnog_data_dir = directory("resources/"),
    conda:
        "eggnog"
    threads:
        config["threads"]
    shell:
        """emapper.py --data_dir {params.eggnog_data_dir} -m diamond --cpu {threads} --output_dir {params.eggnog_dir} --output {wildcards.run_id} -i {input.cleaned_aa} --override --dmnd_ignore_warnings --evalue 0.01 --score 60 --pident 40 --query_cover 20 --subject_cover 20 --itype proteins --tax_scope auto --target_orthologs all --go_evidence non-electronic --pfam_realign none --report_orthologs --decorate_gff yes """
    