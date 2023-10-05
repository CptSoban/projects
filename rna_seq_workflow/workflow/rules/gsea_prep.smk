rule gsea_prep:
    input:
        mod_deg_table = "/home/marc/projects/rna_seq_workflow/results/DEG_analysis/{run_id}_mod_deseq2_results.csv",

    output:
        preranked_genes = "/home/marc/projects/rna_seq_workflow/results/GSEA/{run_id}_prerank.csv",

    conda:
        "../envs/deseq2_prep.yaml"
        
    script:
        "../scripts/gsea_prep.py"