rule gsea_prep:
    input:
        mod_deg_table = "results/DEG_analysis/{run_id}/{run_id}_mod_deseq2_results.csv",

    output:
        preranked_genes = "results/DEG_analysis/{run_id}/{run_id}_ranked_positive_genes.csv",

    conda:
        "../envs/deseq2_prep.yaml"
        
    script:
        "../scripts/gsea_prep.py"