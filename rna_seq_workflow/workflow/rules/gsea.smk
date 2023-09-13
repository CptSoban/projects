rule gsea:
    input:
        mod_deg_table = "/home/marc/projects/rna_seq_workflow/results/DEG_analysis/{run_id}_mod_deseq2_results.csv",
    output:
        gsea_results = "/home/marc/projects/rna_seq_workflow/results/GSEA/{run_id}_gsea_gobp.csv",
    params:
        seed = 27
    conda:
        "../envs/gsea.yaml"
    script:
        "../scripts/gsea.py"