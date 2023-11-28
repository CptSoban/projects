rule g_profiler:
    input:
        preranked_genes = "/home/marc/projects/rna_seq_workflow/results/DEG_analysis/{run_id}/{run_id}_ranked_positive_genes.csv",
    output:
        enriched_functions = "/home/marc/projects/rna_seq_workflow/results/GSEA/{run_id}/{run_id}_gsea_table.csv",
    conda:
        "../envs/g_profiler.yaml"
    script:
        "../scripts/g_profiler.R"