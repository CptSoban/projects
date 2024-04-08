rule g_profiler:
    input:
        preranked_genes = "results/DEG_analysis/{run_id}/{run_id}_ranked_positive_genes.csv",
    output:
        enriched_functions = "results/GSEA/{run_id}/{run_id}_gsea_table.csv",
    conda:
        "../envs/g_profiler.yaml"
    script:
        "../scripts/g_profiler.R"