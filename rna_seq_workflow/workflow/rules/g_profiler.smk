rule g_profiler:
    input:
        preranked_genes = "/home/marc/projects/rna_seq_workflow/results/GSEA/{run_id}_prerank.csv",
    output:
        enriched_functions = "/home/marc/projects/rna_seq_workflow/results/GSEA/{run_id}_gsea_table.csv",
    conda:
        "../envs/g_profiler.yaml"
    script:
        "../scripts/g_profiler.R"