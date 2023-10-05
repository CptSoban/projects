rule gsea:
    input:
        preranked_genes = "/home/marc/projects/rna_seq_workflow/results/GSEA/{run_id}_prerank_genes.csv",

    output:
        gsea_results = "/home/marc/projects/rna_seq_workflow/results/GSEA/{run_id}_gsea_{gene_set}.csv",

    params:
        seed = 27,
        gene_set = config["GSEA gene set"]

    conda:
        "../envs/gsea.yaml"

    script:
        "../scripts/gsea.py"