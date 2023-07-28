rule gsea:
    input:
        mod_deg_table = f"{results_dir}/DEG_analysis/{run_id}/{run_id}_mod_deseq2_results.csv",
    output:
        gsea_results = f"{results_dir}/DEG_analysis/gsea/{run_id}/{run_id}_gsea_gobp.csv",
    params:
        seed = 27
    conda:
        "gsea.yaml"
    script:
        "gsea.py"