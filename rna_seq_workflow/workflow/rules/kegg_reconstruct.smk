rule kegg_reconstruct:
    input:
        combined_annotations = "results/{run_id}/functional_annotations/{contrast}_combined_annotations.csv"

    output:
        annotate_id_with_ko = "results/{run_id}/functional_annotations/{contrast}_annotate_id_with_ko.csv",
        annotate_ko_with_log2fc = "results/{run_id}/functional_annotations/{contrast}_annotate_ko_with_log2fc.csv",

    conda:
        "../envs/kegg_request.yaml"
    
    script:
        "../scripts/kegg_reconstruct.py"