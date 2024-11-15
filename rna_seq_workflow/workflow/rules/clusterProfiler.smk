rule clusterProfiler:
    input:
        annotated_results = "results/functional_annotations/{run_id}/{contrast}/{run_id}_combined_annotations.csv"
        
    output:
    conda:
        "../envs/clusterProfiler.yaml"
    script:
        "../scripts/clusterProfiler.R"