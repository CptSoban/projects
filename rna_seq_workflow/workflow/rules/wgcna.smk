rule wgcna:
    input:
        normalized_counts = "/home/marc/projects/rna_seq_workflow/results/DEG_analysis/{run_id}_normalized_counts.csv",    
        
    output:
        histogram = "/home/marc/projects/rna_seq_workflow/results/WGCNA/{run_id}/final-kme-histogram.pdf",
        temp_file = "/home/marc/projects/rna_seq_workflow/results/WGCNA/{run_id}/temp.tsv"

    params:
        outdir = "/home/marc/projects/rna_seq_workflow/results/WGCNA/{run_id}"
    conda:
        "../envs/wgcna.yaml"
    
    shell:  """sed "s/,/\t/g" {input.normalized_counts} > {output.temp_file}| iterativeWGCNA -i {output.temp_file} -o {params.outdir} | rm {output.temp_file}"""