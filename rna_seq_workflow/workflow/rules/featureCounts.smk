rule featureCounts:
    input:
        gtf = config["gtf"],
        bam_files = expand(config["workflow_dir"]+"/results/hisat2_align/{run_id}/{sample}.sortedByCoord.out.bam", run_id=config["Run ID"], sample=set(samples.sample))

    output:
        counts_table = config["workflow_dir"]+"/results/feature_counts/{run_id}/{run_id}_counts.txt",
        
    params:
        strandness = config["Strandness"][1]
    
    threads: config["threads"]
    
    conda:
        "../envs/featureCounts.yaml"

    shell: """featureCounts -T {threads} -s {params.strandness} -p --countReadPairs --extraAttributes 'product' -a {input.gtf} -o {output.counts_table} {input.bam_files}"""