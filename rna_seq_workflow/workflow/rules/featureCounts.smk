rule featureCounts:
    input:
        gtf = config["gtf"],
        bam_files = expand("results/hisat2_align/{run_id}/{sample}.sortedByCoord.out.bam", run_id=RUN_ID, sample=set(SAMPLES.sample))
    
    output:
        counts_table = "results/feature_counts/{run_id}/{run_id}_counts.txt",
        
    params:
        strandness = config["Strandness"][1]
    
    threads: config["threads"]
    
    conda:
        "../envs/featureCounts.yaml"

    shell:  """ featureCounts \
            -T {threads}  \
            -s {params.strandness} \
            -p --countReadPairs \
            -a {input.gtf} \
            -o {output.counts_table} \
            {input.bam_files}"""
