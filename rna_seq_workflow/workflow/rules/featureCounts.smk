if config["Strandness"] == "Forward":
    strandness = 1
elif config["Strandness"] == "Reverse":
    strandness = 2
elif config["Strandness"] == "Unstranded":
    strandness = 0

rule featureCounts:
    input:
        gtf = config["gtf"],
        bam_files = expand("/home/marc/projects/rna_seq_workflow/results/{aligner}_align/{sample}.sortedByCoord.out.bam", aligner=config["Aligner"], sample=set(samples.sample))

    output:
        counts_table = "/home/marc/projects/rna_seq_workflow/results/feature_counts/{run_id}_counts.txt",
        
    params:
        strandness = 2
    
    threads: config["threads"]
    
    conda:
        "../envs/featureCounts.yaml"

    shell: """featureCounts -T {threads} -s {params.strandness} -p -a {input.gtf} -o {output.counts_table} {input.bam_files}"""