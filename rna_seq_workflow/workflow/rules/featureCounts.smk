rule featureCounts:
    input:
        gtf = config["gtf"],
        bam_files = expand(results_dir+"/{aligner}/alignment/"+run_id+"/{sample}.sortedByCoord.out.bam", aligner=Alignment_tool, sample=set(samples.sample))
    output:
        counts_table = f"{results_dir}/feature_counts/{run_id}/{run_id}_counts.txt"

    threads: config["threads"]
    conda:
        "subread.yaml"

    shell: """featureCounts -T {threads} -p --extraAttributes "gene" -a {input.gtf} -o {output.counts_table} {input.bam_files}"""