rule star_index:
    input:
        ref_genome
    output:
        directory(results_dir+"/STAR/index/"+run_id)
    threads:
        config["threads"]
    conda:
        "star.yaml"
    shell: """STAR --runMode genomeGenerate --genomeDir {output} --genomeFastaFiles {input} 