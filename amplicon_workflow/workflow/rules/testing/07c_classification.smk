rule train_classifier:
    input:
        fasta = config["reference_fasta"],
        tax = config["reference_taxonomy"],
    output:
        model = "resources/naive_bayes_classifier/model.joblib",
    conda:
        "../envs/naive_bayes_classifier.yaml"
    script:
        "../scripts/train_classifier.py"

rule classify_reads:
    input:
        model = "resources/naive_bayes_classifier/model.joblib",
        clean_reads = "results/{run_id}/filtering/{barcode}_clean.fastq.gz",
    output:
        classification = "results/{run_id}/naive_bayes/{barcode}_classification.tsv",
    conda:
        "../envs/naive_bayes_classifier.yaml"
    script:
        "../scripts/classify_reads.py"