import joblib
import numpy as np

CONFIDENCE_THRESHOLD = 0.6

def read_fasta(file_path):
    sequences = []
    labels = []
    with open(file_path, 'r') as f:
        seq = []
        for line in f:
            if line.startswith('>'):
                if seq:
                    sequences.append(''.join(seq))
                    seq = []
                labels.append(line[1:].strip())
            else:
                seq.append(line.strip()).upper()
        if seq:
            sequences.append(''.join(seq))
    return sequences, labels

clf = joblib.load('resources/emu_db/emu_classifier.joblib')
labels, sequences = read_fasta('amplicon_workflow/workflow/scripts/classify_reads_input.fasta')

predictions = clf.predict(sequences)
probabilities = clf.predict_proba(sequences)
conf = probabilities.max(axis=1)

with open("taxonomic_classification.tsv", "w") as out_f:
    out_f.write("Read_ID\tTaxonomy\tConfidence\n")
    for i, t, c in zip(labels, predictions, conf):
        if c < CONFIDENCE_THRESHOLD:
            t = "Unclassified"
        out_f.write(f"{i}\t{t}\t{c:.4f}\n")