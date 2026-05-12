import joblib
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.naive_bayes import MultinomialNB
from sklearn.pipeline import Pipeline

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

def read_taxonomy(file_path):
    taxonomy_dict = {}
    with open(file_path, 'r') as f:
        for line in f:
            k, v = line.rstrip().split('\t')
            taxonomy_dict[k] = v
    return taxonomy_dict

labels, sequences = read_fasta('resources/emu_db/sequences.fasta')
taxonomy_dict = read_taxonomy('resources/emu_db/taxonomy.tsv')
ids = [taxonomy_dict[label] for label in labels]

pipeline = Pipeline([
    ('vectorizer', CountVectorizer(
        analyzer='char', 
        ngram_range=(6, 6),
        lowercase=False
    )),
    ('classifier', MultinomialNB(alpha=0.001))
])

pipeline.fit(sequences, ids)
joblib.dump(pipeline, 'resources/emu_db/emu_classifier.joblib')
                
