#!/usr/bin/env python3
import sys
import re
import pickle
import pandas as pd
import warnings
from Bio import SeqIO
warnings.filterwarnings("ignore", category=UserWarning)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <input.fasta>")
        sys.exit(1)

fasta = sys.argv[1]
entries = []
kmers = []
sites = []

window = 3
motifs = ["KR", "KK", "RR", "K", "R"]
motifs.sort(key=len, reverse=True)  # sort motifs by size, largerst first

for rec in SeqIO.parse(fasta, "fasta"):
    header = str(rec.id)
    seq = str(rec.seq)
    n = 0
    used_positions = set()  # positions already covered by a larger motif

    for motif in motifs:
        for match in re.finditer(motif, seq):
            start = match.start()
            end = match.end()
            positions = set(range(start, end))

            n += 1

            # ignore position if its covered by a larger motif
            if positions & used_positions:
                continue

            site = end - 1
            win_start = max(0, site - window)
            win_end = min(len(seq), site + window + 1)
            kmer = seq[win_start:win_end]
            
            if len(kmer) == 7:
                entries.append(f"{header}_kmer{n}")
                kmers.append(kmer)
                sites.append(site)
                

            # used positions
            used_positions.update(positions)

# Make final dataframe
df = pd.DataFrame({
    "ENTRY": entries,
    "K-MER": kmers,
    "SITE": sites
})


with open("model/encoder.pkl", "rb") as f:
    encoder = pickle.load(f)

with open("model/svm_model.pkl", "rb") as f:
    svm = pickle.load(f)

# Certifique-se de ter a coluna KMER_LIST
df["KMER_LIST"] = df["K-MER"].apply(list)

# One-hot-enconder
X_enc = encoder.transform(df["KMER_LIST"].tolist())

# Prediction
y_pred = svm.predict(X_enc)
probs = svm.predict_proba(X_enc)
scores = svm.decision_function(X_enc)

df["PREDICTED_CLASS"] = y_pred
df["PROBABILITY"] = probs[:,1]
df["SCORE"] = scores


df= df[["ENTRY", "K-MER", "SITE", "PREDICTED_CLASS", "PROBABILITY", "SCORE"]]


# Save results to csv
df.to_csv("results/predictions.csv", index=False)

counts = df["PREDICTED_CLASS"].value_counts()

num_class_0 = counts.get(0, 0)
num_class_1 = counts.get(1, 0)

print(f"N. of putative cleavage sites: {num_class_1}")
