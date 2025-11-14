#!/usr/bin/env python3
import sys
import re
from Bio import SeqIO

# Uso: script.py <input.fasta> <df.csv>
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <input.fasta> <df.csv>")
        sys.exit()

fasta = sys.argv[1]
csv = sys.argv[2]

# Read fasta and save in a dictionary
fasta_dict = {str(rec.id): str(rec.seq) for rec in SeqIO.parse(fasta, "fasta")}
frag_dict = {}

# Read csv and filter lines with predicted_class = 1
with open(csv) as f:
    next(f)
    for row in f:
        row = row.strip()
        (entry, kmer, site, predicted_class, probability, score) = row.split(",")
        entry = entry.split("_kmer")[0]
        if fasta_dict.get(entry) and int(predicted_class) == 1:
            frag_dict.setdefault(entry, []).append({
                "site": int(site),
                "probability": float(probability),
                "score": float(score)
            })

# Process fragments before and after cleavage sites
for entry, fragments in frag_dict.items():
    seq = fasta_dict[entry]
    fragments.sort(key=lambda x: x["site"])  # sort sites
    full_len = len(seq)

    for i, frag in enumerate(fragments):
        start = frag["site"]

        # Check if the first fragment of the site
        if i == 0:
            # antes do primeiro sítio
            frag_before = seq[:start+1]  # inclui o aminoácido do sítio
            print(f">{entry}\tSITE=[{start}]\tPROABILITY=[{frag['probability']}]\tSCORE=[{frag['score']}]")
            print(frag_before)

        # Fragment after site
        end = fragments[i+1]["site"] if i+1 < len(fragments) else full_len
        frag_after = seq[start+1:end]
        print(f">{entry}\tSITE=[{start}]\tPROBABILITY=[{frag['probability']}]\tSCORE=[{frag['score']}]")
        print(frag_after)
