#!/usr/bin/env python3
import re
import sys
import pandas as pd
from Bio import SeqIO
#from sklearn.preprocessing import OneHotEnconder
#import ast

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} <input.fasta> <input.tsv> <input.clstr>")
        sys.exit(1)

fasta = sys.argv[1]
tsv = sys.argv[2]
clstr = sys.argv[3]

kmers = []
labels = []
window = 3
motifs = ["KR", "KK", "RR", "K", "R"]
motifs.sort(key=len, reverse=True)  # sort motifs by size, largerst first

site_dict = {}

#ENTRY   SIGNALP FRAGMENTS
#A0A088MIT0      MAFLKKSLFLVLFLGVVSLSFC, 22      [('EEEKREEHEEEKRDEEDAESLGKR', 23, 46), ('YGGLSPLRISKR', 47, 58), ('VPPGFTPFRSPAR', 59, 71), ('SISGLTPIRLSKR', 72, 84), ('VPPGFTPFRSPARR', 85, 98), ('ISEADPGFTPSFVVIKGLSPLRGKRR', 99, 124), ('PPGFSPFRVD', 125, 134)]

print(f"ENTRY\tCLUSTER\tK-MER_ID\tK-MER\tSITE\tSTATUS")

with open(tsv) as tsv:
    for line in tsv:
        if re.search("ENTRY", line): continue
        entry, signalp, fragments = line.split("\t")
        
        if entry not in site_dict:
            site_dict[entry] = {}
    
        # encontra cada tupla
        tuples = re.findall(r"\('([^']*)'.*?(\d+)\)", fragments)
        # 'tuples' agora é lista de (sequencia, ultimo numero)
        
        for seq, site in tuples:
            site = int(site)
            site = site - 1
            if seq.endswith("K") or seq.endswith("R"):
                site_dict[entry][site] = {seq}

#>Cluster 0
#0       303aa, >sp|P11006|MAGA_XENL... *
#1       23aa, >sp|C0HKN6|MAGR1_XEN... at 100.00%
#>Cluster 1
#0       294aa, >sp|Q7T3L1|BRK1A_BOM... *
#1       208aa, >sp|Q90WZ1|BRK1B_BOM... at 97.60%
#2       152aa, >sp|Q90W88|BRK1C_BOM... at 96.05%

cluster_id_dict = {}
with open(clstr) as clstr:
    cluster = 0
    for line in clstr:
        if re.search(">Cluster", line):
            cluster = line.split(">Cluster ")[1].strip()
        else:
            seq_id = line.split("|")[1].strip()
            cluster_id_dict[seq_id]=cluster

                
for rec in SeqIO.parse(fasta, "fasta"):
    header = str(rec.id).split("|")[1]
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
            # Get Cluster
            cluster = cluster_id_dict.get(header, "NA")
            # Get sites
            if site in site_dict.get(header, {}):
                status = True
            else:
                status = False
            # Check kmer's length and print    
            if len(kmer) == 7:
                print(f"{header}\t{cluster}\t{header}_kmer{n}\t{kmer}\t{site}\t{status}")

# used positions
used_positions.update(positions)

#            if site in site_dict.get(header, {}):
#                if len(kmer) == 7:
#                    print(f"{header}_kmer{n}\t{kmer}\t{site}\t{True}")
#            else:
#                if len(kmer) == 7:
#                    print(f"{header}_kmer{n}\t{kmer}\t{site}\t{False}")
            
            # used positions
#            used_positions.update(positions)
