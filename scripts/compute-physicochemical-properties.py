#!/usr/bin/env python3

from Bio.SeqUtils.ProtParam import ProteinAnalysis
from modlamp.descriptors import GlobalDescriptor
from Bio import SeqIO
import sys
import re

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <input.fasta>")
        sys.exit(1)
    else:
        fasta_file = sys.argv[1]

print("ID\tSequence\tLen\tMolecular_Weight\tNet_Charge\tIsoelectric_Point\tHydrophobicity\tAmphipathicity\tBoman_index\tHelix\tSheet\tCoil")
for rec in SeqIO.parse(fasta_file, "fasta"):
    seq = str(rec.seq)
    header = rec.id
    description = rec.description
    if re.search(r"\[SIGNALP\]", description): continue
    
    # ProteinAnalysis
    analysis = ProteinAnalysis(seq)

    # modlamp descriptors
    desc = GlobalDescriptor([seq])
    desc.calculate_all()

    hydrophobicity = desc.descriptor[0][0]
    amphipathicity = desc.descriptor[0][4]
    mol_weight = desc.descriptor[0][1]
    net_charge = desc.descriptor[0][2]
    boman_index = desc.descriptor[0][3]

    ss = analysis.secondary_structure_fraction()

    # Ponto isoelétrico
    pI = analysis.isoelectric_point()

    print(f"{header}\t{seq}\t{len(seq)}\t{mol_weight:.2f}\t{net_charge:.2f}\t{pI:.2f}\t"
          f"{hydrophobicity:.2f}\t{amphipathicity:.2f}\t{boman_index:.2f}\t{ss[0]:.2f}\t{ss[2]:.2f}\t{ss[1]:.2f}")
