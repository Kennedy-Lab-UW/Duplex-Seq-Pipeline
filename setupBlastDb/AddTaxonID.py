from Bio import SeqIO
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("inFasta")
parser.add_argument("inTaxonID")
parser.add_argument("outFasta")
o = parser.parse_args()
# open output fasta
Fout = open(o.outFasta, 'w')
# load fasta and parse it
for seqRec in SeqIO.parse(o.inFasta, "fasta"):
    Fout.write(
        f">{o.inTaxonID}|{seqRec.id}\n"
        )
    for i in range(0, len(seqRec.seq), 50):
        Fout.write(
            f"{str(seqRec.seq[i:i+50])}\n"
            )
Fout.close()
