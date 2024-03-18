from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-m", "--match", required=True)
parser.add_argument("-o", "--out", required=True)
parser.add_argument("-d", "--database", required=True)
args = parser.parse_args()

match = args.match
out = args.out
database = args.database

outlist  = []

for seq in SeqIO.parse(database,"fasta"):

    if match in seq.id:
        outlist.append(seq)

SeqIO.write(outlist, out,"fasta")
