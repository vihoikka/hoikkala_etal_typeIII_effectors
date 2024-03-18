import argparse
import os
import subprocess
import sys
import pandas as pd

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Parse command-line arguments for input, output and msa and faa. Also use_existing_alignment and existing_alignment_path
parser = argparse.ArgumentParser(
    description='This script takes a Cas10 GGDD domain sequence, extracts a region around it and makes a HMM profile for it.')
parser.add_argument('-i', '--input_table', type=str, metavar='input_table', required=True,
                    help='Table containing info on each locus and its Cas10')
parser.add_argument('-c', '--cas10_sequences', type=str, metavar='cas10s', required=True,
                    help='All Cas10 sequences in fasta format.')
parser.add_argument('-o', '--output', type=str, metavar='output', required=True,
                    help='The output file to write the HMM profile to.')
parser.add_argument('-m', '--msa', type=str, metavar='msa', required=True,
                    help='The output file to write the multiple sequence alignment to.')
parser.add_argument('-f', '--faa', type=str, metavar='faa', required=True,
                    help='The output file to write the amino acid sequences to.')
parser.add_argument('-mo', '--motif', type=str, metavar='faa', required=True,
                    help='The sequence to filter proteins by, e.g. GGDD or GGDE')
parser.add_argument('-e', '--existing_alignment_path', type=str, metavar='existing_alignment_path', required=False)
parser.add_argument('-x', '--use_existing_alignment', action='store_true', required=False)
args = parser.parse_args()

motif = str(args.motif)
cyclase_column = "Cas10_" + motif + "_coord"

if not args.use_existing_alignment:
  # Read in the table
  table = pd.read_csv(args.input_table, sep='\t', header=0)
  print(table)
  #only keep the locus and GGDD coordinates, removing entries that do not have a GGDD
  table = table[['Locus', cyclase_column]]
  table = table.dropna()
  print(table)
  #read the fasta file using seqIO
  cas10_objects = SeqIO.to_dict(SeqIO.parse(args.cas10_sequences, 'fasta'))
  print("Opened ggdd cas objects:")
  print(cas10_objects)
  cyclase_motifs = []

  #for each locus, extract the sequence defined by the GGDD coordinates
  cas10s_GGDD = {}
  for index, row in table.iterrows():
      print(row["Locus"])
      if row[cyclase_column] != "-":
        start = int(row[cyclase_column])
        if start > 100:
          print("adding " + motif + " motif for locus " + row["Locus"])
          cas10_object = cas10_objects[row["Locus"]]
          cyclase_motifs.append(SeqRecord(cas10_object.seq[start-50:start+50], id=row["Locus"], description=""))

  #Write the trimmed sequences to a fasta file
  SeqIO.write(cyclase_motifs, args.faa, 'fasta')
  print("starting muscle for ggdd")
  # Run muscle to align the sequences
  subprocess.run(['muscle', '-super5', args.faa, '-output', args.msa, '-threads', '40'])

  # Run hmmbuild to make the HMM profile
  subprocess.run(['hmmbuild', args.output, args.msa])
  subprocess.run(['hmmpress', args.output])

else:
  # Run hmmbuild to make the HMM profile using a premade alignment
  subprocess.run(['hmmbuild', args.output, args.existing_alignment_path])
  subprocess.run(['hmmpress', args.output])
#run hmmpress on the file
