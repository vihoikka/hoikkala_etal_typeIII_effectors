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
    description='This script takes a Cas10 HD domain sequence and makes a HMM profile for it.')
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
parser.add_argument('-e', '--existing_alignment_path', type=str, metavar='existing_alignment_path', required=False)
parser.add_argument('-x', '--use_existing_alignment', action='store_true', required=False)

args = parser.parse_args()

if not args.use_existing_alignment:
  # Read in the table
  table = pd.read_csv(args.input_table, sep='\t', header=0)

  #Extract locus IDs when Cas10_HD is True
  cas10s = table.loc[table['Cas10_HD'] == True, 'Locus'].tolist()

  # Extract the Cas10 sequences from the fasta file
  cas10_sequences = []
  for record in SeqIO.parse(args.cas10_sequences, 'fasta'):
      if record.id in cas10s:
          cas10_sequences.append(record)

  #Extract the amino acids between coordinates 10 to 35 from each sequence
  cas10_sequences_HD = [
    SeqRecord(r.seq[6:40], id=r.id, description="") for r in cas10_sequences  # offset index to match biological indexing
  ]

  #Write the trimmed sequences to a fasta file
  print("writing to " + args.faa)
  SeqIO.write(cas10_sequences_HD, args.faa, 'fasta')

  # Run muscle to align the sequences
  print("starting muscle")
  subprocess.run(['muscle', '-super5', args.faa, '-output', args.msa, '-threads', '40'])

  # Run hmmbuild to make the HMM profile
  print("starting hmmbuild")
  subprocess.run(['hmmbuild', args.output, args.msa])

  # Run hmmpress to make the HMM profile searchable
  print("starting hmmpress")
  subprocess.run(['hmmpress', args.output])

else:
  # Run hmmbuild to make the HMM profile using a premade alignment
  print("starting hmmbuild after else")
  subprocess.run(['hmmbuild', args.output, args.existing_alignment_path])
  print("starting hmmpress after else")
  subprocess.run(['hmmpress', args.output])