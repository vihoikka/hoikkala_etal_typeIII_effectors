import argparse
import pandas as pd
import subprocess
import sys
import os
from Bio import SeqIO, SearchIO
from hmmer_to_pandas import parse_hmmer_domtblout

argparser = argparse.ArgumentParser()
#args for input and output
argparser.add_argument('--input', type=str, metavar='input', required=True,
                    help='Multifasta of an effector')
argparser.add_argument('--output_basepath', type=str, required=True, metavar='output_basepath')
argparser.add_argument('--effector', type=str, required=True, metavar='effector name')


args = argparser.parse_args()

input_effector_fasta = args.input
output_basepath = args.output_basepath
effector = args.effector

E_value_pfam = str(1e-2) #this is for the full sequence, not domain-specific
pfams_hmm_db = "/media/volume/st_andrews/databases/pfam/oct23/Pfam-A.hmm"

#read in the fasta file
effector_sequences = SeqIO.parse(open(input_effector_fasta),'fasta')

effector_df = pd.DataFrame(columns=["ID",
                                    "length",
                                    ])


#each row should be a domain that is found in an effector

#run hmmscan against pfam database
hmmout = output_basepath + "/" + effector + ".hmm"
with open(os.devnull, 'w') as DEVNULL:
    subprocess.run(['hmmscan', '--domtblout', hmmout, '--cpu', '40', '-E', E_value_pfam, pfams_hmm_db, input_effector_fasta], stdout=DEVNULL, stderr=subprocess.STDOUT)