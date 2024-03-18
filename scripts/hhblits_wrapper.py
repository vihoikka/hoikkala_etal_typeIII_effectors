import os
import argparse
import pandas as pd
import requests
import json
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import subprocess
import concurrent.futures

parser = argparse.ArgumentParser(description='Handle arguments for script')
parser.add_argument('-i', '--input', help='Input file multifasta', required=True)
parser.add_argument('-o', '--output_basepath', help='Output dir', default='None')
parser.add_argument('-d', '--database', help='Database', default='None')

args = parser.parse_args()

print("Starting hhblits_wrapper.py with arguments: " + str(args))

#read the multifasta using biopython
records = list(SeqIO.parse(args.input, "fasta"))

individual_fastas_dir = args.output_basepath + "/fastas"
#create directory for individual fastas
if not os.path.exists(individual_fastas_dir):
    os.makedirs(individual_fastas_dir)

hhblits_dir = args.output_basepath + "/hhblits"
if not os.path.exists(hhblits_dir):
    os.makedirs(hhblits_dir)

#extract proteins from multifasta
for protein in records:
    protein_id = protein.id
    protein_seq = protein.seq
    seq_record = SeqRecord(Seq(str(protein_seq)), id=protein_id, description="")
    SeqIO.write(seq_record, individual_fastas_dir + "/" + protein_id + ".faa", "fasta")

faa_files = [os.path.join(individual_fastas_dir, file) for file in os.listdir(individual_fastas_dir) if file.endswith('.faa')]

print(faa_files)

def run_hhblits(file):
    print("Running hhblits on " + str(file))
    
    #get basename of file
    basename = os.path.basename(file)
    print("Base name " + str(file))
    outpath = args.output_basepath + "/hhblits/" + basename + ".tsv"
    cmd = ["hhblits", "-cpu", "1", "-i", file, "-blasttab", outpath, "-d", args.database]
    print("Running " + str(file) + " with command " + str(cmd))
   
    # Run the command
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    # Check for errors
    if result.returncode != 0:
        print(f"Error processing file {file}: {result.stderr}")
    else:
        print(f"Processed file {file}")

# You can adjust the max_workers parameter to match the number of cores you want to use.
with concurrent.futures.ThreadPoolExecutor(max_workers=40) as executor:
    executor.map(run_hhblits, faa_files)
