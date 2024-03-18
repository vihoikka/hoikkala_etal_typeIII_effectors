import argparse
import pandas as pd
import subprocess
import sys
import os
from Bio import SeqIO

# scripts/cora_neighbour_analysis.py --dedd {input.dedd} --nrn {input.nrn} --samlyase {input.samlyase} --outdir {params.outdir}

parser = argparse.ArgumentParser()
#args for input and output
parser.add_argument('--dedd', type=str, metavar='input', required=True)
parser.add_argument('--dedd_neo', type=str, metavar='input', required=True)
parser.add_argument('--dedd_clost', type=str, metavar='input', required=True)
parser.add_argument('--nrn', type=str, metavar='input', required=True)
parser.add_argument('--samlyase', type=str, required=True)
parser.add_argument('--locus', type=str, required=True)
parser.add_argument('--outdir', type=str, metavar='outpath', required=True)

args = parser.parse_args()
dedd = args.dedd
dedd_neo = args.dedd_neo
dedd_clost = args.dedd_clost
nrn = args.nrn
samlyase = args.samlyase
outdir = args.outdir
locus = args.locus

output_df = pd.DataFrame(columns = ["locus", "dedd", "dedd_neo", "dedd_clost", "nrn", "samlyase"])
output_df.loc[0] = [locus, "False", "False", "False", "False", "False"]

if os.path.getsize(dedd) > 0:
    output_df.loc[0, "dedd"] = "True"

if os.path.getsize(dedd_neo) > 0:
    output_df.loc[0, "dedd_neo"] = "True"

if os.path.getsize(dedd_clost) > 0:
    output_df.loc[0, "dedd_clost"] = "True"

if os.path.getsize(nrn) > 0:
    output_df.loc[0, "nrn"] = "True"

if os.path.getsize(samlyase) > 0:
    output_df.loc[0, "samlyase"] = "True"

out_file = os.path.join(outdir, "neighbourhood_results.tsv")

output_df.to_csv(out_file, sep = "\t", index = False)