import pandas as pd
import os
import argparse
import re
from pprint import pprint
import sys

#This script concises cA hmming information from a single proteome
print("starting typer")
parser = argparse.ArgumentParser()
parser.add_argument("-hmr", "--hmm_result", required=True)
parser.add_argument("-o", "--out", required=True)
parser.add_argument("-n", "--name", required=True)
#parser.add_argument("-o", "--output", required = True)
args = parser.parse_args()
hmm_result = args.hmm_result
name = args.name
out = args.out


#check if hmm result exists
if (os.stat(hmm_result).st_size == 0):
    print("No Cas10 info found. Exiting...")
    sys.exit()

results = {'Cas10': "False", 'Cas10_E-value': [], "Cas10_accession": []}


hmm = pd.read_csv(hmm_result, delim_whitespace=True, header = 0) #header refers to row number

counter = 0
for index, row in hmm.iterrows(): #check each row in hmm result per cA type
   # print(row)
    evalue = row["qlen"] #due to a weird offset, we need to refer to these columns to get the desired data
    query = row["tlen"] #due to a weird offset, we need to refer to these columns to get the desired data
    results["Cas10"] = "True"
    results["Cas10_E-value"].append(evalue)
    results["Cas10_accession"].append(query)
    counter += 1


results = pd.DataFrame([results])
results["host"] = name

filename = out
print("Writing file " + filename)
results.to_csv(filename, index=True, sep = '\t', header = True)