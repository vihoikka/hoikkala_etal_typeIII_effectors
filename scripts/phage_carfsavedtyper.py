import pandas as pd
import os
import argparse
import re
from pprint import pprint
import sys

#This script concises phage carfsaved information
print("starting phage carfsaved typer")
parser = argparse.ArgumentParser()
parser.add_argument("-hmr", "--hmm_result", required=True)
parser.add_argument("-o", "--out", required=True)
parser.add_argument("-n", "--name", required=True)
parser.add_argument("-hmt", "--hmm_targets", required=True)
#parser.add_argument("-o", "--output", required = True)
args = parser.parse_args()
hmm_result = args.hmm_result
hmm_targets = args.hmm_targets #keeps track of what was in the database
name = args.name
out = args.out


effector_list = ["carf", "saved"]
effector_dict = {"carf":0, "saved":0}
effector_precise = {}
effector_dict_binary = {"carf":False, "saved":False}

#construct the precise effector dictionary using the original msa files
hmm_target_file = open(hmm_targets, "r").read().splitlines()
for line in hmm_target_file:
    filename = os.path.basename(line)
    effector = re.split("\.", filename)[0] + "." + re.split("\.", filename)[1] #ca3_nucc.fa -> ca3_nucc
    effector_precise[effector] = 0

#check if hmm result exists
if (os.stat(hmm_result).st_size == 0):
    print("No carfsaved info found. Exiting...")
    sys.exit()

hmm = pd.read_csv(hmm_result, delim_whitespace=True, header = 0) #header refers to row number
print("Shape at start: " + str(hmm.shape))
#print(hmm.loc["r1"]["Subject_accession_ID_version"])
print(hmm)
for index, row in hmm.iterrows(): #check each row in hmm result per carfsaved type
   # print(row)
    print("Shape of DF: " + str(hmm.shape))
    hit = row["target_name"]
    hit_lowercase = hit.lower()
    try:
        if "saved" in hit_lowercase:
            effector_dict["saved"] += 1
            effector_precise[hit] += 1
        elif "carf" in hit_lowercase:
            effector_dict["carf"] += 1
            effector_precise[hit] += 1
    except:
        print(row)
        print(hit)
        sys.exit()
    

#Transform numbers to boolean
for key, value in effector_dict_binary.items():
    if effector_dict[key] > 0: #if there are any with this cA effector
        effector_dict_binary[key] = True #then it's positive

results = {}
results.update(effector_dict_binary)
results.update(effector_precise)

results = pd.DataFrame([results])
results["phage"] = name

filename = out
print("Writing file " + filename)
results.to_csv(filename, index=True, sep = '\t', header = True)
