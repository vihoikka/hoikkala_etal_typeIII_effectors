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
parser.add_argument("-hmt", "--hmm_targets", required=True)
#parser.add_argument("-o", "--output", required = True)
args = parser.parse_args()
hmm_result = args.hmm_result
hmm_targets = args.hmm_targets #keeps track of what was in the database
name = args.name
out = args.out


effector_list = ["ca3", "ca4", "ca5", "ca6", "corA"]
effector_dict = {"ca3":0, "ca4":0, "ca5":0, "ca6":0, "corA": 0}
effector_precise = {}
effector_dict_binary = {"ca3":False, "ca4":False, "ca5":False, "ca6":False, "corA":False}

#construct the precise effector dictionary using the original msa files
#hmm_target_file = open(hmm_targets, "r")
hmm_target_file = open(hmm_targets, "r").read().splitlines()
for line in hmm_target_file:
    filename = os.path.basename(line)
    effector = re.split("\.", filename)[0] #ca3_nucc.fa -> ca3_nucc
    effector = re.split("_", effector)[1] #ca3_nucc -> nucc
    effector_precise[effector] = 0

print(effector_precise)

#check if hmm result exists
if (os.stat(hmm_result).st_size == 0):
    print("No cA info found. Exiting...")
    sys.exit()

hmm = pd.read_csv(hmm_result, delim_whitespace=True, header = 0, index_col=False) #header refers to row number
print("Shape at start: " + str(hmm.shape))
#print(hmm.loc["r1"]["Subject_accession_ID_version"])
print(hmm)
for index, row in hmm.iterrows(): #check each row in hmm result per cA type
    print(row)
    hit = re.split('_', row["target_name"])
    print(hit)
    hit_type = hit[0] #returns ca3, ca4, ca5, ca6 or corA
    precise_hit = hit[1] #returns the hmm target, e.g. nucc, can1, can2, cora...
    print(hit_type)
    print(precise_hit)
    print("...Discovered cA type in hmm: " + str(hit_type) + " (" + precise_hit + ")")
    effector_dict[hit_type] += 1 #add to corresponding cA in the dictionary
    effector_precise[precise_hit] += 1

#Transform numbers to boolean
for key, value in effector_dict_binary.items():
    if effector_dict[key] > 0: #if there are any with this cA effector
        effector_dict_binary[key] = True #then it's positive

results = {}
results.update(effector_dict_binary)
results.update(effector_precise)

results = pd.DataFrame([results])
results["host"] = name

print("Results:")
print(results)

filename = out
print("Writing file " + filename)
results.to_csv(filename, index=True, sep = '\t', header = True)
