import pandas as pd
import os
import argparse
import re
from pprint import pprint
import sys

#This script concises cA blasting information from a single proteome

parser = argparse.ArgumentParser()
parser.add_argument("-b", "--blast_result", required=True)
parser.add_argument("-o", "--out", required=True)
parser.add_argument("-n", "--name", required=True)
#parser.add_argument("-o", "--output", required = True)
args = parser.parse_args()
blast_result = args.blast_result
name = args.name
out = args.out

effector_list = ["cA3", "cA4", "cA5", "cA6"]
effector_dict = {"cA3":0, "cA4":0, "cA5":0, "cA6":0}
effector_dict_binary = {"cA3":False, "cA4":False, "cA5":False, "cA6":False}


#check if blast result exists
if (os.stat(blast_result).st_size == 0):
    print("No cA info found. Exiting...")
    sys.exit()

blast = pd.read_csv(blast_result, sep = '\t', header = 0) #header refers to row number
print("Shape at start: " + str(blast.shape))
#print(blast.loc["r1"]["Subject_accession_ID_version"])
#print(blast)
for index, row in blast.iterrows(): #check each row in blast result per cA type
   # print(row)
    hit_type = re.split('@', row["sseqid"])[0] #match before @ (resulting in e.g. cA4)
    print("...Discovered cA type in blast: " + str(hit_type))
    effector_dict[hit_type] += 1 #add to corresponding cA in the dictionary

#Transform numbers to boolean
for key, value in effector_dict_binary.items():
    if effector_dict[key] > 0: #if there are any with this cA effector
        effector_dict_binary[key] = True #then it's positive

condensed = pd.DataFrame([effector_dict_binary])
condensed["host"] = name

print("Results:")
print(condensed)

filename = out
print("Writing file " + filename)
condensed.to_csv(filename, index=True, sep = '\t', header = True)
