import argparse
import os
 
'''
Takes a cd-hit cluster file and outputs names of cluster representatives
'''


parser = argparse.ArgumentParser(description='Process some inputs.')
parser.add_argument('--input', type=str, help='Input', required=True)
parser.add_argument('--output', type=str, help='Output (list of accessions)', required=True)

args = parser.parse_args()

in_file = args.input
out_file = args.output

acc_list = []

with open(in_file, "r") as f:
    for line in f:
        acc = line.split(">")[1].split("...")[0]
        acc_list.append(acc)

with open(out_file, "w") as outfile:
    for line in acc_list:
        outfile.write(line + "\n")