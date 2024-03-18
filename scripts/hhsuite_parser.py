import os
import argparse
import pandas as pd
import requests
import json


def fetch_PDB_description(ID):
    pdb_id = ID.upper().split("_")[0]
    url = f'https://data.rcsb.org/rest/v1/core/entry/{pdb_id}'
    response = requests.get(url)
    data = response.json()
    try:
        data = data["struct"]["title"]
    except:
        data = "No description found for " + pdb_id
    #print(pdb_id)
    #print(data)
    return data

import chardet

def guessEncoding(file):
    with open(file, 'rb') as rawdata:
        result = chardet.detect(rawdata.read(100000))["encoding"]
    print("Encoding: " + str(result))
    return result


parser = argparse.ArgumentParser(description='Handle arguments for script')
parser.add_argument('-i', '--infile', help='Input file (*.tsv)', required=True)
parser.add_argument('-o', '--outfile', help='Output file', default='None')
parser.add_argument('-d', '--database', help='Output file', default='None')
parser.add_argument('-m', '--mapping', help="Tab separated file for mapping descriptions to IDs", required = False)

args = parser.parse_args()

df = pd.read_csv(args.infile, sep='\t', header = None)

headers = ['query', 'target', 'matches_divided_by_tLen', 'alignment_length', 'mismatch', 'gapOpen', 'qstart', 'qend', 'tstart', 'tend', 'eval', 'score']

#add the headers to df
df.columns = headers
print(df)

#remove entries where eval is less than 0.001
df = df[df['eval'] < 0.001]

if args.database == "PDB":
    #create new column in df "description"
    df['description'] = "-"
    #iterate over df and fetch description for each target from PDB
    for index, row in df.iterrows():
        df.loc[index, 'description'] = fetch_PDB_description(row['target'])
    df.to_csv(args.outfile, sep="\t", index=False)
    exit()

if args.mapping:
    print(args.mapping)
    encoding_guess = guessEncoding(args.mapping)
    map_df = pd.read_csv(str(args.mapping), sep = "\t", header = None, encoding=encoding_guess)
    map_df.columns = ["id", "functional_category", "description", "associated_gene", "associated_pathway", "pubmed_id", "pdb_id"]
    merged = pd.merge(df, map_df, left_on="target", right_on="id")
    merged.to_csv(args.outfile, sep="\t", index=False)
    exit()