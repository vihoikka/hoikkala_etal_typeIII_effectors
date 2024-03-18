import argparse
import pandas as pd
import subprocess
import sys
import os
from Bio import SeqIO

parser = argparse.ArgumentParser()
#args for input and output
parser.add_argument('-ib', '--input_blast', type=str, metavar='input', required=True,
                    help='Blast results of unknown proteins against type III loci proteins')
parser.add_argument('-ip', '--input_pfam_annotations', type=str, metavar='input', required=True,
                    help='input_pfam_annotations')
parser.add_argument('-o', '--output', type=str, required=True,
                    help='Results')
parser.add_argument('-op', '--outpath', type=str, metavar='outpath', required=True,
                    help='plots etc')

args = parser.parse_args()
input_blast = args.input_blast
output = args.output
outpath = args.outpath
input_pfam_annotations = args.input_pfam_annotations

#read in the blast results
blast_results = pd.read_csv(input_blast, sep='\t', header=0)
print(blast_results)

#read in the pfam annotations. This is difficult because of the space in the last column
# Read the file as a single string column skipping the first row, using a fake separator
pfam_annotations = pd.read_csv(input_pfam_annotations, sep='\a', skiprows=1, names=['raw'], engine='python')

# Split only the first 18 columns, leaving the rest as a single column. 
# This works under the assumption that there are 18 fields before 'description_of_target'
pfam_annotations = pfam_annotations['raw'].str.split('\s+', n=18, expand=True)

# Assign column names manually
pfam_annotations.columns = ['target_name', 'accession', 'query_name', 'accession1', 'E-value1', 'score1', 'bias1', 'E-value2', 'score2', \
              'bias2', 'exp', 'reg', 'clu', 'ov', 'env', 'dom', 'rep', 'inc', 'description_of_target']

print(pfam_annotations)

#get unique values in query_name column
unique_pfam_proteins = pfam_annotations['query_name'].unique()
print(unique_pfam_proteins)

#for each unique pfam protein, concatenate values from column target_name in the df pfam_annotations
#into a dictionary where key is the pfam protein and value is the list of target_names as string where each value is separated by a comma
pfam_protein_functions_dict = {}
for pfam_protein in unique_pfam_proteins:
    print(pfam_protein)
    pfam_protein_functions_dict[pfam_protein] = pfam_annotations[pfam_annotations['query_name'] == pfam_protein]['description_of_target'].str.cat(sep=',')

print(pfam_protein_functions_dict)

#transform the dicitinary into a df where column 0 is the key and column 1 is the list of hits separated by comma
pfam_protein_functions_df = pd.DataFrame.from_dict(pfam_protein_functions_dict, orient='index')
print(pfam_protein_functions_df)

#obtain a list of unique values in the locus column
loci = blast_results['locus'].unique()

#obtain a list of unique values in the query column
queries = blast_results['qseqid'].unique()
print(queries)

#turn the queries list into a df
queries_df = pd.DataFrame(queries, columns=['qseqid'])

#merge that with pfam_protein_functions_df
pfam_protein_functions_df = pd.merge(queries_df, pfam_protein_functions_df, left_on='qseqid', right_index=True, how='left')
#rename the column 0 to "pfam_prediction"
pfam_protein_functions_df.rename(columns={0: "pfam_prediction"}, inplace=True)

print(pfam_protein_functions_df)

#calculate the number of times each query appears in the blast results
query_counts = blast_results['qseqid'].value_counts()

#order the queries by the number of times they appear in the blast results
ordered_queries = query_counts.index.tolist()

print("There are " + str(len(loci)) + " loci in the blast results.")
print("There are " + str(len(queries)) + " unique queries in the blast results.")
print("The number of hits each query has is as follows:")
print(query_counts)

query_counts.to_csv(output)

#create a graph that on the x axis shows each query and on the y axis the number of times they appear in the blast results
#this is to show the distribution of hits across the queries
import matplotlib.pyplot as plt
plt.bar(range(len(query_counts)), query_counts)
plt.xticks(range(len(query_counts)), ordered_queries, rotation='vertical')
#tilt the axis labels to make them easier to read
plt.tight_layout()
plt.xlabel("Effector candidate")
plt.ylabel("Number of hits")
plt.title("Number of hits across loci per effector candidate")
#save the plot
plt.savefig(outpath + "/query_counts.png", bbox_inches='tight')

