import sys
import os
from Bio import SeqIO, SearchIO
from hmmer_to_pandas import parse_hmmer_domtblout
import argparse
import pandas as pd

argparser = argparse.ArgumentParser()
#args for input and output
argparser.add_argument('--input', type=str, metavar='input', required=True,
                    help='original raw hmm file of hmmscan results')
argparser.add_argument('--multifasta', type=str, metavar='input', required=True,
                    help='multifasta file of given effector. Used to insert dummy hmm entries for effector proteins with no hits.')
argparser.add_argument('--output_basepath', type=str, required=True, metavar='output_basepath')
argparser.add_argument('--effector', type=str, required=True, metavar='effector name')
argparser.add_argument('--effector_locus_map', type=str, required=True, metavar='effector name')


args = argparser.parse_args()

input = args.input
output_basepath = args.output_basepath
effector = args.effector
effector_locus_map = args.effector_locus_map

pfam_domain_i_Evalue = "1e-03" #cutoff for domain specific hits with independent E-value


#read in the hmmscan results
hmmscan_results = parse_hmmer_domtblout(input)
print(hmmscan_results)

#manually add all entries in the multifasta file, even if there is already a hmm hit. For these hits, made target_name "protein" and fill up other values with dummy values
#read in the fasta file
effector_sequences = SeqIO.parse(open(args.multifasta),'fasta')

for record in effector_sequences:
    print(record.id)
    hmmscan_results = hmmscan_results.append({'target_name': "protein", 'target_accession': "protein", 'tlen': len(record.seq), 'query_name': record.id, 'accession': "protein", 'qlen': len(record.seq), 'E-value': 0.0, 'score': 0.0, 'bias': 0.0, '#': 0, 'of': 0, 'c-Evalue': 0.0, 'i-Evalue': 0.0, 'domain_score': 0.0, 'domain_bias': 0.0, 'hmm_start': 0, 'hmm_end': 0, 'ali_start': 1, 'ali_end': len(record.seq), 'env_from': 0, 'env_to': 0, 'acc': 0.0, 'description of target': "protein"}, ignore_index=True)

# Sort rows by i-evalue in descending order
df_sorted = hmmscan_results.sort_values(by=['target_name', 'i-Evalue'], ascending=[True, True])
print(df_sorted)

#write the sorted df to a file
sorted_out = output_basepath + "/" + effector + "/" + effector +  "_sorted.tsv"
df_sorted.to_csv(sorted_out, sep='\t', header=True, index=False)

#change the type of column i-Evalue to float
df_sorted['i-Evalue'] = df_sorted['i-Evalue'].astype(float)

#remove any rows where the i-Evalue is greater than the cutoff
df_sorted_filtered = df_sorted[df_sorted['i-Evalue'] < float(pfam_domain_i_Evalue)]

#write the sorted and filtered df to a file
sorted_filtered_out = output_basepath + "/" + effector + "/" + effector + "_sorted_filtered.tsv"
df_sorted_filtered.to_csv(sorted_filtered_out, sep='\t', header=True, index=False)

# Function to retain only the best non-overlapping domains 
def filter_domains(data):
    data.reset_index(drop=True, inplace=True)
    filtered_data = data.loc[[0]]
    for i in range(1,len(data)):
        overlap = False
        for j in range(len(filtered_data)):
            if max(0, min(data['ali_start'][i], filtered_data['ali_end'].iloc[j]) - max(data['ali_start'][i], filtered_data['ali_start'].iloc[j])) > 0:
                overlap = True
                break
        if not overlap:
            filtered_data = pd.concat([filtered_data, data.loc[[i]]])
    return filtered_data

# Apply the function to filter domains by overlapping
#df_best_domains = df_sorted.groupby('target_name').apply(filter_domains).reset_index(drop=True)

protein_ID_map = pd.read_csv(effector_locus_map, sep='\t', header=None)
print(protein_ID_map)

#drop the first column
protein_ID_map = protein_ID_map.drop([0], axis=1)
protein_ID_map.columns = ["protein", "effector", "locus"]
print(protein_ID_map)

mapped_hmms_sorted_filtered = pd.merge(df_sorted_filtered, protein_ID_map, left_on = "query_name", right_on = "locus", how = "left")
print(mapped_hmms_sorted_filtered)

mapped_hmms_sorted_filtered.to_csv(output_basepath + "/" + effector + "/" + effector + "_sorted_filtered_mapped.tsv", sep='\t', index=False)