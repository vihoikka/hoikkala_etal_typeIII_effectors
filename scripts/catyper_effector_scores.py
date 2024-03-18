import argparse
import pandas as pd
import subprocess
import sys
import os

import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches
from matplotlib.ticker import FuncFormatter


parser = argparse.ArgumentParser()
#args for --pte {input.pte} --etp {input.etp} --output {output.effector_scores} --summary {output.effector_scores_summary} --outdir
parser.add_argument('-pte', '--pte', type=str, metavar='pte', required=False,
                    help='pte')
parser.add_argument('-etp', '--etp', type=str, metavar='etp', required=True,
                    help='etp')
parser.add_argument('-o', '--output', type=str, required=True,
                    help='Results')
parser.add_argument('-os', '--output_summary', type=str, required=True,
                    help='Results summary')
parser.add_argument('-od', '--outdir', type=str, required=True,
                    help='Results directory')
                    

args = parser.parse_args()

#pte = args.pte
etp = args.etp
output = args.output
output_summary = args.output_summary
outdir = args.outdir

#read in the pte file
#pte_df = pd.read_csv(pte, sep='\t', header=None)
#pte_df.columns = ['Number', 'Protein', 'Effector', 'Locus']
#pte_df.drop('Number', axis=1, inplace=True)
#print(pte_df)

#read in the etp file
etp_df = pd.read_csv(etp, sep='\t', header=None)
etp_df.columns = ['Number', 'Effector','Protein', 'Locus']
etp_df.drop('Number', axis=1, inplace=True)
#print(etp_df)

effectors_df = etp_df.groupby('Effector')['Protein'].apply(list).reset_index(name='Proteins') #group by effector and get a list of proteins for each effector
print(effectors_df)

#effectors_df["Number of Proteins"] = effectors_df["Proteins"].str.len()
#similarly calculate unique proteins
effectors_df["Number of Proteins"] = effectors_df["Proteins"].apply(set).str.len()
print(effectors_df)

def cross_hits(protein_list):
    # First, we remove the current protein list from the master list
    master_list = [protein for sublist in effectors_df["Proteins"].tolist() if sublist!=protein_list for protein in sublist]
    # Then we calculate the cross hits as before
    cross_hits = len(set(protein_list) & set(master_list))
    return cross_hits

# Now we use our new function in place of the lambda function
effectors_df["Number of cross-hits"] = effectors_df["Proteins"].apply(cross_hits)

#create a new column with the percentage of cross-hits
effectors_df["Percentage of cross-hits"] = effectors_df["Number of cross-hits"] / effectors_df["Number of Proteins"] * 100

effectors_df.to_csv(os.path.join(outdir,"effectors.tsv"), sep='\t', index=False)

# Cross matrix for effector-effector pairs
def shared_protein_counter(list1, list2):
    return len(set(list1) & set(list2)) #gets intersection of the two sets and returns its length, i.e. number of shared protein IDs

#the pairwise matrix
df_shared = pd.DataFrame(index=effectors_df["Effector"], columns=effectors_df["Effector"])

for effector1 in effectors_df["Effector"]:
    for effector2 in effectors_df["Effector"]:
        shared_proteins = shared_protein_counter(effectors_df.loc[effectors_df["Effector"]==effector1, "Proteins"].values[0],
                                                 effectors_df.loc[effectors_df["Effector"]==effector2, "Proteins"].values[0])
        df_shared.loc[effector1, effector2] = shared_proteins

# P
plt.figure(figsize=(10, 8))
integer_formatter = FuncFormatter(lambda x, pos: '%d' % x) #formats values to integers
heatmap = sns.heatmap(df_shared.astype(int), annot=True, cmap='viridis', fmt='d', cbar_kws={'format': integer_formatter})
plt.title('Heatmap of shared proteins between effectors')
plt.savefig(os.path.join(outdir,"heatmap.png"), bbox_inches='tight')

plt.clf()

# Generate a scatter plot where 'Total hits' is used for the x-axis values and 'Number of cross-hits' is used for the y-axis values
plt.scatter(x=effectors_df['Number of Proteins'], y=effectors_df['Percentage of cross-hits'])

#annotate the points with the effector names
for i, txt in enumerate(effectors_df['Effector']):
    plt.annotate(txt, (effectors_df['Number of Proteins'][i], effectors_df['Percentage of cross-hits'][i]))

#y axis limits 0 - 100
plt.ylim(0, 100)

plt.xlabel('Number of hits')  # Set the label for the x-axis
plt.ylabel('Percentage of cross-hits')  # Set the label for the y-axis
plt.title('Total hits vs Percentage of cross-hits per effector')  # Set the title of the plot

#save the plot
plt.savefig(os.path.join(outdir,"scatterplot.png"), bbox_inches='tight')