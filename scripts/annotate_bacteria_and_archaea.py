import argparse
import os
import pandas as pd

#create arguments for archaea, bacteria and output
parser = argparse.ArgumentParser(description='Annotate bacteria and archaea')
parser.add_argument('-a', '--archaea', type=str, metavar='', required=True, help='Path to archaea file')
parser.add_argument('-b', '--bacteria', type=str, metavar='', required=True, help='Path to bacteria file')
parser.add_argument('-o', '--output', type=str, metavar='', required=True, help='Path to output file')

args = parser.parse_args()

#read in archaea and bacteria files
archaea = pd.read_csv(args.archaea, sep='\t')
bacteria = pd.read_csv(args.bacteria, sep='\t')

#the files currently do not contain a column header. Create this header "sample"
archaea.columns = ['sample']
bacteria.columns = ['sample']

#add column domain to both archaea and bacteria files and insert the domain
archaea.insert(1, 'domain', 'Archaea')
bacteria.insert(1, 'domain', 'Bacteria')

#the archaea and bacteria files contain two columns: sample name and taxonomy (either bacteria or archaea).
#the bacteria file contains both bacteria and archaea (all annotated as bacteria) and the archaea file contains only archaea.
#merge the files so that the merged file contains both archaea and bacteria and the "misannotated" archaea are corrected by the archaea file
merged = pd.merge(bacteria, archaea, on='sample', how='outer') #outer join to keep all samples
print(merged)
merged['domain_y'] = merged['domain_y'].fillna(merged['domain_x']) # fill in the missing archaea column with the bacteria taxonomy
merged = merged.drop(columns=['domain_x']) #drop the bacteria taxonomy column
merged = merged.rename(columns={'domain_y': 'domain'}) #rename the archaea taxonomy column to taxonomy

#write the merged file to the output file
merged.to_csv(args.output, sep='\t', index=False)
