import pandas as pd
import os
import argparse
import re
from pprint import pprint
import sys

#This script concises CRISPR_info into one row = one bacterium
#Information is retained on the number of CRISPR loci and their types.
#For more in-depth info, see CRISPR_info.tab, where each row = CRISPR locus

#Updated version (Jan 2023) that uses only CRISPR-Cas operon data as input

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--species", required=True)
parser.add_argument("-bp", "--base_path", required=True)
#parser.add_argument("-o", "--output", required = True)
args = parser.parse_args()
species = args.species
base_path = args.base_path
#output = args.output

if (os.stat(base_path + "/" + species + "/02_cctyper/CRISPR_info.tab").st_size == 0):
    print("No CRISPR info found. Exiting...")
    sys.exit()

crispr = pd.read_csv(base_path + "/" + species + "/02_cctyper/CRISPR_info.tab", sep = '\t', header = 0) #header refers to row number
print("Shape at start: " + str(crispr.shape))

# 3. Count the number of remaining entries and make a df out of it
print("Counting multilocus entries")
counts = crispr.groupby('assembly').size()
counts = counts.to_frame()

# 4. Isolate subtypes for each
print("Isolating subtypes")
subtypes = crispr.groupby('assembly')['Prediction'].apply(list)
subtypes = subtypes.to_frame()

# 5. Isolate N_repeats for each locus
#print("Isolating and summing repeats")
#repeats = crispr.groupby('assembly')['N_repeats'].apply(list)
#repeats = repeats.to_frame()

# 6. Sum the repeats across loci
#repeats['summed_repeats'] = repeats.apply(lambda row: sum(row.N_repeats), axis=1)

# 7. Get locations
locs = crispr.groupby('assembly')['Operon_Pos'].apply(list)

# 8. Get host ID
hosts = crispr.groupby('assembly')['Contig'].apply(list)

#Merge above info
print("Merging")
merged = pd.merge(counts, subtypes, left_on = 'assembly',  right_on = 'assembly')
#merged = pd.merge(merged, repeats, on = 'assembly')
merged = pd.merge(merged, locs, on = 'assembly')
merged = pd.merge(merged, hosts, on = 'assembly')

print(merged)

#merged["assembly"] = merged["assembly"].str[0]



filename = base_path + "/" + species + "/02_cctyper/" + species + "_crispr_concise.csv"
print("Writing file " + filename)
merged = merged.rename(columns = {0:"N_CRISPR_loci"}, inplace = False)

#Write to file
merged.to_csv(filename, index=True, sep = '\t', header = True)
