import pandas as pd
import os
import argparse
import re
from pprint import pprint
import sys

#This script concises CRISPR_info into one row = one bacterium
#Information is retained on the number of CRISPR loci and their types.
#For more in-depth info, see CRISPR_info.tab, where each row = CRISPR locus

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

# 1. Remove CRISPR classifications that are not "trusted"
trusted = crispr[crispr['Trusted'] == True]
print("Shape after removing untrusted locus predictions: " + str(trusted.shape))

#2. Remove locus predictions that are too close to another prediction (i.e. same locus counted twice)
#TODO!

# 3. Count the number of remaining entries and make a df out of it
print("Counting multilocus entries")
counts = trusted.groupby('assembly').size()
counts = counts.to_frame()

# 4. Isolate subtypes for each
print("Isolating subtypes")
subtypes = trusted.groupby('assembly')['Subtype'].apply(list)
subtypes = subtypes.to_frame()

# 5. Isolate N_repeats for each locus
print("Isolating and summing repeats")
repeats = trusted.groupby('assembly')['N_repeats'].apply(list)
repeats = repeats.to_frame()

# 6. Sum the repeats across loci
repeats['summed_repeats'] = repeats.apply(lambda row: sum(row.N_repeats), axis=1)

# 7. Get locations
start_locs = trusted.groupby('assembly')['Start'].apply(list)
end_locs = trusted.groupby('assembly')['End'].apply(list)

# 8. Get host ID
hosts = trusted.groupby('assembly')['Contig'].apply(list)

#Merge above info
print("Merging")
merged = pd.merge(counts, subtypes, left_on = 'assembly',  right_on = 'assembly')
merged = pd.merge(merged, repeats, on = 'assembly')
merged = pd.merge(merged, start_locs, on = 'assembly')
merged = pd.merge(merged, end_locs, on = 'assembly')
merged = pd.merge(merged, hosts, on = 'assembly')

print(merged)

#merged["assembly"] = merged["assembly"].str[0]



filename = base_path + "/" + species + "/02_cctyper/" + species + "_crispr_concise.csv"
print("Writing file " + filename)
merged = merged.rename(columns = {0:"N_CRISPR_loci"}, inplace = False)

#Write to file
merged.to_csv(filename, index=True, sep = '\t', header = True)
