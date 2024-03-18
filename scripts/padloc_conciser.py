import pandas as pd
import os
import argparse
import re
from pprint import pprint
import sys

#This script concises padloc into one row = one bacterium. Most code copied from crispr_conciser_jan23.py

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--species", required=True)
parser.add_argument("-bp", "--base_path", required=True)
#parser.add_argument("-o", "--output", required = True)
args = parser.parse_args()
species = args.species
base_path = args.base_path
#output = args.output

if (os.stat(base_path + "/" + species + "/05_padloc/" + species + "_padloc_info.csv").st_size == 0):
    print("No padloc info found. Exiting...")
    sys.exit()

padloc = pd.read_csv(base_path + "/" + species + "/05_padloc/" + species + "_padloc_info.csv", sep = ',', header = 0) #header refers to row number
print("Shape at start: " + str(padloc.shape))

# 3. Count the number of entries and make a df out of it
print("Counting multilocus entries")
counts = padloc.groupby('host').size()
counts = counts.to_frame()

# 4. Isolate subtypes for each
print("Isolating defense systems")
systems = padloc.groupby('host')['system'].apply(list)
systems = systems.to_frame()

# 5. Isolate N_repeats for each locus
#print("Isolating and summing repeats")
#repeats = padloc.groupby('assembly')['N_repeats'].apply(list)
#repeats = repeats.to_frame()

# 6. Sum the repeats across loci
#repeats['summed_repeats'] = repeats.apply(lambda row: sum(row.N_repeats), axis=1)

# 7. Get locations
start_locs = padloc.groupby('host')['start'].apply(list)
end_locs = padloc.groupby('host')['end'].apply(list)

# 8. Get host ID
hosts = padloc.groupby('host')['seqid'].apply(list)

# 9. Get protein type (e.g. cyclase or effector)
protein_names = padloc.groupby('host')['protein.name'].apply(list)

# 10. target descr.
targets = padloc.groupby('host')['target.description'].apply(list)

# system ID
system_IDs = padloc.groupby('host')['system.number'].apply(list)


#Merge above info
print("Merging")
merged = pd.merge(counts, systems, left_on = 'host',  right_on = 'host')
#merged = pd.merge(merged, repeats, on = 'assembly')
merged = pd.merge(merged, start_locs, on = 'host')
merged = pd.merge(merged, end_locs, on = 'host')
merged = pd.merge(merged, protein_names, on = 'host')
merged = pd.merge(merged, targets, on = 'host')
merged = pd.merge(merged, system_IDs, on = 'host')
merged = pd.merge(merged, hosts, on = 'host')

print(merged)

#merged["assembly"] = merged["assembly"].str[0]



filename = base_path + "/" + species + "/05_padloc/" + species + "_padloc_summarized.csv"
print("Writing file " + filename)
merged = merged.rename(columns = {0:"no_of_systems"}, inplace = False)

#Write to file
merged.to_csv(filename, index=True, sep = '\t', header = True)
