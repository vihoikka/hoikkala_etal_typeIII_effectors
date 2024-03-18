import argparse
import pandas as pd
import subprocess
import sys
import os

#create class effector that houses various interesting information about the effectors
class effector:
    def __init__(self, name, mean_temperature, std_temperature, GGDD_prop, count, unknown_proteins, most_common_genus, most_common_genus_prop, bacteria_prop):
        self.name = name
        self.mean_temperature = mean_temperature
        self.std_temperature = std_temperature
        self.GGDD_prop = GGDD_prop
        self.count = count
        self.unknown_proteins = unknown_proteins
        self.most_common_genus = most_common_genus
        self.most_common_genus_prop = most_common_genus_prop
        self.domain_prop = bacteria_prop


parser = argparse.ArgumentParser()
#args for input and output
parser.add_argument('--input', type=str, metavar='input', required=True,
                    help='A mastertable with all the information')
parser.add_argument('--output_basepath', type=str, required=True, metavar='output_basepath')

args = parser.parse_args()
input = args.input
output_basepath = args.output_basepath

mastertable = pd.read_csv(input, sep='\t', header=0)

#count the number of occurrences for each effector. The columns are nucc	can1	can2	csm6	csx1	csx23	cora
#the rows are the different loci
effectors = ['can1-2', 'calpl', 'csm6-ca6', 'csx1', 'csx23', 'cora', 'cami1', 'nucc']
effector_objects = []
for e in effectors:
    #create new effector object
    effector_objects.append(effector(e, 0, 0, 0, 0, '', 0, 0, 0))

#add attributes to the effector objects
for e_object in effector_objects:
    name = e_object.name
    print(name)
    e_object.count = mastertable[name].sum()
    print(e_object.count)
    e_object.mean_temperature = mastertable[mastertable[name] == True]['mean_temperature'].mean()
    e_object.std_temperature = mastertable[mastertable[name] == True]['mean_temperature'].std()
    e_object.GGDD_prop = mastertable[mastertable[name] == True]['GGDD_hmm_boolean'].mean()
    e_object.unknown_proteins = mastertable[mastertable[name] == True]['unknown_proteins'].mean()
    print(mastertable[mastertable[name] == True]['genus'])
    try:
        e_object.most_common_genus = mastertable[mastertable[name] == True]['genus'].value_counts().idxmax()
        e_object.most_common_genus_prop = mastertable[mastertable[name] == True]['genus'].value_counts(normalize=True)[0]
        e_object.bacteria_prop = mastertable[mastertable[name] == True]['domain'].value_counts(normalize=True)[0]
    except:
        e_object.most_common_genus = "None"
        e_object.most_common_genus_prop = float(0)
        e_object.bacteria_prop = float(0)
    print(e_object.name)
    print(e_object.count)
    print(e_object.mean_temperature)
    print(e_object.std_temperature)
    print(e_object.GGDD_prop)
    print(e_object.unknown_proteins)
    print(e_object.most_common_genus)
    print(e_object.most_common_genus_prop)
    print(e_object.bacteria_prop)


#transform all the effector objects into a df
effector_objects_df = pd.DataFrame([vars(e) for e in effector_objects])

#write the df to a file
effector_objects_df.to_csv(output_basepath + '_master.tsv', sep='\t', header=True, index=False)