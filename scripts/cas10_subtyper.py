import argparse
import os
from Bio import SeqIO
import pandas as pd

# create args inputs for metadata and outputfolder
parser = argparse.ArgumentParser(description='Subtype cas10')
parser.add_argument('-i', '--metadata', help='input file', required=True)
parser.add_argument('-o', '--outputfolder', help='output folder', required=True)
parser.add_argument('-f', '--fasta', help='output folder', required=True)
args = parser.parse_args()

#assign args to similarly named variables
metadata = args.metadata
outputfolder = args.outputfolder
fasta = args.fasta

#Using pandas from the input .csv file, read the column "sample" to get list of sample names and
#column "Subtype" to get corresponding list of subtypes. Then using biopython,
#scan through the fasta file and write out a new fasta file for each subtype.
#The new fasta files are named after the subtype and are written to the output folder.

#read in metadata
df = pd.read_csv(metadata, sep ="\t", header = 0)

#remove any rows in which the Subtype column contains (does not have to fully match) the string hybrid
print(df)
df = df[df['Subtype'].str.contains('Hybrid') == False]

#make a list of sample names
samples = df['Sample'].tolist()

#make a list of subtypes
subtypes = df['Subtype'].tolist()

#make a list of CorA association using the column "CorA"
cora = df['CorA'].tolist()

#make a dictionary of sample names and subtypes
sample_subtype = dict(zip(samples, subtypes))

#make a dictionary of sample names and cora association
sample_cora = dict(zip(samples, cora))

#make a list of unique subtypes
unique_subtypes = list(set(subtypes))


#make a dictionary of subtypes and empty lists
subtype_dict = {}
for subtype in unique_subtypes:
    subtype_dict[subtype] = []

#make a similar dictionary, but for samples that are cora-positive
subtype_dict_cora = {}
for subtype in unique_subtypes:
    subtype_dict_cora[subtype] = []

#load fasta file using biopython
for record in SeqIO.parse(fasta, "fasta"):
    #get sample name from fasta header
    sample = record.id
    #get the sample's subtype from the dictionary, if such a subtype exists
    if sample in sample_subtype:
        subtype = sample_subtype[sample]
        subtype_dict[subtype].append(record)

    #if sample is cora-positive, add it to the cora-positive dictionary under its corresponding subtype
    if sample in sample_cora:
        if sample_cora[sample] == True:
            subtype_dict_cora[subtype].append(record)


#write function for writing the fasta files
def write_fasta(dictionary, outputfolder, cora):
    #write each subtype to a fasta file
    for subtype in unique_subtypes:
        #make output file name
        output_file = os.path.join(outputfolder, "cas10_" + str(subtype) + cora + ".faa")
        #write out fasta file
        SeqIO.write(dictionary[subtype], output_file, "fasta")

#call the function twice, once for each dictionary
write_fasta(subtype_dict, outputfolder, "")
write_fasta(subtype_dict_cora, outputfolder, "_cora")
