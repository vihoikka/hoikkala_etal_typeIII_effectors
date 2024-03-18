'''
Once CCTyper has run for all samples, this script takes each sample and creates a subfolder for each locus_id.
This way new Snakemake rules can be created for each locus_id.
'''

#add argparse and as arguments input and outputfolder
import argparse
import pandas as pd
import os

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_folder", help="input folder")
parser.add_argument("-o", "--output_folder", help="output folder")
parser.add_argument("-in", "--interference_cutoff", help="interference cutoff")
args = parser.parse_args()

#assign input and output folder to variables
inputfolder = args.input_folder
outputfolder = args.output_folder
interference_cutoff = args.interference_cutoff

#create a list of samples by examining the subfolders in the inputfolder
samples = [f for f in os.listdir(inputfolder) if os.path.isdir(os.path.join(inputfolder, f))]

#loop through each sample and extract CRISPR-Cas information into separate folders and files. These will be used as new wildcards in the Snakemake pipeline

#create a multithreaded loop that runs for each sample
import multiprocessing
from joblib import Parallel, delayed

def processInput(sample):
    #open cas_operons file if it exists
    if os.path.isfile(inputfolder + "/" + sample + "/cas_operons.tab"):
        cas_operons = inputfolder + "/" + sample + "/cas_operons.tab"
    else:
        print("No cas_operons file found for " + sample)
        #exit the function
        return

    #using panda, open CRISPR_Cas and cas_operons into dataframes
    cas_operons_df = pd.read_csv(cas_operons, sep='\t', header = 0)

    #in cas_operons, only retain rows where column "Prediction" includes "III-"
    cas_operons_df = cas_operons_df[cas_operons_df['Prediction'].str.contains("III-")]

    #then only retain rows where the list of genes in the "Genes" column includes "cas10" exactly once, case insensitive
    cas_operons_df = cas_operons_df[cas_operons_df['Genes'].str.contains("cas10", case=False)]

    #discard any rows where the list of genes in the "Genes" column includes "cas10" more than once, case insensitive
    cas_operons_df = cas_operons_df[~cas_operons_df['Genes'].str.contains("cas10.*cas10", case=False)]
    print("Discarded " + str(len(cas_operons_df)) + " rows where cas10 was not found exactly once")

    #then only retain rows where the list of g

    #filter out any hybrid loci by matching hybrid (case insensitive) in the "Prediction" column
    #cas_operons_df = cas_operons_df[~cas_operons_df['Prediction'].str.contains("hybrid", case=False)]

    #transform the column Complete_Interference in cas_operons_df from percentage to int and remove % symbol
    #cas_operons_df['Complete_Interference'] = cas_operons_df['Complete_Interference'].str.replace('%', '').astype(int)

    #filter cas_operons on the basis of Complete_Interference column (0-100), only keep rows with > interference_cutoff
    #cas_operons_df = cas_operons_df[cas_operons_df['Complete_Interference'] > int(interference_cutoff)]

    print("No of complete cas_operons found: " + str(len(cas_operons_df)))



    #create a new column in cas_operons_df called "locus_id" which consists of "sample" + "$" + a running integer
    cas_operons_df["locus_id"] = sample + "_" + cas_operons_df.index.astype(str)

    #create a new column that contains the sample
    cas_operons_df['sample'] = sample

    #create subfolders in outputfolder for each locus_id
    for index, row in cas_operons_df.iterrows():
        print(row["locus_id"])
        os.makedirs(outputfolder + "/" + row["locus_id"])

        #create a transposed version of the row
        row_transposed = row.to_frame().T #T refers to transpose

        #save the same file but as a regular csv file with header on top and the locus as one row below that
        row_transposed.to_csv(str(outputfolder + "/" + row["locus_id"] + "/cas_operons.tsv"), sep='\t', header=True, index=True)

#run the loop
Parallel(n_jobs=multiprocessing.cpu_count())(delayed(processInput)(sample) for sample in samples)



# for sample in samples:
#     #using panda, open CRISPR_Cas and cas_operons into dataframes
#     CRISPR_Cas_df = pd.read_csv(CRISPR_Cas, sep='\t', header = 0)
#     cas_operons_df = pd.read_csv(cas_operons, sep='\t', header = 0)
#     tax_info_df = pd.read_csv(tax_info, sep='\t', header = 0)

#     #in cas_operons, remove rows where column "Prediction" does not include "III-"
#     cas_operons_df = cas_operons_df[cas_operons_df['Prediction'].str.contains("III-")]

#     #transform the column Complete_Interference in cas_operons_df from percentage to int and remove % symbol
#     cas_operons_df['Complete_Interference'] = cas_operons_df['Complete_Interference'].str.replace('%', '').astype(int)

#     #filter cas_operons on the basis of Complete_Interference column (0-100), only keep rows with > 80
#     cas_operons_df = cas_operons_df[cas_operons_df['Complete_Interference'] > int(interference_cutoff)]

#     print("No of complete cas_operons found: " + str(len(cas_operons_df)))

#     #create a new column in cas_operons_df called "locus_id" which consists of "sample" + "$" + a running integer
#     cas_operons_df['locus_id'] = sample + "$" + cas_operons_df.index.astype(str)

#     #create a new column that contains the sample
#     cas_operons_df['sample'] = sample
#     cas_operons_df['species'] = tax_info_df[0]["Species"]
#     cas_operons_df['genus'] = tax_info_df[0]["Genus"]

#     #create subfolders in outputfolder for each locus_id
#     for index, row in cas_operons_df.iterrows():
#         os.makedirs(outputfolder + "/" + row['locus_id'])

#         #save current row of df (with header) as tab separated file in the subfolder with the name of the locus_id
#         cas_operons_df.loc[index].to_csv(outputfolder + "/" + row['locus_id'] + "/" + row['locus_id'] + "_cas_operons.tsv", sep='\t')