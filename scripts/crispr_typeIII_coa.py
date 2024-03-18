import pandas as pd
import os
import argparse
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor
import ast
'''
1. Is a host CRISPR III positive? If so..
    1.1) extract CRISPR-Cas region plus n kb around it
    1.2) use gff file to obtain protein headers in the this region
    1.3) use protein headers to obtain sequences from faa file
    1.4) make a note whether the host is cA4 positive or negative
    1.5) extract type III proteins 
    And if it's not CRISPR III positive...
    1.1) Mark the output file accordingly
'''

parser = argparse.ArgumentParser()
parser.add_argument("-g", "--gff", required = True)
parser.add_argument("-p", "--proteins", required = True)
parser.add_argument("-m", "--merge", required = True)
parser.add_argument("-hi", "--host_id", required = True)
#parser.add_argument("-op", "--out_proteins", required = True)
parser.add_argument("-om", "--out_metadata", required = True)
#parser.add_argument("-bp", "--base_path", required=False)
#parser.add_argument("-s", "--species", required=False)

args = parser.parse_args()
gff = args.gff
proteins = args.proteins
merge = args.merge
#out_proteins = args.out_proteins
out_metadata = args.out_metadata
#base_path = args.base_path
#species = args.species
host_id = args.host_id

#Load host metadata from "merged"
merged = pd.read_csv(merge, sep = '\t', header=0)
merged = merged[merged["Id"] == host_id] #get current host

#print(merged["Id"])
#print(merged["CRISPR_prediction"])
#print(merged.head(1))

host_dict = merged.to_dict(orient="records")
#print(host_dict)
#print(host_dict[0]["CRISPR_prediction"])


CRISPR_prediction = host_dict[0]["CRISPR_prediction"]

if (pd.isna(CRISPR_prediction) == False): #if host has CRISPR
    CRISPR_prediction = CRISPR_prediction.strip('][').split(", ") #turn this into an actual list (was string)
    #If there are multiple CRISPR-Cas loci, we need the list indeces that refer to type III loci
    indices_iii = [i for i, s in enumerate(CRISPR_prediction) if "III-" in s]
    if len(indices_iii) > 0: #this indicates there are type III loci
        #get corresponding start and stop coordinates
        loci_count = len(indices_iii)
        if loci_count == 1:
            operon_list = ast.literal_eval(host_dict[0]["CRISPR_Operon_Pos"])
           # print(operon_list)
           # print(operon_list[0])
            for o in operon_list:
                print(o)
                coordinates = o.strip('][')
                list = coordinates.split(", ")
                print(list[0])
          #  for o in operon_list:
          #      nested 
            #the coordinates are stored like this
            #['[34178, 40518]', '[16239, 23656]']
          #  print(host_dict[0]["CRISPR_Operon_Pos"])
            #CRISPR_Operon_Pos = pos_string[1:len(pos_string)-1:1]
          #  print(pos_string)
            #CRISPR_Operon_Pos = host_dict[0]["CRISPR_Operon_Pos"].split(", ") #convert positions into list
            
          #  print(test[0])
            #print(type(CRISPR_Operon_Pos))
#            start = host_dict[0]["CRISPR_"]
 #           stop = 
       # print(CRISPR_prediction)
       # print(indices_iii)
   # if len(indices_iii) > 0:
   #     print("CRISPR III indices in dataframe: " + str(indices_iii))

    merged.to_csv(out_metadata, sep="\t", index = False)
else:
    merged.to_csv(out_metadata, sep="\t", index = False)

# if "III-" in host_dict["CRISPR_prediction"]:
#     III_list_position = 
# else:
#     print("Host is type III CRISPR-Cas negative. Writing blank result file")


#1. Subset into III positive and cA4 negative hosts
#hosts_IIItrue_ca4neg = merged[(merged["CRISPR_prediction"].str.contains('III-')) & (merged["ca4True"] == "FALSE")] #subset those with CRISPR III and no ca4


#"host   |   start   |   stop    |   list_of_gene_acc    |   list_of_AA_seq      |"
#"esa    |   23      |   120     |   ['asdasd','ofgnsd'] |   ['POJASOD','AOISDF']|"

