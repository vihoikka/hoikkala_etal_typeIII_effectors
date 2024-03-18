import pandas as pd
import os
import argparse
import re
from pprint import pprint
import numpy as np
import matplotlib.pyplot as plt
import mpld3


parser = argparse.ArgumentParser()
parser.add_argument("-s", "--species", required=True)
parser.add_argument("-on", "--out_nodes", required = True)
parser.add_argument("-oe", "--out_edges", required = True)
parser.add_argument("-om", "--out_merged_all", required = True)
parser.add_argument("-p", "--phage_data", required = True)
parser.add_argument("-cr", "--crispr_info", required = True)
parser.add_argument("-bp", "--base_path", required=True)
parser.add_argument("-ah", "--all_hosts", required=True)
args = parser.parse_args()
species = args.species
crispr_info = args.crispr_info
out_nodes = args.out_nodes
out_edges = args.out_edges
out_merged_all = args.out_merged_all
phage_metadata = args.phage_data
base_path = args.base_path
all_hosts = args.all_hosts

phage_csm6_out = base_path + "/" + species + "/06_merged/" + species + "_phage_csm6.tab"

phages = pd.read_csv(phage_metadata, sep = '\t', header =0) #read Milliard phage metadata into panda df, separator tab
blast = pd.read_csv(base_path + "/" + species + "/04_blast/blast.out", sep = '\t', header = 0) #header refers to header row number. Blast.out contains all hits from blasting the spacers.
crispr = pd.read_csv(crispr_info, sep = '\t', header = 0) #header refers to row number
#csm6 = pd.read_csv(csm6, sep = "\t", header=0)
#ca = pd.read_csv(cATyper, sep ="\t", header = 0)
bacteria = pd.read_csv(all_hosts, sep="\t", header = None) #list of all bacteria in dataset regardless of interactions
bacteria.columns = ["bacterium"]

#TODO: the tables resulting from lines below do not contain bacteria that are not targeting phages. Fix this: we need all hosts.



blast["host_nt_spacer"] = blast['Query_id'].str.extract(r'(.+?(?=@))', expand=False) #everything until @
blast["host_nt"] = blast['host_nt_spacer'].str.extract(r'(.*(?=\_))', expand=False)

blast["host"] = blast['Query_id'].str.extract(r'((?<=\@).*)', expand=False) #everything after @

bacteria = pd.merge(blast, bacteria, left_on = "Subject_accession_ID_version", right_on = "bacterium", how = "outer")

#merge the blast and phage database
merged_1 = pd.merge(bacteria, phages, left_on="Subject_accession_ID_version", right_on = 'Accession',  how="left")
pprint("Blast shape: " + str(blast.shape))
pprint("Phage db shape: " + str(phages.shape))
pprint("Merged 1 shape: " + str(merged_1.shape))

#The merged_all.csv contains one row = one blast hit with phage information
#merged_1 = pd.merge(merged_1, csm6, on="host", how="left")
#merged_1 = pd.merge(merged_1, ca, on="host", how="left")
merged_1.to_csv(out_merged_all, index = False, sep = '\t', header = True)

#creating a phage-only table that contains phage info + whether it is being targeted by a Csm6 containing host
#phage_info = pd.merge(merged_1, csm6, on="host", how="left")
phagecounts = merged_1.groupby(['Subject_id'])['Subject_id'].count().reset_index(name="interactions") #total number of interactions

#get total positives and negatives of csm6-interactions
#phagecounts_csm6 = pd.crosstab(merged_1.Subject_id,merged_1.csm6)
#transform above into a dataframe
#phagecounts_csm6 = phagecounts_csm6.add_prefix('csm6').reset_index().rename_axis(None, axis=1) 
#merge with total interactions
#phagecounts = pd.merge(left=phagecounts, right=phagecounts_csm6, on="Subject_id", how="left")

#sometimes we might have none or all csm6 hosts. In those cases, the other column is not created and that will cause problems later
#if "csm6False" not in phagecounts:
#    phagecounts["csm6False"] = 0
#if "csm6True" not in phagecounts:
#    phagecounts["csm6True"] = 0

#calculate fraction of csm6-positive interactions out of all
#phagecounts["fraction_csm6"] = phagecounts["csm6True"]/phagecounts["interactions"]

#rename subject_id column
phagecounts.rename(columns = {'Subject_id': 'phage'}, inplace = True)

#save csm6 fraction as histogram
#phagecounts["fraction_csm6"].plot(kind = 'hist')
#plt.savefig("hist.png")
#phagecounts.to_csv(phage_csm6_out, index=False, sep = '\t', header = True)

#For Gephi visualisation, we need the interaction data (edges) and bacterium/phage data (nodes).
#Nodes contain metadata on each node. These can be used in Gephi for filtering or aggregating things

#Gephi nodes. A node is either a bacterium or a phage.
gephi_nodes_phages_temp = pd.DataFrame()
gephi_nodes_bacteria_temp = pd.DataFrame()

gephi_nodes_phages_temp["Id"] = merged_1["Subject_accession_ID_version"]
gephi_nodes_bacteria_temp["Id"] = merged_1["bacterium"]
pprint(gephi_nodes_bacteria_temp)
pprint(crispr)
gephi_nodes_phages_temp = gephi_nodes_phages_temp.drop_duplicates()
gephi_nodes_bacteria_temp = gephi_nodes_bacteria_temp.drop_duplicates()

gephi_nodes_phages_temp["type"] = "phage"
gephi_nodes_bacteria_temp["type"] = "host"

gephi_nodes = pd.merge(gephi_nodes_phages_temp, gephi_nodes_bacteria_temp, on = ["type", "Id"], how='outer')
gephi_nodes = pd.merge(gephi_nodes, crispr, left_on="Id", right_on="assembly", how="left") #add CRISPR info to nodes that are hosts
pprint(gephi_nodes)
#gephi_nodes = pd.merge(gephi_nodes, csm6, left_on="Id", right_on="host", how="left") #add CRISPR info to nodes that are hosts
gephi_nodes = pd.merge(gephi_nodes, phages, left_on="Id", right_on="Accession", how="left") #add CRISPR info to nodes that are hosts
gephi_nodes = pd.merge(gephi_nodes, phagecounts, left_on = "Id", right_on="phage", how="left") #add phage csm6 host fractions
#gephi_nodes = pd.merge(gephi_nodes, ca, left_on = "Id", right_on="host", how="left") #add phage csm6 host fractions

gephi_nodes.to_csv(out_nodes, index=False, sep = '\t', header = True)

#Gephi edges require Source and Target columns. These represent the arrows in Gephi. Each source and target must
#have equivalent values in the nodes file
gephi_edges = merged_1
gephi_edges = gephi_edges.rename(columns = {
    'host': 'Source',
    'Subject_accession_ID_version': 'Target'},
    inplace = False)

gephi_edges.to_csv(out_edges, index=False, sep = '\t', header = True)