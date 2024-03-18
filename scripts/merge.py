import pandas as pd
import os
import argparse
import re
from pprint import pprint
import numpy as np
import matplotlib.pyplot as plt
import mpld3
from statistics import mean


parser = argparse.ArgumentParser()
parser.add_argument("-s", "--species", required=True)
parser.add_argument("-on", "--out_nodes", required = True)
parser.add_argument("-oe", "--out_edges", required = True)
parser.add_argument("-bpio", "--blast_with_phage_info_out", required = True)
parser.add_argument("-p", "--phage_data", required = True)
parser.add_argument("-cr", "--crispr_info", required = True)
parser.add_argument("-cs", "--csm6", required = True)
parser.add_argument("-ct", "--cATyper", required = True)
parser.add_argument("-cg", "--customgenes", required = True)
parser.add_argument("-bp", "--base_path", required=True)
parser.add_argument("-ah", "--all_hosts", required=True)
parser.add_argument("-pl", "--padloc", required=True)
parser.add_argument("-ca", "--cas10", required=True)
args = parser.parse_args()
species = args.species
crispr_info = args.crispr_info
out_nodes = args.out_nodes
out_edges = args.out_edges
blast_with_phage_info_out = args.blast_with_phage_info_out
phage_metadata = args.phage_data
csm6 = args.csm6
cATyper = args.cATyper
customgenes = args.customgenes
base_path = args.base_path
all_hosts = args.all_hosts
padloc = args.padloc
cas10 = args.cas10

phages_out = base_path + "/" + species + "/06_merged/" + species + "_phages.tab"

phages = pd.read_csv(phage_metadata, sep = '\t', header =0) #read Milliard phage metadata into panda df, separator tab
blast = pd.read_csv(base_path + "/" + species + "/04_blast/blast.out", sep = '\t', header = 0) #header refers to header row number. Blast.out contains all hits from blasting the spacers.
crispr = pd.read_csv(crispr_info, sep = '\t', header = 0, converters={'Prediction': pd.eval, 'Contig': pd.eval, 'Operon_Pos':pd.eval})#"converters" retains the python list
csm6 = pd.read_csv(csm6, sep = "\t", header = 0)
ca = pd.read_csv(cATyper, sep ="\t", header = 0)
cas10 = pd.read_csv(cas10, sep = "\t", header = 0)
padloc_summarized = pd.read_csv(padloc, sep = "\t", header = 0, converters={'system': pd.eval,'start': pd.eval,'end': pd.eval,'protein.name': pd.eval,'target.description': pd.eval,'system.number': pd.eval,'seqid': pd.eval})
customgenes = pd.read_csv(customgenes, sep ="\t", header = 0)
bacteria = pd.read_csv(all_hosts, sep="\t", header = None) #list of all bacteria in dataset regardless of interactions
bacteria.columns = ["bacterium"]

blast["host_nt_spacer"] = blast['Query_id'].str.extract(r'(.+?(?=@))', expand=False) #everything until @
blast["host_nt"] = blast['host_nt_spacer'].str.extract(r'(.*(?=\_))', expand=False)

blast["host"] = blast['Query_id'].str.extract(r'((?<=\@).*)', expand=False) #everything after @

#merge blast hits with all bacteria with outer join. This retains all entries. Those bacteria without blast hits are only in column "bacterium"
#bacteria = pd.merge(blast, bacteria, left_on = "Subject_accession_ID_version", right_on = "bacterium", how = "outer")


crispr.rename(columns={"Operon_Pos": "CRISPR_Operon_Pos", "Prediction": "CRISPR_prediction", "Contig": "CRISPR_contig"},inplace=True)


#merge the blast and phage database
blast_with_phage_info = pd.merge(blast, phages, left_on="Subject_accession_ID_version", right_on = 'Accession',  how="left")
pprint("Blast shape: " + str(blast.shape))
pprint("Phage db shape: " + str(phages.shape))
pprint("Merged 1 shape: " + str(blast_with_phage_info.shape))

#The blast_with_phage_info.csv contains one row = one blast hit with phage information
#merged_1 = pd.merge(merged_1, csm6, on="host", how="left")
#print(merged_1)

#We are also adding cATyper info to the blast/phage table. TODO: Consider not doing this
blast_with_phage_info = pd.merge(blast_with_phage_info, ca, on="host", how="left")

blast_with_phage_info.to_csv(blast_with_phage_info_out, index = False, sep = '\t', header = True)
#phage_info = pd.merge(merged_1, csm6, on="host", how="left")

phagecounts = blast_with_phage_info.groupby(['Subject_id'])['Subject_id'].count().reset_index(name="interactions") #total number of interactions
print(phagecounts)
interesting_attributes = ["ca6", "ca4", "ca3"]

for attribute in interesting_attributes:
    #creating a phage-only table that contains phage info + whether it is being targeted by a Csm6 containing host

    #get total positives and negatives of csm6-interactions
    phagecounts_interesting = pd.crosstab(blast_with_phage_info["Subject_id"],blast_with_phage_info[attribute])
    #print(phagecounts_interesting)


    #transform above into a dataframe
    phagecounts_interesting = phagecounts_interesting.add_prefix(attribute).reset_index().rename_axis(None, axis=1)
    #print(phagecounts_interesting)
    #rename attribute column to attribute_count
    old_column_name = str(attribute+"0.0")
    column_name = str(attribute+"_count")
    phagecounts_interesting.rename(columns={old_column_name: column_name},inplace=True)
    #print(phagecounts_interesting)
    #merge with total interactions
    phagecounts = pd.merge(left=phagecounts, right=phagecounts_interesting, on="Subject_id", how="left")
    #print(phagecounts)
    #sometimes we might have non or all attributed hosts. In those cases, the other column is not created and that will cause problems later
    #here we create those empty columns if necessary
    f = str(attribute + "False")
    t = str(attribute + "True")

    if f not in phagecounts:
        phagecounts[f] = 0
    if t not in phagecounts:
        phagecounts[t] = 0

    phagecounts["fraction_"+attribute] = phagecounts[t]/phagecounts["interactions"]
    #print(phagecounts)

    phagecounts["fraction_"+attribute].plot(kind = 'hist')
    plt.savefig("fraction_"+attribute+"_hist.png")



#calculate fraction of csm6-positive interactions out of all

#rename subject_id column
phagecounts.rename(columns = {'Subject_id': 'phage'}, inplace = True)

phagecounts.to_csv(phages_out, index=False, sep = '\t', header = True) #note that this is also combined to the gephi nodes file

#For Gephi visualisation, we need the interaction data (edges) and bacterium/phage data (nodes).
#Nodes contain metadata on each node. These can be used in Gephi for filtering or aggregating things

#Gephi nodes. A node is either a bacterium or a phage.
gephi_nodes_phages_temp = pd.DataFrame()
gephi_nodes_bacteria_temp = pd.DataFrame()

gephi_nodes_phages_temp["Id"] = blast_with_phage_info["Subject_accession_ID_version"]
gephi_nodes_bacteria_temp["Id"] = bacteria["bacterium"]

gephi_nodes_phages_temp = gephi_nodes_phages_temp.drop_duplicates()
gephi_nodes_bacteria_temp = gephi_nodes_bacteria_temp.drop_duplicates()

gephi_nodes_phages_temp["type"] = "phage"
gephi_nodes_bacteria_temp["type"] = "host"


gephi_nodes = pd.merge(gephi_nodes_phages_temp, gephi_nodes_bacteria_temp, on = ["type", "Id"], how='outer')
gephi_nodes = pd.merge(gephi_nodes, crispr, left_on="Id", right_on="assembly", how="left") #add CRISPR info to nodes that are hosts

#The above introduces two empty rows. We first mark them as NaN and then remove NaNs
gephi_nodes['Id'].replace('', np.nan, inplace=True)
gephi_nodes.dropna(subset=['Id'], inplace=True)

pprint(gephi_nodes)

CRISPRCas_types = ["I-","II-","III-","IV-","V-","VI-"]
CRISPRCas_types_dict = {"I": ["I-A", "I-B", "I-B2", "I-C", "I-D", "I-E", "I-F1", "I-F2", "I-F3", "I-G"],
                        "II": ["II-A", "II-B", "II-C"],
                        "III": ["III-A", "III-B", "III-C", "III-D", "III-E", "III-F", "III-G", "III-H"],
                        "IV": ["IV-A", "IV-B", "IV-C", "IV-D", "IV-E"],
                        "V": ["V", "V-K"],
                        "VI": ["VI"],
                        "other": ["other"]
                        }

gephi_nodes = pd.merge(gephi_nodes, csm6, left_on="Id", right_on="host", how="left") #add CRISPR info to nodes that are hosts
gephi_nodes = pd.merge(gephi_nodes, phages, left_on="Id", right_on="Accession", how="left")
gephi_nodes = pd.merge(gephi_nodes, phagecounts, left_on = "Id", right_on="phage", how="left")
gephi_nodes = pd.merge(gephi_nodes, ca, left_on = "Id", right_on="host", how="left")
gephi_nodes = pd.merge(gephi_nodes, customgenes, right_on = "host", left_on="Id", how="left")
gephi_nodes = pd.merge(gephi_nodes, padloc_summarized, right_on = "host", left_on="Id", how="left")
gephi_nodes = pd.merge(gephi_nodes, cas10, right_on = "host", left_on="Id", how="left")


#Create new boolean columns for each CRISPR-Cas subtype.
#If any of the subtypes are present in the Subtype column, add True to corresponding type column
gephi_nodes["CRISPR_prediction"] = gephi_nodes["CRISPR_prediction"].fillna('')
for crisprtype, subtypes in CRISPRCas_types_dict.items():
    gephi_nodes[crisprtype]=(gephi_nodes.CRISPR_prediction.apply(set)!=(gephi_nodes.CRISPR_prediction.apply(set)-set(subtypes)))

#similarly, create boolean columns for CBASS types
gephi_nodes["system"] = gephi_nodes["system"].fillna('')
cbass_types = ["cbass_type_I","cbass_type_II","cbass_type_III"]
for cbass_system in cbass_types:
    #we take all systems in a genome. Then we take all systems in a genome _except_ for the system we are interested in.
    #we then compare the whole set to the all-but-the-interesting set using negation (!=)
    #if they are the same, then the interesting system was not there and the function will return False
    #if they are different, that means the interesting system was there and will return True 
    difference = (gephi_nodes.system.apply(set)!=(gephi_nodes.system.apply(set)-set([cbass_system])))
    gephi_nodes[cbass_system] = difference

#newlist =                     [expression for item in iterable if condition == True]
#gephi_nodes["cbass_indeces"] = [system for system in zip(gephi_nodes.system) if "cbass" in system] #enumerate gives indices
#print(gephi_nodes["cbass_indeces"])

#gephi_nodes['cbass_indeces'] = gephi_nodes.apply(system for system in gephi_nodes.system if "cbass" in system)
#gephi_nodes['cbass_indeces'] = np.where("cbass" in gephi_nodes['system'], True, False)
#gephi_nodes.apply(lambda col: col.str[0])

gephi_nodes["cbass_indices"] = np.empty((len(gephi_nodes), 0)).tolist() #these are linked via list indices
gephi_nodes["cbass_contigs"] = np.empty((len(gephi_nodes), 0)).tolist() #these are linked via list indices
gephi_nodes["cbass_locations"] = np.empty((len(gephi_nodes), 0)).tolist() #these are linked via list indices

gephi_nodes["CRISPR_III_CBASS_sameContig"] = np.empty((len(gephi_nodes), 0)).tolist() #indices here refer to cbass indices above
gephi_nodes["CRISPR_III_CBASS_distances"] = np.empty((len(gephi_nodes), 0)).tolist() #indices here refer to cbass indices above

#gephi_nodes['CRISPR_contig'].replace('', np.nan, inplace=True)

for index, row in gephi_nodes.iterrows():
    if row["type"] == "host":
        for element in row["system"]: #loop through all defense systems
            if "cbass" in element: #if current system is cbass
                print("found cbass!")
                cbass_index = int(row["system"].index(element)) #get its index
                print("index: " + str(cbass_index))
                if cbass_index not in row["cbass_indices"]: #if index not already in list, add it
                    row["cbass_indices"].append(cbass_index)
                    contig = row["seqid"][cbass_index]
                    start = row["start"][cbass_index]
                    row["cbass_contigs"].append(contig) #also add corresponding contig...
                    row["cbass_locations"].append(start) #and start position
                    #print(gephi_nodes["seqid"][cbass_index])
                    #print(gephi_nodes["start"][cbass_index])
                    
        if row["N_CRISPR_loci"] > 0: #if we have crispr-cas loci

            for crispr_contig in row["CRISPR_contig"]: #loop through them
                crispr_index = row["CRISPR_contig"].index(crispr_contig) #store crispr list index
                crispr_type = row["CRISPR_prediction"][crispr_index] #get crispr type
                if (("III-" in crispr_type)): #if crispr is type iii
                    cbasscounter = 0
                    for cbass_contig in row["cbass_contigs"]:
                        print("CBASS contig: " + str(cbass_contig))
                        print("CRISPR contig: " + str(crispr_contig))
                        if cbass_contig == crispr_contig: #if a cbass and crispr are colocalised
                            #cbass_index = row["cbass_contigs"].index(cbass_contig)
                            row["CRISPR_III_CBASS_sameContig"].append(True) #mark this cbass index as being on same locus as a crispr iii system
                            print(row["CRISPR_Operon_Pos"][crispr_index])
                            #CRISPR_loc = row["CRISPR_Operon_Pos"][crispr_index].strip("'").split(",") #fetch current crispr location and transform into list
                            CRISPR_loc = row["CRISPR_Operon_Pos"][crispr_index].strip("][").split(",") #fetch current crispr location and transform into list
                            #print(CRISPR_loc)
                            #CRISPR_loc = list(map(int, CRISPR_loc)) #transform locations to int
                            CRISPR_loc = mean([int(CRISPR_loc[0]),int(CRISPR_loc[1])]) #get center of CRISPR locus
                            print(row["cbass_locations"])
                            print(cbass_index)
                            CBASS_loc = int(row["cbass_locations"][cbasscounter]) #get current cbass location
                            distance_crispr_cbass = abs(CRISPR_loc-CBASS_loc) #get difference between their coordinates
                            row["CRISPR_III_CBASS_distances"].append(distance_crispr_cbass) #mark this cbass index as being on same locus as a crispr iii system
                        else:
                            row["CRISPR_III_CBASS_sameContig"].append(False) #mark this cbass index as being on same locus as a crispr iii system
                            row["CRISPR_III_CBASS_distances"].append(0) #zero distance = not on same contig
                        cbasscounter += 1

gephi_nodes.to_csv(out_nodes, index=False, sep = '\t', header = True)


#Gephi edges require Source and Target columns. These represent the arrows in Gephi. Each source and target must
#have equivalent values in the nodes file
gephi_edges = blast_with_phage_info
gephi_edges = gephi_edges.rename(columns = {
    'host': 'Source',
    'Subject_accession_ID_version': 'Target'},
    inplace = False)

gephi_edges.to_csv(out_edges, index=False, sep = '\t', header = True)
