import pandas as pd
import os
import argparse
import gffutils
from Bio import SeqIO
import regex
import traceback
import ast


'''
This script looks at a CCTyper output and fishes out the Cas10, Cas5 and Cas7 proteins.

This version incorporates the HD domain search for Cas10. HD domains are searched for in the Cas10 protein
and within the first 5 to 30 residues.
'''

min_protein_length = 150 #the minimum length of a protein (in bp) when finding the CCTyper annotated proteins from the original .gff
gene_search_range = 250 #this is the range of nucleotides to search for the gene from the cctyper coordinates
cas10_min_length = 500 #minimum length in AA of Cas10

#add argparse and as arguments "cctyper", outputs files for the following (use "out" as prefix): Cas10, Cas5, Cas7, info, folder
parser = argparse.ArgumentParser()
parser.add_argument("-l", "--locus_id", help="unique ID of the CRISPR-Cas locus")
parser.add_argument('--sample_folder', type=str, help='sample folder', required=True)
parser.add_argument('--this_folder', type=str, help='current folder', required=True)
parser.add_argument("-o", "--outputfolder", help="output folder")
parser.add_argument("-c", "--cas_operons", help="cas_operons output file for this locus")
parser.add_argument("-i", "--info_out", help="info output file")
parser.add_argument("-cc", "--cctyper_path", help="path to common cctyper folder")

args = parser.parse_args()

#assign all arguments to variables above
locus_id = args.locus_id
sample_folder = args.sample_folder
this_folder = args.this_folder
outputfolder = args.outputfolder
cas_operons = args.cas_operons
info_out = args.info_out
cctyper_path = args.cctyper_path


#open cas_operons as a pandas df
cas_operons = pd.read_csv(cas_operons, sep="\t", header=0, converters={"Genes":ast.literal_eval, 
                                         "Positions":ast.literal_eval, 
                                         "E-values":ast.literal_eval})

print(cas_operons)

#naming schemes for the different Cas proteins
cas5_names = ["Cas5", "Cmr3", "Csm4"]
cas7_names = ["Cas7", "Csm3", "Cmr4"]
cas10_names = ["Cas10", "Csm1", "Cmr2"]

#This dictionary contains alternative names for the different Cas proteins
dict_of_genes_to_examine = {"Cas10" : cas10_names,
                            "Cas5" : cas5_names,
                            "Cas7" : cas7_names}

gene_search_ranges = {"Cas10" : 300,
                      "Cas5" : 130,
                      "Cas7" : 150}

protein_min_lengths = {"Cas10" : 500,
                    "Cas5" : 150,
                    "Cas7" : 150}

#generate a new dictionary called "info_table" using keys from dict_of_genes_to_examine
info_table = {key: False for key in dict_of_genes_to_examine}

#get the sample name from the info_table column "sample" and index 0
sample = cas_operons["sample"][0]

#add rest of the keys
info_table["Locus"] = locus_id
info_table["Sample"] = sample
info_table["Cas10_GGDD"] = False
info_table["Cas10_GGDD_coord"] = "-"
info_table["Cas10_GGDD_seq"] = ""
info_table["Cas10_GGDE"] = False
info_table["Cas10_GGDE_coord"] = "-"
info_table["Cas10_GGDE_seq"] = ""
info_table["Cas10_GGED"] = False
info_table["Cas10_GGED_coord"] = "-"
info_table["Cas10_GGED_seq"] = ""
info_table["Cas10_HD"] = False
info_table["Cas10_HD_list"] = ""
info_table["Cas10_DH"] = False
info_table["Cas10_HD_coord"] = "-"
info_table["Cas10_DH_coord"] = "-"
info_table["Cas10_coord"] = "-"
info_table["Cas10_length"] = 0
#info_table["Unknown_genes"] = False
info_table["Subtype"] = cas_operons["Prediction"][0]

#convert the dictionary to a pandas dataframe
info_table = pd.DataFrame(info_table, index=[0])

#read gff using gffutils
gff_path = os.path.join(str(sample_folder), str(sample), str(sample) + "_features.gff")
gff_db = gffutils.create_db(gff_path, dbfn=this_folder + "/" + locus_id + '_gff.db', force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)
gff = gffutils.FeatureDB(this_folder + "/" + locus_id + '_gff.db', keep_order=True)

#load proteins with biopython
protein_path = os.path.join(sample_folder, sample, sample + "_proteins.faa")
proteins = SeqIO.to_dict(SeqIO.parse(protein_path, "fasta"))

#load gene position information
cctyper_gene_positions_df = pd.read_csv(os.path.join(cctyper_path, sample, "genes.tab"), sep="\t", header=0)

def getGeneLocationCRISPRCas(pos, contig, position_order):
    positions = pos #df of where gene position translates to coordinates
    print("Determining position of gene " + str(position_order) + " on contig " + contig + "...")
    genepos = positions.loc[(positions["Contig"] == contig) & (positions["Pos"] == position_order)] #find the row that corresponds to current operon
    print(genepos)
    return [genepos["Start"].values[0],genepos["End"].values[0]]


def getCas10Domain(cas10, domain_sequence):
    seq = str(cas10)
    #find the HD domain within the first 5 to 30 residues
    seq = seq[0:50]
    #search for "HD" or "DH" in seq and output the coordinates of the match
    start_index = seq.find(domain_sequence)
    end_index = start_index + 2
    #if the motif is found, return the matched motif
    if start_index != -1:
        return str(start_index)
    #if no match is found, return "-"
    else:
        return "-"

def getCas10Domain_list(cas10, domain_sequence):
    #This new version outputs all occurrences of the domain as list
    seq = str(cas10)  # make sure to cast cas10 to a string and trim
    
    # find the starting indices of all occurrences of domain_sequence in seq
    start_indices = [i for i in range(len(seq)) if seq[i:i+len(domain_sequence)] == domain_sequence]
    
    # if at least one motif was found, return the list of indices as strings
    if start_indices:
        return [str(i) for i in start_indices]
    # if no match was found, return empty list
    else:
        return []

def checkGenePresence(genes, dict_of_genes_to_examine, info_table, positions):
    return_dict = {} #we will return multiple dicts, so we will use a master dictionary to store them
    return_dict["positions"] = {}

    #Loop through the cas_operon genes
    for index, gene in enumerate(genes):
        print("Checking gene: " + gene)
        #Loop through the dictionary of genes to examine
        for key, value in dict_of_genes_to_examine.items():
            print("Checking against: " + key + "(" + str(value) + ")")
            for alternative_name in value:
                if alternative_name in gene:
                    print("Found hit")
                    info_table[key] = True  #if the gene is found, set the corresponding value in the info table to True
                    return_dict["positions"][key] = positions[index] #add the position of the gene to the return dictionary

    return_dict["info_table"] = info_table #add the info table to the return dictionary
    return return_dict #return the return dictionary that contains both the positions and the info table

def getGeneObject(gff, proteins, contig, start, end):
    print("Contig: " + contig)
    print("Start: " + str(start))
    print("End: " + str(end))
    protein_objects = gff.region(seqid=contig, start = start, end = end, featuretype="CDS", completely_within = True)
    #gene = gff.region(region=(contig, start, end), featuretype="CDS", completely_within = True)
    #print("Gene object: " + str(list(protein_object)))
    print(protein_objects)
    protein_objects = list(protein_objects) #converts generator to list (gffutils returns a generator)
    print(protein_objects)
    print(type(protein_objects))
    #create a dictionary to hold the protein objects
    protein_candidates = {}
    #populate the dictionary with the protein objects. In this dictionary, the key will be the position of the protein object in the list and the value will be the length of the protein object
    for index, protein in enumerate(protein_objects):
        length = int(protein.end) - int(protein.start)
        if length > min_protein_length: #consider changing this if the pipeline crashes
            print("Adding protein object to candidates")
            protein_candidates[index] = int(length)
    #print the size of the dictionary with length of protein objects
    print("The number of protein candidates in search window: " + str(len(protein_candidates)))
    #get the index of the protein object with the longest length
    longest_protein_index = max(protein_candidates)
    print("Chose the longest protein candidate at index " + str(longest_protein_index) + " with length: " + str(protein_candidates[longest_protein_index]) + " bp")
    #assign this protein object to the variable "protein_object"
    protein_object = protein_objects[longest_protein_index]
    #
    #print all protein attributes
    print(protein_object.attributes)
    print("Start: " + str(protein_object.start))
    print("End: " + str(protein_object.end))
    name = protein_object.attributes["Name"][0]
    print(name)
    gene_object = proteins[name]
    gene_object.id = locus_id
    gene_object.description = ""
    return gene_object

def getCas10CyclaseDomainDict(cas10_seq, motif_seq):
    seq = str(cas10_seq)
    motif_seq = str(motif_seq)
    motif_fuzzy = regex.findall(f"({motif_seq}){{e<=0}}", seq, overlapped=True) #allow 0 errors
    returnDict = {
        "seq": "",
        "cyclase": False,
        "coordinate": "-"
    }
    if motif_fuzzy: # check if motif_fuzzy is not empty
        returnDict["seq"] = motif_fuzzy[0] # return the sequence
        returnDict["cyclase"] = True # return that it is a cyclase
        for match in regex.finditer(f"({motif_seq}){{e<=0}}", seq, overlapped=True): 
            returnDict["coordinate"] = match.start() # return the coordinate
            break # break the loop after the first match
    return returnDict


#Read genes cas_operons column "Genes"
genes = cas_operons["Genes"].tolist()[0] #split the genes into a list and remove the brackets
cas_operon_positions = cas_operons["Positions"].tolist()[0] #split the positions into a list

#check if the genes list contain any of the genes in the dictionary "dict_of_genes_to_examine"
genePresence_dict = checkGenePresence(genes, dict_of_genes_to_examine, info_table, cas_operon_positions)

#update the info table and positions by extracting them from the returned dictionary that contains both
info_table = genePresence_dict["info_table"] #update the info table
positions_of_found_genes = genePresence_dict["positions"] #update the positions

#For each gene that was found, find their positions in the genome and extract as fasta sequences
#Also do some extra checks for Cas10 (GGDD and HD motifs)
print("Genes found: " + str(len(positions_of_found_genes)))

genes_not_found = [] #list to store the genes that were not found even if cctyper found them

#going through all the genes found by CCTyper
for key, value in positions_of_found_genes.items():
    print(key)
    start_end_tuple = getGeneLocationCRISPRCas(cctyper_gene_positions_df, cas_operons["Contig"][0], value) #value is the order position of the gene, not the actual coordinate
    
    #expand the search range in case the gene is not found in the exact position in the gff. This is calibrated separately for each gene
    start = int(start_end_tuple[0]) - gene_search_ranges[key]
    end = int(start_end_tuple[1]) + gene_search_ranges[key]
    print("**************************************")
    print("Searching for gene " + key + " in range: " + str(start) + "-" + str(end))
    #get corresponding gene from the gff file and turn it into a SeqRecord object
    try:
        gene_object = getGeneObject(gff, proteins, cas_operons["Contig"][0], start, end)
        min_protein_length = protein_min_lengths[key]

        #check if protein length is above threshold
        if len(gene_object.seq) >= min_protein_length:

        #check if the gene is Cas10 and if it contains the HD or GGDD motif
            if key == "Cas10" and len(gene_object.seq):
                #check if the GGDD domain is present in the Cas10 protein
                cas10CyclaseDict = getCas10CyclaseDomainDict(gene_object.seq, "GGDD") #store Cas10 GGDD info in a dictionary
                info_table["Cas10_GGDD"] = cas10CyclaseDict["cyclase"] #using the dictionary, update the info table
                info_table["Cas10_GGDD_seq"] = cas10CyclaseDict["seq"]
                info_table["Cas10_GGDD_coord"] = cas10CyclaseDict["coordinate"]

                cas10CyclaseDict = getCas10CyclaseDomainDict(gene_object.seq, "GGDE") #store Cas10 GGDD info in a dictionary
                info_table["Cas10_GGDE"] = cas10CyclaseDict["cyclase"] #using the dictionary, update the info table
                info_table["Cas10_GGDE_seq"] = cas10CyclaseDict["seq"]
                info_table["Cas10_GGDE_coord"] = cas10CyclaseDict["coordinate"]

                cas10CyclaseDict = getCas10CyclaseDomainDict(gene_object.seq, "GGED") #store Cas10 GGDD info in a dictionary
                info_table["Cas10_GGED"] = cas10CyclaseDict["cyclase"] #using the dictionary, update the info table
                info_table["Cas10_GGED_seq"] = cas10CyclaseDict["seq"]
                info_table["Cas10_GGED_coord"] = cas10CyclaseDict["coordinate"]

                #presume cyclase_literal is False
                info_table["cyclase_literal"] = False
                #make cyclase_literal true if any of the cyclase motifs were found (GGDD, GGDE or GGED)
                info_table["cyclase_literal"] = info_table["Cas10_GGDD"] | info_table["Cas10_GGDE"] | info_table["Cas10_GGED"]                
                #check if the HD domain is present in the Cas10 protein
                info_table["Cas10_HD_coord"] = getCas10Domain(gene_object.seq, "HD")
                info_table["Cas10_HD_list"] = str(getCas10Domain_list(gene_object.seq, "HD"))
                if info_table["Cas10_HD_coord"].iloc[0] != "-": #if HD was found
                    info_table["Cas10_HD"] = True

                info_table["Cas10_DH_coord"] = getCas10Domain(gene_object.seq, "DH") #sometimes HD is flipped as DH
                if info_table["Cas10_DH_coord"].iloc[0] != "-": #if DH was found
                    info_table["Cas10_DH"] = True

                info_table["Cas10_length"] = len(gene_object.seq)

            #Write the gene to a fasta file
            filepath = outputfolder + "/" + locus_id + "_" + key + ".faa"
            print("Writing file: " + filepath)
            SeqIO.write(gene_object, filepath, "fasta")
        #except with error reporting

        else:
            #break the for loop and continue to the next gene
            print("Protein " + key + " is too short " + str(len(gene_object.seq)) + ". Skipping")
            continue

    except Exception:
        print("Error when mapping gene coordinates: " + str(traceback.print_exc()))
        print("Could not find gene " + key + " in range: " + str(start) + "-" + str(end))
        genes_not_found.append(key)
        continue



#write the info table to a tsv file
info_table.to_csv(info_out, sep="\t", index=False, header = True)

#write the genes that were not found to a file
with open(outputfolder + "/" + locus_id + "_genes_not_found_eventhough_cctyper_found_them.txt", "w") as f:
    for gene in genes_not_found:
        f.write(gene + "\n")
