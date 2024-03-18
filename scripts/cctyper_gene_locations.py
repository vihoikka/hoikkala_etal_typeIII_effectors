import pandas as pd
import os
import argparse
import gffutils
from Bio import SeqIO
import regex
import traceback
import ast


#add argparse and as arguments "cctyper", outputs files for the following (use "out" as prefix): Cas10, Cas5, Cas7, CorA, info, folder
parser = argparse.ArgumentParser()
parser.add_argument("-l", "--locus_id", help="unique ID of the CRISPR-Cas locus")
parser.add_argument('--sample_folder', type=str, help='sample folder', required=True)
parser.add_argument('--this_folder', type=str, help='current folder', required=True)
parser.add_argument("-o", "--outputfolder", help="output folder")
parser.add_argument("-c", "--cas_operons", help="cas_operons output file for this locus")
parser.add_argument("-cc", "--cctyper_path", help="path to common cctyper folder")
parser.add_argument("-p", "--protein_fasta", help="protein fasta file for this locus")

args = parser.parse_args()

#assign all arguments to variables above
locus_id = args.locus_id
sample_folder = args.sample_folder
this_folder = args.this_folder
outputfolder = args.outputfolder
cas_operons = args.cas_operons
cctyper_path = args.cctyper_path
protein_fasta = args.protein_fasta


#open cas_operons as a pandas df
cas_operons = pd.read_csv(cas_operons, sep="\t", header=0, converters={"Genes":ast.literal_eval, 
                                         "Positions":ast.literal_eval, 
                                         "E-values":ast.literal_eval})

print(cas_operons)

#open protein fasta as a dictionary
protein_sequences = SeqIO.to_dict(SeqIO.parse(protein_fasta, "fasta"))

sample = cas_operons["sample"][0]
min_protein_length = 40 #minimum length of protein to be considered

#read gff using gffutils
gff_path = os.path.join(str(sample_folder), str(sample), str(sample) + "_features.gff")
gff_db = gffutils.create_db(gff_path, dbfn=this_folder + "/" + locus_id + '_gff.db', force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)
gff = gffutils.FeatureDB(this_folder + "/" + locus_id + '_gff.db', keep_order=True)

#load proteins with biopython
protein_path = os.path.join(sample_folder, sample, sample + "_proteins.faa")
proteins = SeqIO.to_dict(SeqIO.parse(protein_path, "fasta"))

#load gene position information
genes_tab_df = pd.read_csv(os.path.join(cctyper_path, sample, "genes.tab"), sep="\t", header=0)

def getGeneLocationCRISPRCas(genes_tab_df, contig, position_order):
    genes_tab = genes_tab_df #df of where gene position translates to coordinates (genes.tab)
    print("Determining position of gene " + str(position_order) + " on contig " + contig + "...")
    genepos = genes_tab.loc[(genes_tab["Contig"] == contig) & (genes_tab["Pos"] == position_order)] #find the row that corresponds to current operon
    print(genepos)
    return_dict = {}
    return_dict["Start"] = genepos["Start"].values[0]
    return_dict["End"] = genepos["End"].values[0]
    return_dict["Strand"] = genepos["Strand"].values[0]
    return return_dict

def checkGenePosition(genes, positions):
    return_dict = {}

    #Loop through the genes in the gene file
    for index, gene in enumerate(genes):
        print("Checking gene: " + gene)
        return_dict[gene] = positions[index]

    return return_dict #return the return dictionary that contains both the positions and the info table

def getGeneObject(gff, proteins, contig, start, end, strand):
    '''
    This function uses information from CCTyper to fetch the corresponding gene from the gff file.
    It returns a SeqRecord object of the gene.
    '''
    print("Contig: " + contig)
    print("Start: " + str(start))
    print("End: " + str(end))
    print("Strand: " + str(strand))

    #get all proteins from that match the CCTyper coordinates. Usually returns just one
    protein_objects = gff.region(seqid=contig, start = start, end = end, featuretype="CDS", completely_within = True)
    protein_objects = list(protein_objects) #converts generator to list (gffutils returns a generator)
    protein_candidates = {}

    #in case there are multiple hits in the gff file, we gather them all and choose the longest one
    #populate the dictionary with the protein objects. In this dictionary, the key will be the position of the protein object in the list and the value will be the length of the protein object
    for index, protein in enumerate(protein_objects):
        length = int(protein.end) - int(protein.start)
        if length > min_protein_length: #consider changing this if the pipeline crashes
            print("Adding protein object to candidates")
            protein_candidates[index] = int(length)
    print("The number of protein candidates in search window: " + str(len(protein_candidates)))
    longest_protein_index = max(protein_candidates)
    print("Chose the longest protein candidate at index " + str(longest_protein_index) + " with length: " + str(protein_candidates[longest_protein_index]) + " bp")
    #assign this protein object to the variable "protein_object"
    protein_object = protein_objects[longest_protein_index]

    strand_dict = {"+": "1", "-": "-1"}

    print(protein_object.attributes)
    print("Start: " + str(protein_object.start))
    print("End: " + str(protein_object.end))
    print("Strand: " + str(protein_object.strand))

    #create a dictionary with the gene's start and end coordinates, and the strand, id, product, length, sequence
    gene_properties = {}
    gene_properties["start"] = protein_object.start
    gene_properties["end"] = protein_object.end
    gene_properties["protein_id"] = protein_object.attributes["Name"][0]
    gene_properties["gene_id"] = protein_object.attributes["ID"][0]
    gene_properties["product"] = protein_object.attributes["product"][0]
    gene_properties["strand"] = str(protein_object.strand) #inherited from the cctyper output, not gff

    gene_properties["length"] = int(protein_object.end) - int(protein_object.start)

    return gene_properties

#this dataframe will contain all genes and some of their properties for annotation purposes
output_df = pd.DataFrame(columns=["protein_id", "start", "end", "annotation_gff", "annotation_cctyper", "strand", "evalue", "sequence"])

#Read genes cas_operons column "Genes"
genes_list = cas_operons["Genes"].tolist()[0] #split the genes into a list and remove the brackets
positions_list = cas_operons["Positions"].tolist()[0] #split the positions into a list
evalues = cas_operons["E-values"].tolist()[0] #split the evalues into a list

#create dictionary with genes as keys and positions as values
genes_positions_dict = {}
for index, gene in enumerate(genes_list):
    name = genes_list[index]
    pos = positions_list[index]
    evalue = evalues[index]
    genes_positions_dict[index] = name,pos,evalue #storing the gene, position and evalue as tuple in a dictionary with the index as key

genes_not_found = [] #list to store the genes that were not found even if cctyper found them
strand_dict = {"+": "1", "-": "-1"}

casR = [] #storing all CasR annotated genes here

#going through all the genes found by CCTyper
proteinCounter = 0
for key, value in genes_positions_dict.items():
    gene_name = value[0]
    pos = value[1]
    evalue = value[2]
    print("Gene index: " + str(key))
    print("Gene name: " + str(gene_name))
    print("Gene position: " + str(pos))

    #get start, end and strand values for this gene from the genes.tab file
    start_end_strand_dict = getGeneLocationCRISPRCas(genes_tab_df, cas_operons["Contig"][0], pos) #pos value is the order position of the gene, not the actual coordinate
    
    #the gene is now searched in the .gff file
    #expand the search range in case the gene is not found in the exact position in the gff. Using a fixed value for all genes
    start = int(start_end_strand_dict["Start"]) - 10
    end = int(start_end_strand_dict["End"]) + 10
    strand = start_end_strand_dict["Strand"]
    length = int(end) - int(start)
    protein_id = "cctyper_protein_" + str(proteinCounter)

    print("**************************************")
    print("Searching for gene " + gene_name + " in range: " + str(start) + "-" + str(end))
    #get corresponding gene from the gff file and turn it into a SeqRecord object
    try:
        #check if protein length is above threshold
        if length >= min_protein_length:
            #check if we can find a match for the protein in the original gff
            try:
                gene_object = getGeneObject(gff, proteins, cas_operons["Contig"][0], start, end, strand)
                #find the sequence from the fasta file
                sequence = protein_sequences[gene_object["protein_id"]].seq
                output_df = output_df.append({"protein_id":gene_object["protein_id"], "start_cctyper":gene_object["start"], "end_cctyper":gene_object["end"], "annotation_gff_cctyper":gene_object["product"], "annotation_cctyper_cctyper": gene_name, "strand_cctyper":gene_object["strand"], "evalue_cctyper": evalue,"sequence_cctyper": sequence}, ignore_index=True)
            #if not, then just add it without linking it to any gene
            except:
                print("Error when mapping gene coordinates: " + str(traceback.print_exc()))
                print("Could not find gene " + value[0] + " in range: " + str(start) + "-" + str(end))
                print("Will add CCTyper protein as is")
                strand = strand_dict[strand]
                output_df = output_df.append({"protein_id": protein_id, "start_cctyper":start, "end_cctyper":end, "annotation_gff_cctyper":"-", "annotation_cctyper_cctyper": gene_name, "strand_cctyper":strand, "evalue_cctyper": evalue, "sequence_cctyper": "-"}, ignore_index=True)

        else:
            #break the for loop and continue to the next gene
            print("Protein " + value[0] + " is too short " + str(length) + ". Skipping")
            continue

    except Exception:
        print("Error when mapping gene coordinates: " + str(traceback.print_exc()))
        print("Could not find gene " + value[0] + " in range: " + str(start) + "-" + str(end))
        genes_not_found.append(value[0])
        continue

    if "CASR" in gene_name.upper():
        print("Found a CasR gene")
        casR.append(gene_object)

    proteinCounter += 1

print(output_df)
#write the info table to a tsv file
output_df_path = os.path.join(outputfolder, locus_id + "_cctyper_gene_locations.tsv")
output_df.to_csv(output_df_path, sep="\t", index=False, header = True)

#write the genes that were not found to a file
with open(outputfolder + "/" + locus_id + "_genes_not_found_eventhough_cctyper_found_them.txt", "w") as f:
    for gene in genes_not_found:
        f.write(gene + "\n")

#from the CasR table, output the CasR genes to a fasta file
casR_fasta_path = os.path.join(outputfolder, locus_id + "_casR.faa")
with open(casR_fasta_path, "w") as f:
    for gene in casR:
        f.write(">" + gene["protein_id"] + "\n" + str(protein_sequences[gene_object["protein_id"]].seq) + "\n")