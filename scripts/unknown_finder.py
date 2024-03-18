
import argparse
import pandas as pd
import os
import gffutils
from Bio import SeqIO
import ast
import subprocess

'''
Searches type III loci for unknown effectors. Relies on CCTyper predictions of proteins annotated with "Unk" prefix and on
poor E-values for known genes
'''

evalue_cutoff = 1e-7

#create a new class unknown_protein with properties locus_id, sample, sequence, length, evalue, position
class unknown_protein:
    def __init__(self, locus_id, sample, sequence, length, evalue, position, id, cctyper):
        self.locus_id = locus_id
        self.sample = sample
        self.sequence = sequence
        self.length = length
        self.evalue = evalue
        self.position = position
        self.id = id
        self.cctyper = cctyper

    def to_dict(self):
        return {
            'locus_id': self.locus_id,
            'sample': self.sample,
            'sequence': self.sequence,
            'length': self.length,
            'evalue': self.evalue,
            'position': self.position,
            'id': self.id,
            'cctyper': self.cctyper
        }
    
class locus_unknown_info:
    def __init__(self, locus_id, sample, no_of_unknowns, unknown_proteins):
        self.locus_id = locus_id
        self.sample = sample
        self.no_of_unknowns = no_of_unknowns
        self.unknown_proteins = unknown_proteins

    def to_dict(self):
        return {
            'locus_id': self.locus_id,
            'sample': self.sample,
            'no_of_unknowns': self.no_of_unknowns,
            'unknown_proteins': self.unknown_proteins,
        }

def literal_converter(val):
    try:
        return ast.literal_eval(val)
    except (ValueError, SyntaxError):
        return val  # or substitute with a default value

parser = argparse.ArgumentParser()
#create args for --locus_id --cctyper --hmm --outputfolder --info_out, --unknown_proteins, --gff, --host_genomes_folder, cctyper_path, sample_folder
parser.add_argument("-l", "--locus_id", help="locus_id", required = True)
parser.add_argument("-c", "--cctyper", help="cctyper", required = True)
parser.add_argument("-hm", "--hmm", help="hmm", required = True)
parser.add_argument("-o", "--outputfolder", help="outputfolder", required = True)
parser.add_argument("-i", "--info_out", help="info_out", required = True)
parser.add_argument("-u", "--unknown_proteins_output", help="unknown_proteins", required = True)
parser.add_argument("-hg", "--host_genomes_folder", help="host_genomes_folder", required = True)
parser.add_argument("-cp", "--cctyper_path", help="cctyper_path", required = True)
parser.add_argument("-sa", "--sample_folder", help="sample_folder", required = True)
parser.add_argument("-ol", "--output_locus_info", help="output_locus_info", required = True)
parser.add_argument("-add", "--additional_cas_db", help="diamond db for additional cas genes", required = True)
args = parser.parse_args()

#assign input and output folder to variables
locus_id = args.locus_id
cctyper = args.cctyper
hmm = args.hmm
outputfolder = args.outputfolder
info_out = args.info_out
unknown_proteins_output = args.unknown_proteins_output
host_genomes_folder = args.host_genomes_folder
cctyper_path = args.cctyper_path
sample_folder = args.sample_folder
output_locus_info = args.output_locus_info
additional_cas_db = args.additional_cas_db

#get sample name by splitting the locus name at the underscore and taking the first two parts
sample = locus_id.split("_")[0] + "_" + locus_id.split("_")[1]

hmm_filesize = int(os.path.getsize(hmm))
print("HMM filesize is " + str(hmm_filesize) + " bytes")
if hmm_filesize > 0:
    print("Known effectors found. Writing info file for locus") #the empty output files are made in the snakemake pipeline by touching
    unknown_protein_object = locus_unknown_info(locus_id, sample, 0, False)
    #create a pandas dataframe from the unknown_protein_object
    unknown_protein_df = pd.DataFrame([unknown_protein_object.to_dict()])
    #write the dataframe to a csv file
    unknown_protein_df.to_csv(output_locus_info, sep="\t", index=False, header = False)
    exit()

#read the cctyper into pandas
converters = {"Genes": literal_converter, 
              "E-values": literal_converter, 
              "Positions": literal_converter}


cctyper_df = pd.read_csv(cctyper, sep='\t', header = 0, converters=converters)

#find the gff file for the sample from the host_genomes_folder/samples folder
gff_file = os.path.join(host_genomes_folder, sample, sample + "_features.gff")

#Using gffutils, create a database of the gff file
db_path = os.path.join(outputfolder,locus_id)

#load gene position information from cctyper
cctyper_gene_positions_df = pd.read_csv(os.path.join(cctyper_path, sample, "genes.tab"), sep="\t", header=0)

#Read genes cas_operons column "Genes"
genes = cctyper_df["Genes"].tolist()[0] #split the genes into a list and remove the brackets
cas_operon_positions = cctyper_df["Positions"].tolist()[0] #split the positions into a list

#load proteins with biopython
protein_path = os.path.join(sample_folder, sample, sample + "_proteins.faa")
proteins = SeqIO.to_dict(SeqIO.parse(protein_path, "fasta"))

#create folder for the above path if it does not exist
if not os.path.exists(os.path.dirname(db_path)):
    print("Creating folder " + str(os.path.dirname(db_path)))
    os.makedirs(os.path.dirname(db_path))
db = gffutils.create_db(gff_file, 
                        dbfn=outputfolder + "/" + locus_id + ".db", 
                        force=True, 
                        keep_order=True, 
                        merge_strategy='merge', 
                        sort_attribute_values=True,
                        id_spec=["protein_id", "Name", "ID"])

#read the database
print("Reading database from " + str(db_path) + ".db")
db = gffutils.FeatureDB(outputfolder + "/" + locus_id + ".db", keep_order=True)

def getGeneLocationCRISPRCas(pos, contig, position_order):
    positions = pos #df of where gene position translates to coordinates
    print("Determining position of gene " + str(position_order) + " on contig " + contig + "...")
    genepos = positions.loc[(positions["Contig"] == contig) & (positions["Pos"] == position_order)] #find the row that corresponds to current operon
    print(genepos)
    return [genepos["Start"].values[0],genepos["End"].values[0]]

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
        if length > 50: #consider changing this if the pipeline crashes
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

def blast_against_additional_cas_genes(sequence, id):
    '''
    Runs diamond as subprocess against the additional cas genes / effector database
    '''
    #create a temporary file for the sequence
    id = str(id)
    temp_path = outputfolder + "/" + id + ".faa"
    out_path = outputfolder + "/" + id + ".out"
    with open(temp_path, "w") as f:
        f.write(">" + id + "\n")
        f.write(sequence)
    #run diamond
    subprocess.run(["diamond", "blastp", "-d", additional_cas_db, "-q", temp_path, "-o", out_path, "-f", "6"])
    #read the size of the output file
    diamond_output_size = int(os.path.getsize(out_path))
    #if the output file is empty, return False
    if diamond_output_size == 0:
        print("No matches found in the additional cas genes database for " + id)
        return False
    else:
        print("Found match(es) in the additional cas genes database for " + id)
        return True

gene_search_range = 200 #how far from the cctyper coordinate are we looking for the gene in the gff file in bp

#find proteins that contain the string "unk" (case insensitive) in the column "Genes" and create a new unknown_protein object for each of these
#note that the column "Genes" is a list in string format
unknown_proteins = []
for index, row in cctyper_df.iterrows():
    #loop through the list in the "Genes" column using enumerate to get the index of each element
    for i, gene in enumerate(row["Genes"]):
        #if the string "unk" is found in the gene name or if an E-value is very low, create a new unknown_protein object
        evalue = float(row["E-values"][i])
        if ("Unk" in gene) or (evalue > evalue_cutoff):
            print("Found an unknown gene or a gene with poor E-value: " + gene + " with E-value: " + str(evalue))
            #get gene location
            start_end_tuple = getGeneLocationCRISPRCas(cctyper_gene_positions_df, row["Contig"], row["Positions"][i])
            
            #expand the search range in case the gene is not found in the exact position
            start = int(start_end_tuple[0]) - gene_search_range
            end = int(start_end_tuple[1]) + gene_search_range
            
            length = (end - start)/3 #length in AA

            #load the gene object from the gff and protein file
            try:
                gene_object = getGeneObject(db, proteins, row["Contig"], start, end)
            except:
                print("Corresponding gene from the GFF file could not be found. Skipping...")
                break
            #get sequence from the gene object
            seq = gene_object.seq
            #check if this gene is present in the additional cas/effector database
            additional_cas_check = blast_against_additional_cas_genes(str(seq), gene_object.id) #returns True or False
            if additional_cas_check == False:
                unknown_protein_object = unknown_protein(row["locus_id"], row["sample"], seq, length, row["E-values"][i], row["Positions"][i], row["locus_id"] + "_" + str(row["Positions"][i]), gene)
                unknown_proteins.append(unknown_protein_object)

            #append the object to the list unknown_proteins


#write the unknown_proteins to a fasta file
with open(unknown_proteins_output, "w") as f:
    for protein in unknown_proteins:
        f.write(">" + protein.id + "\n")
        f.write(str(protein.sequence) + "\n")

#create a pandas database from the unknown_proteins list. Pandas columns are created from the attributes of the unknown_protein object
unknown_proteins_df = pd.DataFrame([protein.to_dict() for protein in unknown_proteins])
print(unknown_proteins_df)

#output this database to a file with header
unknown_proteins_df.to_csv(info_out, index=False, sep = "\t", header=False)

#get number of unknown proteins
unknown_proteins_count = len(unknown_proteins_df.index)

if unknown_proteins_count == 0:
    unknown_proteins_boolean = False
else:
    unknown_proteins_boolean = True

#also output a locus-specific info file. If we made it this far in the script, we assume unknown proteins were found
unknown_protein_object = locus_unknown_info(locus_id, sample, unknown_proteins_count, unknown_proteins_boolean)
#create a pandas dataframe from the unknown_protein_object
unknown_protein_df = pd.DataFrame([unknown_protein_object.to_dict()])
#write the dataframe to a csv file
unknown_protein_df.to_csv(output_locus_info, index=False, sep = "\t", header=False)