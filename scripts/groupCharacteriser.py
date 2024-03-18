'''
The current script does the following:
1. Extracts all "unknown proteins" that have the locus name in the header
2. Clusters said proteins with CD-HIT
3. Subjects each such protein to a SAVED/CARF search
4. Extracts tidied file, then reopens it to run a hmmscan against local pfam
5. Outputs everything in one neat file
'''
import argparse
import pandas as pd
import subprocess
import sys
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastpCommandline

#parse command line arguments for the following: input genome list, path to unknown proteins fasta file, path to output file, project
parser = argparse.ArgumentParser(description='This script takes a list of loci and extracts the unknown proteins for each locus, then subjects them to a blastp search against the NCBI nr database.')
parser.add_argument('-i', '--input', type=str, metavar='input', required=True,
                    help='The input file containing the list of loci to extract unknown proteins for.')
parser.add_argument('-u', '--unknown_proteins', type=str, metavar='unknown_proteins', required=True,
                    help='The fasta file containing all unknown proteins.')
parser.add_argument('-o', '--output_basename', type=str, metavar='output', required=True,
                    help='The output file to write the results to.')
parser.add_argument('-cu', '--clustered_unknowns', type=str, metavar='clustered_unknowns', required=True,
                    help='Output file containing clustered unknown proteins from the group')
parser.add_argument('-if', '--individual_fastas_dir', type=str, metavar='individual_fastas_dir', required=True,
                    help='A folder where individual protein fastas are outputted to.')
parser.add_argument('-of', '--outputfolder', type=str, metavar='outputfolder', required=True)
parser.add_argument('-tm', '--tmhmm_model', type=str, metavar='tmhmm_model', required=True)

#parser.add_argument('-p', '--project', type=str, metavar='project', required=True)

args = parser.parse_args()

CARFSAVED_hmm_db_path = "/media/volume/st_andrews/databases/carfsaved/02_hmm_profiles/carfsaved.hmm"
pfams_hmm_db = "/media/volume/st_andrews/databases/pfam/Pfam-A.hmm"
project_root = "/media/volume/st_andrews/new_effectors"

#output files
all_against_pfam = args.output_basename + ".all.pfam.rawresults.txt.temp"
carfsaved_raw_results = args.output_basename + ".CARFSAVED.rawresults.tsv"
carfsaved_info = args.output_basename + ".CARFSAVED.info.tsv"
carfsaved_fasta = args.output_basename + ".CARFSAVED.fasta"
carfsaved_against_pfam = args.output_basename + ".CARFSAVED.pfam.rawresults.tsv"
carfsaved_against_pfam_info = args.output_basename + ".CARFSAVED.pfam.info.tsv"
outputfolder = args.outputfolder
tmhmm_model = args.tmhmm_model

pfam_fasta = args.output_basename + ".pfam.fasta"
pfam_info = args.output_basename + ".pfam.info.tsv"

E_value_pfam = "1e-05"
E_value_CARFSAVED = "1e-05"

def patchTMHMM(environment):
    '''
    Patches the tmhmm api.py file to ignore the RuntimeWarning
    '''
    path_to_apipy = os.path.join(".snakemake/conda", environment, "lib/python3.6/site-packages/tmhmm/api.py")
    with open(path_to_apipy, 'r') as f:
        first_line = f.readline()
        if first_line == "import warnings\n":
            print("TMHMM api.py already patched")
        else:
            print("Patching a warning bug in TMHMM apy.py")
            rows_to_add = ["import warnings", "warnings.filterwarnings('ignore', category=RuntimeWarning)"]
            #insert rows to start of file
            with open(path_to_apipy, 'r+') as f:
                content = f.read()
                f.seek(0, 0)
                for row in rows_to_add:
                    f.write(row.rstrip('\r\n') + '\n' + content)

#check if tmhmm in installed in currently active environment. If not, install with pip and patch the warning bug
try:
    print("Checking TMHMM installation")
    subprocess.run(["tmhmm", "-h"])
    print("TMHMM already installed")
    print("Patching a warning bug in TMHMM apy.py")
    patchTMHMM(os.environ['CONDA_DEFAULT_ENV'])

except FileNotFoundError:
    print("TMHMM not installed. Installing with pip")
    subprocess.run(["pip", "install", "tmhmm.py"])
    print("TMHMM installed")
    #check if the tmhmm api.py file is patched to ignore the RuntimeWarning
    current_env = os.environ['CONDA_DEFAULT_ENV']
    print("Current environment: " + current_env)
    path_to_apipy = os.path.join(".snakemake/conda", current_env, "lib/python3.6/site-packages/tmhmm/api.py")
    with open(path_to_apipy, 'r') as f:
        first_line = f.readline()
        if first_line == "import warnings\n":
            print("TMHMM api.py already patched")
        else:
            patchTMHMM(os.environ['CONDA_DEFAULT_ENV'])

#create output folder for TMHMM
subprocess.run(["mkdir", "-p", outputfolder + "/TMHMM"])

#read in the list of loci
loci = []
with open(args.input, 'r') as f:
    for line in f:
        loci.append(line.strip())

#read in the unknown proteins fasta file
unknown_proteins = SeqIO.to_dict(SeqIO.parse(args.unknown_proteins, 'fasta'))

#group_loci
group_loci = [] #a list of loci with unknown proteins

unknown_proteins_df = pd.DataFrame(columns=['Locus', 'Protein'])

#extract any proteins from the fasta file that have the locus name in the header
unknown_proteins_to_blast = []
for locus in loci:
    for protein in unknown_proteins:
        if locus in protein:
            unknown_proteins_to_blast.append(unknown_proteins[protein])
            locus_name = "_".join(protein.split("_")[:3]) #extracts the locus name from the protein name
            protein_name = str(protein)
            #add the protein with the locus to the dataframe
            unknown_proteins_df = unknown_proteins_df.append({'Locus': locus_name, 'Protein': protein_name}, ignore_index=True)
print("unknown_proteins_df")
print(unknown_proteins_df)

#write the extracted proteins to a fasta file
SeqIO.write(unknown_proteins_to_blast, "unknown_proteins_to_blast.fasta", 'fasta')

#cluster the extracted proteins with CD-HIT
subprocess.run(['cd-hit', '-i', 'unknown_proteins_to_blast.fasta', '-o', args.clustered_unknowns, '-c', '0.4', '-n', '2', '-d', '0', '-M', '0', '-T', '40'])

#open the clustered unknown proteins fasta file and extract the sequences into individual fasta files
for i, record in enumerate(SeqIO.parse(args.clustered_unknowns, "fasta")):
    # Writing each individual sequence to separate file
    dir_path = args.individual_fastas_dir
    with open(dir_path + f"/group4_{i+1}.fasta", "w") as output_handle:
        SeqIO.write(record, output_handle, "fasta")

#we need to add each protein to a pandas df manually first, since some of them won't have hits in all databases
group4_sequences = SeqIO.parse(open(args.clustered_unknowns),'fasta')

##### TM PREDICTION #####
def run_TMHMM(sequence, id):
    '''
    Runs TMHMM as a subprocess. Create temporary .faa files for the sequence and the output in their own folders
    We need to add 'desc' as placeholder to the sequence header for TMHMM to work
    '''
    #create a temporary file for the sequence
    id = str(id)
    temp_path = outputfolder + "/TMHMM/" + id + "/" + id + ".faa"
    out_path = outputfolder + "/TMHMM/" + id + "/" + id + ".out"
    out_folder = outputfolder + "/TMHMM/" + id
    header = id + " desc"
    #create the folder
    subprocess.run(["mkdir", "-p", out_folder])
    #convert the Seq into string
    sequence = str(sequence)
    with open(temp_path, "w") as f:
        f.write(">" + header + "\n")
        f.write(sequence)
    #run TMHMM
    print("Running TMHMM for " + id)
    #print the subprocess command before running it
    print(["tmhmm", "-m", tmhmm_model, "-f" , temp_path])
    subprocess.run(["tmhmm", "-m", tmhmm_model, "-f" , temp_path])

    #tmhmm creates the output files into the folder where tmhmm is run. We need to move it to the output folder
    outputfilesnames = [id + ".summary", id + ".plot", id + ".annotation"]
    for file in outputfilesnames:
        subprocess.run(["mv", file, out_folder])

    summary_filepath = out_folder + "/" + id + ".summary"

    has_TM = False
    tmhmm_count = 0

    #check if this file has lines
    if os.stat(summary_filepath).st_size != 0:
        #read the file and print it
        with open(summary_filepath, 'r') as f:
            for line in f:
                if "trans" in line:
                    tmhmm_count += 1
        if tmhmm_count > 0:
            has_TM = True
            
    return has_TM, tmhmm_count #return tuple

#for each of the group4_sequences proteins, run TMHMM using the run_TMHMM function
#then add the results to the dataframe
import multiprocessing

TM_df = pd.DataFrame(columns=['Locus', 'Protein', 'has_TM', 'tmhmm_count'])

def process_record(record):
    has_TM, tmhmm_count = run_TMHMM(record.seq, record.id)
    return record.id, has_TM, tmhmm_count

if __name__ == '__main__':
    with multiprocessing.Pool() as pool:
        results = pool.map(process_record, group4_sequences) #uses multiple cores to run the function using the list of Seq objects 

    for result in results:
        record_id, has_TM, tmhmm_count = result
        #extract locus from protein name. A protein may be called GCA_000968375.1_1_663, and the locus is then GCA_000968375.1_1
        locus_name = "_".join(record_id.split("_")[:3])
        TM_df = TM_df.append({'Locus': locus_name, 'Protein': record_id, 'has_TM': has_TM, 'tmhmm_count': tmhmm_count}, ignore_index=True)
    
    print(TM_df)
    TM_df.to_csv(outputfolder + "/TMHMM/TM_results.tsv", sep='\t', index=False)

##### CARFSAVED CHARACTERISATION #####
def process_hmm(hmm_raw, hmm_tidy):
    for line in hmm_raw:
        if not line.startswith("#"):  # skip the header
            line = line.strip().split()
            #print(line)
            target = line[0]
            query = line[3]
            query_length = line[5]
            evalue = line[6]
            score_overall = line[7]
            bias_score_overall = line[8]
            this_domain = line[9]
            total_domains = line[10]
            c_evalue = line[11]
            i_evalue = line[12]
            score_domain = line[13]
            bias_score_domain = line[14]
            start = line[17]
            end = line[18]
            env_start = line[19]
            env_end = line[20]
            description = line[22]
            hmm_tidy.write(target + "\t" + 
                           query + "\t" +
                           query_length + "\t" +
                           evalue + "\t" +
                           score_overall + "\t" +
                           bias_score_overall + "\t" +
                           this_domain + "\t" +
                           total_domains + "\t" +
                           c_evalue + "\t" +
                           i_evalue + "\t" +
                           score_domain + "\t" +
                           bias_score_domain + "\t" +
                           start + "\t" +
                           end + "\t" +
                           env_start + "\t" +
                           env_end + "\t" +
                           env_start + "\t" +
                           description + "\n")

def insert_sequence(hmm_tidy, protein_list): #inserts the sequence into the dataframe and to a protein list
    for index, row in hmm_tidy.iterrows():
        query = row["Query"]
        seq = unknown_proteins[query].seq
        hmm_tidy.at[index, "Seq"] = seq
        protein_list.append(SeqRecord(seq, id=query, description=""))


#run hmmscan subprocess against our SAVED/CARF database
with open(os.devnull, 'w') as DEVNULL:
    subprocess.run(['hmmscan', '--domtblout', carfsaved_raw_results, '--cpu', '40', '-E', E_value_CARFSAVED, CARFSAVED_hmm_db_path, args.clustered_unknowns], stdout=DEVNULL, stderr=subprocess.STDOUT)

#open the hmmscan result file (hmmscan_results.txt) and extract the locus name, the protein name, the e-value and the score
#write these to a new file
with open(carfsaved_raw_results, 'r') as hmm_raw: #open the hmmscan results file
    with open(carfsaved_info, 'w') as hmm_tidy: #open a new file to write the results to
        process_hmm(hmm_raw, hmm_tidy) #processes line by line and creates a new file

##open the file created above and extract the locus name, the protein name, the e-value and the score and create a pandas dataframe
hmm_tidy_carfsaved = pd.read_csv(carfsaved_info, sep='\t', header=None)
print(hmm_tidy_carfsaved)
#define headers
hmm_tidy_carfsaved.columns = ["Target_Name", "Query", "Query_Length", "Evalue", "Overall_Score", 
                             "Overall_Bias_Score", "This_Domain", "Total_Domains", "C_Evalue", "I_Evalue",
                             "Domain_Score", "Domain_Bias_Score", "Start", "End", "Env_Start",
                             "Env_End", "Another_Env_Start", "Description"]
hmm_tidy_carfsaved["Seq"] = "" #create a new column called "Seq" and leave it blank for now
print("Reopened tidied carfsaved:")
print(hmm_tidy_carfsaved)
carfsavedproteins = []

#inject sequence into df and to protein list
insert_sequence(hmm_tidy_carfsaved, carfsavedproteins) #modification done in place within function for df and protein list

#add all protein sequences to the pandas df to enable plotting even the ones with no hmm hits
for record in group4_sequences:
    hmm_tidy_carfsaved = hmm_tidy_carfsaved.append({'Target_Name': "protein", 'Query': record.id, 'Query_Length': len(record.seq), 'Evalue': 0, 'Overall_Score': 0, "Overall_Bias_Score": 0, 'This_Domain': 0, 'Total_Domains': 0, 'C_Evalue': 0.0, 'I_Evalue': 0.0, 'Domain_Score': 0.0, 'Domain_Bias_Score': 0.0, 'Start': 0, 'End': 0, 'Env_Start': 0, 'Env_End': 0, 'Another_Env_Start': 0, 'Description': "protein", 'Seq': record.seq}, ignore_index=True)

#write the tidied dataframe to a file, overwriting the previous file
hmm_tidy_carfsaved.to_csv(carfsaved_info, sep='\t', index=False)

#write the carfsaved proteins to a fasta file
SeqIO.write(carfsavedproteins, carfsaved_fasta, 'fasta')

print("CARFSAVED done")
print(hmm_tidy_carfsaved)


##### PFAM CHARACTERISATION #####

#run all unknown proteins against the pfam database. Silence outputs using DEVNULL
with open(os.devnull, 'w') as DEVNULL:
    subprocess.run(['hmmscan', '--domtblout', all_against_pfam, '--cpu', '40', '-E', E_value_pfam, pfams_hmm_db, args.clustered_unknowns], stdout=DEVNULL, stderr=subprocess.STDOUT)

with open(all_against_pfam, 'r') as hmm_raw: #open the hmmscan results file
    with open(pfam_info, 'w') as hmm_tidy: #open a new file to write the results to
        process_hmm(hmm_raw, hmm_tidy) #processes line by line and creates a new file

# open the file created above and extract the locus name, the protein name, the e-value and the score and create a pandas dataframe
hmm_tidy_pfam = pd.read_csv(pfam_info, sep='\t', header=None)
hmm_tidy_pfam.columns = ["Target_Name", "Query", "Query_Length", "Evalue", "Overall_Score", 
                             "Overall_Bias_Score", "This_Domain", "Total_Domains", "C_Evalue", "I_Evalue",
                             "Domain_Score", "Domain_Bias_Score", "Start", "End", "Env_Start",
                             "Env_End", "Another_Env_Start", "Description"]
hmm_tidy_pfam["Seq"] = "" #create a new column called "Seq" and leave it blank for now
pfamproteins = []

#inject sequence into df and to protein list
insert_sequence(hmm_tidy_pfam, pfamproteins) #modifies the dataframe and protein list in place

#add all protein sequences to the pandas df to enable plotting even the ones with no hmm hits
group4_sequences = SeqIO.parse(open(args.clustered_unknowns),'fasta')
for record in group4_sequences:
    print("Adding record " + str(record.id) + " to pfam df")
    hmm_tidy_pfam = hmm_tidy_pfam.append({'Target_Name': "protein", 'Query': record.id, 'Query_Length': len(record.seq), 'Evalue': 0, 'Overall_Score': 0, "Overall_Bias_Score": 0, 'This_Domain': 0, 'Total_Domains': 0, 'C_Evalue': 0.0, 'I_Evalue': 0.0, 'Domain_Score': 0.0, 'Domain_Bias_Score': 0.0, 'Start': 0, 'End': 0, 'Env_Start': 0, 'Env_End': 0, 'Another_Env_Start': 0, 'Description': "protein", 'Seq': record.seq}, ignore_index=True)

#write the tidied dataframe to a file, overwriting the previous file
hmm_tidy_pfam.to_csv(pfam_info, sep='\t', index=False)

#write the carfsaved proteins to a fasta file
SeqIO.write(pfamproteins, pfam_fasta, 'fasta')

print("PFAM done")
print(hmm_tidy_pfam)

# The code below merges the two dataframes, but this is not really used so commented out for now


# ##### MERGE CARFSAVED AND PFAM #####
# #before merging, rename columns in the carfsaved and the pfam dataframes to avoid confusion
# hmm_tidy_carfsaved = pd.read_csv(carfsaved_info, sep='\t', header=0)
# hmm_tidy_pfam = pd.read_csv(pfam_info, sep='\t', header=0)

# print("CARFSAVED reloaded")
# print(hmm_tidy_carfsaved)

# hmm_tidy_carfsaved.rename(columns={"Target_name": "CARFSAVED_target_name", "E-value": "CARFSAVED_E-value", "Score": "CARFSAVED_score", "Seq": "CARFSAVED_seq"}, inplace=True)
# hmm_tidy_carfsaved.drop(columns=["CARFSAVED_seq"], inplace=True) #drop the seq column because it's the same as the seq column in the pfam dataframe

# hmm_tidy_pfam.rename(columns={"Target_name": "pfam_target_name", "E-value": "pfam_E-value", "Score": "pfam_score"}, inplace=True)

# print(hmm_tidy_carfsaved)
# print(hmm_tidy_carfsaved.columns)
# #print just the Query column
# print(hmm_tidy_carfsaved["Query"])

# hmm_tidy_carfsaved['Query'] = hmm_tidy_carfsaved['Query'].str.strip()
# hmm_tidy_pfam['Query'] = hmm_tidy_pfam['Query'].str.strip()
# unknown_proteins_df['Protein'] = unknown_proteins_df['Protein'].str.strip()

# #merge the dataframes with all proteins
# merged = pd.merge(unknown_proteins_df, hmm_tidy_carfsaved, how='left', left_on='Protein', right_on='Query')
# merged = pd.merge(merged, hmm_tidy_pfam, how='left', left_on='Protein', right_on='Query')

# merged.fillna("-", inplace=True)

# merged = merged.groupby(['Locus', 'Protein']).agg({
#     'pfam_target_name': lambda x: list(x),
#     'pfam_E-value': lambda x: list(x),
#     'pfam_score': lambda x: list(x),
#     'CARFSAVED_target_name': lambda x: list(x),
#     'CARFSAVED_E-value': lambda x: list(x),
#     'CARFSAVED_score': lambda x: list(x),
#     'Seq' : 'first'
# }).reset_index()

# #combine rows with the same protein name by concatenating the values in the pfam_target_name and CARFSAVED_target_name to lists. Do the same for the corresponding e-values and scores etc
# #merged = merged.groupby(['Locus', 'Protein']).agg({'pfam_target_name': lambda x: list(x), 'pfam_E-value': lambda x: list(x), 'pfam_score': lambda x: list(x), 'CARFSAVED_target_name': lambda x: list(x), 'CARFSAVED_E-value': lambda x: list(x), 'CARFSAVED_score': lambda x: list(x)})

# #write the dataframe to a file
# merged.to_csv(args.output_basename + ".all.tsv", sep='\t', index=False)

# def get_best_hit(group, target, evalues):
#     min_index = group[evalues].index(min(group[evalues]))
#     return group[target][min_index]

# #create a simplified dataframe in which the targets with best e-value are chosen
# merged_simplified = pd.DataFrame(columns=['Locus', 'Protein', 'pfam_target_name', 'pfam_E-value', 'pfam_score', 'CARFSAVED_target_name', 'CARFSAVED_E-value', 'CARFSAVED_score', 'Seq'])
# merged_simplified['Locus'] = merged['Locus']
# merged_simplified['Protein'] = merged['Protein']
# merged_simplified['best_pfam_hit'] = merged.apply(lambda row: get_best_hit(row, 'pfam_target_name', 'pfam_E-value'), axis=1)
# merged_simplified['best_savedcarf_hit'] = merged.apply(lambda row: get_best_hit(row, 'CARFSAVED_target_name', 'CARFSAVED_E-value'), axis=1)
# merged_simplified['Seq'] = merged['Seq']


# #write the simplified dataframe to a file
# merged_simplified.to_csv(args.output_basename + ".all.simplified.tsv", sep='\t', index=False)