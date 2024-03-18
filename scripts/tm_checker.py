'''
This script contains the function to run TMHMM and check if a sequence has transmembrane domains.
Callable from other scripts. Requires installation of tmhmm.py via pip and modifiation of api.py
to bypass a warning using these lines:
    import warnings
    warnings.filterwarnings('ignore', category=RuntimeWarning)
'''
import subprocess
import os

ambiguous_amino_acids = ['X', 'B', 'Z', 'J']

def run_TMHMM(sequence, id, outputfolder, tmhmm_model_path):
    '''
    Runs TMHMM as a subprocess. Create temporary .faa files for the sequence and the output in their own folders
    We need to add 'desc' as placeholder to the sequence header for TMHMM to work
    '''
    print("Current environment: " + os.environ['CONDA_DEFAULT_ENV'])
    tmhmm_model = tmhmm_model_path
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
    
    print("Original sequence had these ambigious characters:")
    for letter in sequence:
        if letter in ambiguous_amino_acids:
            print("Ambigious letter found: " + letter)
    sequence = ''.join([i for i in sequence if i not in ambiguous_amino_acids])
    #save the dequence as a new fasta file. Make the unambigious dir first
    subprocess.run(["mkdir", "-p", outputfolder + "/TMHMM/" + id + "/unambigious"])
    temp_path_unambigious = outputfolder + "/TMHMM/" + id + "/unambigious/" + id + ".faa"
    print(sequence)
    #search for letter X in sequence and count them
    print("Updated sequence has these ambigious characters:")
    for letter in sequence:
        if letter in ambiguous_amino_acids:
            print("Ambigious letter found: " + letter)
    with open(temp_path_unambigious, "w") as f:
        f.write(">" + header + "\n")
        f.write(sequence)
    #run TMHMM again
    print("Running TMHMM for " + id + " without ambigious characters")
    print(["tmhmm", "-m", tmhmm_model, "-f" , temp_path_unambigious])
    result2 = subprocess.run(["tmhmm", "-m", tmhmm_model, "-f" , temp_path_unambigious], cwd=out_folder)
    if result2.returncode != 0:
        print("TMHMM failed again")
        #print the fasta file contens from temp_path_unambigious
        with open(temp_path_unambigious, "r") as f:
            print(f.read())

    summary_filepath = out_folder + "/" + id + ".summary"

    #list files in the new folder
    print("Files in " + out_folder + ":")
    subprocess.run(["ls", "-l", out_folder])

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