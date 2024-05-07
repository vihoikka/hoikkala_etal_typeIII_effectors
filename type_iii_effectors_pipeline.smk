'''
Snakemake pipeline for characterising CRISPR-Cas type III effectors and other proteins.
Example command:
snakemake --snakefile new_effectors.smk --use-conda --cores 40 --config cas10_tree_cluster="False" getGenomesBy="local" genome_mode="all" cas10_anchor="True" --rerun-triggers mtime

Options:
- getGenomesBy: "local" or "remote". If local, the program will look for genomes in a local folder. If remote, it will download genomes from NCBI.

NAR tidied down version (excluding redundant code)
'''

project = "test_170923" #defines the output folder

base_path = "/media/volume/st_andrews/new_effectors" + "/" + project #base path for storing the project folder
program_root = "/home/ubuntu/st_andrews/new_effectors" #root folder for the program. Some 3rd party programs are run directly from here

cas10_cluster_threshold = 0.9 #initial Cas10 clustering threshold
crispr_locus_interference_cutoff = 0 #cutoff for CRISPR loci interference completeness. Loci with less than this percentage of interference genes present are discarded. Default 0 (no filtering)

cas10_tree_cluster = str(config["cas10_tree_cluster"])
getGenomesBy = str(config["getGenomesBy"])
cas10_anchor = config["cas10_anchor"]
catyper_hmm_evalue = "1e-10"

hmm_msa_folder = "/media/volume/st_andrews/databases/cATyper/hmm/sep23/profiles" #folder names within this folder are used for creating effector dictionary
hmm_database_folder = "/media/volume/st_andrews/databases/cATyper/hmm/sep23/concatenated_profiles" #contains the concatenated and hmmpressed hmm profiles
hmm_database_file = "all.hmm" #filename for the concatenated hmmpressed hmm profiles

validated_effectors_hmm_db_path = "/media/volume/st_andrews/databases/validated_new_effectors/concatenated_profiles/all.hmm"

modified_cas10_hd = "/media/volume/st_andrews/new_effectors/manual_inputs/HD_HMM.msa" #modified cas10 hmm profile for hmmsearch

temperature_data = "/media/volume/st_andrews/new_effectors/manual_inputs/200617_TEMPURA.csv"

if cas10_anchor == False:
    prefiltering_wildcards = "05_host_genomes"
    prefiltering_host_genomes = "06_host_genomes"
elif cas10_anchor == True:
    prefiltering_wildcards = "02_host_wildcards"
    prefiltering_host_genomes = "03_host_genomes"

def aggregate_crisprcas(wildcards):
    if getGenomesBy == "local":
        checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{i}.txt")).i
    elif getGenomesBy == "remote":
        checkpoint_output = checkpoints.getHosts.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{sample}/{i}_genome.fna")).i
    return expand(base_path + "/07_cctyper/{i}/{i}_renaming.done", i=ivals)

def aggregate_host_genomes(wildcards):
    if getGenomesBy == "local":
        checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{i}.txt")).i
    elif getGenomesBy == "remote":
        checkpoint_output = checkpoints.getHosts.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{sample}/{i}_genome.fna")).i
    return expand(base_path + "/02_genome_wildcards/{i}.txt", i=ivals)

def aggregate_cas10_sequences(wildcards):
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/09_crispr_iii_CorA/loci/{c}/{c}_Cas10.faa", c=cvals)


def aggregate_renamed_crisprs(wildcards):
    checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards.get(**wildcards).output[0]
    ivals = glob_wildcards(os.path.join(checkpoint_output,"{j}.txt")).j
    return expand(base_path + "/07_cctyper/{j}/{j}_renaming.done", j=ivals)

def aggregate_download_genomes(wildcards):
    if getGenomesBy == "local":
        checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{i}.txt")).i
    elif getGenomesBy == "remote":
        checkpoint_output = checkpoints.getHosts.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{sample}/{i}_genome.fna")).i
    return expand(base_path + "/06_host_genomes/{i}/{i}_genome.fna", i=ivals)


    if getGenomesBy == "local":
        if cas10_anchor == True:
            checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards_postcas10.get(**wildcards).output[0]
        else:
            checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{j}.txt")).j
    elif getGenomesBy == "remote":
        checkpoint_output = checkpoints.getHosts.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{sample}/{i}_genome.fna")).i
    print(expand(base_path + "/19_16S_seq/{j}_16S.fna", j=ivals))
    return expand(base_path + "/19_16S_seq/{j}_16S.fna", j=ivals)

def aggregate_cas10_booleans(wildcards):
    if getGenomesBy == "local":
        checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{i}.txt")).i
    elif getGenomesBy == "remote":
        checkpoint_output = checkpoints.getHosts.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{sample}/{i}_genome.fna")).i
    return expand(base_path + "/062_genomes_cas10/{i}/{i}_cas10_boolean.tsv", i=ivals)

def aggregate_cas10_seq_postboolean(wildcards):
    if getGenomesBy == "local":
        checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{i}.txt")).i
    elif getGenomesBy == "remote":
        checkpoint_output = checkpoints.getHosts.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{sample}/{i}_genome.fna")).i
    return expand(base_path + "/063_cas10_seq/{i}/{i}_cas10.faa", i=ivals)


def aggregate_cas10_sequences_prior_to_1st_clustering(wildcards):
    if getGenomesBy == "local":
        checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{i}.txt")).i
    elif getGenomesBy == "remote":
        checkpoint_output = checkpoints.getHosts.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{sample}/{i}_genome.fna")).i
    return expand(base_path + "/063_cas10_seq/{i}/{i}_cas10.faa", i=ivals)


def aggregate_cas10_genomes(wildcards):
    if getGenomesBy == "local":
        checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards_postcas10.get(**wildcards).output[0]
        jvals = glob_wildcards(os.path.join(checkpoint_output,"{j}.txt")).j
    elif getGenomesBy == "remote":
        checkpoint_output = checkpoints.getHosts.get(**wildcards).output[0]
        jvals = glob_wildcards(os.path.join(checkpoint_output,"{sample}/{j}_genome.fna")).j
    return expand(base_path + "/03_postfiltering_genome_wildcards/{j}.txt", j=jvals)


    if getGenomesBy == "local":
        checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards_postcas10.get(**wildcards).output[0]
        jvals = glob_wildcards(os.path.join(checkpoint_output,"{j}.txt")).j
    elif getGenomesBy == "remote":
        checkpoint_output = checkpoints.getHosts.get(**wildcards).output[0]
        jvals = glob_wildcards(os.path.join(checkpoint_output,"{sample}/{j}_genome.fna")).j
    return expand(base_path + "/30_clustertable/{j}/{j}.tsv", j=jvals)

def aggregate_typeIII_info(wildcards):
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/09_crispr_iii_CorA/loci/{c}/{c}_crispr_iii_info.tsv", c=cvals)

def aggregate_taxInfo(wildcards):
    if getGenomesBy == "local":
        if cas10_anchor == True:
            checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards_postcas10.get(**wildcards).output[0]
        else:
            checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{j}.txt")).j
    elif getGenomesBy == "remote":
        checkpoint_output = checkpoints.getHosts.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{sample}/{i}_genome.fna")).i
    return expand(base_path + "/06_host_genomes/{j}/{j}_taxon.txt", j=ivals)

def aggregate_renamed(wildcards):
    if getGenomesBy == "local":
        if cas10_anchor == True:
            checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards_postcas10.get(**wildcards).output[0]
        else:
            checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{j}.txt")).j
    elif getGenomesBy == "remote":
        checkpoint_output = checkpoints.getHosts.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{sample}/{j}_genome.fna")).j
    return expand(base_path + "/07_cctyper/{j}/{j}_renaming.done", j=ivals)

def aggregate_cATyper_hmm(wildcards):
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/10_cATyper_hmm/{c}/{c}_cATyper_hmm_rows.tsv", c=cvals)

def aggregate_validated_new_effectors_hmm(wildcards):
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/60_validated_new_effectors/{c}/{c}_validated_new_effectors_hmm.tsv", c=cvals)

def aggregate_validated_new_effectors_analysis(wildcards):
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/61_validated_new_effectors_analysis/{c}/{c}_cATyper_results.tsv", c=cvals)

def aggregate_cATyper_analysis(wildcards):
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/11_cATyper_analysis/{c}/{c}_cATyper_results.tsv", c=cvals)

def aggregate_cATyper_etp(wildcards):
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/11_cATyper_analysis/{c}/{c}_effector_to_protein.tsv", c=cvals)

def aggregate_cATyper_pte(wildcards):
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/11_cATyper_analysis/{c}/{c}_protein_to_effector.tsv", c=cvals)

def aggregate_validated_new_effectors_etp(wildcards):
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/61_validated_new_effectors_analysis/{c}/{c}_effector_to_protein.tsv", c=cvals)

def aggregate_validated_new_effectors_pte(wildcards):
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/61_validated_new_effectors_analysis/{c}/{c}_protein_to_effector.tsv", c=cvals)

def aggregate_CorA_sequences_catyper(wildcards):
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    list_of_coras = []
    for c in cvals:
        cas10_dir_path = base_path + "/09_crispr_iii_CorA"
        cora_dir_path = base_path + "/11_cATyper_analysis/" + c
        if (os.path.exists(cas10_dir_path)) & (os.path.exists(cora_dir_path)):
            if os.path.getsize(base_path + "/09_crispr_iii_CorA/loci/" + c + "/" + c + "_Cas10.faa") > 0:
                list_of_coras.append(c)
                if not os.path.exists(base_path + "/11_cATyper_analysis/" + c + "/cora.faa"):
                    print("Touching an empty CorA")
                    open(base_path + "/11_cATyper_analysis/" + c + "/cora.faa", 'a').close()
    print("Length of list of coras: " + str(len(list_of_coras)))
    #return the list of coras after adding the filepaths to them
    fullpath_coras = []
    for locus in list_of_coras:
        fullpath_coras.append(base_path + "/11_cATyper_analysis/" + locus + "/cora.faa")
    print(fullpath_coras)
    return fullpath_coras

def aggregate_CorA_sequences(wildcards):
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/09_crispr_iii_CorA/loci/{c}/CorA.faa", c=cvals)


def aggregate_unknowns(wildcards):
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/30_unknown_effectors/{c}/{c}_unknown_proteins.faa", c=cvals)

def aggregate_unknowns_locus_info(wildcards):
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/30_unknown_effectors/{c}/{c}_locus_unknown_info.tsv", c=cvals)

def aggregate_crispr_locus_proteins(wildcards):
    '''
    Aggregates the outputs of the crispr_locus_proteins rule.
    '''
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/072_crispr_locus_proteins/{c}/{c}_crispr_locus_proteins.faa", c=cvals)

def aggregate_locus_viz(wildcards):
    '''
    Aggregates the outputs of the crispr_locus_proteins rule.
    '''
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/90_locus_viz/{c}/{c}_viz.png", c=cvals)

def aggregate_group4_blasts(wildcards):
    '''
    Aggregates the outputs of the crispr_locus_proteins rule.
    '''
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/group4_prober/{c}/{c}.blast", c=cvals)

def aggregate_known_effector_wildcards(wildcards):
    checkpoint_output = checkpoints.effector_wildcarder.get(**wildcards).output[0]
    effector_vals = glob_wildcards(os.path.join(checkpoint_output,"{effector}.eff")).effector
    return expand(base_path + "/45_effector_tree/{effector}_tree.txt", effector=effector_vals)

def aggregate_known_effector_wildcarder(wildcards):
    checkpoint_output = checkpoints.effector_wildcarder.get(**wildcards).output[0]
    effector_vals = glob_wildcards(os.path.join(checkpoint_output,"{effector}.eff")).effector
    return expand(base_path + "/40_known_effector_wildcards/{effector}.eff", effector=effector_vals)

def aggregate_known_effector_annotations(wildcards):
    checkpoint_output = checkpoints.effector_wildcarder.get(**wildcards).output[0]
    effector_vals = glob_wildcards(os.path.join(checkpoint_output,"{effector}.eff")).effector
    return expand(base_path + "/43_effector_hmmer_analysis/{effector}/{effector}_sorted_filtered_mapped.tsv", effector=effector_vals)

def aggregate_cATyper_hhsuite(wildcards):
    checkpoint_output = checkpoints.effector_wildcarder.get(**wildcards).output[0]
    effector_vals = glob_wildcards(os.path.join(checkpoint_output,"{effector}.eff")).effector
    return expand(base_path + "/42_effector_hhsuite/{effector}/{effector}_all.tsv", effector=effector_vals)

def aggregate_cATyper_hhsuite_parser(wildcards):
    checkpoint_output = checkpoints.effector_wildcarder.get(**wildcards).output[0]
    effector_vals = glob_wildcards(os.path.join(checkpoint_output,"{effector}.eff")).effector
    return expand(base_path + "/42_effector_hhsuite/{effector}/{effector}_hhsuite_parsed.tsv", effector=effector_vals)

def aggregate_cATyper_hhsuite_parser_cogs(wildcards):
    checkpoint_output = checkpoints.effector_wildcarder.get(**wildcards).output[0]
    effector_vals = glob_wildcards(os.path.join(checkpoint_output,"{effector}.eff")).effector
    return expand(base_path + "/42_effector_hhsuite_cogs/{effector}/{effector}_hhsuite_parsed_cogs.tsv", effector=effector_vals)

def aggregate_cora_neighbourhoods(wildcards):
    checkpoint_output = checkpoints.cora_neighbourhood_preparation.get(**wildcards).output[0]
    cora_loci = glob_wildcards(os.path.join(checkpoint_output,"{cora_locus}_crispr_locus_proteins.faa")).cora_locus
    return expand(base_path + "/52_cora_neighbour_analysis/{cora_locus}/neighbourhood_results.tsv", cora_locus=cora_loci)


#Cas10 blast parameters
blast_headers = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle sacc"
blast_headers_group4 = "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tlocus"
cas10_e_value = "1e-20"
corA_hmm_value = "1e-20"

group4_pfam_headers = "target_name\taccession\tquery_name\taccession\tE-value\tscore\tbias\tE-value\tscore\tbias\texp\treg\tclu\tov\tenv\tdom\trep\tinc\tdescription_of_target"

print("Starting CRISPR-Cas type III characterisation pipeline")

rule all: 
    input: base_path + "/done"

#global filepaths
genomes_folder = "/media/volume/st_andrews/databases/ncbi_genomes/ncbi_dataset/data"
archaea_folder = "/media/volume/st_andrews/cas10_corA/genomes/archaea/ncbi_dataset/data"
genome_count = 1000
subsampling_seed = 666

genomes_json = os.path.join(genomes_folder,"assembly_data_report.jsonl")

rule write_down_genome_info:
    '''
    Outputs a text file genome names based on folders in genomes_folders
    '''
    output: base_path + "/01_genomelist/genomelist.txt"
    shell:
        """
        cd {genomes_folder}
        find * -maxdepth 0 -type d > {output}
        """

rule annotate_bacteria_and_archaea_domains:
    """
    Marks down whether a sample is an archaeon or a bacterium
    1. Create a .csv file in the **archaea folder**. This two-column file shows that each sample ID there is, is an archaeon
	2. Create another .csv file in the **bacteria folder**. Note that bacteria folder also contains the archaea (they have been copied there to simplify the fetching of genomes).
        This csv will claim that all samples, including the archaea, are bacteria.
	3. Use the archaea-csv file to **reannotate** the .csv file created in step two, making all archaea appear archaea
	4. The final output file will then be the one created in step three

    This table can then be used in subsequent rules by merging it to any output table using sample name as common identifier.

    *** Note that this solution is not ideal and will break down easily if the input genomes are in a different format.
    So make sure that the bacteria folder contains **all** samples and that the archaea folder only contains archaea ***
    """
    output: base_path + "/01_genomelist/annotated_genomelist.csv"
    params: 
        archaea_csv = base_path + "/01_genomelist/archaea.csv",
        bacteria_csv = base_path + "/01_genomelist/bacteria.csv"
    shell:
        """
        cd {archaea_folder}
        find * -maxdepth 0 -type d > {params.archaea_csv}

        cd {genomes_folder}
        find * -maxdepth 0 -type d > {params.bacteria_csv}

        cd {program_root}
        python3 scripts/annotate_bacteria_and_archaea.py --archaea {params.archaea_csv} --bacteria {params.bacteria_csv} --out {output}
        """


if getGenomesBy == "local":
    checkpoint expand_host_genomelist_to_wildcards:
        '''
        Takes as input a list of genome accessions separated by newline. These genomes are then turned into wildcards.
        In random mode, uses subsampling to randomly sample a set of genomes.
        '''
        input:
            genome_info = rules.write_down_genome_info.output,
            domain_annotations = rules.annotate_bacteria_and_archaea_domains.output
        output: directory(base_path + "/" + prefiltering_wildcards)
        run:
            import pandas as pd
            if not os.path.exists(str(output)):
                os.makedirs(str(output))
            genomefiles = pd.read_csv(str(input.genome_info), sep = "\t", header = None)
            #print("HELLOO " + genomefiles)
            genomes = [] #final list of accessions

            #Random mode. Subsamples n genomes from all genomes randomly with seed.
            if config["genome_mode"] == "random":
                subsample_genomes = genomefiles.sample(n = genome_count, random_state=subsampling_seed)
                genomes = subsample_genomes[0]

            elif config["genome_mode"] == "all":
                genomes = genomefiles[0]

            #Taxid id mode. NCBI datasets json is probed to find representatives of a taxid.
            elif config["genome_mode"] == "taxid":
                json_df = pd.read_json(genomes_json, lines=True)
                json_df = json_df.join(pd.json_normalize(json_df['organism'])).drop('organism', axis='columns') #expand pandas dictionary column into new columns in the actual df
                taxidlistlist = []

                with open(config["taxidlistfile"]) as f: #read taxid list file
                    taxids_comma = f.read()
                    idlist = taxids_comma.split(",")

                taxidlist = pd.DataFrame(idlist) #convert list to df
                #taxidlist = pd.read_csv(config["taxidlistfile"], sep=",", header=None)
                taxidlist.columns = ["taxId"] #add column to the df
                json_df['taxId'] = json_df['taxId'].astype("string")
                taxidlist['taxId'] = taxidlist['taxId'].astype("string")
                chosen = pd.merge(json_df, taxidlist, on="taxId", how="inner")
                #chosen = json_df.loc[json_df["taxId"] == int(config["taxid"])]
                print("Shape of subsampled dataframe using TaxIDs from " + str(config["taxidlistfile"]) + ": " + str(chosen.shape))
                genomes = chosen["accession"].values.tolist()


            #add genome to list if .faa, .gff and .fna files are present
            for i in genomes:
                if (os.path.isfile(os.path.join(genomes_folder,i,"protein.faa")) and (os.path.isfile(os.path.join(genomes_folder,i,"genomic.gff")))): #make sure it's annotated
                    for fname in os.listdir(os.path.join(genomes_folder,i)): #checking fna is more complicated because the filename is not consistent
                        if fname.endswith('.fna'):
                            with open(str(output) + "/" + str(i) + '.txt', 'w') as f:
                                f.write(str(output))
        


    rule download_genomes:
        '''
        Based on the wildcards from expand_host_genomelist_to_wildcards,
        makes symlinks for the desired genomes from the genomes_folder
        into the working folder (prefiltering_host_genomes).
        '''
        input:
            ivalue = base_path + "/" + prefiltering_wildcards + "/{i}.txt",
        output: 
            fna = base_path + "/" + prefiltering_host_genomes + "/{i}/{i}_genome.fna",
            faa = base_path + "/" + prefiltering_host_genomes + "/{i}/{i}_proteins.faa",
            gff = base_path + "/" + prefiltering_host_genomes + "/{i}/{i}_features.gff"
        params: #these are mostly deprecated
            folder = base_path + "/" + prefiltering_host_genomes,
        conda: "envs/ncbidownload.yaml"
        #threads: 4
        log:
            out = base_path + "/logs/" + prefiltering_host_genomes + "/{i}.log",
            err = base_path + "/logs/" + prefiltering_host_genomes + "/{i}.err"
        shell:
            '''
            echo "Creating symlinks for chosen genomes"
            ln -s {genomes_folder}/{wildcards.i}/*.fna {output.fna}
            ln -s {genomes_folder}/{wildcards.i}/*.faa {output.faa}
            ln -s {genomes_folder}/{wildcards.i}/*.gff {output.gff}
            '''


elif getGenomesBy == "remote":
    species_space = str(config["species"])
    checkpoint getHosts:
        '''
        Downloads desired host genomes from NCBI using the wonderful python program ncbi-genome-download.
        Files are then gunzipped and named after their parent folder
        '''
        output: directory(base_path + "/" + prefiltering_host_genomes)
        conda: "envs/ncbidownload.yaml"
        params:
            parallels = "40",
            folder = base_path + "/06_host_genomes"
        shell:
            '''
            mkdir -p {params.folder}
            cd {params.folder}
            count=$(ncbi-genome-download --dry-run --genera "{species_space}" bacteria | wc -l)
            echo "Will download $count genomes from {species_space}. Note that the download will likely fail if the dataset is large (>1000 genomes). In such a case, just restart the script."
            printf 'Do you want to continue? (y/n)? '
            read -p "Do you want to continue? (y/n) " -n 1 -r
            if [[ $REPLY =~ ^[Yy]$ ]] ;then
                ncbi-genome-download --genera "{species_space}" bacteria --parallel {params.parallels} --formats fasta,protein-fasta,gff,assembly-report
                mv refseq/bacteria/* {params.folder}
                find {params.folder} -type f -name '*.gz' | tqdm | xargs gunzip
                echo renaming
                rename 's!(.*)/(.*)genomic\.fna!$1/$1_genome.fna!' */*genomic.fna
                rename 's!(.*)/(.*)protein\.faa!$1/$1_proteins.faa!' */*protein.faa
                rename 's!(.*)/(.*)genomic\.gff!$1/$1_features.gff!' */*genomic.gff
                rename 's!(.*)/(.*)assembly_report\.txt!$1/$1_report.txt!' */*assembly_report.txt
            else
                exit 1
            fi

            '''


if cas10_anchor == True: #if we filter analyzable genomes by the presence of Cas10. Uses Cas10 HMM profiles from CCtyper.
    print("Running Cas10 stuff")
    rule Cas10_genomes:
        '''
        First filter all sequences to get only >500 AA long proteins. This prevents the inclusion of
        truncated Cas10s or those that are cut artifically due to contig ending. Also reduces HMM searches.
        
        Then:
        Searches the filtered proteins of host against a user-specified preprepared hmm profile.
        Strips file of "#" -rows and checks if there is anything left. If there is, the gene has been found.
        Boolean output format:
        i,True/False,cas10_accession
        '''
        input:
            proteins = base_path + "/" + prefiltering_host_genomes + "/{i}/{i}_proteins.faa"
        output:
            hmm = base_path + "/062_genomes_cas10/{i}/{i}_cas10.tsv",
            boolean = base_path + "/062_genomes_cas10/{i}/{i}_cas10_boolean.tsv"
        params:
            proteins_filt = base_path + "/062_genomes_cas10/{i}/{i}_proteins_lengthfiltered.faa",
            out = base_path + "/062_genomes_cas10/{i}/{i}_temp.out",
            rows1 = base_path + "/062_genomes_cas10/{i}/{i}_temp_rows_1.out",
            cas10_db = "/media/volume/st_andrews/cas10_corA_2/profiles/Profiles/all_cas10s.hmm",
            rows = base_path + "/062_genomes_cas10/{i}/{i}_temp_rows.out",
            headers = base_path + "/062_genomes_cas10/{i}/{i}_headers.out",
            all_data = base_path + "/062_genomes_cas10/{i}/{i}_all_data.out",
        conda: "envs/hmmer.yaml"
        shell:
            '''
            echo "Running rule Cas10_genomes"
            echo -e 'target_name\taccession\ttlen\tquery_name\taccession\tqlen\tE-value_fullseq\tscore_fullseq\tbias_fullseq\t#_domain\tof_domain\tc-Evalue_domain\tscore_domain\tbias_domain\t_hmm_from\thmm_to\t_ali_from\tali_to\tenv_from\tenv_to\tacc\tdescription' > {output.hmm}
            if [ -s "{input.proteins}" ]; then
                cat {input.proteins} | seqkit seq -m 500 > {params.proteins_filt}
                if [ -s "{params.proteins_filt}" ]; then
                    hmmscan --domtblout {params.out} --cpu 1 -E {corA_hmm_value} {params.cas10_db} {params.proteins_filt}  &> /dev/null
                    grep -v "#" {params.out} > {params.rows}||:
                    head -1 {params.rows} > {params.rows1}
                    echo "id,cas10_boolean,cas10_acc" > {output.boolean}
                    if [ -s {params.rows1} ]; then
                        cat {params.rows1} >> {output.hmm}
                        ACC=$(awk -F ' ' '{{print $4}}' {params.rows1})
                        echo "{wildcards.i},True","${{ACC}}" > {output.boolean}
                    else
                        echo "{wildcards.i},False,-" > {output.boolean}
                    fi
                    touch {output.hmm}
                else
                    echo "{wildcards.i},False,-" > {output.boolean}
                fi
            else
                echo "{wildcards.i},False,-" > {output.boolean}
            fi
            '''

        
    rule extract_cas10_sequence:
        input:
            boolean = rules.Cas10_genomes.output.boolean,
            proteins = base_path + "/" + prefiltering_host_genomes + "/{i}/{i}_proteins.faa"
        output: base_path + "/063_cas10_seq/{i}/{i}_cas10.faa"
        shell: 
            '''
            CAS10_ID=$(awk -F ',' '{{print $3}}' {input.boolean})
            if [ "$CAS10_ID" != "-" ]; then
                echo ">"{wildcards.i} > {output}
                seqkit grep -r -p ${{CAS10_ID}} {input.proteins} | seqkit seq -w 0 | tail -n +2 >> {output}
            else
                touch {output}
            fi
            '''

    rule concat_cluster_cas10_proteins:
        '''
        Takes all Cas10 sequences, concatenates them into one fasta and
        clusters this fasta at cas10_cluster_threshold (0-1) similarity.
        Greps only representative lines "marked by *"
        A python script extracts the accession from these lines into output.reps
        '''
        input: aggregate_cas10_seq_postboolean
        output:
            all_cas10 = base_path + "/064_cas10_clusters/cas10_all.faa",
            clusters = base_path + "/064_cas10_clusters/cas10_clust.faa.clstr",
            proteins = base_path + "/064_cas10_clusters/cas10_clust.faa",
            reps = base_path + "/064_cas10_clusters/cas10_unique_genomes.txt",
        params:
            clusterlines = base_path + "/064_cas10_clusters/clusterlines.txt"
        threads: 40
        shell:
            '''
            echo "concatenating cas10 sequences for clustering"
            find '{base_path}/063_cas10_seq' -maxdepth 2 -type f -wholename '*/*_cas10.faa' -print0 | xargs -0 cat > {output.all_cas10}
            echo clustering
            cd-hit -i {output.all_cas10} -o {output.proteins} -c {cas10_cluster_threshold} -n 5 -d 0 -M 16000 -T {threads}
            grep "*" {output.clusters} > {params.clusterlines}
            python scripts/getClusterRep.py --input {params.clusterlines} --output {output.reps}
            '''

    checkpoint expand_host_genomelist_to_wildcards_postcas10:
        '''
        '''
        input:
            genomes = rules.concat_cluster_cas10_proteins.output.reps
        output: directory(base_path + "/03_postfiltering_genome_wildcards")
        run:
            print("Expanding host list to wildcards...")

            list_of_genomes_to_remove = ["GCA_019977735.1"] #use this list to manually exclude any unwanted genomes

            if not os.path.exists(str(output)):
                os.makedirs(str(output))

            with open(str(input.genomes)) as filehandle:
                genomes = [line.rstrip('\n') for line in filehandle]

            for i in genomes:
                if i not in list_of_genomes_to_remove:
                    sample = str(str(i).split(",")[0])
                    with open(str(output) + "/" + sample + '.txt', 'w') as f:
                        f.write(sample)



    rule download_genomes_postCas10:
        '''
        After filtering by presence of Cas10, create new symlinks for such genomes.
        '''
        input: 
            genomes = base_path + "/03_postfiltering_genome_wildcards/{j}.txt",
        output: 
            fna = base_path + "/06_host_genomes/{j}/{j}_genome.fna",
            faa = base_path + "/06_host_genomes/{j}/{j}_proteins.faa",
            gff = base_path + "/06_host_genomes/{j}/{j}_features.gff"
        params: #these are mostly deprecated
            folder = base_path + "/06_host_genomes",
            zipped = base_path + "/06_host_genomes/{j}.zip",
            original_fna = base_path + "/06_host_genomes/{j}/*.fna", #for checking successful download
            original_prot = base_path + "/06_host_genomes/{j}/protein.faa",#for checking successful download
            original_gff = base_path + "/06_host_genomes/{j}/genomic.gff", #for checking successful download
        conda: "envs/ncbidownload.yaml"
        #threads: 4
        log:
            out = base_path + "/logs/06_host_genomes/{j}.log",
            err = base_path + "/logs/06_host_genomes/{j}.err"
        shell:
            '''
            ln -s {genomes_folder}/{wildcards.j}/*.fna {output.fna}
            ln -s {genomes_folder}/{wildcards.j}/*.faa {output.faa}
            ln -s {genomes_folder}/{wildcards.j}/*.gff {output.gff}
            '''

rule getTaxInfo:
    input: base_path + "/03_postfiltering_genome_wildcards/{j}.txt"
    output:
        taxon = base_path + "/06_host_genomes/{j}/{j}_taxon.txt"
    run:
        import pandas as pd
        json_df = pd.read_json(genomes_json, lines=True)
        json_df = json_df.join(pd.json_normalize(json_df['organism'])).drop('organism', axis='columns') #expand pandas dictionary column into new columns in the actual df
        json_df = json_df.join(pd.json_normalize(json_df['averageNucleotideIdentity'])).drop('averageNucleotideIdentity', axis='columns') #expand pandas dictionary column into new columns in the actual df 
        try:
            row = json_df.loc[json_df["accession"] == wildcards.j]
            species = row['submittedSpecies'].iloc[0].strip("][")
            genus = species.split(" ")[0]
            with open(output.taxon, "w") as f:
                f.write("Sample\tGenus\tSpecies\n")
                f.write(wildcards.j + "\t" + genus + "\t" + species + "\n")
        except:
            with open(output.taxon, "w") as f:
                            f.write(wildcards.j + "\tUnknown\tUnknown\n")
            pass


rule concat_taxInfo:
    input: aggregate_taxInfo
    output: base_path + "/06_host_genomes/taxInfo.txt"
    shell:
        """
        echo "Sample\tgenus\tspecies\n" > {output}
        find '{base_path}/06_host_genomes/' -maxdepth 2 -type f -wholename '*_taxon.txt' -print0 | xargs -0 cat >> {output}
        """

rule CRISPRCasTyper:
    '''
    Types CRISPR-Cas loci in the genomes based on fasta.
    Takes output of the checkpoint. Automatically scales to the number of files
    produced by the checkpoint, so no need to call for the aggregator function here.

    '''
    input: base_path + "/06_host_genomes/{j}/{j}_genome.fna"
    output: base_path + "/07_cctyper/{j}/{j}.done"
    conda: "envs/CRISPRCasTyper.yaml"
    log: base_path + "/logs/07_cctyper/{j}/{j}_cctyper.log"
    params:
        outdir = base_path + "/07_cctyper/{j}"
    shell:
        '''
        rm -rf {params.outdir}
        cctyper '{input}' '{base_path}/07_cctyper/{wildcards.j}' --prodigal single 2>{log}
        touch '{output}'
        '''

rule CRISPRCasTyper_rename:
    '''
    UPDATED VERSION. This uses the CRISPR-Cas.tab file to include only fully functional CRISPR-Cas loci
    The sed adds assembly information to the results file, if it exists (which is originally based on just the nucleotide accessions)
    '''
    input: rules.CRISPRCasTyper.output
    output: base_path + "/07_cctyper/{j}/{j}_renaming.done"
    params:
        outdir = base_path + "/07_cctyper/{j}"
    shell:
        '''
        if test -e "{params.outdir}/CRISPR_Cas.tab";then
            mv {params.outdir}/CRISPR_Cas.tab {params.outdir}/CRISPR_Cas_old.tab
            sed '1s/$/\tassembly/; 2,$s/$/\t{wildcards.j}/' {params.outdir}/CRISPR_Cas_old.tab > {params.outdir}/CRISPR_Cas.tab
            rm {params.outdir}/CRISPR_Cas_old.tab
        fi
        touch {output}
        '''

rule concat_renamed_crisprs:
    '''
    This rule aggregates the outputs of the CRISPRCasTyper_rename rule.
    It is necessary to do this because the CRISPRCasTyper_rename rule is called by the checkpoint,
    and the checkpoint does not automatically aggregate the outputs of the rule.
    '''
    input: aggregate_renamed
    output: base_path + "/07_cctyper/renaming.done"
    shell:
        '''
        touch {output}
        ''' 

checkpoint type_iii_wildcarder:
    '''
    Extracts only type III loci from the CCTyper output
    This checkpoint examines outputs from CCTyper and creates new wildcards based on complete loci.
    This marks a point where the pipeline no longer looks at individual strains, but rather at the CRISPR-Cas loci individually.
    As outputs, this checkpoint extracts locus-specific information from CCTyper outputs (CRISPR-Cas.tab and cas_operons.tab)
    NOTE: now excludes all loci with >1 Cas10s
    '''
    input:
        cctyper_done = aggregate_renamed,
    output: directory(base_path + "/071_cctyper_loci")
    params:
        interference_cutoff = crispr_locus_interference_cutoff, #0-100. If the interference score is lower than this, the locus is discarded
        cctyper_folder = base_path + "/07_cctyper",
    shell:
        '''
        python3 scripts/loci_wildcarder.py --input_folder {params.cctyper_folder} --output_folder {output} --interference_cutoff {params.interference_cutoff}
        '''


rule cATyper_hmm_search:
    '''
    Runs HMMscan for all proteins in each locus. As database we use the concatenated hmm database of all known effectors.
    The main output file is hmm_rows, which contains all hmm hits for all proteins in the locus.
    The column target_name is the hmm hit and query_name is the protein accession
    '''
    input:
        crispr_positive_samples = base_path + "/071_cctyper_loci/{c}/cas_operons.tsv",
    output:
        contig_proteins = base_path + "/10_cATyper_hmm/{c}/{c}_contig_proteins.faa",
        temp_rows = base_path + "/10_cATyper_hmm/{c}/{c}_temp_rows.out",
        hmm_rows = base_path + "/10_cATyper_hmm/{c}/{c}_cATyper_hmm_rows.tsv",
        cora = base_path + "/10_cATyper_hmm/{c}/CorA.faa",
    params:
        outdir = base_path + "/10_cATyper_hmm/{c}",
        hmm_msa_folder = hmm_msa_folder,
        host_genomes_folder = base_path + "/06_host_genomes",
        hmm_profile = hmm_database_folder + "/" + hmm_database_file,
        temp_hmm = base_path + "/10_cATyper_hmm/{c}/{c}_temp.out",
    conda: "envs/hmmer.yaml"
    log:
        out = base_path + "/10_cATyper_hmm/logs/{c}.out",
        err = base_path + "/10_cATyper_hmm/logs/{c}.err",
    shell:
        '''
        python scripts/catyper_prepper_10.py --locus {wildcards.c} --cas_operons_file {input.crispr_positive_samples} --output_folder {params.outdir} --host_genomes_folder {params.host_genomes_folder} --mode pre_hmm 2> {log.err} 1> {log.out}
        echo "Running hmmscan" >> {log.out}
        hmmscan --domtblout {params.temp_hmm} --cpu 8 -E {catyper_hmm_evalue} {params.hmm_profile} {output.contig_proteins}  &> /dev/null
        echo "Removing commented rows" >> {log.out}
        grep -v "^#" {params.temp_hmm} > {output.temp_rows} ||:
        echo "Writing header" >> {log.out}
        echo -e 'target_name\taccession\ttlen\tquery_name\taccession\tqlen\tE-value_fullseq\tscore_fullseq\tbias_fullseq\t#_domain\tof_domain\tc-Evalue_domain\ti-Evalue_domain\tscore_domain\tbias_domain\t_hmm_from\thmm_to\t_ali_from\tali_to\tenv_from\tenv_to\tacc\tdescription' > {output.hmm_rows}
        echo "Checking if hits were found" >> {log.out}
        if [ -s {output.temp_rows} ]; then 
            echo "Hits found for {wildcards.c}" >> {log.out}
            cat {output.temp_rows} >> {output.hmm_rows}
        else
            echo "No hits found for {wildcards.c}" >> {log.out}
            touch {output.hmm_rows}
        fi
        touch {output.cora}
        '''

rule concatenate_cATyper_hmm:
    input: aggregate_cATyper_hmm
    output: base_path + "/10_cATyper_hmm/cATyper_all.tsv"
    shell:
        """
        awk '(NR == 1) || (FNR > 1)' {input} > {output}
        """

rule cATyper_analysis:
    '''
    This is the analysis step of cATyper. Note that we are still using the same python script as in the previous rule,
    but the mode is now set to "post_hmm".
    '''
    input:
        crispr_positive_samples = base_path + "/071_cctyper_loci/{c}/cas_operons.tsv",
        contig_proteins = base_path + "/10_cATyper_hmm/{c}/{c}_contig_proteins.faa",
        hmm_rows = base_path + "/10_cATyper_hmm/{c}/{c}_cATyper_hmm_rows.tsv",
    output:
        catyper = base_path + "/11_cATyper_analysis/{c}/{c}_cATyper_results.tsv",
        hmm_targets = base_path + "/11_cATyper_analysis/{c}/{c}_cATyper_hmm_targets.tsv",
        effector_to_protein = base_path + "/11_cATyper_analysis/{c}/{c}_effector_to_protein.tsv",
        protein_to_effector = base_path + "/11_cATyper_analysis/{c}/{c}_protein_to_effector.tsv",
        plottable_effector_positions = base_path + "/11_cATyper_analysis/{c}/{c}_plottable_effector_positions.tsv",
    params:
        outdir = base_path + "/11_cATyper_analysis/{c}",
        host_genomes_folder = base_path + "/06_host_genomes",
        hmm_msa_folder = hmm_msa_folder,
        tmhmm_model_path = "/media/volume/st_andrews/new_effectors/TM/TMHMM2.0.model",
    conda: "envs/hmmer.yaml"
    log:
        out = base_path + "/11_cATyper_analysis/logs/{c}.out",
        err = base_path + "/11_cATyper_analysis/logs/{c}.err",
    shell:
        '''
        echo "Running cATyper analysis" >> {log.out}
        echo "Listing targets" >> {log.out}
        ls {params.hmm_msa_folder}/*/*.hmm > {output.hmm_targets}
        python scripts/catyper_prepper_10.py --locus {wildcards.c} --cas_operons_file {input.crispr_positive_samples} --output_folder {params.outdir} --host_genomes_folder {params.host_genomes_folder} --mode post_hmm --hmm_targets {output.hmm_targets} --hmm_rows {input.hmm_rows} --catyper_out {output.catyper} --catyper_type "known_effector" --effector_plot_data {output.plottable_effector_positions} --tmhmm_model_path {params.tmhmm_model_path} 2> {log.err} 1> {log.out}
        '''

rule concatenate_cATyper_analysis:
    input: aggregate_cATyper_analysis
    output: base_path + "/11_cATyper_analysis/cATyper_all.tsv"
    shell:
        """
        awk '(NR == 1) || (FNR > 1)' {base_path}/11_cATyper_analysis/*/*_cATyper_results.tsv > {output}
        touch {output}
        """

rule concatenate_cATyper_analysis_effector_scores:
    '''
    Concatenates info on the effector hits from cATyper analysis
    '''
    input:
        protein_to_effector = aggregate_cATyper_pte,
        effector_to_protein = aggregate_cATyper_etp
    output:
        protein_to_effector_concatenated = base_path + "/11_cATyper_analysis/cATyper_protein_to_effector.tsv",
        effector_to_protein_concatenated = base_path + "/11_cATyper_analysis/cATyper_effector_to_protein.tsv",
    shell:
        '''
        cat {input.protein_to_effector} > {output.protein_to_effector_concatenated}
        cat {input.effector_to_protein} > {output.effector_to_protein_concatenated}
        '''

rule analyse_cATyper_effector_scores:
    '''
    The algorithm that analyses how well each hmm profile performs in cATyper.
    '''
    input:
        pte = rules.concatenate_cATyper_analysis_effector_scores.output.protein_to_effector_concatenated,
        etp = rules.concatenate_cATyper_analysis_effector_scores.output.effector_to_protein_concatenated,
    output:
        effector_scores = base_path + "/11_cATyper_analysis/cATyper_effector_scores.tsv",
        effector_scores_summary = base_path + "/11_cATyper_analysis/cATyper_effector_scores_summary.tsv",
        #effector_scores_plot1 = base_path + "/11_cATyper_analysis/cATyper_effector_scores_plot1.png",
    params:
        outdir = base_path + "/11_cATyper_analysis"
    shell:
        '''
        python scripts/catyper_effector_scores.py --pte {input.pte} --etp {input.etp} --output {output.effector_scores} --output_summary {output.effector_scores_summary} --outdir {params.outdir}
        touch {output.effector_scores}
        touch {output.effector_scores_summary}
        '''

checkpoint effector_wildcarder:
    '''
    This rule creates wildcards for each previosly known effector.
    Wildcards are used downstream when creating phylogenetic trees etc for each effector.
    '''
    output:
        directory(base_path + "/40_known_effector_wildcards")
        #hmm_targets = base_path + "/40_known_effector_wildcards/{effector}/{effector}.txt"
    params:
        hmm_targets = base_path + "/40_known_effector_wildcards/alltargets.txt",
        hmm_msa_folder = hmm_msa_folder,
        outdir = base_path + "/40_known_effector_wildcards",
    run:
        import glob
        #Equivalent of `mkdir`
        Path(params.outdir).mkdir(parents=True, exist_ok=True)
        
        #Equivalent of `ls {hmm_msa_folder}/*/*.hmm > {hmm_targets}`
        with open(params.hmm_targets, 'w') as f:
            files = glob.glob(f"{params.hmm_msa_folder}/*/*.hmm") #get all files in the hmm_msa_folder
            f.write('\n'.join(files))

        with open(params.hmm_targets, 'r') as f:
            lines = f.readlines()
            for line in lines:
                filename = os.path.basename(line) # get filename from path
                effector = re.split("\.", filename)[0] #remove extension
                if "#" not in line:
                    effector = re.split("_", effector)[1]
                elif "#" in line: #in case a single effector is characterised by multiple hmm profiles, we mark the filenames with # followed by the effector, e.g ca6_csm6-ca6_italicus#csm6-ca6
                    effector = re.split("#", effector)[1]
                print(effector)
                with open(f'{params.outdir}/{effector}.eff', 'w') as f:
                    f.write(effector)

rule effector_fetcher:
    '''
    Fetches protein sequences for an effector from catyper_analysis outputs
    '''
    input:
        aggregate_known_effector_wildcarder,
        catyper_done = rules.concatenate_cATyper_analysis.output, #this indicates we can safely probe the folders for effector sequences
    output:
        multifasta = base_path + "/41_known_effector_mf/{effector}.faa",
    shell:
        '''
        effector_name={wildcards.effector}
        echo "Fetching effector sequences for {wildcards.effector}"
        touch {output.multifasta}
        find {base_path}/11_cATyper_analysis/ -name \"*${{effector_name}}*.faa\" -exec cat {{}} \; >> {output.multifasta}
        '''

rule effector_hmmer:
    '''
    Uses hmmscan against pfam to annotate the discovered effectors' domains.
    Output used when plotting tree with R
    '''
    input:
        multifasta = rules.effector_fetcher.output.multifasta
    output:
        raw_hmm = base_path + "/42_effector_hmmer/{effector}.hmm",
        #filtered_hmm_tsv = base_path + "/42_effector_hmmer/{effector}.tsv",
    params:
        base_path = base_path + "/42_effector_hmmer",
    shell:
        '''
        python3 scripts/effector_hmmer.py --input {input.multifasta} --output_basepath {params.base_path} --effector {wildcards.effector}
        '''

rule cATyper_hmm_search_hhsuite_pdb:
    '''
    HMM search against the PDB database using HHsuite.
    Using a custom made wrapper to first divide the multifasta into separate fastas for each protein as required by hhblits.
    The wrapper divides a multifasta (wildcard {c}) into its constituent proteins and then runs hhblits on each protein.
    The resulting files are .tsv files, each corresponding to a single protein within the effector wildcard (e.g. RelE protein X)
    These individual .tsv files are then concatenated into a single .tsv file for each effector and fed to the parse_hhsuite rule.
    '''
    input:
        multifasta = rules.effector_fetcher.output.multifasta
    output:
        hhsuite_concat = base_path + "/42_effector_hhsuite/{effector}/{effector}_all.tsv",
    params:
        outdir = base_path + "/42_effector_hhsuite/{effector}",
        pdb70 = "/media/volume/st_andrews/databases/pdb/pdb70",
        pdb30 = "/media/volume/st_andrews/databases/pdb30/pdb30",
    conda: "envs/hhsuite.yaml"
    threads: 40
    log:
        out = base_path + "/42_effector_hhsuite/logs/{effector}.out",
        err = base_path + "/42_effector_hhsuite/logs/{effector}.err",
    shell:
        '''
        python3 scripts/hhblits_wrapper.py --input {input.multifasta} --output_basepath {params.outdir} --database {params.pdb30}
        cat {params.outdir}/hhblits/*.tsv > {output.hhsuite_concat}
        '''

rule parse_hhsuite:
    '''
    Parser output from HHSuite (PDB)
    '''
    input: rules.cATyper_hmm_search_hhsuite_pdb.output.hhsuite_concat
    output: base_path + "/42_effector_hhsuite/{effector}/{effector}_hhsuite_parsed.tsv"
    conda: "envs/hhsuite.yaml"
    params:
        database = "PDB"
    shell:
        '''
        python3 scripts/hhsuite_parser.py --infile {input} --outfile {output}  --database {params.database}
        '''

rule concatenate_cATyper_hmm_hhsuite:
    input: aggregate_cATyper_hhsuite_parser
    output: base_path + "/42_effector_hhsuite/cATyper_all.tsv"
    shell:
        """
        awk '(NR == 1) || (FNR > 1)' {input} > {output}
        """

rule cATyper_hmm_search_hhsuite_cog:
    '''
    HMM search against the COG database using HHsuite.
    Using a custom made wrapper to first divide the multifasta into separate fastas for each protein as required by hhblits.
    The wrapper divides a multifasta (wildcard {c}) into its constituent proteins and then runs hhblits on each protein.
    The resulting files are .tsv files, each corresponding to a single protein within the effector wildcard (e.g. RelE protein X)
    These individual .tsv files are then concatenated into a single .tsv file for each effector and fed to the parse_hhsuite rule.
    '''
    input:
        multifasta = rules.effector_fetcher.output.multifasta
    output:
        hhsuite_concat = base_path + "/42_effector_hhsuite_cogs/{effector}/{effector}_all.tsv",
    params:
        outdir = base_path + "/42_effector_hhsuite_cogs/{effector}",
        cogs = "/media/volume/st_andrews/databases/cog/COG_KOG/COG_KOG",
    conda: "envs/hhsuite.yaml"
    threads: 40
    log:
        out = base_path + "/42_effector_hhsuite_cogs/logs/{effector}.out",
        err = base_path + "/42_effector_hhsuite_cogs/logs/{effector}.err",
    shell:
        '''
        python3 scripts/hhblits_wrapper.py --input {input.multifasta} --output_basepath {params.outdir} --database {params.cogs}
        cat {params.outdir}/hhblits/*.tsv > {output.hhsuite_concat}
        '''

rule parse_hhsuite_cogs:
    input: rules.cATyper_hmm_search_hhsuite_cog.output.hhsuite_concat
    output: base_path + "/42_effector_hhsuite_cogs/{effector}/{effector}_hhsuite_parsed_cogs.tsv"
    conda: "envs/hhsuite.yaml"
    params:
        database = "COGs",
        mapping = "/media/volume/st_andrews/databases/cog/COG_KOG/cog-20.def.tab"
    shell:
        '''
        python3 scripts/hhsuite_parser.py --infile {input} --outfile {output} --database {params.database} --mapping {params.mapping}
        '''

rule concatenate_cATyper_hmm_hhsuite_cogs:
    input: aggregate_cATyper_hhsuite_parser_cogs
    output: base_path + "/42_effector_hhsuite_cogs/cATyper_all_hmm_cogs.tsv"
    shell:
        """
        awk '(NR == 1) || (FNR > 1)' {input} > {output}
        """

rule effector_analyse_domains:
    '''
    Takes the HMM output from effector_hmmer and creates a simplified table of domains for each effector.
    This is only for the Pfam analysis.
    LEGACY: Links protein ID to each effector (the original effector fastas previously contained locus IDs, not protein IDs)
    '''
    input:
        raw_hmm = rules.effector_hmmer.output.raw_hmm,
        multifasta = rules.effector_fetcher.output.multifasta
    output:
        domains = base_path + "/43_effector_hmmer_analysis/{effector}/{effector}_sorted_filtered_mapped.tsv",
        filtered_hmm_tsv = base_path + "/43_effector_hmmer_analysis/{effector}/{effector}_sorted_filtered.tsv",
        sorted_hmm_tsv = base_path + "/43_effector_hmmer_analysis/{effector}/{effector}_sorted.tsv"
    params:
        base_path = base_path + "/43_effector_hmmer_analysis",
        effector_locus_map = base_path + "/43_effector_hmmer_analysis/{effector}/{effector}.locus_map.tsv",
    shell:
        '''
        cat {base_path}/11_cATyper_analysis/*/*_protein_to_effector.tsv | grep "{wildcards.effector}" > {params.effector_locus_map}
        python3 scripts/effector_hmmer_analyse.py --input {input.raw_hmm} --multifasta {input.multifasta} --output_basepath {params.base_path} --effector {wildcards.effector} --effector_locus_map {params.effector_locus_map}
        '''

rule effector_align:
    input: rules.effector_fetcher.output.multifasta
    output: base_path + "/44_effector_alignments/{effector}.afa",
    conda: "envs/trees.yaml"
    log:
        out = base_path + "/44_effector_alignments/logs/{effector}.out",
        err = base_path + "/44_effector_alignments/logs/{effector}.err"
    threads: 40
    shell:
        '''
        muscle -super5 "{input}" -output "{output}" -threads {threads} 2> {log.err} 1> {log.out} 
        '''


rule effector_tree:
    input: rules.effector_align.output
    output: base_path + "/45_effector_tree/{effector}_tree.txt",
    conda: "envs/trees.yaml"
    shell:
        '''
        FastTree -wag -gamma {input} > {output}
        '''


rule concatenate_effector_wildcards:
    input:
        trees = aggregate_known_effector_wildcards,
        annotations = aggregate_known_effector_annotations,
        checkpoint = aggregate_known_effector_wildcarder
    output: base_path + "/known_effectors_finished.txt"
    shell:
        '''
        cat {input} > {output}
        '''


checkpoint cora_neighbourhood_preparation:
    '''
    Once the basic analyses for known effectors are finished, we look in detail into CorA containing loci.
    For each CorA, we find out if it is associated with NrN and/or SAM-lyase
    These data along with the tree and CRISPR-Cas subtype/host information are passed onto R for visualisation.
    Get CorA loci from the mapping file and use path base_path + "/072_crispr_locus_proteins/{c}/{c}_crispr_locus_proteins.faa"

    This checkpoint output can be used for wildcard generation for each CorA associated locus (see next rule)
    '''
    input:   
        effectors_finished = rules.concatenate_effector_wildcards.output, #signifies all known effector batch analyses, incl. CorA, are finished
    output: directory(base_path + "/50_cora_proteomes")
    params:
        cora_multifasta = base_path + "/41_known_effector_mf/cora.faa", #multifasta for cora
        cora_tree = base_path + "/45_effector_tree/cora_tree.txt", #tree for cora
        protein_locus_map = base_path + "/43_effector_hmmer_analysis/cora/cora.locus_map.tsv" #tsv file with protein id linking it to its original locus
    run:
        import pandas as pd
        import shutil
        if not os.path.exists(str(output)):
            os.makedirs(str(output))
        locus_map = pd.read_csv(params.protein_locus_map, sep = "\t", header = None)
        print(locus_map)
        locus_map.columns = ["something", "protein_id", "effector", "locus_id"]
        print(locus_map)
        #find proteomes for CorA loci and copy the .faas into a new folder
        proteome_root_path = base_path + "/072_crispr_locus_proteins"
        for index, row in locus_map.iterrows():
            print(row)
            print(row["locus_id"])
            locus_path = proteome_root_path + "/" + row["locus_id"] + "/" + row["locus_id"] + "_crispr_locus_proteins.faa"
            print("Copying proteome from " + locus_path)
            shutil.copy(locus_path, str(output))

rule blast_cora_associated_proteins_against_dbs:
    '''
    This rule blasts proteomes from CorA-positive CRISPR-Cas III loci
    against NrN/DEDD and SAM-lyase databases. The results are parsed in another rule
    '''
    input: base_path + "/50_cora_proteomes/{cora_locus}_crispr_locus_proteins.faa"
    output:
        dedd = base_path + "/51_cora_neighbour_blast/{cora_locus}/dedd.tsv",
        dedd_neo = base_path + "/51_cora_neighbour_blast/{cora_locus}/dedd_neo.tsv",
        dedd_clost = base_path + "/51_cora_neighbour_blast/{cora_locus}/dedd_clost.tsv",
        nrn = base_path + "/51_cora_neighbour_blast/{cora_locus}/nrn.tsv",
        samlyase = base_path + "/51_cora_neighbour_blast/{cora_locus}/samlyase.tsv",
    params:
        dedd_db = "/media/volume/st_andrews/databases/cora_neighbours/dedd.dmnd",
        dedd_neo_db = "/media/volume/st_andrews/databases/cora_neighbours/dedd_neowinella.dmnd", #additional dedd family
        dedd_clost_db = "/media/volume/st_andrews/databases/cora_neighbours/dedd_clostridiisali.dmnd", #additional dedd family
        nrn_db = "/media/volume/st_andrews/databases/cora_neighbours/nrn.dmnd",
        samlyase_db = "/media/volume/st_andrews/databases/cora_neighbours/samlyase.dmnd",
    conda: "envs/diamond.yaml"
    shell:
        '''
        diamond blastp --db {params.dedd_db} --query {input} --out {output.dedd} --outfmt 6 --quiet
        diamond blastp --db {params.dedd_neo_db} --query {input} --out {output.dedd_neo} --outfmt 6 --quiet
        diamond blastp --db {params.dedd_clost_db} --query {input} --out {output.dedd_clost} --outfmt 6 --quiet
        diamond blastp --db {params.nrn_db} --query {input} --out {output.nrn} --outfmt 6 --quiet
        diamond blastp --db {params.samlyase_db} --query {input} --out {output.samlyase} --outfmt 6 --quiet
        '''

rule analyse_cora_associated_proteins_against_dbs:
    '''
    This rule analyses blast hits for NrN/DEDD and SAM-lyase hits,
    outputting a mostly boolean .tsv file with columns:
    locus | genome | NrN_DEDD | SAM_lyase | NrN_DEDD_sequence | SAM_lyase_sequence | CorA_sequence
    '''
    input:
        dedd = base_path + "/51_cora_neighbour_blast/{cora_locus}/dedd.tsv",
        dedd_neo = base_path + "/51_cora_neighbour_blast/{cora_locus}/dedd_neo.tsv",
        dedd_clost = base_path + "/51_cora_neighbour_blast/{cora_locus}/dedd_clost.tsv",
        nrn = base_path + "/51_cora_neighbour_blast/{cora_locus}/nrn.tsv",
        samlyase = base_path + "/51_cora_neighbour_blast/{cora_locus}/samlyase.tsv",
    output: base_path + "/52_cora_neighbour_analysis/{cora_locus}/neighbourhood_results.tsv"
    params:
        outdir = base_path + "/52_cora_neighbour_analysis/{cora_locus}"
    shell:
        '''
        python scripts/cora_neighbour_analysis.py --dedd {input.dedd} --dedd_neo {input.dedd_neo} --dedd_clost {input.dedd_clost} --nrn {input.nrn} --samlyase {input.samlyase} --outdir {params.outdir} --locus {wildcards.cora_locus}
        '''

rule concatenate_cora_neighbourhoods:
    input: aggregate_cora_neighbourhoods
    output: base_path + "/52_cora_neighbour_analysis/cora_neighbours.tsv"
    shell:
        '''
        awk '(NR == 1) || (FNR > 1)' {base_path}/52_cora_neighbour_analysis/*/neighbourhood_results.tsv > {output}
        '''

rule crispr_locus_proteins:
    '''
    Extracts all proteins from a CRISPR-Cas locus +- 4000 bp
    '''
    input:
        cas_operons_tsv = base_path + "/071_cctyper_loci/{c}/cas_operons.tsv",
    output:
        crispr_locus_proteins = base_path + "/072_crispr_locus_proteins/{c}/{c}_crispr_locus_proteins.faa",
    params:
        output_folder = base_path + "/072_crispr_locus_proteins/{c}",
        locus = "{c}",
        host_genomes_folder = base_path + "/06_host_genomes",
    run:
        import os
        import pandas as pd
        from Bio import SeqIO
        import gffutils
        import re
        import sys

        #when looking for effectors in the locus, the effector search range is the number of bases up or downstream of the cctyper defined cas operon boundaries
        effector_search_range = 4000

        #get sample name by splitting the locus name at the underscore and taking the first two parts
        sample = params.locus.split("_")[0] + "_" + params.locus.split("_")[1]

        #read in the cas operons file
        cas_operons_df = pd.read_csv(input.cas_operons_tsv, sep ="\t", header = 0)

        #find the gff file for the sample from the host_genomes_folder/samples folder
        gff_file = os.path.join(params.host_genomes_folder, sample, sample + "_features.gff")

        #Using gffutils, create a database of the gff file
        db_path = os.path.join(params.output_folder, params.locus)

        #get the path to the proteins fasta file for the current sample
        proteins_fasta = os.path.join(params.host_genomes_folder, sample, sample + "_proteins.faa")

        #create folder for the above path if it does not exist
        #if not os.path.exists(os.path.dirname(db_path)):
        #    print("Creating folder " + str(os.path.dirname(db_path)))
        #    os.makedirs(os.path.dirname(params.db_path))
        db = gffutils.create_db(gff_file, 
                                dbfn=db_path + ".db", 
                                force=True, 
                                keep_order=True, 
                                merge_strategy='merge', 
                                sort_attribute_values=True,
                                id_spec=["protein_id", "Name", "ID"])

        #read in the db
        db = gffutils.FeatureDB(db_path + ".db", keep_order=True)

        #get contig from row 0, column Contig
        contig = cas_operons_df.iloc[0]["Contig"]

        cas_operon_start = int(cas_operons_df['Start'][0]) - effector_search_range
        cas_operon_end = int(cas_operons_df['End'][0]) + effector_search_range

        #from the gff file, extract all the features that are on the contig. The featuretype must be "CDS"
        protein_ids = []
        proteins_on_contig = db.region(featuretype='CDS', seqid=contig)

        #then, extract all proteins whose coordinates are between the start and end of the cas operon +- the effector search range
        proteins_on_crispr_locus = db.region(featuretype='CDS', seqid=contig, start = cas_operon_start, end = cas_operon_end)

        #convert the returned generator to a list
        print("Extracting protein IDs from gff file")
        for i in proteins_on_crispr_locus:
            #check if the attributes of the feature has a key called 'protein_id'
            if "protein_id" in i.attributes:
                id_in_gff = str(i.attributes['protein_id'][0])
                #id = str(i.attributes['ID'][0]).split("-")[1]
                protein_ids.append(id_in_gff)

        #using biopython, extract the protein sequences using the list protein_ids from the proteins fasta file
        #the proteins fasta file is in the same folder as the gff file
        protein_seqs = []
        print("Reading protein fasta file from " + str(proteins_fasta))
        for record in SeqIO.parse(proteins_fasta, "fasta"):
            if record.id in protein_ids:
                protein_seqs.append(record)

        #write the protein sequences to a multifasta file
        print("Writing protein sequences to " + str(output.crispr_locus_proteins))
        SeqIO.write(protein_seqs, output.crispr_locus_proteins, "fasta")



rule aggregate_crispr_locus_proteins:
    '''
    Aggregates the outputs of the crispr_locus_proteins rule.
    '''
    input: aggregate_crispr_locus_proteins
    output: base_path + "/072_crispr_locus_proteins/crispr_locus_proteins_all.faa"
    shell:
        '''
        cat {input} > {output}
        '''



rule typeIII_characterizer:
    '''
    Characterizes previously known genes from type III loci (currently Cas10, Cas7, Cas5, CorA).
    This version incorporates the HD domain search. HD domains are searched for in the Cas10 protein
    and within the first 5 to 30 residues. Also adds Cas10 length.
    Now also looks for GGDD domain in Cas10.
    '''
    input:
        cctyper = base_path + "/071_cctyper_loci/{c}/cas_operons.tsv", #the cas_operons file is trimmed down to only include this specific {c} locus
        #gff = base_path + "/06_host_genomes/{j}/{j}_features.gff",
        #proteins = base_path + "/06_host_genomes/{j}/{j}_proteins.faa",
    output:
        Cas10_fasta = base_path + "/09_crispr_iii_CorA/loci/{c}/{c}_Cas10.faa",
        Cas5_fasta = base_path + "/09_crispr_iii_CorA/loci/{c}/{c}_Cas5.faa",
        Cas7_fasta = base_path + "/09_crispr_iii_CorA/loci/{c}/{c}_Cas7.faa",
        #CorA_fasta = base_path + "/09_crispr_iii_CorA/loci/{c}/{c}_CorA.faa", deprecated as of 17.1.2024
        info = base_path + "/09_crispr_iii_CorA/loci/{c}/{c}_crispr_iii_info.tsv",
    params:
        this_folder = base_path + "/09_crispr_iii_CorA/loci/{c}",
        outputfolder = base_path + "/09_crispr_iii_CorA/loci/{c}",
        sample_folder = base_path + "/06_host_genomes/", #this is used to find the gff and proteins files that are strain-, not locus-specific
        cctyper_folder = base_path + "/07_cctyper"
    conda: "envs/gff_utils.yaml"
    log:
        out = base_path + "/09_crispr_iii_CorA/logs/{c}.out",
        err = base_path + "/09_crispr_iii_CorA/logs/{c}.err"
    shell:
        '''
        python3 scripts/type_iii_effector_finder_2.0_HD.py --locus_id {wildcards.c} --sample_folder {params.sample_folder} --this_folder {params.this_folder} --outputfolder {params.outputfolder} --cas_operons {input.cctyper} --info_out {output.info} --cctyper_path {params.cctyper_folder} 2> {log.err} 1> {log.out}
        touch {output.Cas10_fasta}
        touch {output.Cas5_fasta}
        touch {output.Cas7_fasta}
        '''

rule unknown_finder:
    '''
    Probes the CCTyper outputs for unknown genes and/or genes with low E-values in the HMM search.
    In the future, the info file could contain the following columns:
        - locus_id
        - protein id
        - protein length
        - protein sequence
        - presence of CARF/SAVED domain
        - presence of HEPN domain
        - transmembrane properties
    Also outputs a fasta file with the unknown proteins. The output files are completely empty if no unknown proteins are found.
    
    For transmembrane prediction, we use TMHMM. This is not a conda package, so is installed
    using pip. The pip freeze command + grep is used to check if the package is already installed.
    NOTE: TMHMM will fail out of the box. A file called /home/ubuntu/miniconda3/envs/environmentname/lib/python3.9/site-packages/tmhmm/api.py in the conda installation will need this 
    to bypass a non-crucial warning that causes the crash:
    import warnings
        warnings.filterwarnings('ignore', category=RuntimeWarning)
    '''
    input:
        cctyper = base_path + "/071_cctyper_loci/{c}/cas_operons.tsv",
        hmm = rules.cATyper_hmm_search.output.temp_rows
    output:
        unknown_proteins = base_path + "/30_unknown_effectors/{c}/{c}_unknown_proteins.faa",
        info = base_path + "/30_unknown_effectors/{c}/{c}_unknown_proteins_info.tsv", #each row is a protein
        locus_info = base_path + "/30_unknown_effectors/{c}/{c}_locus_unknown_info.tsv", #each row is a locus
    params:
        outputfolder = base_path + "/30_unknown_effectors/{c}",
        host_genomes_folder = base_path + "/06_host_genomes",
        cctyper_folder = base_path + "/07_cctyper",
        sample_folder = base_path + "/06_host_genomes/", #this is used to find the gff and proteins files that are strain-, not locus-specific
        extra_cas_db = "/media/volume/st_andrews/databases/new_effectors_additional_cas/fastas/new_effectors_additional_cas.dmnd",
    conda: "envs/gff_utils.yaml"
    log:
        out = base_path + "/30_unknown_effectors/logs/{c}.out",
        err = base_path + "/30_unknown_effectors/logs/{c}.err"
    shell:
        '''
        pip freeze | grep tmhmm.py || pip install tmhmm.py
        python3 scripts/unknown_finder.py --locus_id {wildcards.c} --cctyper {input.cctyper} --hmm {input.hmm} --additional_cas_db {params.extra_cas_db} --outputfolder {params.outputfolder} --info_out {output.info} --unknown_proteins_output {output.unknown_proteins} --host_genomes_folder {params.host_genomes_folder} --cctyper_path {params.cctyper_folder} --sample_folder {params.sample_folder} --output_locus_info {output.locus_info} 2> {log.err} 1> {log.out}
        touch {output.unknown_proteins}
        touch {output.info}
        '''

rule concatenate_unknowns:
    input: aggregate_unknowns
    output:
        proteins = base_path + "/30_unknown_effectors/unknowns.faa",
        info = base_path + "/30_unknown_effectors/unknowns_info.tsv"
    shell:
        '''
        find '{base_path}/30_unknown_effectors' -maxdepth 2 -type f -wholename '*/*_unknown_proteins.faa' -print0 | xargs -0 cat >> {output.proteins}
        echo "locus_id\tsample\tsequence\tlength\tevalue\tposition\tid\tcctyper" > {output.info}
        find '{base_path}/30_unknown_effectors' -maxdepth 2 -type f -wholename '*/*_unknown_proteins_info.tsv' -print0 | xargs -0 cat >> {output.info}
        '''

rule concatenate_unknowns_locus_info:
    '''
    Aggregates the locus specific regarding unknown effectors info files
    '''
    input: aggregate_unknowns_locus_info
    output:
        info = base_path + "/30_unknown_effectors/unknowns_info_loci.tsv"
    shell:
        '''
        echo "locus_id\tsample\tno_of_unknowns\tunknown_proteins" > {output.info}
        find '{base_path}/30_unknown_effectors' -maxdepth 2 -type f -wholename '*/*_locus_unknown_info.tsv' -print0 | xargs -0 cat >> {output.info}
        '''

rule cluster_unknowns:
    '''
    Makes clusters of the unknown proteins using CD-HIT.
    Choose of word size:
        -n 5 for thresholds 0.7 ~ 1.0
        -n 4 for thresholds 0.6 ~ 0.7
        -n 3 for thresholds 0.5 ~ 0.6
        -n 2 for thresholds 0.4 ~ 0.5
    '''
    input: rules.concatenate_unknowns.output.proteins
    output:
        proteins = base_path + "/31_unknowns_cluster/unknowns_cluster.faa",
        clusterinfo = base_path + "/31_unknowns_cluster/unknowns_cluster.faa.clstr",
    conda: "envs/trees.yaml"
    log:
        out = base_path + "/31_unknowns_cluster/logs/unknowns_cluster.out",
        err = base_path + "/31_unknowns_cluster/logs/unknowns_cluster.err"
    params:
        cluster_cutoff = 0.4
    threads: 40
    shell:
        '''
        cd-hit -i {input} -o {output.proteins} -c {params.cluster_cutoff} -n 2 -d 0 -M 16000 -T {threads}
        '''


rule concatenate_type_iii_info:
    '''
    1. Creates header in output
    2. Concatenates individual loci info files
    3. Removes the headers from the concatenated file using inverse grep
    4. Joins this modified file with the header file
    '''
    input: aggregate_typeIII_info
    output: base_path + "/09_crispr_iii_CorA/loci/type_iii_info.tsv"
    params:
        temp_out = base_path + "/09_crispr_iii_CorA/loci/temp.tsv" 
    shell:
        '''
        echo "Cas10\tCas5\tCas7\tLocus\tSample\tCas10_GGDD\tCas10_GGDD_coord\tCas10_GGDD_seq\tCas10_GGDE\tCas10_GGDE_coord\tCas10_GGDE_seq\tCas10_GGED\tCas10_GGED_coord\tCas10_GGED_seq\tCas10_HD\tCas10_HD_list\tCas10_DH\tCas10_HD_coord\tCas10_DH_coord\tCas10_coord\tCas10_length\tSubtype\tcyclase_literal" > {params.temp_out}
        find '{base_path}/09_crispr_iii_CorA/loci' -maxdepth 2 -type f -wholename '*/*_crispr_iii_info.tsv' -print0 | xargs -0 tail -n +2 >> {params.temp_out}
        grep -v "==>" {params.temp_out} > {output}
        '''

rule Cas10_concatenate:
    '''
    Concatenates Cas10 sequences output by rule concatenate_type_iii_info
    '''
    input: aggregate_cas10_sequences
    output: base_path + "/10_Cas10_cluster/cas10s.faa"
    shell:
        '''
        cat {input} > {output} 
        '''


rule Cas10_HD_hmm_maker:
    '''
    1. Extracts the first 10-35 AA of Cas10s that are marked as having an HD domain.
    2. Aligns the Cas10s using Muscle
    3. Creates a new HD-HMM profile from this alignment

    Also has the option to use a premade alignment. Use --use_existing_alignment to do this (does not need argument)
    '''
    input:
        info_table = rules.concatenate_type_iii_info.output,
        cas10_sequences = rules.Cas10_concatenate.output,
    output:
        hmm = base_path + "/Cas10_HD_hmm_profiles/HD_HMM.hmm",
        msa = base_path + "/Cas10_HD_hmm_profiles/HD_HMM.msa",
        faa = base_path + "/Cas10_HD_hmm_profiles/HD_HMM.faa",
    conda: "envs/hmmer.yaml"
    log:
        out = base_path + "/Cas10_HD_hmm_profiles/logs/HD_HMM_maker.log",
        err = base_path + "/Cas10_HD_hmm_profiles/logs/HD_HMM_maker.err"
    shell:
        '''
        python scripts/Cas10_HD_hmm_maker.py --input_table {input.info_table} --cas10_sequences {input.cas10_sequences} --output {output.hmm} --msa {output.msa} --faa {output.faa}  --existing_alignment_path {modified_cas10_hd} 2> {log.err} 1> {log.out}
        touch {output.msa}
        touch {output.faa}
        '''

rule Cas10_HD_hmmer:
    '''
    Using the HD-HMM profiles made in rule Cas10_HD_hmm_maker, searches for HD domains in all Cas10s
    '''
    input:
        hmm_db = rules.Cas10_HD_hmm_maker.output.hmm,
        cas10_sequences = rules.Cas10_concatenate.output
    output:
        raw_table = base_path + "/Cas10_HD_hmm_profiles/HD_HMM_hits.out"
    params:
        out = base_path + "/Cas10_HD_hmm_profiles/HD_HMM_hits_temp.out",
        rows = base_path + "/Cas10_HD_hmm_profiles/HD_HMM_hits_rows.out",
        rows1 = base_path + "/Cas10_HD_hmm_profiles/HD_HMM_hits_rows1.out",
        E = "1e-01"
    conda: "envs/hmmer.yaml"
    shell:
        '''
        hmmscan --tblout {params.out} --cpu 8 -E {params.E} {input.hmm_db} {input.cas10_sequences} &> /dev/null
        grep -v "#" {params.out} > {params.rows}||:
        echo "target_name accession query_name accession E-value score bias E-value_best_domain score bias exp reg clu over env dom rep inc description_of_target" > {output.raw_table}
        cat {params.rows} >> {output.raw_table}
        '''

rule merge_HD_hmm_with_info:
    '''
    Merges the HD-HMM hits with the type III info table
    '''
    input:
        info_table = rules.concatenate_type_iii_info.output,
        HD_hmm_hits = rules.Cas10_HD_hmmer.output.raw_table
    output:
        merged_table = base_path + "/Cas10_HD_hmm_profiles/HD_HMM_hits_merged.tsv"
    run:
        import pandas as pd
        info = pd.read_csv(str(input.info_table), sep = "\t")
        print(info)
        hmm_hits = pd.read_csv(str(input.HD_hmm_hits), sep = '\s+')
        print(hmm_hits)

        #reduce the hmm_hits to only query_name and E-value_best_domain
        hmm_hits = hmm_hits[["query_name", "E-value_best_domain"]]

        #rename the columns
        hmm_hits.columns = ["Locus", "HD_E-value"]
        #add column "HD_hmm_boolean" and set it to True
        hmm_hits["HD_hmm_boolean"] = True

        #merge the info table with the hmm_hits table
        merged = pd.merge(info, hmm_hits, on = "Locus", how = "left")
        #if some rows were left without a hit, set the HD_hmm_boolean to False
        merged["HD_hmm_boolean"] = merged["HD_hmm_boolean"].fillna(False)

        #from the final file, remove all columns but Locus, HD_E-value and HD_hmm_boolean
        merged = merged[["Locus", "HD_E-value", "HD_hmm_boolean"]]

        #save the merged table
        merged.to_csv(str(output.merged_table), sep = "\t", index = False)


rule R_HD:
    '''
    R script that produces a histogram of Cas10s and their HD domains
    '''
    input: rules.concatenate_type_iii_info.output
    output:
        HD_histogram = base_path + "/4_R_HD/cas10_HD_lengths.png"
    log:
        out = base_path + "/4_R_HD/logs/HD_hist.out",
        err = base_path + "/4_R_HD/logs/HD_hist.err"
    params:
        outputfolder = base_path + "/4_R_HD"
    conda:
        "envs/R.yaml"
    shell:
        '''
        Rscript R/HD_R.R --input {input} --output {params.outputfolder} 2> {log.out} 1> {log.err}
        '''


rule Cas10_GGDD_hmm_maker:
    '''
    1. Extracts the X AA of sequence around GGDDs in Cas10s.
    2. Aligns the Cas10s using Muscle
    3. Creates a new HD-HMM profile from this alignment

    As of 23.8.2023, running the script two times: once for GGDD and once for GGDE.
    '''
    input:
        info_table = rules.concatenate_type_iii_info.output,
        cas10_sequences = rules.Cas10_concatenate.output,
    output:
        hmm = base_path + "/Cas10_GGDD_hmm_profiles/GGDD_HMM.hmm",
        msa = base_path + "/Cas10_GGDD_hmm_profiles/GGDD_HMM.msa",
        faa = base_path + "/Cas10_GGDD_hmm_profiles/GGDD_HMM.faa",
        hmm_GGDE = base_path + "/Cas10_GGDD_hmm_profiles/GGDE_HMM.hmm",
        msa_GGDE = base_path + "/Cas10_GGDD_hmm_profiles/GGDE_HMM.msa",
        faa_GGDE = base_path + "/Cas10_GGDD_hmm_profiles/GGDE_HMM.faa",
    conda: "envs/hmmer.yaml"
    shell:
        '''
        python scripts/Cas10_GGDD_hmm_maker.py --input_table {input.info_table} --cas10_sequences {input.cas10_sequences} --output {output.hmm} --msa {output.msa} --faa {output.faa} --motif "GGDD"
        python scripts/Cas10_GGDD_hmm_maker.py --input_table {input.info_table} --cas10_sequences {input.cas10_sequences} --output {output.hmm_GGDE} --msa {output.msa_GGDE} --faa {output.faa_GGDE} --motif "GGDE"
        touch {output.msa}
        touch {output.faa}
        '''

rule Cas10_GGDD_hmmer:
    '''
    Using the GGDD-HMM profiles made in rule Cas10_GGDD_hmm_maker, searches for GGDD domains in all Cas10s
    '''
    input:
        hmm_db_GGDD = rules.Cas10_GGDD_hmm_maker.output.hmm,
        hmm_db_GGDE = rules.Cas10_GGDD_hmm_maker.output.hmm_GGDE,
        cas10_sequences = rules.Cas10_concatenate.output
    output:
        raw_table_GGDD = base_path + "/Cas10_GGDD_hmm_profiles/GGDD_HMM_hits.out",
        raw_table_GGDE = base_path + "/Cas10_GGDD_hmm_profiles/GGDE_HMM_hits.out"
    params:
        out_GGDD = base_path + "/Cas10_GGDD_hmm_profiles/GGDD_HMM_hits_temp.out",
        out_GGDE = base_path + "/Cas10_GGDD_hmm_profiles/GGDE_HMM_hits_temp.out",
        rows_GGDD = base_path + "/Cas10_GGDD_hmm_profiles/GGDD_HMM_hits_rows.out",
        rows_GGDE = base_path + "/Cas10_GGDD_hmm_profiles/GGDE_HMM_hits_rows.out",
        rows1 = base_path + "/Cas10_GGDD_hmm_profiles/GGDD_HMM_hits_rows1.out",
        E = "1e-02"
    conda: "envs/hmmer.yaml"
    shell:
        '''
        hmmscan --tblout {params.out_GGDD} --cpu 8 -E {params.E} {input.hmm_db_GGDD} {input.cas10_sequences}
        grep -v "#" {params.out_GGDD} > {params.rows_GGDD}||:
        echo "target_name accession query_name accession E-value score bias E-value_best_domain score bias exp reg clu over env dom rep inc description_of_target" > {output.raw_table_GGDD}
        cat {params.rows_GGDD} >> {output.raw_table_GGDD}

        hmmscan --tblout {params.out_GGDE} --cpu 8 -E {params.E} {input.hmm_db_GGDE} {input.cas10_sequences}
        grep -v "#" {params.out_GGDE} > {params.rows_GGDE}||:
        echo "target_name accession query_name accession E-value score bias E-value_best_domain score bias exp reg clu over env dom rep inc description_of_target" > {output.raw_table_GGDE}
        cat {params.rows_GGDE} >> {output.raw_table_GGDE}
        '''

rule merge_GGDD_hmm_with_info:
    '''
    Merges the GGDD-HMM hits with the type III info table
    '''
    input:
        info_table = rules.concatenate_type_iii_info.output,
        GGDD_hmm_hits = rules.Cas10_GGDD_hmmer.output.raw_table_GGDD,
        GGDE_hmm_hits = rules.Cas10_GGDD_hmmer.output.raw_table_GGDE,
        cas10s_faa = rules.Cas10_concatenate.output
    output:
        merged_table = base_path + "/Cas10_GGDD_hmm_profiles/GGDD_HMM_hits_merged.tsv"
    params:
        output_folder = base_path + "/Cas10_GGDD_hmm_profiles"
    run:
        import pandas as pd
        from Bio import SeqIO

        info = pd.read_csv(str(input.info_table), sep = "\t")
        hmm_hits_GGDD = pd.read_csv(str(input.GGDD_hmm_hits), sep = '\s+')
        hmm_hits_GGDE = pd.read_csv(str(input.GGDE_hmm_hits), sep = '\s+')

        #reduce the hmm_hits to only query_name and E-value_best_domain
        hmm_hits_GGDD = hmm_hits_GGDD[["query_name", "E-value_best_domain"]]
        hmm_hits_GGDE = hmm_hits_GGDE[["query_name", "E-value_best_domain"]]

        #rename the columns
        hmm_hits_GGDD.columns = ["Locus", "GGDD_E-value"]
        hmm_hits_GGDE.columns = ["Locus", "GGDE_E-value"]
        #add column "GGDD_hmm_boolean" and set it to True
        hmm_hits_GGDD["GGDD_hmm_boolean"] = True
        hmm_hits_GGDE["GGDE_hmm_boolean"] = True

        #merge the info table with the hmm_hits_GGDD table
        merged = pd.merge(info, hmm_hits_GGDD, on = "Locus", how = "left")
        merged = pd.merge(merged, hmm_hits_GGDE, on = "Locus", how = "left")

        #if some rows were left without a hit, set the GGDD_hmm_boolean to False
        merged["GGDD_hmm_boolean"] = merged["GGDD_hmm_boolean"].fillna(False)
        merged["GGDE_hmm_boolean"] = merged["GGDE_hmm_boolean"].fillna(False)


        #read in the Cas10 sequences using biopython
        cas10s = SeqIO.to_dict(SeqIO.parse(str(input.cas10s_faa), "fasta"))

        all_cyclase_literals = ["SGDD", "AGDD", "GGED", "GGDD", "GGED", "GGDE", "EGDD", "GEDD", "DGDD", "AGDE", "KGDD", "AGDE"] #updated list of possible cyclase motifs

        #create a dictionary of the literals and add count as value
        all_cyclase_literals_dic = {literal: 0 for literal in all_cyclase_literals}

        #create a column cyclase_literal_sequence in the merged table
        merged["cyclase_literal_sequence"] = ""

        #the merged contains a column "cyclase_literal". This column is true if the sequence contains one of the motifs in all_cyclase_literals. Go through each cas10 sequence
        for index, row in merged.iterrows():
            locus = row["Locus"]
            try:
                cas10 = cas10s[locus].seq
                #check if the sequence contains any of the motifs
                for motif in all_cyclase_literals:
                    if motif in cas10:
                        merged.at[index, "cyclase_literal"] = True
                        print("Cyclase found in " + locus + " with motif " + motif)
                        merged.at[index, "cyclase_literal_sequence"] = motif
                        break
                else:
                    print("No cyclase found for " + locus)
                    merged.at[index, "cyclase_literal"] = False
            except:
                print("No Cas10 found in fasta for " + locus)

        #create new column cyclase which is true only if column GGDD_hmm_boolean or column GGDE_hmm_boolean is true and cyclase_literal is true
        merged["cyclase"] = merged["cyclase_literal"] & (merged["GGDD_hmm_boolean"] | merged["GGDE_hmm_boolean"])

        #in the all_cyclase_literals_dic, add the count of each motif
        for index, row in merged.iterrows():
            motif = row["cyclase_literal_sequence"]
            if motif in all_cyclase_literals_dic:
                all_cyclase_literals_dic[motif] += 1

        #create a new dataframe from the dictionary
        cyclase_literal_df = pd.DataFrame(list(all_cyclase_literals_dic.items()), columns = ["motif", "count"])

        #count percentages of each motif
        cyclase_literal_df["percentage"] = round(cyclase_literal_df["count"] / len(merged) * 100, 1)

        #save the dataframe to a file
        cyclase_literal_df.to_csv(params.output_folder + "/cyclase_literal_count.tsv", sep = "\t", index = False)

        #from the final file, remove all columns but Locus, GGDD_E-value and GGDD_hmm_boolean and cyclase_literal_sequence
        merged = merged[["Locus", "GGDD_E-value", "GGDD_hmm_boolean", "GGDE_E-value", "GGDE_hmm_boolean", "cyclase", "cyclase_literal", "cyclase_literal_sequence"]]

        #save the merged table
        merged.to_csv(str(output.merged_table), sep = "\t", index = False)


rule combine_GGDD_HMM_to_mastertable:
    '''
    Merges the HD and GGDD-HMM data with the master table. Also merges domains info.
    Note that this is not the final mastertable (see rule mastercombiner)
    '''
    input:
        info = rules.concatenate_type_iii_info.output,
        HD_hmm_hits = rules.merge_HD_hmm_with_info.output.merged_table,
        GGDD_hmm_hits = rules.merge_GGDD_hmm_with_info.output.merged_table,
        domains = rules.annotate_bacteria_and_archaea_domains.output,
    output:
        final_info_table = base_path + "/mastertable.tsv"
    run:
        import pandas as pd
        info = pd.read_csv(str(input.info), sep = "\t")
        domains = pd.read_csv(str(input.domains), sep = "\t")
        hmm_hits_HD = pd.read_csv(str(input.HD_hmm_hits), sep = "\t")
        hmm_hits_GGDD = pd.read_csv(str(input.GGDD_hmm_hits), sep = "\t")
        print(hmm_hits_GGDD)
        merged = pd.merge(info, hmm_hits_HD, on = "Locus", how = "left")
        merged = pd.merge(merged, hmm_hits_GGDD, on = "Locus", how = "left")
        merged = pd.merge(merged, domains, left_on = "Sample", right_on = "sample", how = "left")
        merged.to_csv(str(output.final_info_table), sep = "\t", index = False)


#a rule that uses the output file cora_type_iii_info from the rule above filter samples with CorA and further divide them into CRISPR-Cas subtypes
rule cora_plot_extractor:
    '''
    Gets cctyper generated plots for each sample with a CorA and a type III CRISPR-Cas system.
    '''
    input: 
        info = rules.concatenate_type_iii_info.output[0]
    output:
        done = base_path + "/xtra1_cora_iii_loci_plots/done.done",
    run:
        import pandas as pd
        import shutil
        #using pandas, get the info from the cora_type_iii_info file and filter out samples that do not have a CorA
        cora_type_iii_info = pd.read_csv(input.info, sep = "\t")
        cora_type_iii_info = cora_type_iii_info[cora_type_iii_info["CorA"] == True]

        #create output folder if it does not exist
        if not os.path.exists(base_path + "/xtra1_cora_iii_loci_plots"):
            os.makedirs(base_path + "/xtra1_cora_iii_loci_plots")

        #For each sample in the pandas dataframe, extract the corresponding plot from the cctyper folder (path is base_path + "/07_cctyper/{j}/plot.png) and copy it to the output folder in a CRISPR-Cas subtype subfolder (e.g. III-A or III-B) depending on the Subtype column in the cora_type_iii_info file
        for index, row in cora_type_iii_info.iterrows():
            sample = row["Sample"]
            subtype = row["Subtype"]
            print(sample + "," + subtype)
            #check if subtype does not contain the substring "Hybrid" (this is because some samples have a hybrid subtype, e.g. III-A/B, and we don't want to create a III-A/B folder)
            if "Hybrid" not in subtype:
                #if the subtype folder does not exist, create it
                if not os.path.exists(base_path + "/xtra1_cora_iii_loci_plots/" + subtype):
                    os.makedirs(base_path + "/xtra1_cora_iii_loci_plots/" + subtype)
                shutil.copyfile(base_path + "/07_cctyper/" + sample + "/plot.png", base_path + "/xtra1_cora_iii_loci_plots/" + subtype + "/" + sample + "_plot.png")

        #create the done.done file to indicate that the rule has finished running
        open(output.done, "w").close()


rule CorA_concatenate:
    '''
    This version works on the cATyper outputs instead of CCTyper outputs
    '''
    input:
        coras = aggregate_CorA_sequences_catyper,
        type_iii_wildcarder_finished = rules.concatenate_type_iii_info.output
    output: base_path + "/09_crispr_iii_CorA/CorAs.faa"
    shell:
        '''
        cat {input.coras} > {output}
        '''

rule CorA_cluster:
    '''

    '''
    input: rules.CorA_concatenate.output
    output:
        proteins = base_path + "/15_CorA_cluster/CorA_cluster.faa",
        clusterinfo = base_path + "/15_CorA_cluster/CorA_cluster.faa.clstr",
    conda: "envs/trees.yaml"
    log:
        out = base_path + "/15_CorA_cluster/logs/CorA_align.out",
        err = base_path + "/15_CorA_cluster/logs/CorA_align.err"
    threads: 40
    shell:
        '''
        cd-hit -i {input} -o {output.proteins} -c 0.90 -n 5 -d 0 -M 16000 -T {threads}
        '''


rule CorA_align:
    input: rules.CorA_cluster.output.proteins
    output: base_path + "/16_CorA_align/CorA_alignment.afa",
    conda: "envs/trees.yaml"
    log:
        out = base_path + "/16_CorA_align/logs/CorA_align.out",
        err = base_path + "/16_CorA_align/logs/CorA_align.err"
    threads: 40
    shell:
        '''
        muscle -super5 "{input}" -output "{output}" -threads {threads} 2> {log.err} 1> {log.out}
        '''

rule CorA_align_unclustered:
    input: rules.CorA_concatenate.output
    output: base_path + "/16_CorA_align_unclustered/CorA_alignment_unclustered.afa",
    conda: "envs/trees.yaml"
    log:
        out = base_path + "/16_CorA_align/logs/CorA_align.out",
        err = base_path + "/16_CorA_align/logs/CorA_align.err"
    threads: 40
    shell:
        '''
        muscle -super5 "{input}" -output "{output}" -threads {threads} 2> {log.err} 1> {log.out}
        '''

rule CorA_tree:
    '''

    '''
    input: rules.CorA_align.output
    output: base_path + "/17_CorA_tree/CorA_tree.txt",
    conda: "envs/trees.yaml"
    shell:
        '''
        FastTree -wag -gamma {input} > {output}
        '''

rule CorA_tree_unclustered:
    '''

    '''
    input: rules.CorA_align_unclustered.output
    output: base_path + "/17_CorA_tree_unclustered/CorA_tree_unclustered.txt",
    conda: "envs/trees.yaml"
    shell:
        '''
        FastTree -wag -gamma {input} > {output}
        '''


rule Cas10_cluster:
    '''
    '''
    input: 
        cas10 = rules.Cas10_concatenate.output
    output:
        proteins = base_path + "/10_Cas10_cluster/cas10_cluster.faa",
        clusterinfo = base_path + "/10_Cas10_cluster/cas10_cluster.faa.clstr",
    conda: "envs/trees.yaml"
    log:
        out = base_path + "/10_Cas10_cluster/logs/cas10_align.out",
        err = base_path + "/10_Cas10_cluster/logs/cas10_align.err"
    threads: 40
    shell:
        '''
        echo {cas10_tree_cluster}
        if [ {cas10_tree_cluster} = "True" ]; then
            cd-hit -i {input.cas10} -o {output.proteins} -c 0.99 -n 5 -d 0 -M 16000 -T {threads}
        else
            cp {input.cas10} {output.proteins}
            touch {output.clusterinfo}
        fi
        '''


rule Cas10_align:
    '''
    
    '''
    input: rules.Cas10_cluster.output.proteins
    output: base_path + "/11_Cas10_align/cas10_alignment.afa",
    conda: "envs/trees.yaml"
    log:
        out = base_path + "/11_Cas10_align/logs/cas10_align.out",
        err = base_path + "/11_Cas10_align/logs/cas10_align.err"
    threads: 40
    shell:
        '''
        muscle -super5 "{input}" -output "{output}" -threads {threads} 2> {log.err} 1> {log.out}
        '''

rule Cas10_tree:
    '''

    '''
    input: rules.Cas10_align.output
    output: base_path + "/12_Cas10_tree/cas10_tree.txt",
    conda: "envs/trees.yaml"
    shell:
        '''
        FastTree -wag -gamma {input} > {output}
        '''

rule validated_new_effectors:
    '''
    This rule is updated as potential new effectors from group4 emerge.
    Creates tables for each locus showing whether a validated effector is
    found in that locus. Does not interfere with running the group4 analysis
    and the validated effectors found here will also be listed in group4 proteins,
    as that is where they were originally found.
    '''
    input:
        proteins = rules.cATyper_hmm_search.output.contig_proteins
    output:
        #validated_new_effectors = base_path + "/60_validated_new_effectors/{c}/{c}_validated_new_effectors_orig.tsv",
        temp_rows = base_path + "/60_validated_new_effectors/{c}/{c}_validated_new_effectors_temp.tsv",
        hmm_rows = base_path + "/60_validated_new_effectors/{c}/{c}_validated_new_effectors_hmm.tsv"
    conda: "envs/hmmer.yaml"
    params:
        validated_effectors_hmm_db = validated_effectors_hmm_db_path,
        temp_hmm = base_path + "/60_validated_new_effectors/{c}/{c}_validated_new_effectors.temp",
        evalue = "1e-10"
    log:
        out = base_path + "/60_validated_new_effectors/logs/{c}/{c}_validated_new_effectors.out",
        err = base_path + "/60_validated_new_effectors/logs/{c}/{c}_validated_new_effectors.err"
    shell:
        '''
        hmmscan --domtblout {params.temp_hmm} --cpu 8 -E {params.evalue} {params.validated_effectors_hmm_db} {input.proteins} &> /dev/null
        echo "Removing commented rows" >> {log.out}
        grep -v "^#" {params.temp_hmm} > {output.temp_rows} ||:
        echo "Writing header" >> {log.out}
        echo -e 'target_name\taccession\ttlen\tquery_name\taccession\tqlen\tE-value_fullseq\tscore_fullseq\tbias_fullseq\t#_domain\tof_domain\tc-Evalue_domain\ti-Evalue_domain\tscore_domain\tbias_domain\t_hmm_from\thmm_to\t_ali_from\tali_to\tenv_from\tenv_to\tacc\tdescription' > {output.hmm_rows}
        echo "Checking if hits were found" >> {log.out}
        if [ -s {output.temp_rows} ]; then 
            echo "Hits found for {wildcards.c}" >> {log.out}
            cat {output.temp_rows} >> {output.hmm_rows}
        else
            echo "No hits found for {wildcards.c}" >> {log.out}
            touch {output.hmm_rows}
        fi
        '''

rule concatenate_validate_new_effectors_hmm:
    input: aggregate_validated_new_effectors_hmm
    output: base_path + "/60_validated_new_effectors/validated_new_effectors_all_hmm.tsv"
    shell:
        """
        awk '(NR == 1) || (FNR > 1)' {input} > {output}
        """

rule validated_new_effectors_analysis:
    input:
        crispr_positive_samples = base_path + "/071_cctyper_loci/{c}/cas_operons.tsv",
        contig_proteins = base_path + "/10_cATyper_hmm/{c}/{c}_contig_proteins.faa",
        hmm_rows = base_path + "/60_validated_new_effectors/{c}/{c}_validated_new_effectors_hmm.tsv",
    output:
        catyper = base_path + "/61_validated_new_effectors_analysis/{c}/{c}_cATyper_results.tsv",
        hmm_targets = base_path + "/61_validated_new_effectors_analysis/{c}/{c}_cATyper_hmm_targets.tsv",
        effector_to_protein = base_path + "/61_validated_new_effectors_analysis/{c}/{c}_effector_to_protein.tsv",
        protein_to_effector = base_path + "/61_validated_new_effectors_analysis/{c}/{c}_protein_to_effector.tsv",
        plottable_effector_positions = base_path + "/61_validated_new_effectors_analysis/{c}/{c}_plottable_effector_positions.tsv",
    params:
        outdir = base_path + "/61_validated_new_effectors_analysis/{c}",
        host_genomes_folder = base_path + "/06_host_genomes",
        hmm_msa_folder = "/media/volume/st_andrews/databases/validated_new_effectors/profiles",
        tmhmm_model_path = "/media/volume/st_andrews/new_effectors/TM/TMHMM2.0.model",
    conda: "envs/groupChar.yaml"
    log:
        out = base_path + "/61_validated_new_effectors_analysis/logs/{c}.out",
        err = base_path + "/61_validated_new_effectors_analysis/logs/{c}.err",
    shell:
        '''
        echo "Running validated effector analysis" >> {log.out}
        echo "Listing targets" >> {log.out}
        ls {params.hmm_msa_folder}/*/*.hmm > {output.hmm_targets}
        python scripts/catyper_prepper_10.py --locus {wildcards.c} --cas_operons_file {input.crispr_positive_samples} --output_folder {params.outdir} --host_genomes_folder {params.host_genomes_folder} --mode post_hmm --hmm_targets {output.hmm_targets} --hmm_rows {input.hmm_rows} --catyper_out {output.catyper} --catyper_type "validated_new_effector" --effector_plot_data {output.plottable_effector_positions} --tmhmm_model_path {params.tmhmm_model_path} 2> {log.err} 1> {log.out}
        '''

rule concatenate_validated_new_effectors_analysis:
    input: aggregate_validated_new_effectors_analysis
    output: base_path + "/61_validated_new_effectors_analysis/validated_new_effectors_all_hmm.tsv"
    shell:
        """
        awk '(NR == 1) || (FNR > 1)' {input} > {output}
        """


rule concatenate_validated_new_effectors_scores:
    '''
    Concatenates info on the effector hits from cATyper analysis
    '''
    input:
        protein_to_effector = aggregate_validated_new_effectors_pte,
        effector_to_protein = aggregate_validated_new_effectors_etp
    output:
        protein_to_effector_concatenated = base_path + "/61_validated_new_effectors_analysis/val_new_effectors_protein_to_effector.tsv",
        effector_to_protein_concatenated = base_path + "/61_validated_new_effectors_analysis/val_new_effectors_to_protein.tsv",
    shell:
        '''
        cat {input.protein_to_effector} > {output.protein_to_effector_concatenated}
        cat {input.effector_to_protein} > {output.effector_to_protein_concatenated}
        '''

rule analyse_validated_new_effectors_scores:
    '''
    The algorithm that analyses how well each hmm profile performs in validated_new_effectors.
    '''
    input:
        pte = rules.concatenate_validated_new_effectors_scores.output.protein_to_effector_concatenated,
        etp = rules.concatenate_validated_new_effectors_scores.output.effector_to_protein_concatenated,
    output:
        effector_scores = base_path + "/61_validated_new_effectors_analysis/val_new_effectors_scores.tsv",
        effector_scores_summary = base_path + "/61_validated_new_effectors_analysis/val_new_effectors_scores_summary.tsv",
        #effector_scores_plot1 = base_path + "/11_cATyper_analysis/cATyper_effector_scores_plot1.png",
    params:
        outdir = base_path + "/61_validated_new_effectors_analysis"
    shell:
        '''
        python scripts/catyper_effector_scores.py --pte {input.pte} --etp {input.etp} --output {output.effector_scores} --output_summary {output.effector_scores_summary} --outdir {params.outdir}
        touch {output.effector_scores}
        touch {output.effector_scores_summary}
        '''

rule heatmap_known_validated_effectors:
    input:
        etp_validated = rules.concatenate_validated_new_effectors_scores.output.effector_to_protein_concatenated,
        etp_known = rules.concatenate_cATyper_analysis_effector_scores.output.effector_to_protein_concatenated
    output:
        effectors_combined = base_path + "/70_validated_and_known_effectors_heatmap/etp_combined.tsv",
        effector_scores_summary = base_path + "/70_validated_and_known_effectors_heatmap/effectors.tsv",
    params:
        outdir = base_path + "/70_validated_and_known_effectors_heatmap",
        dummyout = base_path + "/70_validated_and_known_effectors_heatmap/dummyout.tsv"
    shell:
        '''
        cat {input.etp_validated} {input.etp_known} > {output.effectors_combined}
        python scripts/catyper_effector_scores.py --etp {output.effectors_combined} --output {params.dummyout} --output_summary {output.effector_scores_summary} --outdir {params.outdir}
        '''




rule validated_effectors_cas10_fusions:
    '''
    Runs Cas10 sequences as query against the validated effectors database to detect any fusions
    '''
    input:
        cas10 = rules.Cas10_cluster.output.proteins
    output:
        temp_rows = base_path + "/101_validated_effectors_cas10/validated_effectors_cas10_temp.tsv",
        hmm_rows = base_path + "/101_validated_effectors_cas10/validated_effectors_cas10_hmm.tsv"
    conda: "envs/hmmer.yaml"
    params:
        validated_effectors_hmm_db = "/media/volume/st_andrews/databases/validated_new_effectors/concatenated_profiles/all.hmm",
        temp_hmm = base_path + "/101_validated_effectors_cas10/validated_effectors_cas10_temp.temp",
        evalue = "1e-10"
    log:
        out = base_path + "/101_validated_effectors_cas10/logs/validated_effectors_cas10.out",
        err = base_path + "/101_validated_effectors_cas10/logs/validated_effectors_cas10.err"
    shell:
        '''        
        hmmscan --domtblout {params.temp_hmm} --cpu 40 -E {params.evalue} {params.validated_effectors_hmm_db} {input.cas10}
        echo "Removing commented rows" >> {log.out}
        grep -v "^#" {params.temp_hmm} > {output.temp_rows} ||:
        echo -e 'target_name\taccession\ttlen\tquery_name\taccession\tqlen\tE-value_fullseq\tscore_fullseq\tbias_fullseq\t#_domain\tof_domain\tc-Evalue_domain\ti-Evalue_domain\tscore_domain\tbias_domain\t_hmm_from\thmm_to\t_ali_from\tali_to\tenv_from\tenv_to\tacc\tdescription' > {output.hmm_rows}
        cat {output.temp_rows} >> {output.hmm_rows}
        '''

rule known_effectors_cas10_fusions:
    '''
    Runs Cas10 sequences as query against the known effectors database to detect any fusions
    '''
    input:
        cas10 = rules.Cas10_cluster.output.proteins
    output:
        temp_rows = base_path + "/101_known_effectors_cas10/known_effectors_cas10_temp.tsv",
        hmm_rows = base_path + "/101_known_effectors_cas10/known_effectors_cas10_hmm.tsv"
    conda: "envs/hmmer.yaml"
    params:
        known_effectors_db = hmm_database_folder + "/" + hmm_database_file,
        temp_hmm = base_path + "/101_known_effectors_cas10/known_effectors_cas10_temp.temp",
        evalue = "1e-10"
    log:
        out = base_path + "/101_known_effectors_cas10/logs/known_effectors_cas10.out",
        err = base_path + "/101_known_effectors_cas10/logs/known_effectors_cas10.err"
    shell:
        '''        
        hmmscan --domtblout {params.temp_hmm} --cpu 40 -E {params.evalue} {params.known_effectors_db} {input.cas10}
        echo "Removing commented rows" >> {log.out}
        grep -v "^#" {params.temp_hmm} > {output.temp_rows} ||:
        echo -e 'target_name\taccession\ttlen\tquery_name\taccession\tqlen\tE-value_fullseq\tscore_fullseq\tbias_fullseq\t#_domain\tof_domain\tc-Evalue_domain\ti-Evalue_domain\tscore_domain\tbias_domain\t_hmm_from\thmm_to\t_ali_from\tali_to\tenv_from\tenv_to\tacc\tdescription' > {output.hmm_rows}
        cat {output.temp_rows} >> {output.hmm_rows}
        '''

rule mastercombiner:
    '''
    Added 4.9.2023 to incorporate merging and data handling previously done in R in here
    '''
    input:
        concat_taxInfo = rules.concat_taxInfo.output,
        catyper = rules.concatenate_cATyper_analysis.output,
        type_iii_info = rules.combine_GGDD_HMM_to_mastertable.output.final_info_table,
        clustered_unknowns_info = rules.concatenate_unknowns_locus_info.output.info,
        temperature_data = temperature_data,
        validated_new_effectors = rules.concatenate_validated_new_effectors_analysis.output,
    output:
        final_info_table = base_path + "/mastertable_v2.tsv",
        group4_IDs = base_path + "/group4/group4_IDs.txt"
    run:
        import pandas as pd
        mastertable = pd.read_csv(str(input.type_iii_info), sep = "\t")
        taxInfo = pd.read_csv(str(input.concat_taxInfo), sep = "\t")
        catyper = pd.read_csv(str(input.catyper), sep = "\t")
        clustered_unknowns_info = pd.read_csv(str(input.clustered_unknowns_info), sep = "\t")
        validated_new_effectors_df = pd.read_csv(str(input.validated_new_effectors), sep = "\t")


        #### KNOWN EFFECTORS PREP ####
        #in catyper, if any of the following columns are True, set has_known_effector column to True: ca3	ca4	ca5	ca6	sam-amp
        catyper["has_known_effector"] = False
        catyper.loc[catyper[["ca3", "ca4", "ca5", "ca6", "sam-amp"]].any(axis = 1), "has_known_effector"] = True

        #### VALIDATED NEW EFFECTOR PREP ####
#       validated_new_effectors_df = validated_new_effectors_df[["can3", "tirsaved", "cam2", "mem-01", "mem2", "locus"]] #drop cOA columns
        validated_new_effectors_df["has_validated_new_effector"] = False
        #if any column except locus is True, set has_validated_new_effector to True
        validated_new_effectors_df.loc[validated_new_effectors_df[["can3", "tirsaved", "cam3", "cam2"]].any(axis = 1), "has_validated_new_effector"] = True

        #rename column no_effectors to no_validated_new_effectors
        validated_new_effectors_df = validated_new_effectors_df.rename(columns = {"no_effectors": "no_validated_new_effectors"})

        #rename columns ca3	ca4	ca5	ca6 to NE_ca3 NE_ca4 NE_ca5 NE_ca6
        validated_new_effectors_df = validated_new_effectors_df.rename(columns = {"ca3": "NE_ca3", "ca4": "NE_ca4", "ca5": "NE_ca5", "ca6": "NE_ca6", "sam-amp": "NE_sam-amp"})

        #### MERGE ####
        mastertable = pd.merge(mastertable, taxInfo, on = "Sample", how = "left")
        mastertable = pd.merge(mastertable, catyper, left_on = "Locus", right_on = "locus", how = "left")
        mastertable = pd.merge(mastertable, clustered_unknowns_info, left_on = "Locus", right_on = "locus_id", how = "left")
        mastertable = pd.merge(mastertable, validated_new_effectors_df, left_on = "Locus", right_on = "locus", how = "left")

        #create columns con_ca3, con_ca4, con_ca5, con_ca6, con_sam_amp. Con stands for consensus
        mastertable["con_ca3"] = False
        mastertable["con_ca4"] = False
        mastertable["con_ca5"] = False
        mastertable["con_ca6"] = False
        mastertable["con_sam-amp"] = False
        #print columsn from mastertable
        #print(mastertable.columns)
        #if ca3 or NE_ca3 is True, set con_ca3 to True. Same for all signal molecules
        mastertable.loc[(mastertable["ca3"] == True) | (mastertable["NE_ca3"] == True), "con_ca3"] = True
        mastertable.loc[(mastertable["ca4"] == True) | (mastertable["NE_ca4"] == True), "con_ca4"] = True
        mastertable.loc[(mastertable["ca5"] == True) | (mastertable["NE_ca5"] == True), "con_ca5"] = True
        mastertable.loc[(mastertable["ca6"] == True) | (mastertable["NE_ca6"] == True), "con_ca6"] = True
        mastertable.loc[(mastertable["sam-amp"] == True) |(mastertable["NE_sam-amp"]), "con_sam-amp"] = True

        mastertable["multiple_signals"] = False
        #if more than one of the con columns is True, set multiple_signals to True
        mastertable.loc[(mastertable[["con_ca3", "con_ca4", "con_ca5", "con_ca6", "con_sam-amp"]].sum(axis = 1) > 1), "multiple_signals"] = True

        mastertable["effector_count_known_new_sum"] = 0
        #sum the number of effectors in the known and validated new effectors columns. Effectors to sum are can3, tirsaved, cam2, cam3, cam2, saved-chat, nucc, calpl, cami1, can1-2, csx1, csx23 csm6-ca6, cora
        mastertable["effector_count_known_new_sum"] = mastertable[["can3", "tirsaved", "cam3", "cam2", "saved-chat", "nucc", "calpl", "cami1", "can1-2", "csx1", "csx23", "csm6-ca6", "cora", "cam1"]].sum(axis = 1)

        #open the temperature data file
        temperature_data = pd.read_csv(str(input.temperature_data), sep = ",")
        #group by column "genus" and calculate mean temperature for each genus
        mean_temperatures = temperature_data.groupby("genus")["Topt_ave"].mean().dropna()
        print(mean_temperatures)
        #merge mean_temperature with mastertable
        mastertable = pd.merge(mastertable, mean_temperatures, on = "genus", how = "left")
        #rename topt_ave to "temperature"
        mastertable = mastertable.rename(columns = {"Topt_ave": "mean_temperature"})

        #Tidy the mastertable
        #drop the following columns (list)
        droppable_columns = ["sample_x", "Unnamed: 0", "val_x", "locus_x", "sample_y", "val_y", "locus_y", "Unnamed: 0_y", "mem_y", "Unnamed: 0_x"]

        #drop the columns in a for loop with try/catch for each column in case the mastertable changes in the future
        for column in droppable_columns:
            try:
                mastertable = mastertable.drop(columns = [column])
            except:
                pass

        #create a new column multilocus_sample. This a boolean value that reports whether the same sample has other loci
        mastertable["multilocus_sample"] = False

        #create a new column multisignal_sample. This is a boolean value that reports whether the same sample has multiple signals (ca3, ca4, ca6 or sam-amp)
        mastertable["multisignal_sample"] = False

        #create a new column all_signals_sample. This is a string that contains comma separated list of all signals present in the sample
        mastertable["all_signals_sample"] = ""

        #for each sample, go through loci and list which signal molecules are true and update the multisingal_sample and all_signals_sample columns accordingly
        for sample in mastertable["Sample"].unique():
            #get all loci for the sample
            sample_loci = mastertable[mastertable["Sample"] == sample]["Locus"].unique()
            #if the number of loci is greater than 1, set multilocus_sample to True
            if len(sample_loci) > 1:
                mastertable.loc[mastertable["Sample"] == sample, "multilocus_sample"] = True

            #get all signals for the sample
            sample_signals = mastertable[mastertable["Sample"] == sample][["con_ca3", "con_ca4", "con_ca5", "con_ca6", "con_sam-amp"]].sum(axis = 1)
            #if the number of signals is greater than 1, set multisignal_sample to True
            if sample_signals.sum() > 1:
                mastertable.loc[mastertable["Sample"] == sample, "multisignal_sample"] = True

            if sample_signals.sum() > 0:
                mastertable.loc[mastertable["Sample"] == sample, "signal_sample"] = True

            #using columns ca3, ca4, ca6 and sam-amp, create a string of all signals present in the sample by checking if their boolean value is True
            all_signals = ""
            if mastertable.loc[mastertable["Sample"] == sample, "con_ca3"].any() == True:
                all_signals += "ca3,"
            if mastertable.loc[mastertable["Sample"] == sample, "con_ca4"].any() == True:
                all_signals += "ca4,"
            if mastertable.loc[mastertable["Sample"] == sample, "con_ca6"].any() == True:
                all_signals += "ca6,"
            if mastertable.loc[mastertable["Sample"] == sample, "con_sam-amp"].any() == True:
                all_signals += "sam-amp,"
            #remove the last comma from the string
            all_signals = all_signals[:-1]
            #update the all_signals_sample column with the string
            mastertable.loc[mastertable["Sample"] == sample, "all_signals_sample"] = all_signals

        #create column effector_elsewhere_in_sample
        mastertable["effector_elsewhere_in_sample"] = False
        mastertable["effector_elsewhere_in_sample_but_not_here"] = False

        #for each locus, check if the sample has another locus that generates a signal and if so, set cas10_effector_elsewhere_in_sample to True
        for locus in mastertable["Locus"].unique():
            #get the sample for the locus
            sample = mastertable[mastertable["Locus"] == locus]["Sample"].unique()[0]
            #get all loci for the sample
            sample_loci = mastertable[mastertable["Sample"] == sample]["Locus"].unique()
            #if any of the loci in the sample have a signal, set cas10_effector_elsewhere_in_sample to True
            if mastertable[mastertable["Locus"].isin(sample_loci)][["con_ca3", "con_ca4", "con_ca6", "con_sam-amp"]].sum(axis = 1).any() == True:
                mastertable.loc[mastertable["Locus"] == locus, "effector_elsewhere_in_sample"] = True

            #if any of the loci in the sample have a signal but this locus does not have a signal, set cas10_effector_elsewhere_in_sample_but_not_here to True
            #first check if current locus has signal
            if mastertable.loc[mastertable["Locus"] == locus][["con_ca3", "con_ca4", "con_ca6", "con_sam-amp"]].sum(axis = 1).any() == False:
                #if it does not have a signal, check if any other loci in the sample have a signal
                if mastertable[mastertable["Locus"].isin(sample_loci)][["con_ca3", "con_ca4", "con_ca6", "con_sam-amp"]].sum(axis = 1).any() == True:
                    mastertable.loc[mastertable["Locus"] == locus, "effector_elsewhere_in_sample_but_not_here"] = True
            

        #subset the mastertable to include samples where (GGDD_hmm_boolean is True or GGDE_hmm_boolean is True) AND HD_hmm_boolean is False
        mastertable["Cas10_HD_coord"] = pd.to_numeric(mastertable["Cas10_HD_coord"], errors='coerce')
        mastertable_group4 = mastertable[
            ((mastertable["GGDD_hmm_boolean"] == True) | (mastertable["GGDE_hmm_boolean"] == True)) & 
            #(mastertable["HD_hmm_boolean"] == False) & 
            #((mastertable["Cas10_HD_coord"] > 100) | (mastertable["Cas10_HD"] == False)) & 
            (mastertable["unknown_proteins"] == True)
        ]
        mastertable_group4["Locus"].to_csv(str(output.group4_IDs), sep = "\t", index = False)

        #some loci are wrongly subtyped. Use this dictionary to manually correct them. Key shows locus, value shows corrected subtype
        subtype_corrections = {
            "GCA_003201765.2_3": "III-D",
            "GCF_003201765.2_2": "III-D",
            "GCA_022845955.1_1": "III-D",
            "GCA_001548075.1_0": "III-D",
            "GCF_014495845.1_1": "III-D",
            "GCA_000143965.1_0": "III-D",
            "GCA_000010725.1_0": "III-D",
            "GCF_019900845.1_1": "III-D",
            "GCF_000227745.2_1": "III-B",
            "GCA_000970265.1_31": "III-F",
            "GCF_000022425.1_4": "III-F",
            "GCF_000022405.1_3": "III-F",
            "GCF_009602405.1_0": "III-F",
            "GCA_003201675.2_2": "III-F",
            "GCF_003201675.2_3": "III-F",
            "GCF_028472365.1_2": "III-F",
            "GCA_000338775.1_1": "III-F",
            "GCA_000011205.1_3": "III-F",
            "GCA_003967175.1_0": "III-F",
        }

        #correct the subtypes for the given loci (column Locus)
        for locus, subtype in subtype_corrections.items():
            mastertable.loc[mastertable["Locus"] == locus, "Subtype"] = subtype

        #drop any duplicate rows
        mastertable = mastertable.drop_duplicates()
        mastertable.to_csv(str(output.final_info_table), sep = "\t", index = False)


rule node_graph:
    '''
    This rule produces files for Gephi to view effector colocalisation
    using a node graph. It also generates upset plots from the same data.

    For Gephi, it works by first creating a table with all effector combinations (edges) with the headers
    source, target, type and weight. It then looks at each locus in the mastertable
    and upon finding co-occurrence of an effector pair, it adds 1 to the weight column.
    '''
    input:
        mastertable = rules.mastercombiner.output.final_info_table,
        effector_list = rules.heatmap_known_validated_effectors.output.effectors_combined
    output:
        edges =  base_path + "/80_node_graph/edges_all.tsv",
        nodes = base_path + "/80_node_graph/nodes_all.tsv",
    params:
        outdir = base_path + "/80_node_graph",
        edges_basename = base_path + "/80_node_graph/edges.tsv",
        nodes_basename = base_path + "/80_node_graph/nodes.tsv"
    shell:
        '''
        python scripts/effector_nodes_and_other_viz.py -i {input.mastertable} -e {input.effector_list} -o {params.outdir} -n {params.nodes_basename} -d {params.edges_basename}
        '''

rule transmembrane_prediction:
    '''
    Using TMHMM predicts transmembrane regions for all proteins in all type III loci
    '''
    input:
        contig_proteins = rules.cATyper_hmm_search.output.contig_proteins
    output:
        tmhmm_results = base_path + "/90_transmembrane_prediction/{c}/{c}_tmhmm_results.tsv"
    

rule groupCharacteriser:
    '''
    HMMER based databases (Pfam and CARFSAVED)
    This takes as input a list of locus names (from rule mastercombiner) and performs characterisation of all unknown proteins in those loci.
    Characterisation involves:
        - HMM search against Pfam
        - HMM search against CARF/SAVED
    '''
    input:
        group4_list = rules.mastercombiner.output.group4_IDs,
        unknown_proteins = rules.concatenate_unknowns.output.proteins
    output:
        #group4_CARFSAVED = base_path + "/group4/group4.CARFSAVED.rawresults.tsv",
        group4_pfam = base_path + "/group4/group4.all.pfam.rawresults.tsv",
        clustered_unknowns = base_path + "/group4/group4_clustered_unknowns.faa",
        individual_fastas_created = base_path + "/group4/individual_fastas/fastas_created.done"
    params:
        outdir = base_path + "/group4",
        output_base = base_path + "/group4/group4",
        output_individual_fastas = base_path + "/group4/individual_fastas",
        pfams_temp = base_path + "/group4/group4.all.pfam.rawresults.txt.temp",
        tmhmm_model = "/media/volume/st_andrews/new_effectors/TM/TMHMM2.0.model"
    log:
        out = base_path + "/group4/logs/group4_characteriser.out",
        err = base_path + "/group4/logs/group4_characteriser.err"
    conda:
        "envs/groupChar.yaml"
    shell:
        '''
        python3 scripts/groupCharacteriser.py --input {input.group4_list} --unknown_proteins {input.unknown_proteins} --output_basename {params.output_base} --clustered_unknowns {output.clustered_unknowns} --individual_fastas_dir {params.output_individual_fastas} --tmhmm_model {params.tmhmm_model} --outputfolder {params.outdir} 2> {log.err} 1> {log.out}
        touch {output.individual_fastas_created}
        echo -e "{group4_pfam_headers}" > {output.group4_pfam}
        grep -v "#" {params.pfams_temp} >> {output.group4_pfam}
        '''

rule group4_PDB:
    '''
    This rule runs HHBlits using PDB as database for the group4 unknown proteins
    '''
    input:
        multifasta = rules.groupCharacteriser.output.clustered_unknowns
    output:
        hhsuite_concat = base_path + "/group4/pdb/group4_pdb.tsv",
    params:
        outdir = base_path + "/group4/pdb",
        pdb70 = "/media/volume/st_andrews/databases/pdb/pdb70",
        pdb30 = "/media/volume/st_andrews/databases/pdb30/pdb30",
    conda: "envs/hhsuite.yaml"
    threads: 40
    log:
        out = base_path + "/group4/pdb/logs/pdb.out",
        err = base_path + "/group4/pdb/logs/pdb.err",
    shell:
        '''
        python3 scripts/hhblits_wrapper.py --input {input.multifasta} --output_basepath {params.outdir} --database {params.pdb30}
        cat {params.outdir}/hhblits/*.tsv > {output.hhsuite_concat}
        '''

rule parse_hhsuite_group4_pdb:
    '''
    Parser output from HHSuite (PDB)
    '''
    input: rules.group4_PDB.output.hhsuite_concat
    output: base_path + "/group4/pdb/group4_pdb_parsed.tsv"
    conda: "envs/hhsuite.yaml"
    params:
        database = "PDB"
    shell:
        '''
        python3 scripts/hhsuite_parser.py --infile {input} --outfile {output}  --database {params.database}
        '''

rule group4_COGs:
    '''
    This rule runs HHBlits using COGs as database for the group4 unknown proteins
    '''
    input:
        multifasta = rules.groupCharacteriser.output.clustered_unknowns
    output:
        hhsuite_concat = base_path + "/group4/cog/group4_cog.tsv",
    params:
        outdir = base_path + "/group4/cog",
        cogs = "/media/volume/st_andrews/databases/cog/COG_KOG/COG_KOG"
    conda: "envs/hhsuite.yaml"
    threads: 40
    log:
        out = base_path + "/group4/cog/logs/cog.out",
        err = base_path + "/group4/cog/logs/cog.err",
    shell:
        '''
        python3 scripts/hhblits_wrapper.py --input {input.multifasta} --output_basepath {params.outdir} --database {params.cogs}
        cat {params.outdir}/hhblits/*.tsv > {output.hhsuite_concat}
        '''

rule parse_hhsuite_group4_cog:
    '''
    Parser output from HHSuite (PDB)
    '''
    input: rules.group4_COGs.output.hhsuite_concat
    output: base_path + "/group4/cog/group4_cog_parsed.tsv"
    conda: "envs/hhsuite.yaml"
    params:
        database = "COGs",
        mapping = "/media/volume/st_andrews/databases/cog/COG_KOG/cog-20.def.tab"
    shell:
        '''
        python3 scripts/hhsuite_parser.py --infile {input} --outfile {output} --database {params.database} --mapping {params.mapping}
        '''

rule group4_commonness:
    '''
    Blasts all group4 proteins against the proteomes of a CRISPR-Cas locus.
    The idea is to measure how common these proteins are and what kind of type III loci they are associated with.
    The rule first creates a diamond database for the proteome of each locus.
    Then, it blasts each group4 protein against each locus proteome using Diamond
    '''
    input:
        group4_proteins = rules.groupCharacteriser.output.clustered_unknowns,
        locus_proteins = rules.crispr_locus_proteins.output.crispr_locus_proteins
    output:
        blast_result = base_path + "/group4_prober/{c}/{c}.blast",
        diamond_db = base_path + "/group4_prober/{c}/{c}.dmnd"
    params:
        blast_result_temp = base_path + "/group4_prober/{c}/{c}.blast.temp"
    conda:
        "envs/diamond.yaml"
    shell:
        '''
        diamond makedb --in {input.locus_proteins} -d {output.diamond_db} --quiet
        diamond blastp --query {input.group4_proteins} --db {output.diamond_db} --outfmt 6 --out {params.blast_result_temp} --quiet
        awk -v locus={wildcards.c} '{{print $0,"\t",locus}}' {params.blast_result_temp} > {output.blast_result}
        rm {params.blast_result_temp}
        '''

rule concatenate_group4_commonness:
    '''
    Concatenates all blast results from group4_commonness into one file
    '''
    input: aggregate_group4_blasts
    output:
        concatenated_blast_results = base_path + "/group4_prober/group4_prober.blast"
    shell:
        '''
        echo -e "{blast_headers_group4}" > {output}
        cat {input} >> {output}
        '''

rule analyseGroup4Hits:
    input:
        blast_results = rules.concatenate_group4_commonness.output.concatenated_blast_results,
        pfam_hits = rules.groupCharacteriser.output.group4_pfam
    output:
        results_txt = base_path + "/group4_prober/group4_prober_analysis.txt"
    params:
        outpath = base_path + "/group4_prober"
    shell:
        '''
        python3 scripts/group4HitsAnalyser.py --input_blast {input.blast_results} --input_pfam_annotations {input.pfam_hits} --output {output.results_txt} --outpath {params.outpath}
        '''

rule effector_commonness:
    '''
    This rule counts the number of different effectors across the loci
    '''
    input:
        crispr_loci = rules.mastercombiner.output.final_info_table
    output:
        effector_commonness_tsv = base_path + "/13_effector_commonness/effector_commonness_master.tsv",
        #effector_commonness_plot = base_path + "/13_effector_commonness/effector_commonness.png"
    params:
        base_path = base_path + "/13_effector_commonness/effector_commonness"
    shell:
        '''
        python3 scripts/effector_commonness.py --input {input.crispr_loci} --output_basepath {params.base_path}
        '''

rule cctyper_gene_locations:
    '''
    This uses CCTyper outputs to generate a table that shows each cas gene's
    coordinate on the contig and other attributes. All attributes are:
    | protein_id     | start   | end     | effector | locus             | sample           | strand | type
    '''
    input:
        cas_operons_cctyper = base_path + "/071_cctyper_loci/{c}/cas_operons.tsv",
        cas_proteins = rules.crispr_locus_proteins.output.crispr_locus_proteins,
    output:
        cctyper_gene_locations_plottable = base_path + "/14_cctyper_gene_locations/{c}/{c}_cctyper_gene_locations.tsv"
    params:
        sample_folder = base_path + "/06_host_genomes/", #this is used to find the gff and proteins files that are strain-, not locus-specific
        this_folder = base_path + "/14_cctyper_gene_locations/{c}",
        outputfolder = base_path + "/14_cctyper_gene_locations/{c}",
        cctyper_folder = base_path + "/07_cctyper"
    log:
        out = base_path + "/14_cctyper_gene_locations/logs/{c}.out",
        err = base_path + "/14_cctyper_gene_locations/logs/{c}.err"
    shell:
        '''
        python scripts/cctyper_gene_locations.py --locus_id {wildcards.c} --sample_folder {params.sample_folder} --this_folder {params.this_folder} --outputfolder {params.outputfolder} --cas_operons {input.cas_operons_cctyper} --cctyper_path {params.cctyper_folder} --protein_fasta {input.cas_proteins} 2> {log.err} 1> {log.out}
        '''

rule locus_visualiser:
    '''
    This rule visualises a given range of a gff file.
    Using locus as wildcard.
    GFF comes from genome.
    Coordinates for CRISPR locus (and beyond) comes from 
    Takes CRISPR coordinates from 
    '''
    input:
        #genome_gff = rules.crispr_locus_proteins.output.crispr_locus_gff,
        crispr_positive_samples = base_path + "/071_cctyper_loci/{c}/cas_operons.tsv",
        cctyper_plottable = rules.cctyper_gene_locations.output.cctyper_gene_locations_plottable,
        plottable_effector_positions = rules.validated_new_effectors_analysis.output.plottable_effector_positions,
        known_effectors = rules.cATyper_analysis.output.plottable_effector_positions,
    output:
        visualisation = base_path + "/90_locus_viz/{c}/{c}_viz.png"
    conda:
        "envs/locus_visualiser.yaml"
    params:
        outdir = base_path + "/90_locus_viz/{c}",
        host_genomes_folder = base_path + "/06_host_genomes",
        cctyper_folder = base_path + "/07_cctyper",
    log:
        out = base_path + "/90_locus_viz/logs/{c}.out",
        err = base_path + "/90_locus_viz/logs/{c}.err",
    shell:
        '''
        python scripts/locus_visualiser.py --locus {wildcards.c} --cas_operons_file {input.crispr_positive_samples} --output_folder {params.outdir} --host_genomes_folder {params.host_genomes_folder} --cctyper_folder {params.cctyper_folder} --validated_effectors {input.plottable_effector_positions} --cctyper_protein_table {input.cctyper_plottable} --known_effector_table {input.known_effectors}  2> {log.err} 1> {log.out}
        '''

rule concatenate_locus_viz:
    input: aggregate_locus_viz
    output: base_path + "/90_locus_viz/locus_viz.done"
    shell:
        '''
        touch {output}
        '''
        

rule create_html_file:
    '''
    This creates an html file out of the mastertable and adds hyperlinks
    to the loci figures and the associated .tsv files
    '''
    input:
        mastertable = rules.mastercombiner.output.final_info_table,
    output:
        html_file = base_path + "/type_iii_mastertable.html"
    params:
        viz_dir_full = base_path + "/data",
        viz_dir_relative = "data",
    conda:
        "envs/excel_output.yaml"
    shell:
        '''
        python scripts/html_writer.py --mastertable {input.mastertable} --viz_dir_full {params.viz_dir_full} --viz_dir_relative {params.viz_dir_relative} --output_html {output.html_file}
        '''

rule casR_clustering:
    '''
    Creating clusters of everything annotated as CasR by CCTyper.
    Here we concatenate all CasRs, so we are not using loci-based wildcards.
    Instead we use output from concatenate_locus_viz as an indirect signal that all CasRs have been generated
    '''
    input:
        viz_done = rules.concatenate_locus_viz.output,
    output:
        clustered_CasR = base_path + "/casR_clustering/casR_clustered.faa",
        concatenated_casR = base_path + "/casR_clustering/casR_concatenated.faa",
    conda:
        "envs/trees.yaml"
    params:
        casR_wildcarded = base_path + "/14_cctyper_gene_locations/*/*casR.faa",
        cutoff = 0.4,
        n = 2,
    shell:
        '''
        cat {params.casR_wildcarded} > {output.concatenated_casR}
        cd-hit -i {output.concatenated_casR} -o {output.clustered_CasR} -c {params.cutoff} -n {params.n} -d 0 -M 16000 -T {threads}
        '''

rule final:
    input:
        tax = base_path + "/06_host_genomes/taxInfo.txt",
        tree_Cas10 = rules.Cas10_tree.output,
        type_iii_info = rules.combine_GGDD_HMM_to_mastertable.output.final_info_table,
        mastercombiner_final_info_table = rules.mastercombiner.output.final_info_table,
        concat_taxInfo = rules.concat_taxInfo.output,
        catyper_hmm = rules.concatenate_cATyper_hmm.output,
        catyper_analysis = rules.concatenate_cATyper_analysis.output,
        tree_CorA = rules.CorA_tree.output,
        tree_CorA_unclustered = rules.CorA_tree_unclustered.output,
        clustered_unknowns = rules.cluster_unknowns.output.proteins, #each row is an unknown protein
        clustered_unknowns_info = rules.concatenate_unknowns_locus_info.output.info, #contains locus-specific info on unknowns
        cas10_HD_faa = rules.Cas10_HD_hmm_maker.output.faa,
        cas10_hd_hmm = rules.Cas10_HD_hmmer.output,
        cas10_HD_hmm_merged = rules.merge_HD_hmm_with_info.output.merged_table,
        cas10_GGDD_faa = rules.Cas10_GGDD_hmm_maker.output.faa,
        cas10_GGDD_hmm = rules.Cas10_GGDD_hmmer.output,
        cas10_GGDD_hmm_merged = rules.merge_GGDD_hmm_with_info.output.merged_table,
        groupCharacteriser = rules.groupCharacteriser.output.group4_pfam,
        aggregate_crispr_locus_proteins = rules.aggregate_crispr_locus_proteins.output,
        concatenated_blast_results = rules.concatenate_group4_commonness.output,
        analyseGroup4Hits = rules.analyseGroup4Hits.output.results_txt,
        effector_commonness = rules.effector_commonness.output.effector_commonness_tsv,
        effector_scores_summary = rules.analyse_cATyper_effector_scores.output.effector_scores_summary,
        validated_effectors_scores_summary = rules.analyse_validated_new_effectors_scores.output.effector_scores_summary,
        concatenate_effector_wilcards = rules.concatenate_effector_wildcards.output,
        concatenate_cATyper_hmm_hhsuite = rules.concatenate_cATyper_hmm_hhsuite.output,
        parse_hhsuite = rules.concatenate_cATyper_hmm_hhsuite.output,
        parse_hhsuite_cogs = rules.concatenate_cATyper_hmm_hhsuite_cogs.output,
        group4_pdb = rules.parse_hhsuite_group4_pdb.output,
        group4_cog = rules.parse_hhsuite_group4_cog.output,
        cora_neighbourhood = rules.concatenate_cora_neighbourhoods.output,
        concatenate_validated_new_effectors_analysis = rules.concatenate_validated_new_effectors_analysis.output,
        heatmap_known_validated_effectors = rules.heatmap_known_validated_effectors.output.effector_scores_summary,
        node_graph = rules.node_graph.output.edges,
        concatenate_locus_viz = rules.concatenate_locus_viz.output,
        html = rules.create_html_file.output.html_file,
        casR_cluster = rules.casR_clustering.output.clustered_CasR,
        validated_effectors_cas10_fusions = rules.validated_effectors_cas10_fusions.output.hmm_rows,
        known_effectors_cas10_fusions = rules.known_effectors_cas10_fusions.output.hmm_rows,
    output: base_path + "/done"
    shell:
        '''
        touch {output}
        mkdir -p /media/volume/st_andrews/new_effectors/{project}/R
        mkdir -p /media/volume/st_andrews/new_effectors/{project}/R/effector_trees

        mkdir -p /media/volume/st_andrews/new_effectors/{project}/R/effector_hmm
        mkdir -p /media/volume/st_andrews/new_effectors/{project}/R/effector_hmm/pfam
        mkdir -p /media/volume/st_andrews/new_effectors/{project}/R/effector_hmm/pdb
        mkdir -p /media/volume/st_andrews/new_effectors/{project}/R/effector_hmm/cog

        mkdir -p /media/volume/st_andrews/new_effectors/{project}/R/group4
        mkdir -p /media/volume/st_andrews/new_effectors/{project}/R/group4/pfam
        mkdir -p /media/volume/st_andrews/new_effectors/{project}/R/group4/pdb
        mkdir -p /media/volume/st_andrews/new_effectors/{project}/R/group4/cog
        mkdir -p /media/volume/st_andrews/new_effectors/{project}/R/group4/carfsaved
        
        cp /media/volume/st_andrews/new_effectors/{project}/12_Cas10_tree/cas10_tree.txt /media/volume/st_andrews/new_effectors/{project}/R
        cp /media/volume/st_andrews/new_effectors/{project}/06_host_genomes/taxInfo.txt /media/volume/st_andrews/new_effectors/{project}/R
        cp /media/volume/st_andrews/new_effectors/{project}/09_crispr_iii_CorA/loci/type_iii_info.tsv /media/volume/st_andrews/new_effectors/{project}/R
        cp /media/volume/st_andrews/new_effectors/{project}/30_unknown_effectors/unknowns_info_loci.tsv /media/volume/st_andrews/new_effectors/{project}/R
        cp /media/volume/st_andrews/new_effectors/{project}/13_effector_commonness/effector_commonness_master.tsv /media/volume/st_andrews/new_effectors/{project}/R
        cp /media/volume/st_andrews/new_effectors/{project}/mastertable_v2.tsv /media/volume/st_andrews/new_effectors/{project}/R

        cp -r /media/volume/st_andrews/new_effectors/{project}/45_effector_tree/* /media/volume/st_andrews/new_effectors/{project}/R/effector_trees

        cp -r /media/volume/st_andrews/new_effectors/{project}/43_effector_hmmer_analysis/*/*_sorted_filtered_mapped.tsv /media/volume/st_andrews/new_effectors/{project}/R/effector_hmm/pfam
        cp -r /media/volume/st_andrews/new_effectors/{project}/42_effector_hhsuite/*/*_parsed.tsv /media/volume/st_andrews/new_effectors/{project}/R/effector_hmm/pdb
        cp -r /media/volume/st_andrews/new_effectors/{project}/42_effector_hhsuite_cogs/*/*_parsed_cogs.tsv /media/volume/st_andrews/new_effectors/{project}/R/effector_hmm/cog
        
        cp -r /media/volume/st_andrews/new_effectors/{project}/group4/group4.CARFSAVED.info.tsv /media/volume/st_andrews/new_effectors/{project}/R/group4/carfsaved
        cp -r /media/volume/st_andrews/new_effectors/{project}/group4/group4.pfam.info.tsv /media/volume/st_andrews/new_effectors/{project}/R/group4/pfam
        cp -r /media/volume/st_andrews/new_effectors/{project}/group4/cog/group4_cog_parsed.tsv /media/volume/st_andrews/new_effectors/{project}/R/group4/cog
        cp -r /media/volume/st_andrews/new_effectors/{project}/group4/pdb/group4_pdb_parsed.tsv /media/volume/st_andrews/new_effectors/{project}/R/group4/pdb
        cp -r /media/volume/st_andrews/new_effectors/{project}/group4_prober/group4_prober_analysis.txt /media/volume/st_andrews/new_effectors/{project}/R/group4

        mkdir -p /media/volume/st_andrews/new_effectors/{project}/R/cora_neighbourhood
        cp -r /media/volume/st_andrews/new_effectors/{project}/52_cora_neighbour_analysis/cora_neighbours.tsv /media/volume/st_andrews/new_effectors/{project}/R/cora_neighbourhood
        cp -r /media/volume/st_andrews/new_effectors/{project}/45_effector_tree/cora_tree.txt /media/volume/st_andrews/new_effectors/{project}/R/cora_neighbourhood
        cp -r /media/volume/st_andrews/new_effectors/{project}/44_effector_alignments/cora.afa /media/volume/st_andrews/new_effectors/{project}/R/cora_neighbourhood
        cp -r /media/volume/st_andrews/new_effectors/{project}/17_CorA_tree/CorA_tree.txt /media/volume/st_andrews/new_effectors/{project}/R/cora_neighbourhood
        cp -r /media/volume/st_andrews/new_effectors/{project}/16_CorA_align/CorA_alignment.afa /media/volume/st_andrews/new_effectors/{project}/R/cora_neighbourhood
        '''
        