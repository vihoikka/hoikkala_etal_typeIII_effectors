# hoikkala_etal_typeIII_effectors
A Snakemake pipeline that characterises type III CRISPR-Cas loci and aims to discover new cOA-dependent effectors. Scripts used to generate the data on Hoikkala, Graham and White 2024 paper on new type III CRISPR-Cas effectors.

The run command that reproduces the results of the paper:

```snakemake --snakefile new_effectors.smk --use-conda --cores 40 --config getGenomesBy="local" genome_mode="all" cas10_anchor="True" cas10_tree_cluster="False"```

## Requirements
- A conda environment with Snakemake 7 and Pandas 1.5.2 installed
- All other dependencies are installed by the pipeline when running with the --use-conda option
- Genome and HMM databases (see below)
- Ubuntu (not tested on other Linux distributions)

### Expected runtime after downloading and preparing databases is about 24 hours on 40-core HPC.

## Databases
- NCBI complete prokaryotic genomes. Note that these are around 400 gb in total size. Use the NCBI datasets tool to download the genomes and prepare the files using these steps:
  1. Install the datasets tool with Conda using ```conda install -c conda-forge ncbi-datasets-cli```
  2. Download Bacterial and Archaeal genomes to separate folders in steps iii and iv:
  3. Go to your desired folder and download Bacterial genomes ```datasets download genome taxon 2 --filename bacteria_genomes.zip --include gff3,genome,protein --dehydrated --annotated --assembly-level complete --assembly-source RefSeq```
  4. Also download Archaeal genomes ```datasets download genome taxon 2157 --filename archaea_genomes.zip --include gff3,genome,protein --annotated --assembly-level complete```
  5. Unzip both downloads
  6.  "Rehydrate" both downloaded genome sets: ```datasets rehydrate --directory <directory_name>```, where directory points (relative to current location) to the unzipped data folders
  7. Change ```genomes_folder``` path parameter in the snakemake script to match the bacterial folder and ```archaea_folder``` to match the Archaeal folder
  8. Copy all Archaeal genomes to the Bacterial genomes folder. The Bacterial genomes folder now contains *all* prokaryotic genomes.
  9. When the pipeline is run, list of genome IDs with domain information are automatically generated based on contents of the Bacterial and Archaeal folders. This is a crude approach, but since the taxonomy files from NCBI do not clearly state the domain of the organism, we need to do it this way. You could also try symlinks for the Archaeal genomes instead of copying them to save space.

- HMM databases:
  - Known effectors: download from XX and put in path defined by variable "hmm_database_folder" in the Snakemake pipeline
  - New effectors: download from XX and put in path defined by variable "validated_effectors_hmm_db_path" in the Snakemake pipeline
  - Cas10: the nuclease and cyclase domain HMM profiles are automatically generated by the pipeline. You can alternatively provide your own HD HMM profiles using the argument ```modified_cas10_hd [path-to-hmm-file]``` in the rule Cas10_HD_hmm_maker. This option has not been implemented for the cyclase domain.

## Configuration arguments
 - getGenomesBy: "local" or "remote". Currently only "local" works in most downstream rules, but this option is left here in case of further development.
 - genome_mode: "all", "random", "taxid". Applicable when getGenomesBy == "local". "All" uses all given genomes, "random" subsamples the number of genomes denoted by the global parameter genome_count, "taxid" only takes genomes whose taxonomic ID matches an ID in a file provided by the taxidlistfile argument (see below)
 - taxidlistfile: a text file that contains a comma-separated list of taxonomic IDs (e.g. "105,104,103,102,101" without quotation marks). Applicable when genome_mode is "taxid".
 - cas10_anchor: filter genomes by the presence of cas10
 - cas10_tree_cluster: whether to cluster Cas10s further for the phylogenetic tree (after initial clustering)

## Main steps of the pipeline
1. Obtains bacterial and archaeal genomes from local files (there's also a remote option if needed, see "Configuration arguments"). Also creates annotation file to hold information on which acc number corresponds to bacteria and which to archaea
2. Filters all genomes by presence of Cas10 using an HMM profile (if cas10_anchor enabled)
3. Extracts Cas10 protein sequences and creates symlinks to the corresponding host genomes
4. Extracts host taxonomy infromation from the genome files (rule getTaxInfo)
5. Runs CRISPRCasTyper on Cas10-assocation genomes
6. Proceeds only with type III loci
7. Runs an in-house script cATyper on all loci, looking for pre-established effector proteins using custom HMM profiles
8. Extract Cas10, Cas7, Cas5 and CorA fastas from type III loci. Also mark down if Cas10 contains HD or GGDD on a sequence comparison level (rule typeIII_characterizer)
9. Create HMM profiles based on pre-flagged HD-domain or GGDD-domain containing Cas10 (rules rule Cas10_HD_hmm_maker and rule Cas10_GGDD_hmm_maker). Note that different profiles are used for GGDD and GGDE seeded cyclase domains.
10. Looks for unexpected or poor E-value scoring proteins in CRISPR-Cas loci. These proteins are flagged as being potentially interesting (rule unknown_finder)
11. Various rules characterizing different proteins such as CorA and Cas7
12. Create Cas10 tree using muscle -based alignment and fasttree
13. Combine all results together into a mastertable
14. Extract loci that fulfill specific criteria (e.g. no HD, has GGDD, has potential unknown effectors in locus) and examine the unknown proteins in these loci by HMM searches against local SAVED/CARF and PFAM databases

## Plotting and models
A separate R project is used to plot data and create models. The R project and scripts are in the folder "R_scripts". One of the folders the pipeline creates is called R. Copy this folder inside R_scripts/data and then modify the R script to point to this folder (variable "project" must match the folder name, so if not renamed, project = "R").
