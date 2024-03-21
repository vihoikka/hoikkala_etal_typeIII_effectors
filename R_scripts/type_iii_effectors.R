# R scripts for reproducing figures in the manuscript by Hoikkala, Graham and White (2024) on CRISPR-Cas type III effectors

#### LIBRARIES ####

devtools::install_github("quinten-goens/plotly.R@fix/kaleido-export-bug") #fixes saving bug, see https://github.com/plotly/plotly.R/pull/2228

install.packages("wakefield")
install.packages("viridis")
install.packages("phytools")
install.packages("shadowtext")
install.packages("ggnewscale")
install.packages("shiny")
install.packages("xlsx")
install.packages("ggrepel")
install.packages("ggnewscale")
install.packages("ragg")
install.packages("ggtreeExtra")
install.packages("ggrepel")
install.packages("ggforce")
install.packages("tidyverse")
install.packages("ggplot2")
install.packages("UpSetR")
if (!require("processx")) install.packages("processx")
if(!require(devtools)) install.packages("devtools")
devtools::install_github("krassowski/complex-upset")
install.packages('ComplexUpset')
install.packages('ggbeeswarm')

BiocManager::install("ggbio")
BiocManager::install("GenomicRanges")


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ggtree")
BiocManager::install("ggtreeExtra", force = TRUE)

library(tidyverse)
library(RColorBrewer)
library(ggbio)
library(GenomicRanges)
library(S4Vectors)
library(gggenes)
library("treeio")
library("wakefield")
library(stringr)
library("viridis")
library(shadowtext)
library("phytools")
library(dplyr)
library(ggplot2)
library(ggnewscale)
library("xlsx")
library(ggrepel)
library("scales")
library("patchwork")
library(ragg)
library("ggtreeExtra")
library("ggplot2")
library(ggforce)
library("ggtree")
library(tidyverse)
library("gggenes")
library(stringr)
library(dplyr)
library(forcats)
library(RColorBrewer)
library(Biostrings)
library(tidyverse)
library(grid)
library(tidyr)
library(UpSetR)
library(plotly)
library(reticulate)
library(ComplexUpset)
library(ggbeeswarm)

#### Load data and tidy ####
project = "190224" #modify this to the project name. The folder must be in the data folder (relative to current R script)

#define paths
info_path <- paste("data/",project,"/mastertable_v2.tsv", sep = "")
cas10_tree_path <- paste("data/",project,"/cas10_tree.txt", sep = "")
effectors <- paste("data/", project, "/effector_commonness_master.tsv", sep = "")
tax_info_path <- paste("data/", project, "/taxInfo.txt", sep = "")

#open files
cas10_tree <- read.newick(cas10_tree_path)
cas10_tree <- midpoint.root(cas10_tree)
effectors <- read.csv(effectors, sep="\t")
tax_info <- read.csv(tax_info_path, sep="\t") #produced by a separate rule in the Snakemake pipeline
info <- read.csv(info_path, sep="\t")

info <- subset(info, info$Cas10_length > 499) #this was previously set to 500, resulting in a mismatch between tree and info df

#create folder plots
if (!dir.exists(paste0("data/",project,"/plots"))) {
  dir.create(paste0("data/",project,"/plots"))
}

#create folder plots/effectors
if (!dir.exists(paste0("data/",project,"/plots/effectors"))) {
  dir.create(paste0("data/",project,"/plots/effectors"))
}

#Reduce hybrid naming to the type III component
info <- info %>%
  mutate(Subtype = str_extract_all(Subtype, "III-\\w"), # extract all "III-letter"'s and get a list of character vectors
         Subtype = sapply(Subtype, paste0, collapse=", ")) # concatenate multiple III sequences if they exist in one string

#In cases where the subtype is "III-A, III-D" rename to III-D
info$Subtype[info$Subtype == "III-A, III-D"] <- "III-D"

#mark thermophiles based on temperature
info$thermophile <- FALSE
info$thermophile[info$mean_temperature > 45] <- TRUE

#Mark hyperthermophiles
info$hyperthermophile <- FALSE
info$hyperthermophile[info$genus == "Thermus" | info$genus == "Thermotoga" | info$genus == "Aquifex"] <- TRUE

#change variable types
info$has_known_effector <- as.logical(info$has_known_effector)
info$has_validated_new_effector <- as.logical(info$has_validated_new_effector)
info$multiple_signals <- as.logical(info$multiple_signals)

#Rename hybrid loci
info$Hybrid <- FALSE #no hybrids is default
info$Hybrid[str_detect(info$Subtype, 'Hybrid') == TRUE] <- TRUE #mark some as hybrids
info$Subtype[str_detect(info$Subtype, 'Hybrid') == TRUE & str_detect(info$Subtype, 'III-A') == TRUE] <- "III-A"
info$Subtype[str_detect(info$Subtype, 'Hybrid') == TRUE & str_detect(info$Subtype, 'III-B') == TRUE] <- "III-B"
info$Subtype[str_detect(info$Subtype, 'Hybrid') == TRUE & str_detect(info$Subtype, 'III-C') == TRUE] <- "III-C"
info$Subtype[str_detect(info$Subtype, 'Hybrid') == TRUE & str_detect(info$Subtype, 'III-D') == TRUE] <- "III-D"

info$effector_count_known_new_sum <- as.numeric(info$effector_count_known_new_sum)
#create a boolean column version of the above if > 0
info$has_effector <- FALSE
info$has_effector[info$effector_count_known_new_sum > 0] <- TRUE

#transform to some strings to boolean
info$Cas10_HD <- as.logical(info$Cas10_HD)

info$Cas10_GGDD <- as.logical(info$Cas10_GGDD)
info$unknown_proteins <- as.logical(info$unknown_proteins)

info$Cas10_HD_coord <- as.numeric(info$Cas10_HD_coord)

info$effector_elsewhere_in_sample_but_not_here <- as.logical(info$effector_elsewhere_in_sample_but_not_here)
info$effector_elsewhere_in_sample <- as.logical(info$effector_elsewhere_in_sample)
info$multisignal_sample <- as.logical(info$multisignal_sample)
info$GGDD_hmm_boolean <- as.logical(info$GGDD_hmm_boolean)
info$cyclase <- as.logical(info$cyclase)

# in info$signal_sample, fill empty cells with FALSE
info$signal_sample <- as.logical(info$signal_sample)
info$signal_sample[is.na(info$signal_sample)] <- FALSE
#change all these columns to logical: can3, tirsaved, tmcarf, cramp2, cramp1, saved.chat, nucc, calpl, cami1, can12, csx1, csx23, csm6.ca6, cora
info[,c("can3", "tirsaved", "cam1", "cam2", "cam3", "saved.chat", "nucc", "calpl", "cami1", "can1.2", "csx1", "csx23", "csm6.ca6", "cora")] <- lapply(info[,c("can3", "tirsaved", "cam1", "cam2", "cam3", "saved.chat", "nucc", "calpl", "cami1", "can1.2", "csx1", "csx23", "csm6.ca6", "cora")], as.logical)
#same for con_ca3, con_ca4, con_ca6, con_sam_amp
info[,c("con_ca3", "con_ca4", "con_ca6", "con_sam.amp")] <- lapply(info[,c("con_ca3", "con_ca4", "con_ca6", "con_sam.amp")], as.logical)

info$has_multiple_effectors <- ifelse(info$effector_count_known_new_sum > 1, TRUE, FALSE)

#set rownames by locus
rownames(info) <- info$Locus

#setup Cas10 length plot for Cas10 tree
info_for_cas10length <- select(info, Locus, Cas10_length, Cas10_HD, mean_temperature)
#info_for_cas10length <- na.omit(info_for_cas10length)

#this renaming must be done. Otherwise will get an error "Error: object 'Cas10_length' not found"
info_for_cas10length <- dplyr::rename(info_for_cas10length, Cas10_length_plot = Cas10_length)
info_for_cas10length <- dplyr::rename(info_for_cas10length, Cas10_HD_plot = Cas10_HD)
info_for_cas10length <- dplyr::rename(info_for_cas10length, mean_temperature_plot = mean_temperature)

#### Colours and other parameters ####

colorBlindBlack8  <- c("#404040", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
colorBlindBlack8_2  <- c("#000000", "#D55E00")
label_fontsize = 2

#### Known effector trees ####
# Effector trees
tree_dir <- paste("data/", project, "/effector_trees", sep = "")
# Get a list of all the .txt files in that directory
tree_files <- list.files(path = tree_dir, pattern = "*.txt", full.names = TRUE)
trees <- list()
ggtrees <- list()
effector_analyses <- list()
patchwork_list <- list()
patch_index <- 1 # index to control when to start new row
trees_length = as.numeric(length(seq_along(tree_files)))

for(i in seq_along(tree_files)) {
  print(i)
  #i <- 1
  basename <- sub("_tree.txt", "", basename(tree_files[i]))
  trees[[i]] <- read.tree(tree_files[i])
  
  #read effector analysis for Pfam
  effector_analysis <- paste0("data/",project,"/effector_hmm/pfam/", "/", basename, "_sorted_filtered_mapped.tsv")
  effector_analysis_raw <- read.csv(effector_analysis, sep="\t") #this will be unchanged data
  
  effector_analysis <- read.csv(effector_analysis, sep="\t")
  # Set "locus" as the row names
  
  #create locus column
  effector_analysis$locus <- sapply(strsplit(effector_analysis$query_name, "__"), `[`, 1)
  effector_analysis$locus_protein_domain <- paste0(effector_analysis$query_name,"_",effector_analysis$X.)
  
  #group by query and combine other columns into comma separated strings
  effector_analysis <- effector_analysis %>%
    group_by(query_name) %>%
    summarise_all(~paste(., collapse = ","))
  
  #extract protein length
  effector_analysis$protein_length <- sapply(strsplit(effector_analysis$qlen, ","), function(x) as.numeric(x[1]))
  effector_analysis$locus <- sapply(strsplit(effector_analysis$locus, ","), function(x) as.character(x[1]))
  
  # Add the protein_length column from df2 to df1
  effector_analysis_raw <- merge(x = effector_analysis_raw, y = effector_analysis[,c("query_name","protein_length")], by = "query_name", all.x = TRUE)
  #add subtype to effector analysis
  effector_analysis <- left_join(effector_analysis, info[c('Locus', 'Subtype', 'domain')], by = c('locus' = 'Locus'))
  
  #convert from tibble to df
  effector_analysis <- as.data.frame(effector_analysis)
  rownames(effector_analysis) <- effector_analysis$query_name #set rownames
  
  effector_analyses[[i]] <- effector_analysis
  
  effector_analysis_raw_clean <- data.frame(
    query_name_plot = as.character(effector_analysis_raw$query_name),
    ali_start_plot = as.numeric(effector_analysis_raw$ali_start),
    ali_end_plot = as.numeric(effector_analysis_raw$ali_end),
    protein_length = as.numeric(effector_analysis_raw$protein_length),
    target_name = as.character(effector_analysis_raw$target_name)
  )
  
  effector_analysis$Subtype[is.na(effector_analysis$Subtype)] <- "Unknown"
  effector_analysis$Subtype <- as.factor(effector_analysis$Subtype)
  
  #get some more info by merging from info table. In effector_analysis the column is "locus" and in info it is "Locus"
  effector_analysis <- left_join(effector_analysis, dplyr::select(info, Locus, cora), by = c('locus' = 'Locus'))
  
  if (!dir.exists(paste0("data/",project,"/plots/effectors"))) {
    dir.create(paste0("data/",project,"/plots/effectors"))
  }
  
  tsv_dir <- paste0("data/",project,"/plots/effectors/tables")
  
  if (!dir.exists(tsv_dir)) {
    dir.create(tsv_dir)
  }
  
  tsv_filename <- paste0("data/",project,"/plots/effectors/tables/", sub("_tree.txt", "_data.tsv", basename(tree_files[i])))
  write.table(effector_analysis, tsv_filename, sep = "\t", row.names = FALSE, col.names = TRUE)
  
  # create a ggtree object
  gtree <- ggtree(trees[[i]], 
                  aes()) %<+% effector_analysis
  
  gtree <- gtree + geom_tiplab(aes(label=""), 
                               align=TRUE, family='mono',
                               linetype = "dotted", linesize = .7, color = "grey")
  
  gtree <- gtree + ggtitle(basename) + geom_point2(aes(subset = (isTip & domain == "Archaea")), color = "red", size = 1.5, alpha = 0.6)
  
  if (i == 1) { #if Ae2, will mark associated CorA
    gtree <- gtree + ggtitle(basename) + geom_point2(aes(subset = (isTip & cora == "True")), color = "darkgreen", size = 4, shape = 18, alpha = 1)
  }
  
  ggtrees[[i]] <- gtree
  
  patchwork_list[[i]] <- gtree
  
  #Access tree node information
  nodes <- subset(gtree$data, isTip == TRUE) #pick only tips (leafs)
  leaf_nodes <- nodes %>% dplyr::select(label, y) #select the label and y columns
  
  #gets leaf order from tree
  leaf_order <- gtree$data %>% 
    filter(isTip) %>%  
    arrange(y) %>% 
    pull(label)
  
  # Extract unique proteins
  unique_proteins <- unique(effector_analysis_raw$query_name)
  
  # For each protein, construct a GRangesList
  gr_list <- lapply(unique_proteins, function(protein) {
    # Extract data for this protein
    protein_data <- effector_analysis[effector_analysis$query_name == protein, ]
    domains_data <- effector_analysis_raw[effector_analysis_raw$query_name == protein, ]
    
    # Create a GRanges object for the protein length
    gr_protein <- GRanges("yes", IRanges(start = 1, end = protein_data$protein_length, names = protein))
    #seqinfo = Seqinfo(protein, seqlengths = protein_data$protein_length)
    gr_protein$type <- "protein"
    gr_protein$type <- as.factor(gr_protein$type)
    
    # Create a GRanges object for the domains
    gr_domains <- GRanges("yes", IRanges(start = domains_data$ali_start, end = domains_data$ali_end, names = domains_data$query_name))
    gr_domains$type <- domains_data$target_name
    gr_domains$type <- as.factor(gr_domains$type)
    
    # return Combined GRanges
    return(c(gr_protein, gr_domains))
  })
  
  # Convert list to GRangesList
  grl <- GRangesList(gr_list)
  
  #export data from GRangesList
  names <- grl@unlistData@ranges@NAMES
  starts <- grl@unlistData@ranges@start
  ends <- grl@unlistData@ranges@width
  types  <- grl@unlistData@elementMetadata@listData$type
  grl_df <- data.frame(label = names, start=starts, end=ends, type=types)
  
  #pfam y positions
  grl_df <- leaf_nodes %>%
    left_join(grl_df, by = "label")
  names(grl_df)[names(grl_df) == "y"] <- "y_position"
  
  
  #PDB30/70 results
  pdb_path <- paste0("data/",project,"/effector_hmm/pdb/", basename, "_hhsuite_parsed.tsv")
  pdb_df <- read.csv(pdb_path, sep="\t")
  names(pdb_df)[names(pdb_df) == "query"] <- "label"
  
  #pdb30 y positions
  if (nrow(pdb_df) > 0) {
    pdb_df <- leaf_nodes %>%
      left_join(pdb_df, by = "label")
    names(pdb_df)[names(pdb_df) == "y"] <- "y_position"
  }
  
  #retain only best score in PDB30/70
  pdb_df_best_score <- pdb_df %>%
    group_by(label) %>%
    arrange(desc(score)) %>%
    dplyr::slice(1)
  
  
  #COGs results
  cog_path <- paste0("data/",project,"/effector_hmm/cog/", basename, "_hhsuite_parsed_cogs.tsv")
  cog_df <- read.csv(cog_path, sep="\t")
  names(cog_df)[names(cog_df) == "query"] <- "label"
  
  #Remove redundant strings in the PDB description (such as "Crystal structure of...")
  pdb_df_best_score$description <- sub("Crystal structure of a", "", pdb_df_best_score$description)
  pdb_df_best_score$description <- sub("Crystal structure of ", "", pdb_df_best_score$description)
  pdb_df_best_score$description <- sub("Structure of ", "", pdb_df_best_score$description)
  
  pdb_df_best_score <- subset(pdb_df_best_score, !is.na(description))
  
  
  
  backbone_color = c("protein" = "lightgrey")
  other_levels = setdiff(unique(grl_df$type), "protein")
  if(length(other_levels) > 0){
    other_colors = setNames(brewer.pal(length(other_levels)+2, "Set3"), other_levels)
  } else {
    other_colors = setNames(character(0), character(0))
  }
  
  final_colors = c(backbone_color, other_colors)
  
  #arrow height calculation based on number of rows, constraining between 1 and 4
  min_labels = 5
  max_labels = 30
  
  # Calculate proportion of data size
  prop = (length(unique(grl_df$label)) - min_labels) / (max_labels - min_labels)
  
  # Scale arrow height between 1 and 4, with linear approach
  arrow_height = 4 - 3 * min(max(prop, 0), 1)  
  
  
  # Reorder factor levels so that proteins are drawn before domains
  grl_df <- grl_df %>% arrange(end - start)
  
  grl_df <- grl_df %>%
    group_by(label) %>% # Group by each unique sample
    arrange(type == "protein", end - start) %>%  # Sort by type and size within each group
    ungroup() # Remove grouping
  
  grl_df$length <- abs(grl_df$end - grl_df$start)
  
  # First draw only "protein", then draw everything else
  tree_facet <- gtree
  
  composite <- tree_facet +
    geom_facet(
      panel = "Pfam",
      data = grl_df,
      geom = geom_hline,
      alpha = 0.4,
      linetype = "dashed",
      show.legend = FALSE,
      mapping = aes(
        yintercept = y_position,
      )
    )
  
  
  # data for first subplot
  data1 <- grl_df[grl_df$type == 'protein',]
  if (nrow(data1) > 0) {
    composite <- composite +
      geom_facet(
        panel = "Pfam",
        data = data1,
        geom = geom_gene_arrow,
        fill = "grey",
        alpha = 0.8,
        arrowhead_width = grid::unit(1, "mm"),
        arrowhead_height = grid::unit(arrow_height, "mm"),
        arrow_body_height = grid::unit(arrow_height, "mm"),
        inherit.aes = FALSE,
        stat = "identity",
        show.legend = FALSE,
        mapping = aes(
          xmin = start,
          xmax = end,
          y = y_position,
          group = label,
        )
      )
  }
  
  # data for second subplot
  data2 <- subset(grl_df, type != 'protein' & length > 49)
  if (nrow(data2) > 0) {
    composite <- composite +
      geom_facet(
        panel = "Pfam",
        data = data2,
        geom = geom_gene_arrow,
        alpha = 0.8,
        arrowhead_width = grid::unit(0, "mm"),
        arrowhead_height = grid::unit(arrow_height, "mm"),
        arrow_body_height = grid::unit(arrow_height, "mm"),
        inherit.aes = FALSE,
        stat = "identity",
        show.legend = FALSE,
        mapping = aes(
          xmin = start,
          xmax = end,
          y = y_position,
          group = label,
          fill = type
        )
      )
  }
  
  # data for third subplot
  data3 <- subset(grl_df, type != 'protein' & length < 50)
  if (nrow(data3) > 0) {
    composite <- composite +
      geom_facet(
        panel = "Pfam",
        data = data3,
        geom = geom_gene_arrow,
        alpha = 0.8,
        arrowhead_width = grid::unit(0, "mm"),
        arrowhead_height = grid::unit(arrow_height, "mm"),
        arrow_body_height = grid::unit(arrow_height, "mm"),
        inherit.aes = FALSE,
        stat = "identity",
        show.legend = FALSE,
        mapping = aes(
          xmin = start,
          xmax = end,
          y = y_position,
          group = label,
          fill = type
        )
      )
  }
  
  grl_df$protein_id <- sapply(strsplit(grl_df$label, split = "__"), '[', 2)
  grl_df <- grl_df %>%
    left_join(effector_analysis[c('query_name','target_name')], by = c('label' = 'query_name'))
  
  grl_df$target_name <- sub(",protein", "", grl_df$target_name)
  grl_df$target_name <- sub("protein", "", grl_df$target_name)
  grl_df$plot_annotation <- paste(grl_df$protein_id, " (", grl_df$target_name, ")", sep = "")
  
  composite <- composite + geom_facet(
    panel = "Pfam acc. + descr.",
    data = grl_df,
    geom = geom_text2,
    hjust = 0,
    size = 2,
    inherit.aes = FALSE,
    stat = "identity",
    show.legend = FALSE,
    mapping = aes(subset=isTip, label=plot_annotation, y= y_position, x = 0)) +
    xlim_expand(c(0,0.03), panel = "Pfam acc. + descr.")
  
  composite <- composite +
    geom_facet(
      panel = "PDB30",
      data = grl_df,
      geom = geom_hline,
      alpha = 0.4,
      linetype = "dashed",
      show.legend = FALSE,
      mapping = aes(
        yintercept = y_position,
      )
    )
  
  
  if (nrow(data1) > 0) {
    composite <- composite +
      geom_facet(
        panel = "PDB30",
        data = data1,
        geom = geom_gene_arrow,
        alpha = 0.8,
        fill = "grey",
        arrowhead_width = grid::unit(1, "mm"),
        arrowhead_height = grid::unit(arrow_height, "mm"),
        arrow_body_height = grid::unit(arrow_height, "mm"),
        inherit.aes = FALSE,
        stat = "identity",
        show.legend = FALSE,
        mapping = aes(
          xmin = start,
          xmax = end,
          y = y_position,
          group = label
        )
      )
  }
  
  if (nrow(pdb_df_best_score) > 0) {
    composite <- composite +
      geom_facet(
        panel = "PDB30",
        data = pdb_df_best_score,
        geom = geom_gene_arrow,
        alpha = 0.8,
        arrowhead_width = grid::unit(1, "mm"),
        arrowhead_height = grid::unit(arrow_height, "mm"),
        arrow_body_height = grid::unit(arrow_height, "mm"),
        inherit.aes = FALSE,
        stat = "identity",
        show.legend = FALSE,
        mapping = aes(
          xmin = qstart,
          xmax = qend,
          y = y_position,
          group = label,
          fill = description
        )
      )
  }
  
  if (nrow(pdb_df_best_score) > 0) {
    pdb_df_best_score$plot_annotation <- paste(pdb_df_best_score$target, " (", pdb_df_best_score$description, ")", sep = "")
    composite <- composite + geom_facet(
      panel = "PDB30 acc. + descr.",
      data = pdb_df_best_score,
      geom = geom_text2,
      hjust = 0,
      size = 2,
      inherit.aes = FALSE,
      stat = "identity",
      show.legend = FALSE,
      mapping = aes(label=plot_annotation, y= y_position, x = 0)) +
      xlim_expand(c(0,0.05), panel = "PDB30 acc. + descr.")
  }
  
  #plot cog_df grey protein backbones
  composite <- composite +
    geom_facet(
      panel = "COGs",
      data = data1,
      geom = geom_gene_arrow,
      alpha = 0.8,
      fill = "grey",
      arrowhead_width = grid::unit(1, "mm"),
      arrowhead_height = grid::unit(arrow_height, "mm"),
      arrow_body_height = grid::unit(arrow_height, "mm"),
      inherit.aes = FALSE,
      stat = "identity",
      show.legend = FALSE,
      mapping = aes(
        xmin = start,
        xmax = end,
        y = y_position,
        group = label
      )
    )
  
  if (nrow(cog_df) > 0) {
    
    #create grouped version of cogs df
    cog_df_grouped <- cog_df %>%
      group_by(label) %>%
      summarise_all(~paste(., collapse = ","))
    
    #rename the associated gene column. The plot version will be used for the geom_genearrows
    cog_df <- dplyr::rename(cog_df, associated_gene_plot = associated_gene)
    
    #merge with the grouped df where associated_gene column contains all associated genes separated by comma
    cog_df <- left_join(cog_df, dplyr::select(cog_df_grouped, label, associated_gene), by = 'label')
    
    cog_df$protein_id <- sapply(strsplit(cog_df$label, split = "__"), '[', 2)
    
    #cog y positions
    if (nrow(cog_df) > 0) {
      cog_df <- leaf_nodes %>%
        left_join(cog_df, by = "label")
      names(cog_df)[names(cog_df) == "y"] <- "y_position"
    }
    
    composite <- composite +
      geom_facet(
        panel = "COGs",
        data = cog_df,
        geom = geom_gene_arrow,
        alpha = 0.8,
        arrowhead_width = grid::unit(1, "mm"),
        arrowhead_height = grid::unit(arrow_height, "mm"),
        arrow_body_height = grid::unit(arrow_height, "mm"),
        inherit.aes = FALSE,
        stat = "identity",
        show.legend = FALSE,
        mapping = aes(
          xmin = qstart,
          xmax = qend,
          y = y_position,
          group = label,
          fill = associated_gene_plot
        )
      )
    
    
    composite <- composite + geom_facet(
      panel = "COGs acc. + descr.",
      data = cog_df,
      geom = geom_text2,
      hjust = 0,
      size = 2,
      inherit.aes = FALSE,
      stat = "identity",
      show.legend = FALSE,
      mapping = aes(label=associated_gene, y= y_position, x = 0)) +
      xlim_expand(c(0,0.05), panel = "COGs acc. + descr.")
  }
  
  
  composite <- composite + 
    theme(
      legend.position = "right",
      axis.line.x = element_line(color = "black"),
      axis.text.x = element_text(color = "black"),
      axis.ticks.x = element_line(color = "black")
    )
  
  
  filename <- paste0("data/",project,"/plots/effectors/", sub(".txt", ".png", basename(tree_files[i])),"_domains.png")
  ggsave(filename, composite, width = 15, height = 10)
  
  # create patchwork every 3 plots or when it is the last iteration
  if (i %% 3 == 0 || i == trees_length) {
    # handle less than 3 plots in the last row
    if (i - patch_index < 2) {
      plot_row <- patchwork_list[[patch_index]]
      for (n in (patch_index + 1):i) {
        plot_row <- plot_row | patchwork_list[[n]]
      }
    } else {
      plot_row <- patchwork_list[[patch_index]] | patchwork_list[[patch_index + 1]] | patchwork_list[[patch_index + 2]]
    }
    
    # check if it's first row, if so, assign it to plot_group, else, add new row to plot_group
    if (patch_index == 1) {
      plot_group <- plot_row
    } else {
      plot_group <- plot_group / plot_row
    }
    
    # increment j by 3 for next row
    patch_index <- patch_index + 3
  }
  
  
}

# add annotation to patchwork
plot_group <- plot_group + plot_annotation(title = 'Known effectors',
                                           caption = '26.9.2023')

# or save it to file
filename <- paste0("data/",project,"/plots/effectors/all.png")
ggsave(filename, plot = plot_group, height = 17, width = 50)





#### CorA neighbourhood ####

#create color palette from colorblind friendly palette
colorBlindBlack8_cora  <- c("#ffffff", "#E69F00", "#56B4E9", "#009E73", 
                            "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

cora_tree_path <- paste("data/",project,"/cora_neighbourhood/CorA_tree.txt", sep = "")
cora_tree <- read.newick(cora_tree_path) #cora tree contains 52 nodes

#cora_neighbours contain only CorA positive samples (54 total)
cora_neighbours <- read.csv(paste("data/",project,"/cora_neighbourhood/cora_neighbours.tsv", sep = ""), sep = "\t")

#cora_tree <- midpoint.root(cora_tree)

cora_neighbours$dedd <- as.logical(cora_neighbours$dedd)
cora_neighbours$dedd_clost <- as.logical(cora_neighbours$dedd_clost)
cora_neighbours$dedd_neo <- as.logical(cora_neighbours$dedd_neo)
cora_neighbours$nrn <- as.logical(cora_neighbours$nrn)
cora_neighbours$samlyase <- as.logical(cora_neighbours$samlyase)
#create new column dedd_master and set it to True if any of the dedd columns is True
cora_neighbours$dedd <- ifelse(cora_neighbours$dedd | cora_neighbours$dedd_clost | cora_neighbours$dedd_neo, TRUE, FALSE)

rownames(cora_neighbours) <- cora_neighbours$locus #set rownames

#merge cora_mastertable with info table

#cora positives from info have 53
cora_positives_info <- subset(info, cora == TRUE)
#rename Locus column to locus
cora_positives_info <- dplyr::rename(cora_positives_info, locus = Locus)

cora_mastertable <- merge(cora_neighbours, cora_positives_info, by = "locus")
rownames(cora_mastertable) <- cora_mastertable$locus #set rownames

#load fasta
#remove protein ids from headers
cora_fasta_path <- paste("data/",project,"/cora_neighbourhood/CorA_alignment.afa", sep = "")
cora_fasta <- readAAStringSet(cora_fasta_path)
names(cora_fasta) <- sub("__.*", "", names(cora_fasta))
#save again as fasta without the headers
cora_fasta_fixed_path = paste("data/",project,"/cora_neighbourhood/cora_fixed.afa", sep = "")
writeXStringSet(cora_fasta, filepath = cora_fasta_fixed_path)

#remove protein id from tip labels to only include locus
cora_tree[["tip.label"]] <- sub("__.*", "", cora_tree[["tip.label"]])

cora_associations <- c("dedd", "nrn", "samlyase", "can3", "csx1", "can1.2")
custom_labels_cora <- c("DEDD", "NrN", "SAM-LYASE", "Can3", "Csx1", "Can1-2")

cora_associations <- c("dedd", "nrn", "samlyase", "can3")
custom_labels_cora <- c("DEDD", "NrN", "SAM-LYASE", "Can3")

custom_true_colors <- c("dedd" = "#E69F00", "nrn" = "#56B4E9", "samlyase" = "#009E73", "can3" = "#F0E442", "csx1" = "#0072B2", "can1.2" = "#D55E00")


# Get palette of 6 blue colors from light to dark
blue_palette <- colorRampPalette(brewer.pal(9, "Blues"))
# Create a vector of 6 blue colors
custom_blue_colors <- blue_palette(9)[3:8]
names(custom_blue_colors) <- c("dedd", "nrn", "samlyase", "can3", "csx1", "can1.2")

#change can3 to color to #ff743e
custom_blue_colors["can3"] <- "#ff743e"

names(custom_labels_cora) <- cora_associations
offsets_cora <- c(0, 0.1, 0.2, 0.4, 0.7, 0.9)
lane_width_cora <- 0.03

ca4_color <- "#E69F00"
ca6_color <- "#ff743e"
sampamp_color <- "#009E73"
#make a light red color
sampamp_color <- "#FF9999"

clade_alpha <- 0.4

pCora <- ggtree(cora_tree, layout="rectangular") +
  theme(legend.position = "none") +
  #add bottom margin to lift higher
  theme(plot.margin = margin(0, 0, 1, 0, "cm")) +
  geom_hilight(node=81, fill=ca6_color, alpha = clade_alpha) #note that node=81 does not apply to new runs do to stochasticity of tree generation

pCora <- pCora %<+% cora_mastertable

for(i in 1:length(cora_associations)){
  i
  pCora <- pCora + new_scale_fill()
  pCora <- gheatmap(pCora, cora_mastertable[c(cora_associations[i])], offset=offsets_cora[i]+1.9, width=lane_width_cora, colnames = TRUE,
                    colnames_angle=90, colnames_offset_y = 0, colnames_offset_x = 0, font.size = label_fontsize, color = "grey",
                    hjust = 1, custom_column_labels = custom_labels_cora[[cora_associations[i]]]) +
    scale_fill_manual(values = c("white", custom_blue_colors[[cora_associations[i]]]), name=c(cora_associations[i]), guide = "none")
}

pCora <- pCora + new_scale_fill()
pCora <- gheatmap(pCora, cora_mastertable[c("Subtype")], width = lane_width_cora, offset = 2.5, colnames = TRUE, colnames_angle=90, colnames_offset_y = 0, font.size = label_fontsize, color = "grey", hjust = 1) +
  scale_fill_manual(values = colorBlindBlack8, name="CRISPR\nsubtype")

pCora <- pCora + ggtree::vexpand(.03, -1)

#from cora_mastertable, remove spaces in front of column "species" values
cora_mastertable$species <- gsub("^\\s+", "", cora_mastertable$species)
#remove [ from species names
cora_mastertable$species <- gsub("]", "", cora_mastertable$species)

pCora2 <- pCora + geom_tiplab(aes(label=species), hjust = -0.05, size = 2.2, align = TRUE, alpha = 0.6, fontface = 3) + coord_cartesian(clip="off")
pCora2 <- pCora2 +
  geom_point2(aes(subset = (node == 38)), color = "red", size = 2, shape = 16, alpha = 1) + #note that node=38 does not apply to new runs do to stochasticity of tree generation
  geom_point2(aes(subset = (node == 39)), color = "red", size = 2, shape = 16, alpha = 1) #note that node=39 does not apply to new runs do to stochasticity of tree generation

#pCora2 <- pCora + geom_tiplab(aes(label=label), hjust = -0.5, size = 3) #accession numbers

plot_path = paste("data/",project,"/plots/cora_neighbourhood_species.png", sep="")
ggsave(plot_path, pCora2, height = 6, width = 6)

#### Potential novel effectors ("group4") plotting ####

group4_counts = paste0("data/",project,"/group4/group4_prober_analysis.txt")
group4_counts <- read.csv(group4_counts, sep=",") #this will be unchanged data

group4_counts <- dplyr::rename(group4_counts, Query = X)
group4_counts <- dplyr::rename(group4_counts, Count = qseqid)

#read effector analysis for carfsaved
group4_carfsaved <- paste0("data/",project,"/group4/carfsaved/group4.CARFSAVED.info.tsv")
group4_carfsaved_raw <- read.csv(group4_carfsaved, sep="\t")

carfsaved_protein_data <- group4_carfsaved_raw %>% 
  group_by(Query) %>% 
  reframe(Start = 0, End = Query_Length-1, Length = Query_Length) %>%
  distinct(Query, .keep_all = TRUE)

# Create artificial y coordinates
#group4_carfsaved_raw$Y = as.numeric(as.factor(group4_carfsaved_raw$Query))
#carfsaved_protein_data$Y = as.numeric(as.factor(carfsaved_protein_data$Query))

#retain only best scores from hits
group4_carfsaved_raw_best <- subset(group4_carfsaved_raw, Target_Name != "protein") %>%
  group_by(Query) %>%
  arrange(desc(Domain_Score)) %>%
  dplyr::slice(1)

#group4_carfsaved_raw_pruned <- rbind(group4_carfsaved_raw_best, subset(group4_carfsaved_raw, Target_Name == "protein"))

#read effector analysis for pfam
group4_pfam <- paste0("data/",project,"/group4/pfam/group4.pfam.info.tsv")
group4_pfam_raw <- read.csv(group4_pfam, sep="\t") #this will be unchanged data

pfam_protein_data <- group4_pfam_raw %>% 
  group_by(Query) %>% 
  reframe(Start = 0, End = Query_Length-1, Length = Query_Length) %>%
  distinct(Query, .keep_all = TRUE)

# Create artificial y coordinates
#group4_pfam_raw$Y = as.numeric(as.factor(group4_pfam_raw$Query))
#pfam_protein_data$Y = as.numeric(as.factor(pfam_protein_data$Query))

#retain only best scores from hits
group4_pfam_raw_best <- subset(group4_pfam_raw, Target_Name != "protein") %>%
  group_by(Query) %>%
  arrange(desc(Domain_Score)) %>%
  dplyr::slice(1)

#Insert counts and Y values to a df that only contains the protein entries
proteins_with_counts_and_Y <- subset(group4_pfam_raw, Target_Name == "protein")
proteins_with_counts_and_Y <- merge(x = proteins_with_counts_and_Y, y = group4_counts[,c("Query","Count")], by = "Query", all.x = TRUE)
proteins_with_counts_and_Y$Y <- rank(proteins_with_counts_and_Y$Count, ties.method = "first")

group4_pfam_raw_pruned <- rbind(group4_pfam_raw_best, subset(group4_pfam_raw, Target_Name == "protein"))


#merge with count data
#group4_carfsaved_raw <- merge(x = group4_carfsaved_raw, y = proteins_with_counts_and_Y[,c("Query","Count","Y")], by = "Query", all.x = TRUE)
#carfsaved_protein_data <- merge(x = carfsaved_protein_data, y = proteins_with_counts_and_Y[,c("Query","Count","Y")], by = "Query", all.x = TRUE)

#remove entries mentioning saved and carf from the pfams to reduce redundancy
group4_pfam_raw_pruned <- group4_pfam_raw_pruned %>%
  filter(!str_detect(Target_Name, "SAVED"))
group4_pfam_raw_pruned <- group4_pfam_raw_pruned %>%
  filter(!str_detect(Target_Name, "CARF"))

#combine pfam and carfsaved
# we get protein entries and Y values from the pfam df
bound_pfam_carfsaved <- rbind(group4_carfsaved_raw_best, group4_pfam_raw_pruned)
bound_pfam_carfsaved <- merge(x = bound_pfam_carfsaved, y = proteins_with_counts_and_Y[,c("Query","Count","Y")], by = "Query", all.x = TRUE)

#plot definitions
label_size = 2
arrow_height = 2

protein_counts <- ggplot() +
  geom_segment( # grey protein backbones
    data = subset(bound_pfam_carfsaved, Target_Name == "protein"),
    aes(x = 0, xend = Count+5, y = Y, yend = Y),
    size = 2,
    colour = "lightblue"
  ) +
  geom_label_repel(
    data = subset(bound_pfam_carfsaved, Target_Name == "protein"),
    aes(x = Count, y = Y, label = Count),
    force = 0,
    size = 2,
    label.padding = 0.1
  ) +
  theme_minimal() +
  scale_y_continuous(breaks = subset(bound_pfam_carfsaved, Target_Name == "protein")$Y, labels = subset(bound_pfam_carfsaved, Target_Name == "protein")$Query) +
  labs(y = "Protein ID") +
  ggtitle("Count") +
  scale_fill_discrete(guide = FALSE) +
  theme(legend.position = "none") +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()) +
  scale_x_continuous(breaks = c(0, max(subset(bound_pfam_carfsaved, Target_Name == "protein")$Count)))


carfsaved_group4 <- ggplot() +
  geom_gene_arrow( # grey protein backbones
    data = subset(bound_pfam_carfsaved, Target_Name == "protein"),
    aes(xmin = 0, xmax = Query_Length -1, y = Y, forward = TRUE),
    fill = "grey",
    arrowhead_width = grid::unit(1, "mm"),
    arrowhead_height = grid::unit(arrow_height, "mm"),
    arrow_body_height = grid::unit(arrow_height, "mm"),
  ) +
  geom_gene_arrow(
    data = subset(bound_pfam_carfsaved, Target_Name != "protein"),
    aes(xmin = Start, xmax = End, y = Y, forward = TRUE, fill = as.factor(Query)),
    arrowhead_width = grid::unit(1, "mm"),
    arrowhead_height = grid::unit(arrow_height, "mm"),
    arrow_body_height = grid::unit(arrow_height, "mm"),
  ) +
  geom_label_repel(
    data = subset(bound_pfam_carfsaved, Target_Name != "protein"),
    aes(x = (Start + End) / 2, y = Y, label = Target_Name),
    min.segment.length = unit(0, 'lines'), #forces lines
    nudge_y = .2,
    force = 10,
    alpha = 1,
    size = label_size,
    xlim = c(-10, 1000),
    box.padding = 0.5,
    label.padding = 0.1,
    max.overlaps = 30
  ) +
  theme_minimal() +
  #scale_y_continuous(breaks = subset(bound_pfam_carfsaved, Target_Name == "protein")$Y, labels = subset(bound_pfam_carfsaved, Target_Name == "protein")$Query) +
  labs(x = "Protein position", y = "Protein ID") +
  scale_fill_discrete(guide = FALSE) +
  theme(legend.position = "none") +
  ggtitle("CARF/SAVED + Pfam") +
  ylab(NULL) +
  theme(axis.text.y = element_blank(),               # removes Y-axis numbers
        axis.ticks.y = element_blank()) +
  theme(plot.margin = margin(0, 1, 0, -0.5, "cm")) # adjust left margin to go closer to the left plot
#carfsaved_group4

#same for pdb and cogs
#read effector analysis for carfsaved
group4_pdb <- paste0("data/",project,"/group4/pdb/group4_pdb_parsed.tsv")
group4_pdb_raw <- read.csv(group4_pdb, sep="\t")

pdb_protein_data <- group4_pdb_raw %>% 
  group_by(query) %>% 
  distinct(query, .keep_all = TRUE)

#retain only best scores from hits
group4_pdb_raw_best <- subset(group4_pdb_raw, target != "protein") %>%
  group_by(query) %>%
  arrange(desc(score)) %>%
  dplyr::slice(1)

#read effector analysis for cogs
group4_cogs <- paste0("data/",project,"/group4/cog/group4_cog_parsed.tsv")
group4_cogs_raw <- read.csv(group4_cogs, sep="\t") #this will be unchanged data

cogs_protein_data <- group4_cogs_raw %>% 
  group_by(query) %>% 
  distinct(query, .keep_all = TRUE)

#retain only best scores from hits
group4_cogs_raw_best <- subset(group4_cogs_raw, target != "protein") %>%
  group_by(query) %>%
  arrange(desc(score)) %>%
  dplyr::slice(1)

group4_cogs_raw_best <- group4_cogs_raw_best %>%
  mutate(associated_gene = if_else(associated_gene=="", description, associated_gene))

#Insert counts and Y values to a df that only contains the protein entries

#group4_cogs_raw_pruned <- rbind(group4_cogs_raw_best, subset(group4_pfam_raw, Target_Name == "protein"))


#merge with count data
#group4_carfsaved_raw <- merge(x = group4_carfsaved_raw, y = proteins_with_counts_and_Y[,c("Query","Count","Y")], by = "Query", all.x = TRUE)
#carfsaved_protein_data <- merge(x = carfsaved_protein_data, y = proteins_with_counts_and_Y[,c("Query","Count","Y")], by = "Query", all.x = TRUE)

proteins_with_counts_and_Y_hhsuite <- dplyr::rename(proteins_with_counts_and_Y, query = Query)
proteins_with_counts_and_Y_hhsuite <- dplyr::rename(proteins_with_counts_and_Y_hhsuite, count = Count)

#combine pfam and carfsaved
# we get protein entries and Y values from the pfam df
bound_pdb_cogs <- rbind(group4_cogs_raw_best, group4_pdb_raw_best)
bound_pdb_cogs <- rbind(bound_pdb_cogs, group4_pdb_raw_best)
bound_pdb_cogs <- merge(x = bound_pdb_cogs, y = proteins_with_counts_and_Y_hhsuite[,c("query","count","Y")], by = "query", all.x = TRUE)

pdb_plottable <- merge(x = group4_pdb_raw_best, y = proteins_with_counts_and_Y_hhsuite[,c("query","count","Y")], by = "query", all.x = TRUE)
cogs_plottable <- merge(x = group4_cogs_raw_best, y = proteins_with_counts_and_Y_hhsuite[,c("query","count","Y")], by = "query", all.x = TRUE)

pdb_group4 <- ggplot() +
  geom_gene_arrow( # grey protein backbones. We can reuse the carfsaved protein backbones
    data = subset(bound_pfam_carfsaved, Target_Name == "protein"),
    aes(xmin = 0, xmax = Query_Length -1, y = Y, forward = TRUE),
    fill = "grey",
    arrowhead_width = grid::unit(1, "mm"),
    arrowhead_height = grid::unit(arrow_height, "mm"),
    arrow_body_height = grid::unit(arrow_height, "mm"),
  ) +
  geom_gene_arrow(
    data = subset(pdb_plottable, target != "protein"),
    aes(xmin = qstart, xmax = qend, y = Y, forward = TRUE, fill = as.factor(query)),
    arrowhead_width = grid::unit(1, "mm"),
    arrowhead_height = grid::unit(arrow_height, "mm"),
    arrow_body_height = grid::unit(arrow_height, "mm"),
  ) +
  # geom_label_repel(
  #   data = subset(pdb_plottable, target != "protein"),
  #   aes(x = (qstart + qend) / 2, y = Y, label = description),
  #   min.segment.length = unit(0, 'lines'),
  #   #nudge_y = .2,
  #   #nudge_x = 10,
  #   force = 1,
  #   alpha = 1,
  #   size = label_size,
  #   xlim = c(-10, 1000),
  #   ylim = c(-3, 55),
#   box.padding = 1,
#   label.padding = 0.1,
#   max.overlaps = 30
# ) +
geom_label_repel(
  data = subset(pdb_plottable, target != "protein"),
  aes(x=1500, y=Y, label=description),
  size = 2,
  force = 0,
  alpha = 0.8,
  label.padding = 0.1) +
  theme_minimal() +
  scale_y_continuous(breaks = subset(bound_pfam_carfsaved, Target_Name == "protein")$Y, labels = subset(bound_pfam_carfsaved, Target_Name == "protein")$Query) +
  labs(x = "Protein position") +
  ggtitle("PDB") +
  scale_fill_discrete(guide = FALSE) +
  theme(legend.position = "none") +
  theme_minimal() +
  xlim(0,2000)
pdb_group4

cogs_group4 <- ggplot() +
  geom_gene_arrow( # grey protein backbones
    data = subset(bound_pfam_carfsaved, Target_Name == "protein"),
    aes(xmin = 0, xmax = Query_Length -1, y = Y, forward = TRUE),
    fill = "grey",
    arrowhead_width = grid::unit(1, "mm"),
    arrowhead_height = grid::unit(arrow_height, "mm"),
    arrow_body_height = grid::unit(arrow_height, "mm"),
  ) +
  geom_gene_arrow(
    data = subset(cogs_plottable, target != "protein"),
    aes(xmin = qstart, xmax = qend, y = Y, forward = TRUE, fill = as.factor(query)),
    arrowhead_width = grid::unit(1, "mm"),
    arrowhead_height = grid::unit(arrow_height, "mm"),
    arrow_body_height = grid::unit(arrow_height, "mm")
  ) +
  geom_label_repel(
    data = subset(bound_pdb_cogs, target != "protein"),
    aes(x = (qstart + qend) / 2, y = Y, label = associated_gene),
    min.segment.length = unit(0, 'lines'),
    nudge_y = .2,
    force = 10,
    alpha = 1,
    size = label_size,
    xlim = c(-10, 1000),
    box.padding = 0.5,
    label.padding = 0.1,
    max.overlaps = 30
  ) +
  theme_minimal() +
  ggtitle("COGs") +
  scale_y_continuous(breaks = subset(bound_pfam_carfsaved, Target_Name == "protein")$Y, labels = subset(bound_pfam_carfsaved, Target_Name == "protein")$Query) +
  labs(x = "Protein position", y = "Protein ID") +
  scale_fill_discrete(guide = FALSE) +
  theme(legend.position = "none") +
  theme_minimal()
#cogs_group4

pdb_ids <- ggplot() +
  geom_label_repel(
    data = pdb_plottable,
    aes(x=100, y=Y, label=target),
    size = 2,
    force = 0,
    alpha = 0.8,
    label.padding = 0.1) +
  theme_minimal() +
  scale_y_continuous(breaks = subset(bound_pfam_carfsaved, Target_Name == "protein")$Y, labels = subset(bound_pfam_carfsaved, Target_Name == "protein")$Query) +
  theme_minimal() +
  scale_y_continuous(breaks = subset(bound_pfam_carfsaved, Target_Name == "protein")$Y, labels = subset(bound_pfam_carfsaved, Target_Name == "protein")$Query) +
  labs(x = "Protein position") +
  ggtitle("PDB ID") +
  scale_fill_discrete(guide = FALSE) +
  theme(legend.position = "none") +
  theme_minimal() + 
  theme(axis.text.y = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),  # hides y-axis title
        axis.text.x = element_blank(),   # hides y-axis numbers
        axis.ticks.x = element_blank()   # hides y-axis ticks
  )

pdb_group4 <- pdb_group4 + theme(axis.text.y = element_blank(), 
                                 axis.title.y = element_blank())

cogs_group4 <- cogs_group4 + theme(axis.text.y = element_blank(), 
                                   axis.title.y = element_blank())

group4_collage <- {protein_counts | carfsaved_group4 | cogs_group4 | pdb_group4 | pdb_ids} + plot_layout(widths = c(0.1, 1, 1, 1, 0.2))

plot_path <- paste("data/",project,"/plots/group4.png", sep = "")
ggplot2::ggsave(plot_path, plot = group4_collage, width = 14, height = 8)


#### Known effectors ####
known_temp <- ggplot(data = effectors, aes(x = name, y = mean_temperature)) +
  geom_col() +
  geom_errorbar(aes(ymin = mean_temperature - std_temperature, 
                    ymax = mean_temperature + std_temperature), width = 0.2) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

known_GGDD_prop <- ggplot(data = effectors, aes(x = name, y = GGDD_prop)) +
  geom_col() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

known_counts <- ggplot(data = effectors, aes(x = name, y = count)) +
  geom_col()+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

known_genus <- ggplot(data = effectors, aes(x = name, y = most_common_genus_prop)) +
  geom_col() +
  geom_text(aes(label = most_common_genus), vjust = -0.5) +
  ylim(0,1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

known_bacteriaprop <- ggplot(data = effectors, aes(x = name, y = bacteria_prop)) +
  geom_col() +
  ylim(0,1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

known_collage <- (known_temp | known_GGDD_prop | known_counts) / (known_genus | known_bacteriaprop) + plot_annotation(tag_levels = 'A',
                                                                                                                      title = 'Known effectors',
                                                                                                                      caption = '19.9.2023')
plot_path <- paste("data/",project,"/plots/",project,"_known_effectors.png", sep = "")
ggsave(plot_path, plot = known_collage, width = 12, height = 10)

#### Signal molecules and temperature ####
cOA_long_info <- info %>%
  pivot_longer(cols = c(ca3, ca4, ca5, ca6, sam.amp),
               names_to = "ca",
               values_to = "value")
cOA_long_info$value <- as.logical(cOA_long_info$value)

cOA_long_info <- cOA_long_info %>%
  filter(value)

cOA_long_info$Cas10_length_norm <- (cOA_long_info$Cas10_length - min(cOA_long_info$Cas10_length)) / (max(cOA_long_info$Cas10_length) - min(cOA_long_info$Cas10_length))

p1 <- ggplot(cOA_long_info, aes(x = ca, y = mean_temperature)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(height = 0, width = 0.3, aes(colour = domain, size = (Cas10_length_norm)), alpha = 0.7) +
  scale_color_manual(values = colorBlindBlack8)


p2 <- ggplot(cOA_long_info, aes(x = ca, y = mean_temperature)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(height = 0, width = 0.3, aes(colour = HD_hmm_boolean, size = (Cas10_length_norm)), alpha = 0.7) +
  scale_color_manual(values = colorBlindBlack8)


plot_path <- paste("data/",project,"/plots/",project,"_cOA_temps_domain.png", sep = "")
ggsave(plot_path, plot = p1, width = 6, height = 6)

plot_path <- paste("data/",project,"/plots/",project,"_cOA_temps_HD.png", sep = "")
ggsave(plot_path, plot = p2, width = 6, height = 6)


#### Loop based Cas10 tree ####
# Create associated molecules and colors
associated_effectors <- c("nucc", "csx1", "can1.2", "cami1", "cam1", "calpl", "csm6.ca6", "cora")
associated_colors <- c("#a8d2e1", "#7fb1dc", "#5790d4", "#2d74bb", "#1859a4", "#0c4d8d", "#003166", "#9ee493", "#ff743e")
associated_colors <- c("#a8d2e1", "#6fb0dc", "#388dd1", "#1561a3", "#123675", "#091848", "#9ee493", "#ff743e")
#use colorbrewer to make to shades of pleasant blue
HD_colors <- c("#6fb0dc", "#123675")
label_fontsize <- 2
offsets <- c(1, 2.5, 3.5, 4.5, 5.5, 6.5, 8, 9.5)
# Combine into a named list
molecule_color_list <- setNames(associated_colors, associated_effectors) 
lane_width = 0.03

rownames(info) <- info$Locus
cas10_tree$tip.label <- as.character(cas10_tree$tip.label)

#make the Locus colun the first column in info
info <- info[, c('Locus', setdiff(names(info), 'Locus'))]

pX <- ggtree(cas10_tree, layout="circular", branch.length = "none", size=0.15, open.angle=20, aes(color=Subtype, labels=genus) +
               theme(legend.position = "none"))  %<+% info
#geom_point2(aes(subset = (domain == "Archaea")), color = "red", size = 0.3) +


pX <- gheatmap(pX, info[c("Subtype")], offset=-1.3, width=0.06, colnames = TRUE, colnames_angle=85, colnames_offset_y = 0, font.size = label_fontsize, color = "grey", hjust = 1, legend_title = "Subtype") + 
  scale_fill_manual(values = colorBlindBlack8, name = "Subtype") +
  theme(text = element_text(size = 14)) # + geom_point2(aes(subset = isTip == TRUE & domain == "Archaea")), color = "red", size = 0.1 #uncomment to label archaea


custom_labels <- c("NucC", "Csx1", "Can1-2", "Cami1", "Cam1", "CalpL", "Csm6", "CorA")
names(custom_labels) <- associated_effectors

# Iterate over each effector
for(i in 1:length(associated_effectors)){
  pX <- pX + new_scale_fill()
  pX <- gheatmap(pX, info[c(associated_effectors[i])], offset=offsets[i] + 0.3, width=lane_width, colnames = TRUE,
                 colnames_angle=85, colnames_offset_y = 0, colnames_offset_x = 0, font.size = label_fontsize, color = "grey",
                 hjust = 1, custom_column_labels = custom_labels[[associated_effectors[i]]]) +
    scale_fill_manual(values=c("white", associated_colors[i]), guide = "none") + 
    theme(legend.position = 'none')
}

pX <- pX + new_scale_fill() + theme(legend.position = 'none') +
  geom_fruit(
    data = info_for_cas10length, #if you don't inject the data again here, will get an error: "Missing value where true/false needed". Remember to rename columns so they don't match existing tree data
    geom=geom_col,
    offset=0.36, #use 0.05 when only CorA ring. When all rings, 0.18
    alpha = 0.7,
    mapping=aes(y=Locus, x=Cas10_length_plot, fill=Cas10_HD_plot), #Cas10 HD fill based on literal sequence, not HMM
    pwidth=0.5,
    orientation="y",
    axis.params=list(
      axis="x", # add axis text of the layer.
      text.angle=-90, # the text size of axis.
      hjust=0 , # adjust the horizontal position of text of axis.
      text.size = 1.5
    ),
    grid.params=list() # add the grid line of the external bar plot.
  ) + 
  scale_fill_manual(values = HD_colors, name = "HD domain") +
  theme(legend.position=c(0.95, 0.5),
        legend.title=element_text(size=7),
        legend.text=element_text(size=6),
        legend.spacing.y = unit(0.02, "cm"))

pX_opened <- open_tree(pX, 15)

#save pX
plot_path <- paste("data/",project,"/plots/",project,"_Cas10_effectors_loop.png", sep = "")
ggsave(plot_path, plot = pX_opened, width = 10, height = 10)

#### Cas10 cOA ring, looped ####
# Create associated molecules and colors
signals <- c("ca3", "ca4", "ca6", "sam.amp")
associated_colors_signals <- c("turquoise4", "lightcoral", "forestgreen", "dodgerblue4")

pX <- ggtree(cas10_tree, layout="circular", branch.length = "none", size=0.15, open.angle=20, aes(color=Subtype, labels=genus) +
               theme(legend.position = "none") +
               geom_text2(aes(subset=!isTip, label=node), hjust = 2 ))
pX <- pX %<+% info

# Iterate over each effector
for(i in 1:length(signals)){
  pX <- pX + new_scale_fill()
  pX <- gheatmap(pX, info[c(signals[i])], offset=i-1, width=lane_width, colnames = TRUE,
                 colnames_angle=85, colnames_offset_y = 0, font.size = label_fontsize, color = "grey",
                 hjust = 1, legend_title = signals[i]) +
    scale_fill_manual(values=c("white", associated_colors_signals[i]), name=signals[i])
}

pX <- pX + new_scale_fill() +
  geom_fruit(
    data = info_for_cas10length, #if you don't inject the data again here, will get an error: "Missing value where true/false needed". Remember to rename columns so they don't match existing tree data
    geom=geom_col,
    offset=0.13, #use 0.05 when only CorA ring. When all rings, 0.18
    alpha = 0.7,
    mapping=aes(y=Locus, x=Cas10_length_plot, fill=Cas10_HD_plot),
    pwidth=0.5,
    orientation="y",
    axis.params=list(
      axis="x", # add axis text of the layer.
      text.angle=-90, # the text size of axis.
      hjust=0,  # adjust the horizontal position of text of axis.
      text.size = 1.5
    ),
    grid.params=list() # add the grid line of the external bar plot.
  ) + 
  scale_fill_manual(values = colorBlindBlack8_2, name = "HD motif") +
  theme(legend.position=c(0.95, 0.5),
        legend.title=element_text(size=7),
        legend.text=element_text(size=6),
        legend.spacing.y = unit(0.02, "cm"))

pX_opened <- open_tree(pX, 15)

#save pX
plot_path <- paste("data/",project,"/plots/",project,"_Cas10_effectors_loop_signals.png", sep = "")
ggsave(plot_path, plot = pX_opened, width = 10, height = 10)


#### Upset plot of cooccurrence between different effectors ####
#When subsetting for all, archaea or thermophiles, change the following:
# - The subsetting function at the top (first thing in this section)
# - Signal molecule associations at the top of this section
# - The name of the plot file at the end of this section
# - Title of the plot in the plot script


upset_info <- info
#upset_info <- subset(info, domain == "Archaea") #uncomment this for running only for archaea
#upset_info <- subset(info, mean_temperature >= 60) #uncomment this for running only for thermophiles)

info_multiple_loci_effectors <- upset_info %>%
  pivot_longer(cols = c(can3, tirsaved, cam1, cam2, cam3, saved.chat, tirsaved, nucc, calpl, cami1, can1.2, csx1, csx23, csm6.ca6, cora), 
               names_to = "Effector", 
               values_to = "Presence")

#presence as factor
info_multiple_loci_effectors$Presence <- as.factor(info_multiple_loci_effectors$Presence)

# Sort out duplicates and spread the data
info_multiple_loci_effectors_upset <- info_multiple_loci_effectors %>%
  filter(Presence == TRUE) %>%  # Select only rows where effector is present
  select(Locus, Effector) %>%    # Select necessary columns
  mutate(Presence = 1) %>%      # Convert presence to 1
  group_by(Locus, Effector) %>% 
  summarise(Presence = max(Presence)) %>% # Get maximum presence value for each combination
  pivot_wider(names_from = Effector, values_from = Presence, values_fill = list(Presence = 0)) # Pivot data

#create column total in info_multiple_loci_effectors_upset
info_multiple_loci_effectors_upset$total <- rowSums(info_multiple_loci_effectors_upset[,2:ncol(info_multiple_loci_effectors_upset)])

#backup dataframe with all occurrences of effectors
info_multiple_loci_beforeDuplicateRemoval <- info_multiple_loci_effectors_upset

#remove total column
info_multiple_loci_effectors_upset <- subset(info_multiple_loci_effectors_upset, select = -c(total))
info_multiple_loci_beforeDuplicateRemoval <- subset(info_multiple_loci_beforeDuplicateRemoval, select = -c(total))

#convert tibble to df
info_multiple_loci_effectors_upset <- as.data.frame(info_multiple_loci_effectors_upset)
info_multiple_loci_beforeDuplicateRemoval <- as.data.frame(info_multiple_loci_beforeDuplicateRemoval)

# Convert your dataframe to list
info_multiple_loci_effectors_upset_list <- lapply(names(info_multiple_loci_effectors_upset)[-1], function(x){
  info_multiple_loci_effectors_upset$Locus[info_multiple_loci_effectors_upset[,x] == 1]})
info_multiple_loci_beforeDuplicateRemoval_list <- lapply(names(info_multiple_loci_beforeDuplicateRemoval)[-1], function(x){
  info_multiple_loci_beforeDuplicateRemoval$Locus[info_multiple_loci_beforeDuplicateRemoval[,x] == 1]})

# Make names of listInput as the names of Effectors
names(info_multiple_loci_effectors_upset_list) <- names(info_multiple_loci_effectors_upset)[-1]
names(info_multiple_loci_beforeDuplicateRemoval_list) <- names(info_multiple_loci_beforeDuplicateRemoval)[-1]

names(info_multiple_loci_effectors_upset)[names(info_multiple_loci_effectors_upset) == 'cam1'] <- 'Cam1'
names(info_multiple_loci_effectors_upset)[names(info_multiple_loci_effectors_upset) == 'cam2'] <- 'Cam2'
names(info_multiple_loci_effectors_upset)[names(info_multiple_loci_effectors_upset) == 'cam3'] <- 'Cam3'
names(info_multiple_loci_effectors_upset)[names(info_multiple_loci_effectors_upset) == 'can3'] <- 'Can3'
names(info_multiple_loci_effectors_upset)[names(info_multiple_loci_effectors_upset) == 'csm6.ca6'] <- 'Csm6'
names(info_multiple_loci_effectors_upset)[names(info_multiple_loci_effectors_upset) == 'saved.chat'] <- 'SAVED-CHAT'
names(info_multiple_loci_effectors_upset)[names(info_multiple_loci_effectors_upset) == 'tirsaved'] <- 'TIR-SAVED'
names(info_multiple_loci_effectors_upset)[names(info_multiple_loci_effectors_upset) == 'csx1'] <- 'Csx1'
names(info_multiple_loci_effectors_upset)[names(info_multiple_loci_effectors_upset) == 'csx23'] <- 'Csx23'
names(info_multiple_loci_effectors_upset)[names(info_multiple_loci_effectors_upset) == 'calpl'] <- 'CalpL'
names(info_multiple_loci_effectors_upset)[names(info_multiple_loci_effectors_upset) == 'cami1'] <- 'Cami1'
names(info_multiple_loci_effectors_upset)[names(info_multiple_loci_effectors_upset) == 'can1.2'] <- 'Can1-2'
names(info_multiple_loci_effectors_upset)[names(info_multiple_loci_effectors_upset) == 'cora'] <- 'CorA'
names(info_multiple_loci_effectors_upset)[names(info_multiple_loci_effectors_upset) == 'nucc'] <- 'NucC'

#same for info_multiple_loci_beforeDuplicateRemoval
names(info_multiple_loci_beforeDuplicateRemoval)[names(info_multiple_loci_beforeDuplicateRemoval) == 'cam1'] <- 'Cam1'
names(info_multiple_loci_beforeDuplicateRemoval)[names(info_multiple_loci_beforeDuplicateRemoval) == 'cam2'] <- 'Cam2'
names(info_multiple_loci_beforeDuplicateRemoval)[names(info_multiple_loci_beforeDuplicateRemoval) == 'cam3'] <- 'Cam3'
names(info_multiple_loci_beforeDuplicateRemoval)[names(info_multiple_loci_beforeDuplicateRemoval) == 'can3'] <- 'Can3'
names(info_multiple_loci_beforeDuplicateRemoval)[names(info_multiple_loci_beforeDuplicateRemoval) == 'csm6.ca6'] <- 'Csm6'
names(info_multiple_loci_beforeDuplicateRemoval)[names(info_multiple_loci_beforeDuplicateRemoval) == 'saved.chat'] <- 'SAVED-CHAT'
names(info_multiple_loci_beforeDuplicateRemoval)[names(info_multiple_loci_beforeDuplicateRemoval) == 'tirsaved'] <- 'TIR-SAVED'
names(info_multiple_loci_beforeDuplicateRemoval)[names(info_multiple_loci_beforeDuplicateRemoval) == 'csx1'] <- 'Csx1'
names(info_multiple_loci_beforeDuplicateRemoval)[names(info_multiple_loci_beforeDuplicateRemoval) == 'csx23'] <- 'Csx23'
names(info_multiple_loci_beforeDuplicateRemoval)[names(info_multiple_loci_beforeDuplicateRemoval) == 'calpl'] <- 'CalpL'
names(info_multiple_loci_beforeDuplicateRemoval)[names(info_multiple_loci_beforeDuplicateRemoval) == 'cami1'] <- 'Cami1'
names(info_multiple_loci_beforeDuplicateRemoval)[names(info_multiple_loci_beforeDuplicateRemoval) == 'can1.2'] <- 'Can1-2'
names(info_multiple_loci_beforeDuplicateRemoval)[names(info_multiple_loci_beforeDuplicateRemoval) == 'cora'] <- 'CorA'
names(info_multiple_loci_beforeDuplicateRemoval)[names(info_multiple_loci_beforeDuplicateRemoval) == 'nucc'] <- 'NucC'

#merge the info_multiple_loci_effectors_upset and info dataframes
upset_merged_info <- merge(info_multiple_loci_effectors_upset, info, by.x = "Locus", by.y = "Locus")
upsetAll_merged_info <- merge(info_multiple_loci_beforeDuplicateRemoval, info, by.x = "Locus", by.y = "Locus")

# Signal molecule associations
#All
effector_coa_data = data.frame(
  set=names(info_multiple_loci_effectors_upset)[-1],
  cOAs=c('cA4', 'cA4', 'cA4', 'cA4', 'SAM-AMP', 'cA6', 'cA4', 'cA3', 'cA3', 'cA4', 'cA4', 'cA6', 'cA3', 'cA4')
)

# Archaea 
#effector_coa_data = data.frame(
#  set=names(info_multiple_loci_effectors_upset)[-1],
#  cOAs=c('cA4', 'cA4', 'cA4', 'cA4', 'SAM-AMP', 'cA4', 'cA4')
#)

#Thermophiles
#effector_coa_data = data.frame(
#  set=names(info_multiple_loci_effectors_upset)[-1],
#  cOAs=c('cA4', 'cA4', 'cA4', 'cA4', 'cA4', 'cA4', 'SAM-AMP')
#)

png(paste("data/",project,"/plots/",project,"_upset_plot_effectors_complex.png", sep = ""), width = 3000, height = 3000, res = 400)
upset(upset_merged_info[-1],
      names(info_multiple_loci_effectors_upset)[-1],
      height_ratio=1,
      name = "Effector or combination",
      width_ratio = 0.15,
      min_degree=1,
      set_sizes=(
        upset_set_size(
          geom=geom_bar(
            stat = "count",
            position='fill',
            aes(fill=has_multiple_effectors, x=group),
            width=0.8
          ),
          position='right'
        ) +
          scale_fill_manual(values = c("#D6D6D6", "#B3B3B3")) +
          theme(legend.position = "none")
      ),
      base_annotations=list(
        'Count'=intersection_size(
          bar_number_threshold=1,
          text=list(vjust=1.1),
          counts=FALSE,
          mapping=aes(fill=Subtype)) +
          scale_fill_manual(values = colorBlindBlack8) +
          ggtitle('Effector co-occurrences and subtype distributions') +
          annotate(
            geom='text', x=Inf, y=Inf,
            color='#595959',
            size=3,
            label=paste('Total occurrences:', nrow(info_multiple_loci_effectors_upset)),
            vjust=1, hjust=1
          ) +
          geom_text(stat='count', aes(group=intersection, label=..count.., y=..count..), nudge_y = 10, size = 3) + #use nudge_y = 3 for archaea
          theme(axis.title.y = element_text(vjust = -10),
                plot.margin = margin(0, 0, 0, 0, "cm"))
        
      ),
      stripes=upset_stripes(
        mapping=aes(color=cOAs),
        colors=c(
          'cA3'='#d7fcd8',
          'cA4'='#fcf4e1',
          'cA6'='#e6eeff',
          'SAM-AMP' = '#dbdbdb'
        ),
        data=effector_coa_data
      ),
      guides='over'
) +
  xlab('Number of loci') +
  ylab("Co-occurrence\nproportion (dark)") +
  theme(plot.title = element_text(size = 8),
        axis.text.x = element_blank())

dev.off()

#### Sunburst plotly ####
#Create a sunburst graph showing the distribution of effectors so that the associated signal molecule is in the inner ring
#https://plotly.com/r/sunburst-charts/

info_sunburst <- info[c("nucc", "csx1", "can1.2", "cami1", "cam1", "csx23", "calpl", "csm6.ca6", "saved.chat", "cora")]

primary_df <- data.frame(
  signal_molecule = c("ca3", "ca4", "ca4", "ca4", "ca4", "ca4", "ca4", "ca6","ca3", "sam.amp"),
  effector = c("nucc", "csx1", "can1.2", "cami1", "cam1", "csx23", "calpl", "csm6.ca6", "saved.chat", "cora")
)

effector_count <- sapply(info_sunburst, function(x) sum(x == TRUE))

final_df <- merge(primary_df, effector_count, by.x = "effector", by.y = "row.names")

# Adding signal molecule to the effector for unique ids
final_df$id <- paste(final_df$effector, final_df$signal_molecule, sep="-")

# Adding signal molecule alone as id (adding these into the effector column)
updated_effector <- unique(c(final_df$effector, final_df$signal_molecule))

total_df <- data.frame(id = updated_effector, 
                       labels = updated_effector, 
                       parents = ifelse(updated_effector %in% final_df$signal_molecule, "", 
                                        final_df$signal_molecule),
                       values = ifelse(updated_effector %in% final_df$effector, 
                                       final_df$y, 0)
)

blue_shades <- c("#f7fbff", "#deebf7", "#c6dbef", "#9ecae1", "#6baed6", "#4292c6", "#2171b5", "#08519c", "#08306b", "#000000")
blue_shades <- c("#99BBDD", "#7799CC", "#5577BB", "#3366AA", "#115599", "#004488", "#003377", "#002266", "#001155", "#000044")

# calculate numbers of unique 'ca4' effectors
ca4_num <- length(unique(total_df$id[total_df$parents == "ca4"]))

# repeat the colors so the length of the color array equals to the number of unique 'ca4' effectors
ca4_colors <- rep(blue_shades, length.out = ca4_num)

ca3_colors <- rep(c("#ceecf5", "#a5c9d4"), length(unique(total_df$id[total_df$parents == "ca3"])))


# all others might be another color, say, grey
other_colors <- rep("grey", length(unique(total_df$id)) - ca4_num)

total_df$percentage <- total_df$values/sum(total_df$values)

#order total_df by percentage (decreasing)

total_df_only_ca4 <- total_df[total_df$parents == "ca4",]
total_df_only_ca4 <- total_df_only_ca4[order(total_df_only_ca4$percentage, decreasing = F),]
total_df_only_ca4$color <- ifelse(total_df_only_ca4$parents == "ca4", ca4_colors, other_colors)

#merge the ca4_only with total_df
total_df <- merge(total_df, total_df_only_ca4, by = c("id", "labels", "parents", "values", "percentage"), all.x = TRUE)

# color all except ca4 parented as grey
total_df$color[total_df$id != "ca4" & total_df$parents != "ca4"] <- "grey"

total_df$color[total_df$id == "ca3"] <- "#ceecf5"
total_df$color[total_df$id == "ca4"] <- "#336699" 
total_df$color[total_df$id == "ca6"] <- "#9ee493"
total_df$color[total_df$id == "sam.amp"] <- "#ff743e"

#make nucc like peachpuff but darker
total_df$color[total_df$id == "csm6.ca6"] <- "green4"
total_df$color[total_df$id == "cora"] <- "peachpuff3"

total_df$color[total_df$id == "saved.chat"] <- ca3_colors[1]
total_df$color[total_df$id == "nucc"] <- ca3_colors[2]


fig <- plot_ly(total_df, 
               ids = ~id,
               labels = ~labels,
               parents = ~parents,
               values = ~values,
               #display label and percentage
               #textinfo = "label+percent parent",
               #make the textinfo smaller
               textfont = list(size = 20),
               type = 'sunburst',
               marker = list(colors = ~color),
               textfont = list(size = 40)
)
fig

#draw the same sunburst plot but without labels
fig <- plot_ly(total_df, 
               ids = ~id, 
               labels = "",
               parents = ~parents, 
               values = ~values,
               type = 'sunburst', 
               marker = list(colors = ~color),
               textfont = list(size = 40)
)
fig

#save using orca if it works (uncomment installation scripts below if not installed). If not, just save the sunburst from the RStudio viewer panel -> export -> save as png
plot_path = paste("data/",project,"/plots/",project,"_sunburst.png", sep = "")
#reticulate::install_miniconda()
#reticulate::conda_install('r-reticulate', 'python-kaleido')
#reticulate::conda_install('r-reticulate', 'plotly', channel = 'plotly')
#reticulate::use_miniconda('r-reticulate')
#install kaleido too
#reticulate::py_install('kaleido')
save_image(fig, plot_path, width = 1000, height = 1000, scale = 1)



#### Statistics on crosstalk ####

# examining possible crosstalk in samples with multiple loci through inference

p2 <- ggplot(data = info, aes(x=has_effector, y = after_stat(count), fill = multilocus_sample)) +
  geom_bar() +
  ggtitle("Number of loci with signal molecule in sample") +
  ylab("Number of loci") +
  xlab("Has effector") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  #rename x axis labels
  scale_fill_manual(values = colorBlindBlack8)
p2

# model formulated in a generalised linear model
info$no_effector <- !info$has_effector
glm_no_effector <- glm(no_effector ~ multilocus_sample, data = info, family = "binomial")
summary(glm_no_effector)


#### Summary tables ####

#Using the info df, create some summarised stats. These have Cas10 length cutoff > 500

#First, total number of loci (i.e. number of rows in info df)
total_loci <- nrow(info)

#Number of loci with HD domain (info$HD_hmm_boolean)
HD_loci_hmm <- nrow(subset(info, HD_hmm_boolean == "True")) #hmm search for HD
HD_loci_literal <- nrow(subset(info, Cas10_HD == TRUE)) #literal HD 0-50 AA

#Percentage of loci with HD domain
HD_percentage <- (HD_loci_literal / total_loci) * 100

# List of subtypes
subtypes <- c("III-A", "III-B", "III-D", "III-F", "III-C")

# Empty list to store results
results <- list()

# Loop over subtypes
for (subtype in subtypes) {
  
  # Subset data
  subtype_data <- subset(info, Subtype == subtype)
  
  subtype_without_dash <- gsub("-", "", subtype)
  
  # Calculate count and percentages
  subtype_count <- nrow(subtype_data)
  HD_in_subtype <- nrow(subset(subtype_data, HD_hmm_boolean == "True"))
  HD_subtype_percentage <- (HD_in_subtype / subtype_count) * 100
  cyclase_subtype <- nrow(subset(subtype_data, cyclase == TRUE))
  cyclase_subtype_percentage <- (cyclase_subtype / subtype_count) * 100
  
  # Store the computed percentage values into respective variables
  assign(paste0("HD_subtype_", subtype_without_dash, "_percentage"), HD_subtype_percentage)
  assign(paste0("cyclase_", subtype_without_dash, "_percentage"), cyclase_subtype_percentage)
}

#Number of loci with cyclase domain
cyclase_loci <- nrow(subset(info, cyclase == TRUE)) #cyclase = either GGDD or GGDE HMM found AND literal cyclase is True
cyclase_percentage <- (cyclase_loci / total_loci) * 100

#Number of loci with both HD and cyclase
HD_cyclase_loci <- nrow(subset(info, HD_hmm_boolean == "True" & cyclase == TRUE))
HD_cyclase_percentage <- (HD_cyclase_loci / total_loci) * 100

#calculate how many known effectors there are in the dataset. The known effectors are: can3, cam1, saved.chat, nucc, calpl, cami1, can1.2, csx1, csx23, csm6.ca6, cora
info$known_effector_count <- rowSums(info[,c("cam1", "saved.chat", "nucc", "calpl", "cami1", "can1.2", "csx1", "csx23", "csm6.ca6", "cora")])
known_effectors_total <- sum(info$known_effector_count)

#calculate counts for each known effector
known_effectors <- c("cam1", "saved.chat", "nucc", "calpl", "cami1", "can1.2", "csx1", "csx23", "csm6.ca6", "cora")
known_effectors_count <- c()
for (effector in known_effectors) {
  effector_count <- nrow(subset(info, get(effector) == TRUE))
  known_effectors_count <- c(known_effectors_count, effector_count)
}
#store the counts in a dataframe
known_effectors_count_df <- data.frame(known_effectors, known_effectors_count)

#calculate how many new effectors. They are: can3, cam2, cam3, tirsaved
info$has_validated_new_effector <- ifelse(info$can3 == TRUE | info$cam2 == TRUE | info$cam3 == TRUE | info$tirsaved == TRUE, 1, 0)
new_effectors_total <- sum(info$has_validated_new_effector)

#Make a table of the above
summary_table <- data.frame(
  "Total loci" = total_loci,
  "Loci with HD domain (HMM)" = HD_loci_hmm,
  "Loci with HD domain (literal)" = HD_loci_literal,
  "Loci with cyclase domain" = cyclase_loci,
  "Loci with both HD and cyclase" = HD_cyclase_loci,
  "Percentage of HD and cyclase" = HD_cyclase_percentage,
  "Percentage of loci with HD domain" = HD_percentage,
  "Percentage of loci with cyclase domain" = cyclase_percentage,
  "Percentage of loci with HD domain and subtype III-A" = HD_subtype_IIIA_percentage,
  "Percentage of loci with HD domain and subtype III-B" = HD_subtype_IIIB_percentage,
  "Percentage of loci with HD domain and subtype III-D" = HD_subtype_IIID_percentage,
  "Percentage of loci with HD domain and subtype III-F" = HD_subtype_IIIF_percentage,
  "Percentage of loci with HD domain and subtype III-C" = HD_subtype_IIIC_percentage,
  "Percentage of loci with cyclase domain and subtype III-A" = cyclase_IIIA_percentage,
  "Percentage of loci with cyclase domain and subtype III-B" = cyclase_IIIB_percentage,
  "Percentage of loci with cyclase domain and subtype III-D" = cyclase_IIID_percentage,
  "Percentage of loci with cyclase domain and subtype III-F" = cyclase_IIIF_percentage,
  "Percentage of loci with cyclase domain and subtype III-C" = cyclase_IIIC_percentage,
  "Type III-A count" = nrow(subset(info, Subtype == "III-A")),
  "Type III-B count" = nrow(subset(info, Subtype == "III-B")),
  "Type III-D count" = nrow(subset(info, Subtype == "III-D")),
  "Type III-F count" = nrow(subset(info, Subtype == "III-F")),
  "Type III-C count" = nrow(subset(info, Subtype == "III-C")),
  "Effector count all" = sum(info$effector_count_known_new_sum),
  "Known effector count" = known_effectors_total,
  "New effectors" = new_effectors_total,
  "Cam2 positives" = nrow(subset(info, cam2 == TRUE)),
  "Cam3 positives" = nrow(subset(info, cam3 == TRUE)),
  "Can3 positives" = nrow(subset(info, can3 == TRUE)),
  "TIR-SAVED positives" = nrow(subset(info, tirsaved == TRUE))
)

summary_table
