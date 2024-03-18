# install.packages("wakefield")
# install.packages("viridis")
# install.packages("phytools")
# install.packages("shadowtext")
# install.packages("ggnewscale")
# install.packages("xlsx")
# install.packages("ggrepel")
# install.packages("ggnewscale")

library("treeio")
library("ggtree")
library("wakefield")
library("stringr")
library("viridis")
library("shadowtext")
library("phytools")
library("dplyr")
library("ggplot2")
library("ggnewscale")
library("xlsx")
library("ggrepel")
library("scales")
library("patchwork")

#### Load files and tidy ####

#define project folder
project = "run1"

#define paths
cas10_tree_path <- paste("data/",project,"/cas10_tree.txt", sep = "")
CorA_tree_path <- paste("data/",project,"/CorA_tree.txt", sep = "")
CorA_tree_path_unclustered <- paste("data/",project,"/CorA_tree_unclustered.txt", sep = "")

info_path <- paste("data/",project,"/type_iii_info.tsv", sep = "")
cATyper_path <- paste("data/",project,"/cATyper_all.tsv", sep = "")
tax_info_path <- paste("data/",project,"/taxInfo.txt", sep = "")

#open files
cas10_tree <- read.newick(cas10_tree_path)
cas10_tree <- midpoint.root(cas10_tree)
#CorA_tree <- read.newick(CorA_tree_path)
#CorA_tree <- midpoint.root(CorA_tree)
CorA_tree_unclustered <- read.newick(CorA_tree_path_unclustered)
CorA_tree_unclustered <- midpoint.root(CorA_tree_unclustered)
info <- read.csv(info_path, sep="\t")
cATyper = read.csv(cATyper_path, sep="\t")
tax_info <- read.csv(tax_info_path, sep="\t") #produced by a separate rule in the Snakemake pipeline

#merge and tidy opened files
info <- merge(info,tax_info, by="Sample") #merge info tables with taxonomy information
info <- merge(info, cATyper, by.x = "Locus", by.y = "locus") #merge info with cATyper data

#tidy
info$Hybrid <- FALSE #no hybrids is default
info$Hybrid[str_detect(info$Subtype, 'Hybrid') == TRUE] <- TRUE #mark some as hybrids
info$Subtype[str_detect(info$Subtype, 'Hybrid') == TRUE & str_detect(info$Subtype, 'III-A') == TRUE] <- "III-A"
info$Subtype[str_detect(info$Subtype, 'Hybrid') == TRUE & str_detect(info$Subtype, 'III-B') == TRUE] <- "III-B"
info$Subtype[str_detect(info$Subtype, 'Hybrid') == TRUE & str_detect(info$Subtype, 'III-C') == TRUE] <- "III-C"
info$Subtype[str_detect(info$Subtype, 'Hybrid') == TRUE & str_detect(info$Subtype, 'III-D') == TRUE] <- "III-D"

#merge info with protein specific data
rownames(info) <- info$Locus #set rownames

infosubset <- subset(info, info$CorA.y == "True")
info_unk <- subset(info, info$Unknown_genes == "True")

#set a labeling column where only CorA positives have the sample name
info$CorALabel <- with(info, ifelse(CorA.y == "False", '', Sample))

#subsets
info_cora = subset(info, info$CorA.y=="True")
cora_filepath = "cas10_cora_positives.xlsx"
write.xlsx(info_cora, cora_filepath, sheetName = "CorA positives", 
           col.names = TRUE, row.names = TRUE, append = FALSE)

colorBlindBlack8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#### Cas10 circular with branch lengths ####
p1 <- ggtree(cas10_tree, layout="circular", size=0.15, open.angle=20, aes(color=Subtype, labels=genus) +
               theme(legend.position = "none") +
               geom_text2(aes(subset=!isTip, label=node), hjust = 2 )) %<+% info

p1 <- gheatmap(p1, info[c("Subtype")], offset=-0.1, width=0.05, colnames = FALSE) + 
  scale_fill_manual(values = colorBlindBlack8, name = "Subtype") +
  theme(text = element_text(size = 14))

p1 <- p1 + geom_tiplab2(aes(label="",  subset = CorA.y == "True"), 
                            align=TRUE, family='mono',
                            linetype = "dotted", linesize = .2, alpha = 1.0, color = "red")
p2 <- p1 + new_scale_fill()
p2 <- gheatmap(p2, info[c("CorA.y")], offset=0.3, width=0.04,colnames = FALSE) +
  scale_fill_manual(values=c("grey90", "red"), name="CorA") 
#p2

plot_path <- paste("data/",project,"/plots/",project,"_Cas10_branchslengths.png", sep = "")
ggsave(plot_path, plot = p2, width = 8, height = 8)

#### CorA circular nonclustered ####
p1_nonclustered <- ggtree(CorA_tree_unclustered, layout="circular", aes(color=Subtype, labels=genus) +
               theme(legend.position = "none") +
               geom_text2(aes(subset=!isTip, label=node), hjust = 2 )) %<+% info

p1_nonclustered <- gheatmap(p1_nonclustered, info[c("Subtype")], width=0.1, colnames = FALSE) + 
  scale_fill_manual(values = colorBlindBlack8, name = "Subtype") +
  theme(text = element_text(size = 14))

p2_nonclustered <- p1_nonclustered + geom_tiplab2(aes(label=""), 
                        align=TRUE, family='mono',
                        linetype = "dotted", linesize = .7, alpha = 0.3)

plot_path <- paste("data/",project,"/plots/",project,"_CorA_unclustered_branchlengths.png", sep = "")
ggsave(plot_path, plot = p2_nonclustered, width = 6, height = 6)