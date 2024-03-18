from BCBio import GFF
from dna_features_viewer import GraphicFeature, GraphicRecord, BiopythonTranslator
import matplotlib.pyplot as plt
import argparse
import pandas as pd
import os
import gffutils

# create args inputs for cas_operons_file, locus, output_folder, host_genomes_folder, mode, hmm_rows, hmm_targets, catyper_out, cctyper_protein_table
parser = argparse.ArgumentParser(description='Visualises expanded CRISPR loci')
parser.add_argument('-i', '--cas_operons_file', help='input file', required=True)
parser.add_argument('-c', '--cctyper_folder', help='cctyper folder', required=True)
parser.add_argument('-l', '--locus', help='locus', required=True)
parser.add_argument('-o', '--output_folder', help='output folder', required=True)
parser.add_argument('-hg', '--host_genomes_folder', help='host genomes folder', required=True)
parser.add_argument('-v', '--validated_effectors', help='validated effectors', required=True)
parser.add_argument('-cc', '--cctyper_protein_table', help='cctyper plottable table', required=True)
parser.add_argument('-k', '--known_effector_table', help='known effectors', required=True)
args = parser.parse_args()

#generate bash command for the above
# python locus_visualiser.py -i cas_operons.tsv -l locus -o output_folder -hg host_genomes_folder

#assign args to similarly named variables
cas_operons_file = args.cas_operons_file
locus = args.locus
output_folder = args.output_folder
host_genomes_folder = args.host_genomes_folder
cctyper_folder = args.cctyper_folder
validated_effectors = args.validated_effectors
cctyper_protein_table = args.cctyper_protein_table
known_effector_table = args.known_effector_table

validated_effectors = pd.read_csv(validated_effectors, sep = "\t", header = 0)
cctyper_protein_table = pd.read_csv(cctyper_protein_table, sep = "\t", header = 0)
known_effector_table = pd.read_csv(known_effector_table, sep = "\t", header = 0)

#merge all three tables by id
merged = validated_effectors.merge(cctyper_protein_table, on = "protein_id", how = "outer")
merged = merged.merge(known_effector_table, on = "protein_id", how = "outer")

#when looking for effectors in the locus, the effector search range is the number of bases up or downstream of the cctyper defined cas operon boundaries
effector_search_range = 4000

#get sample name by splitting the locus name at the underscore and taking the first two parts
sample = locus.split("_")[0] + "_" + locus.split("_")[1]

#read in the cas operons file
cas_operons_df = pd.read_csv(cas_operons_file, sep ="\t", header = 0)

#read in the crispr gff file. Not all loci have crisprs
try:
    crispr_gff_path = os.path.join(cctyper_folder, sample, "crisprs.gff")
except:
    print("No crispr gff file found for sample " + sample)

gff_file = os.path.join(host_genomes_folder, sample, sample + "_features.gff")

#get contig ID of the CRISPR positive contig
contig = cas_operons_df.iloc[0]['Contig'] #This is wrong. 
print("Cas operons file: " +str(cas_operons_df))
print("The chosen contig: " + contig)

#define start and end coordinates for plotting
cas_operon_start = int(cas_operons_df['Start'][0]) - effector_search_range
cas_operon_end = int(cas_operons_df['End'][0]) + effector_search_range
features_list = []

def create_gff_iterator(gff_file, contig):
    """
    Creates a GFF iterator for a given GFF file and contig.
    Note that the if the in_handle is closed, then the script will crash later when the iterator is used. This is why we leave it open
    """
    in_handle = open(gff_file)
    limit_info = dict(gff_id = [contig], 
                      gff_type = ["CDS", "repeat_region"])  # Add or remove types as needed
    gff_iterator = GFF.parse(in_handle, limit_info=limit_info)
    return gff_iterator

def add_features_from_gff(gff_iterator, features_list, cas_operon_start, cas_operon_end, source):
    """
    Adds GraphicFeature objects to a list from a given GFF iterator
    """
    color_dict = {"original_gff": "#ffcccc", "crispr_gff": "#8f8f8f", "repeat_region": "#e3e2dc"}
    new_features_list = []
    for rec in gff_iterator:
        for feature in rec.features:
            if int(feature.location.start) >= cas_operon_start and int(feature.location.end) <= cas_operon_end:
                print(feature)
                if source == "crispr_gff":
                    label = "Repeat region"
                    strand = 0
                else:
                    label = feature.qualifiers.get('product', [None])[0]
                    #if label is not NoneType and contains the string "unknown" or "hypthetical", then make label None
                    if label is not None and ("unknown" in label or "hypothetical" in label):
                        label = None
                    strand=1 if feature.strand == 1 else -1
                graphic_feature = GraphicFeature(
                    start=int(feature.location.start)-cas_operon_start, 
                    end=int(feature.location.end)-cas_operon_start,
                    strand=strand,
                    color=color_dict[source], 
                    label=label
                    )
                print(graphic_feature)
                new_features_list.append(graphic_feature)
    #combine new and old feature lists
    features_list = features_list + new_features_list
    return features_list

gff_iterator = create_gff_iterator(gff_file, contig)

#check if the file crispr_gff_path exists
if os.path.isfile(crispr_gff_path):
    crispr_gff_iterator = create_gff_iterator(crispr_gff_path, contig)

# Prepare a list to hold the Graphic Features
features_list = []

# Add features from GFF
features_list = add_features_from_gff(gff_iterator, features_list, cas_operon_start, cas_operon_end, source = "original_gff")

if os.path.isfile(crispr_gff_path):
    features_list = add_features_from_gff(crispr_gff_iterator, features_list, cas_operon_start, cas_operon_end, source = "crispr_gff")

# Add features from the validated_effectors table. Columns are protein_id	start	end	effector	locus	sample	strand	type
strand_dict = {"+": 1, "-": -1}

for index, row in validated_effectors.iterrows():
    if row["locus"] == locus:
        graphic_feature = GraphicFeature(
            start = int(row["start"])-cas_operon_start,
            end = int(row["end"])-cas_operon_start,
            strand = strand_dict[row["strand"]], #need to convert + and - to 1 and -1
            color = "#e6652e",
            label = row["effector"]
        )
        print("graphic feature from validated effectors: " + str(graphic_feature))
        features_list.append(graphic_feature)

# Similary add features from the cctyper_protein_table. Columns are protein_id	start	end	annotation_gff	annotation_cctyper	strand
for index, row in cctyper_protein_table.iterrows():

    #Forgive me but...
    if row["strand_cctyper"] == "+":
        updated_strand = 1
    elif row["strand_cctyper"] == "-":
        updated_strand = -1
    elif row["strand_cctyper"] == "+1":
        updated_strand = 1
    elif row["strand_cctyper"] == "-1":
        updated_strand = -1

    graphic_feature = GraphicFeature(
        start = int(row["start_cctyper"])-cas_operon_start,
        end = int(row["end_cctyper"])-cas_operon_start,
        strand = updated_strand,
        color = "#79c77c",
        label = row["annotation_cctyper_cctyper"].split("_")[0]
    )
    print("graphic feature from cctyper protein table: " + str(graphic_feature))
    features_list.append(graphic_feature)


# similarly add features from known_effector_table
for index, row in known_effector_table.iterrows():
    if row["locus"] == locus:
        graphic_feature = GraphicFeature(
            start = int(row["start"])-cas_operon_start,
            end = int(row["end"])-cas_operon_start,
            strand = strand_dict[row["strand"]], #need to convert + and - to 1 and -1
            color = "#4372b0",
            label = row["effector"]
        )
        print("graphic feature from known effectors: " + str(graphic_feature))
        features_list.append(graphic_feature)

# finally add the CRISPR-Cas range
features_list.append(
    GraphicFeature(
        start = effector_search_range,  # No need to adjust by cas_operon_start
        end = cas_operon_end-cas_operon_start-effector_search_range,  # Adjust end by subtracting cas_operon_start and effector_search_range
        strand = 0,
        color = "#e3e2dc",
        label = "CRISPR-Cas locus",
    )
)

print(features_list)

# Generate plot using DNA Features Viewer
record = GraphicRecord(sequence_length=cas_operon_end-cas_operon_start, features=features_list)
ax, _ = record.plot(figure_width=10)

#save
plt.savefig(output_folder + "/" + locus + "_viz.png")

#save the merged
merged.to_csv(output_folder + "/" + locus + "_merged.csv", sep = "\t", index = False)