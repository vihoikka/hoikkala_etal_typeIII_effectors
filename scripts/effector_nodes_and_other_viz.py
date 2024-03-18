import pandas as pd
import itertools
import matplotlib.pyplot as plt
import seaborn as sns
from upsetplot import UpSet
from matplotlib_venn import venn3
from upsetplot import from_memberships, plot
import argparse
argparser = argparse.ArgumentParser()

argparser.add_argument('-i', '--input_mastertable', type=str, required=True,
                    help='Input mastertable')
argparser.add_argument('-e', '--effector_list', type=str, required=True,
                    help='Input effector list')
argparser.add_argument('-o', '--outdir', type=str, required=True,
                    help='Output directory')
argparser.add_argument('-n', '--out_nodes', type=str, required=True,
                    help='Output nodes file')
argparser.add_argument('-d', '--out_edges', type=str, required=True,
                    help='Output edges file')

#generate shell command based on above args
# python effector_nodes_and_other_viz.py -i mastertable.tsv -e effectors.txt -o viz -n nodes.tsv -d edges.tsv

args  = argparser.parse_args()
mastertable = args.input_mastertable
effector_list = args.effector_list
outdir = args.outdir
output_nodes = args.out_nodes
output_edges = args.out_edges


#this dict is used to add a cOA attribute to the nodes output file
coa_dict_known = {"cami1": "cA4",
            "calpl": "cA4",
            "can1-2": "cA4",
            "cora": "sam-amp",
            #"can3": "cA6",
            "csm6-ca6": "cA6",
            "csx1": "cA4",
            "csx23": "cA4",
            "csx6": "cA4",
            "nucc": "cA3",
            #"tirsaved": "cA4",
            "cam1": "cA4",
            #"cam2": "cA3",
            #"cam3": "cA4",
            "saved-chat": "cA4",
            "nucc": "cA3"
        }

coa_dict_all = {"cami1": "cA4",
            "calpl": "cA4",
            "can1-2": "cA4",
            "cora": "sam-amp",
            "can3": "cA6",
            "csm6-ca6": "cA6",
            "csx1": "cA4",
            "csx23": "cA4",
            "csx6": "cA4",
            "nucc": "cA3",
            "tirsaved": "cA4",
            "cam1": "cA4",
            "cam2": "cA3",
            "cam3": "cA4",
            "saved-chat": "cA4",
            "nucc": "cA3"
        }



def create_node_graph(effector_signal_map, mode):
    effectors_all = pd.read_csv(str(effector_list), sep = "\t")
    #get unique effectors by looking at column 1
    effectors_all = effectors_all.iloc[:,1].unique()
    print(effectors_all)
    #pick only the effectors listed in the coa_dict as keys
    effectors = [effector for effector in effectors_all if effector in effector_signal_map]
    print(effectors)
    #create dictionary with effectors as keys and empty string as values
    effectors_coa_dict = dict.fromkeys(effectors, "")
    #create a dataframe with all possible combinations of effectors
    edges = pd.DataFrame(list(itertools.combinations(effectors, 2)), columns = ["source", "target"])
    edges["type"] = "Undirected"
    edges["weight"] = 0
    print(edges)

    #create nodes df from unique values in edges sources column
    nodes = pd.DataFrame(edges["source"].unique(), columns = ["id"])
    #compare to target colunn to see if any effectors were left out because they were only in source
    nodes2 = pd.DataFrame(edges["target"].unique(), columns = ["id"])
    #concatenate the two dataframes
    nodes = pd.concat([nodes, nodes2])
    #get the coa attribute for each effector
    nodes["signal"] = nodes["id"].map(effector_signal_map)
    nodes["totalCount"] = 0
    nodes['cooccurenceCount'] = 0
    nodes['soloCount'] = 0

    mastertable_df = pd.read_csv(str(mastertable), sep = "\t")
    #iterate through mastertable
    for index, row in mastertable_df.iterrows():
        #for each row, check which effectors are present and make weight += 1 when co-occurrence happens
        for effector in effectors:
            secondEffectorFound = False
            if row[effector] == True: #if effector is present
                #add one to effector's total count in nodes
                nodes.loc[nodes["id"] == effector, "totalCount"] += 1
                #check for presence of other effectors in the same locus
                for effector2 in effectors:
                    if row[effector2] == True: #if effector2 is present
                        #if effector and effector2 are not the same
                        if effector != effector2:
                            secondEffectorFound = True #this means that in the current effector in this locus is not soloing
                            #add 1 to the weight column of the row where source = effector and target = effector2
                            edges.loc[(edges["source"] == effector) & (edges["target"] == effector2), "weight"] += 1
                if secondEffectorFound:
                    #add 1 to the cooccurenceCount column of the row where id = effector
                    nodes.loc[nodes["id"] == effector, "cooccurenceCount"] += 1
                else:
                    #if effector is the only effector in the locus, add 1 to soloCount
                    nodes.loc[nodes["id"] == effector, "soloCount"] += 1

    #calculate some statics in the nodes table
    nodes["proportionSolo"] = nodes["soloCount"] / nodes["totalCount"]
    nodes["proportionCooccurence"] = nodes["cooccurenceCount"] / nodes["totalCount"]

    rename_dict = {"cami1": "Cami1",
                "calpl": "CalpL",
                "can1-2": "Can1-2",
                "cora": "CorA",
                "can3": "can3",
                "csm6-ca6": "Csm6-cA6",
                "csx1": "Csx1",
                "csx23": "Csx23",
                "csx6": "Csx6",
                "nucc": "NucC",
                "tirsaved": "TIR-SAVED",
                "cam1": "Cam1",
                "cam2": "Cam2",
                "cam3": "Cam3",
                "saved-chat": "SAVED-CHAT",
                "nucc": "NucC"
            }

    #rename the ids in the nodes table
    nodes["id"] = nodes["id"].map(rename_dict)
    #rename the source and target columns in the edges table
    edges["source"] = edges["source"].map(rename_dict)
    edges["target"] = edges["target"].map(rename_dict)

    #remove file extension from output_nodes
    output_nodes_path = output_nodes.split(".")[0] + "_" + mode + ".tsv"
    output_edges_path = output_edges.split(".")[0] + "_" + mode + ".tsv"

    #save the edges table
    edges.to_csv(str(output_edges_path), sep = "\t", index = False)
    #save nodes
    nodes.to_csv(str(output_nodes_path), sep = "\t", index = False)

create_node_graph(coa_dict_known, "known")
create_node_graph(coa_dict_all, "all")

# pieGraphs = {}
# sizes_dict = {}
# labels_dict = {}
# # in addition to the edges/nodes for the gephi, also make separate bar graphs for each effectors showing how often it co-occurs with other ones
# for effector in effectors:
#     #create a dataframe with the effector as the source
#     effector_df = edges[edges["source"] == effector]
#     #sort by weight
#     effector_df = effector_df.sort_values(by = "weight", ascending = False)
#     #save the dataframe
#     effector_df.to_csv(outdir + "/" + effector + ".tsv", sep = "\t", index = False)
#     #render as bar graph using matplot lib imported above
#     plt.bar(effector_df["target"], effector_df["weight"])
#     plt.xticks(rotation = 90)
#     plt.title(effector)
#     #rotate x axis labels
#     plt.tight_layout()
#     plt.savefig(outdir + "/" + effector + "_cooccurrence.png")
#     plt.close()

#     #create a pie graph showing the proportion of co-occurrence vs soloing
#     labels = ["Co-occurrence", "Solo"]
#     sizes = [nodes.loc[nodes["id"] == effector, "proportionCooccurence"].values[0], nodes.loc[nodes["id"] == effector, "proportionSolo"].values[0]]
#     plt.pie(sizes, labels = labels, autopct='%1.1f%%', shadow=True, startangle=90)
#     plt.title(effector)
#     plt.tight_layout()
#     #add plot to piegraphs dict for later concatenation
#     pieGraphs[effector] = plt
#     #also save alone
#     plt.savefig(outdir + "/" + effector + "_pie.png")
#     plt.close()
#     sizes_dict[effector] = sizes
#     labels_dict[effector] = labels

# #from all the piegraphs, create a collage
# #create a new figure
# fig = plt.figure(figsize=(20, 20))
# #add a grid to the figure
# grid = plt.GridSpec(4, 4, wspace=0.4, hspace=0.4)

# #add the pie graphs to the grid. Organize it so that the ones with the highest co-occurrence are on the top left
# for i, effector in enumerate(sorted(pieGraphs, key = lambda x: sizes_dict[x][0], reverse = True)):
#     ax = fig.add_subplot(grid[i])
#     ax.pie(sizes_dict[effector], labels = labels_dict[effector], autopct='%1.1f%%', shadow=True, startangle=90)
#     ax.set_title(effector)


# #save the figure
# plt.savefig(outdir + "/pie_graphs_all.png")

# #Create an UpSet plot. First take the mastertable and only retain all effectors, which can be fetched from the effectors list
# # Convert DataFrame to Boolean
# upsetTable = mastertable[effectors].astype(bool)

# # Convert any False to NaN and then assign True to loci where effectors are present
# upsetTable = upsetTable[upsetTable] 

# # Replace NaNs to False
# upsetTable.fillna(False, inplace=True)

# # Create a list of effector combinations per locus (row)
# locus_to_effector = upsetTable.apply(lambda row: row.index[row].tolist(), axis=1)

# # Convert list of lists to list of tuples
# locus_to_effector_tuples = locus_to_effector.apply(lambda x: tuple(x)).tolist()

# # Generate the data for upset plot
# upset_data = from_memberships(locus_to_effector_tuples)

# upset_data.to_csv(outdir + "/upset_data.tsv", sep = "\t", index = False)

# # Create UpSet plot and save
# plot(upset_data, subset_size='count', sort_by = "degree")
# plt.savefig(outdir + "/upset_plot.png")