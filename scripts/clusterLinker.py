import pandas as pd
import os
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("-c10", "--cas10_clusters", required=True)
parser.add_argument("-c7", "--cas7_clusters", required=True)
parser.add_argument("-c5", "--cas5_clusters", required=True)
parser.add_argument("-co", "--cora_clusters", required=True)
parser.add_argument("-s", "--sample", required=True)
parser.add_argument("-o", "--out", required = True)
args = parser.parse_args()

cas10_clusters = args.cas10_clusters
cas7_clusters = args.cas7_clusters
cas5_clusters = args.cas5_clusters
cora_clusters = args.cora_clusters
out = args.out
sample = args.sample

# 1. Separate in clusters, return dic
def separate_clusters(clusterfile):
    print("Converting cluster file to dict")
    all_clusters = {}
    current_cluster_samplelist = []
    notFirstCluster = False

    with open(clusterfile) as f:
        lines = f.readlines()
        last = lines[-1]
    for line in lines:
        print(line)
        if (line[0] == ">") and (notFirstCluster == False): #if first line in file
            currentClusterName = "none" #we define cluster names by representatives, not cluster number
            print("First line")
            print("Current cluster " + str(currentClusterName))
            notFirstCluster = True
        elif (line[0] == ">") and (notFirstCluster == True): #start of new cluster. Wrap up previous samples as key/value of dic and reset
            print("New header found. Wrapping up previous lines for cluster " + str(currentClusterName))
            all_clusters[currentClusterName] = current_cluster_samplelist
            currentClusterName = "none"
            current_cluster_samplelist = []
        elif line is last:
            print("End of file. Wrapping up remaining lines for cluster " + str(currentClusterName))
            line_sample = line.split(">")[1].split("...")[0]
            current_cluster_samplelist.append(line_sample)
            all_clusters[currentClusterName] = current_cluster_samplelist
        else:
            line_sample = line.split(">")[1].split("...")[0]
            if "*" in line: #check if this is the representative
                currentClusterName = line_sample
            print("Writing lines sample " + line_sample)
            current_cluster_samplelist.append(line_sample)
        
    
    print(all_clusters)
    return all_clusters

# 2. Find which clusters contains current sample. To speed up, convert to pandas df.
def find_cluster(clusters, current_sample):
    for key, value in clusters.items():
        for i in value:
            if i == current_sample:
                print("Found match")
                return key
    return "none"

properties = {"Cas10": cas10_clusters, "Cas7": cas7_clusters, "Cas5": cas5_clusters, "CorA": cora_clusters}
infoDic = {}
infoDic["Sample"] = sample


for property, value in properties.items():
    print("Examining " + str(property))
    clusters = separate_clusters(value) #build dictionary of clusters
    assigned_cluster = find_cluster(clusters, sample) #find which cluster has this sample
    infoDic[str(property) + "_cluster_link"] = assigned_cluster #assign the cluster to current property

#convert dictionary with cluster info to df
print(infoDic)
infoDic["16S_link"] = sample #16s is not clustered, so we just refer to the sample itself here.
finalDf = pd.DataFrame.from_dict([infoDic])
print(finalDf)

#write
finalDf.to_csv(out, sep = "\t", header = 0, index=False)