#Import Packages and Set directory

import os
import pandas as pd
import matplotlib.pyplot as plt

os.chdir("/gpfs/projects/bsc08/shared_projects/BALL_RELAPSE/results/GRNs/GRNBoost2")

#Load Network and relevant_genes
net = pd.read_csv("GRN.csv", sep=",")
genes = pd.read_csv("../ball_relapse_bsc/annotation/GRNs/Filtering_GRN/Bcell_relevant_genes.csv", sep ="," )
genes = genes.iloc[: , 1:]

#Step 1: Filter GRN based on set of relevant genes
gene_set = set(genes['Gene'])
filtered_net = net[net['TF'].isin(gene_set) | net['target'].isin(gene_set)]

#Step 2: Filter GRN based on score-distribution
filtered_net = filtered_net[filtered_net["importance"] > 0.3]

#Check relevant TFs in filtered GRN
relevant_TFs = filtered_net.merge(genes, left_on='TF', right_on='Gene')
print(relevant_TFs["TF"].unique())

#Check degree of connectivity for relevant TFs
edge_counts = relevant_TFs.groupby('TF').size().reset_index(name='edge_count')
plt.hist(edge_counts["edge_count"], bins=150, color="purple")
plt.xlabel(" number of edges")
plt.ylabel("relevant TFs")
plt.savefig("Degree_of_edges_relevant_TFs.png")


#Save filtered network
filtered_net.to_csv('filtered_GRN.csv')


