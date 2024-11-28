#Import Packages and Set directory

import os
import pandas as pd
import matplotlib.pyplot as plt

os.chdir("/gpfs/projects/bsc08/shared_projects/BALL_RELAPSE/results/GRNs/GRNBoost2")

#Load Network
net = pd.read_csv("GRN.csv", sep=",")

#Check distribution of score between edges
plt.hist(net["importance"], bins=150, color="purple")
plt.xlabel("score")
plt.ylabel("number of edges")
plt.savefig(fname = "Distribution_of_score_edges.png")

#Check distribution of edges per TFs
edge_counts = net.groupby('TF').size().reset_index(name='edge_count')
plt.hist(edge_counts["edge_count"], bins=150, color="purple")
plt.xlabel("number of edges")
plt.ylabel("TFs")
plt.savefig("Distribution_of_edges_TFs.png")


