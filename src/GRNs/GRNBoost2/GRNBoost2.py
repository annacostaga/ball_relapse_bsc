

#Import packages and Set directory to data
import os
import pandas as pd
from arboreto.algo import grnboost2

os.chdir("../ball_relapse_bsc/annotation/GRNs/GRNBoost2")

#Load and Prepare inputs
expression = pd.read_csv("../ball_relapse_bsc/results/rna_seq_primary/rlog_filtered_Bcounts.csv", sep=',', index_col=0)
expression = expression.astype(float)
expression.columns = expression.columns.astype(str)
expression = expression.T

tf_names = pd.read_csv("allTFs_hg38.txt")
tf_list = tf_names["TFs"].to_list()

#Run GRNBoost2
network = grnboost2(expression_data=expression, tf_names=tf_list)

#Save network in local
network.to_csv("../ball_relapse_bsc/results/GRNs/GRNBoost2/GRN.csv")


