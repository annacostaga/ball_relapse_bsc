############### RUN in MN5 ##########################

#Import packages and Set directory to data

import os
import pandas as pd
import numpy as np
import netZooPy
from netZooPy.panda import Panda
from netZooPy.lioness import Lioness

os.chdir("../gpfs/projects/bsc08/shared_projects/BALL_RELAPSE/annotations/Netzoo/Panda_Lioness")

#Load data (The data is located in BALL_RELAPSE/annotations/Netzoo/Panda_Lioness in MN5)

expression = pd.read_csv("average_filtered_counts.csv", index_col=0, sep=",")
motifs = pd.read_csv("motifs_hg38.csv", index_col=0, sep=",")
PPI = pd.read_csv("PPI.csv", index_col=0, sep=",")


#Run Panda and save results

panda_obj = Panda(expression, motifs, PPI, remove_missing=False,keep_expression_matrix=True, save_memory=False, modeProcess='legacy')
panda_obj.save_panda_results('../ball_relapse_bsc/results/Netzoo/Panda_Lioness/output_panda.txt')

#Run Lioness and save results

lioness_obj = Lioness(panda_obj)
lioness_obj.save_lioness_results('../ball_relapse_bsc/results/Netzoo/Panda_Lioness/lioness.txt')





