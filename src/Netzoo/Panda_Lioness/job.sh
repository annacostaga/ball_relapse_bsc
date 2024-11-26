#!/bin/bash
#SBATCH --job-name=lioness_job
#SBATCH --output=lioness_%j.out
#SBATCH --error=lioness_%j.err
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --qos=gp_bscls

module load python
cd ../ball_relapse_bsc/src/Netzoo/Panda_Lioness
python
python Panda_Lioness.py
