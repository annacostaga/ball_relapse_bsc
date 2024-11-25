#!/bin/bash
#SBATCH --job-name=GRaNIE.R
#SBATCH --chdir=./
#SBATCH --output=GRaNIE%j.out
#SBATCH --error=GRaNIE%j.err
#SBATCH --cpus-per-task=90
#SBATCH --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --qos=gp_bscls

cd ../
module load R/4.3.2
Rscript GRaNIE.R




