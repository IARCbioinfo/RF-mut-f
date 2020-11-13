#!/bin/bash
#SBATCH --job-name="gnomad"
#SBATCH --partition=high_p
#SBATCH --share
#SBATCH -c 24
#SBATCH --mem-per-cpu=3000

## conda path
export PATH=/home/digenovaa/miniconda3/bin/:$PATH
make -j 24 all
