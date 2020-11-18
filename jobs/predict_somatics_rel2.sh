#!/bin/bash
#SBATCH --job-name="P-TONLY"
#SBATCH --partition=high_p
#SBATCH -c 30
#SBATCH --mem-per-cpu=2000

##we use the base code and containers
ROOT=/data/scratch/digenovaa/Somatic-reference-free/SNV-INDELs/RF-mut-f
#make file
CM=${ROOT}/code/makefiles/predict_t_only.mk
#container
CONTAINER=${ROOT}/container/rf-mut-f_v2.0.sif
#release to process
REL2DIR=${ROOT}/mesomics/release2/t-only

cd ${REL2DIR}
#exec
singularity exec $CONTAINER make -f ${CM} all -j 30
