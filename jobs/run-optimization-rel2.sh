#!/bin/bash
#SBATCH --job-name="RFTO"
#SBATCH --partition=high_p
#SBATCH -c 1
#SBATCH --mem-per-cpu=30000


##we use the base code and containers
ROOT=/data/scratch/digenovaa/Somatic-reference-free/SNV-INDELs/RF-mut-f
#make file
RS=${ROOT}/code/Rscripts/RF-SNV-optimization-dist-rel2.R
#container
CONTAINER=${ROOT}/container/rf-mut-f_v2.0.sif
#release to process
REL2DIR=${ROOT}/mesomics/release2/matched-t-only

## Grid search parameters a total of 48 combinations
g=`awk NR==${SLURM_ARRAY_TASK_ID} ${ROOT}/code/Rscripts/grid_search_parameters.txt`

#with MESO_061
singularity exec $CONTAINER Rscript ${RS} $g Somatics.snv.matrix.txt Germline.snv.matrix.txt ${REL2DIR} > ${SLURM_ARRAY_TASK_ID}.meso61.r1.log
#singularity exec $CONTAINER Rscript ${RS} $g Somatics.snv.matrix.txt Germline.snv.matrix.txt ${REL2DIR} > ${SLURM_ARRAY_TASK_ID}.meso61.r2.log
#singularity exec $CONTAINER Rscript ${RS} $g Somatics.snv.matrix.txt Germline.snv.matrix.txt ${REL2DIR} > ${SLURM_ARRAY_TASK_ID}.meso61.r3.log
#Three replicates without MESO_061
singularity exec $CONTAINER Rscript ${RS} $g Somatics.snv.matrix.WMESO_061.txt Germline.snv.matrix.WMESO_061.txt ${REL2DIR} > ${SLURM_ARRAY_TASK_ID}.wmeso61.r1.log
singularity exec $CONTAINER Rscript ${RS} $g Somatics.snv.matrix.WMESO_061.txt Germline.snv.matrix.WMESO_061.txt ${REL2DIR} > ${SLURM_ARRAY_TASK_ID}.wmeso61.r2.log
singularity exec $CONTAINER Rscript ${RS} $g Somatics.snv.matrix.WMESO_061.txt Germline.snv.matrix.WMESO_061.txt ${REL2DIR} > ${SLURM_ARRAY_TASK_ID}.wmeso61.r3.log
