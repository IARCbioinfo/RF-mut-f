#!/bin/bash
#SBATCH --job-name="RFT"
#SBATCH --partition=high_p
#SBATCH -c 1
#SBATCH --mem-per-cpu=30000

## conda path
##we use the base code and containers
ROOT=/data/scratch/digenovaa/Somatic-reference-free/SNV-INDELs/RF-mut-f
#make file
MS=${ROOT}/code/Rscripts/RF-BEST-SNV-REL2.R
MI=${ROOT}/code/Rscripts/RF-BEST-INDEL-REL2.R
#container
CONTAINER=${ROOT}/container/rf-mut-f_v2.0.sif
#release to process
REL2DIR=${ROOT}/mesomics/release2/matched-t-only
#before fixing mpos 12 1500 25
#model for SNVs
#singularity exec $CONTAINER Rscript ${MS} 12 1500 25 Somatics.snv.matrix.txt Germline.snv.matrix.txt rf-12_1500_25_snv_meso61.rds ${REL2DIR}
#model for INDELS, inlcuding SNVs+INDELs
#singularity exec $CONTAINER Rscript ${MI} 12 1500 25 Somatics.snv.matrix.txt Germline.snv.matrix.txt Somatics.indel.matrix.txt Germline.indel.matrix.txt rf-12_1500_25_indel_meso61.rds ${REL2DIR}

#After fixing mpos 8 1000 5
#model for SNVs R1
singularity exec $CONTAINER Rscript ${MS} 8 1000 5 Somatics.snv.matrix.txt Germline.snv.matrix.txt rf-8_1000_5_snv_meso61.r1.rds ${REL2DIR}
#model for INDELS, inlcuding SNVs+INDELs
singularity exec $CONTAINER Rscript ${MI} 8 1000 5 Somatics.snv.matrix.txt Germline.snv.matrix.txt Somatics.indel.matrix.txt Germline.indel.matrix.txt rf-8_1000_5_indel_meso61.r1.rds ${REL2DIR}

#model for SNVs R2
singularity exec $CONTAINER Rscript ${MS} 8 1000 5 Somatics.snv.matrix.txt Germline.snv.matrix.txt rf-8_1000_5_snv_meso61.r2.rds ${REL2DIR}
#model for INDELS, inlcuding SNVs+INDELs
singularity exec $CONTAINER Rscript ${MI} 8 1000 5 Somatics.snv.matrix.txt Germline.snv.matrix.txt Somatics.indel.matrix.txt Germline.indel.matrix.txt rf-8_1000_5_indel_meso61.r2.rds ${REL2DIR}

#model for SNVs R3
singularity exec $CONTAINER Rscript ${MS} 8 1000 5 Somatics.snv.matrix.txt Germline.snv.matrix.txt rf-8_1000_5_snv_meso61.r3.rds ${REL2DIR}
#model for INDELS, inlcuding SNVs+INDELs
singularity exec $CONTAINER Rscript ${MI} 8 1000 5 Somatics.snv.matrix.txt Germline.snv.matrix.txt Somatics.indel.matrix.txt Germline.indel.matrix.txt rf-8_1000_5_indel_meso61.r3.rds ${REL2DIR}
