#!/bin/bash
#SBATCH --job-name="RFT"
#SBATCH --partition=high_p
#SBATCH -c 1
#SBATCH --mem-per-cpu=30000

## conda path
##we use the base code and containers
ROOT=/data/scratch/digenovaa/Somatic-reference-free/SNV-INDELs/RF-mut-f
#container
CONTAINER=${ROOT}/container/rf-mut-f_v2.0.sif

#apply SNV model

# We apply the models to the T-only data from tumors
#snvs
singularity exec $CONTAINER  Rscript ${ROOT}/code/Rscripts/RF-APPLY-MODEL-CUTOFF-VOTES-PROB.R -i Mutations.snv.matrix.txt -o ${PWD} -m rf-8_1000_5_snv_meso61.r1.rds  -s rf-8_1000_5_snv_meso61.r1
singularity exec $CONTAINER  Rscript ${ROOT}/code/Rscripts/RF-APPLY-MODEL-CUTOFF-VOTES-PROB.R -i Mutations.snv.matrix.txt -o ${PWD} -m rf-8_1000_5_snv_meso61.r2.rds  -s rf-8_1000_5_snv_meso61.r2
singularity exec $CONTAINER  Rscript ${ROOT}/code/Rscripts/RF-APPLY-MODEL-CUTOFF-VOTES-PROB.R -i Mutations.snv.matrix.txt -o ${PWD} -m rf-8_1000_5_snv_meso61.r3.rds  -s rf-8_1000_5_snv_meso61.r3
#indels
singularity exec $CONTAINER  Rscript ${ROOT}/code/Rscripts/RF-APPLY-MODEL-CUTOFF-VOTES-PROB.R -i Mutations.indel.matrix.txt -o ${PWD} -m rf-8_1000_5_indel_meso61.r1.rds  -s rf-8_1000_5_indel_meso61.r1
singularity exec $CONTAINER  Rscript ${ROOT}/code/Rscripts/RF-APPLY-MODEL-CUTOFF-VOTES-PROB.R -i Mutations.indel.matrix.txt -o ${PWD} -m rf-8_1000_5_indel_meso61.r2.rds  -s rf-8_1000_5_indel_meso61.r2
singularity exec $CONTAINER  Rscript ${ROOT}/code/Rscripts/RF-APPLY-MODEL-CUTOFF-VOTES-PROB.R -i Mutations.indel.matrix.txt -o ${PWD} -m rf-8_1000_5_indel_meso61.r3.rds  -s rf-8_1000_5_indel_meso61.r3


#first model rel2
singularity exec $CONTAINER  Rscript ${ROOT}/code/Rscripts/RF-APPLY-MODEL-CUTOFF-VOTES-PROB.R -i Mutations.snv.matrix.txt -o ${PWD} -m rel2_rf-12_1500_25_snv_meso61.rds -s rel2_rf_snv_m61_1
#first model rel2
singularity exec $CONTAINER Rscript ${ROOT}/code/Rscripts/RF-APPLY-MODEL-CUTOFF-VOTES-PROB.R -i Mutations.indel.matrix.txt -o ${PWD} -m rel2_rf-12_1500_25_indel_meso61.rds -s rel2_rf_indels_m61_1
