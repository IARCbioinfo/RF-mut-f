# RF-mut-f
A random forest implementation to filter germline mutations from tumor-only samples

## Container
To ensure reproducibility, we created a docker/singularity container:

```
#fetch and create the singularity container
export TMPDIR=/tmp
singularity pull docker://adigenova/rf-mut-f:v1.0
```


## Directory structure

```
├── container						# container
│   └── rf-mut-f_v1.0.sif	
├── databases						# external databases
│   ├── cosmic					
│   ├── genome
│   ├── gnomad
│   └── misc						# BED for centromers...
├── mesomics					   # mesomic data release
│   └── release2
└── README.md					   # main README.md	
```


## Mesomics (rel2)

The release2 directory contains 3 subdirectories:


1. **matched**
 
	links to 46 VCFs files of matched tumors.
	
2. **matched-t-only**
 
   links to 46 VCFs files of matched samples called as t-only.
   	  
3. **t-only**
  
 links to 73 VCFs files of tumor-only samples.
 
 	
The below code create the previous file structure (get_mesomics_rel2.sh):

```
#create mesomics rel2 dataset
mkdir -p mesomics/release2
mkdir -p mesomics/release2/matched
mkdir -p mesomics/release2/matched-t-only
mkdir -p mesomics/release2/t-only
#matched data (46)
RELEASE2=/data/gcs/mesomics/files/WGS/variant_calling/somatic_release2_26032020/intermediate_files/normalized_calling
ln -s ${RELEASE2}/Mutect2_somatic_pon_blood_TN_gatk4150_normalized_20200324/*_norm.vcf.gz $PWD/mesomics/release2/matched        
ln -s ${RELEASE2}/Mutect2_somatic_pon_blood_TN_gatk4150_normalized_20200324/*_norm.vcf.gz.tbi $PWD/mesomics/release2/matched
#matched-t-only data (46)
ln -s ${RELEASE2}/Mutect2_somatic_pon_blood_TN_Tonly_mode_gatk4150_normalized_20200324/*_norm.vcf.gz $PWD/mesomics/release2/matched-t-only
ln -s ${RELEASE2}/Mutect2_somatic_pon_blood_TN_Tonly_mode_gatk4150_normalized_20200324/*_norm.vcf.gz.tbi $PWD/mesomics/release2/matched-t-only
#t-only data (73)
ln -s ${RELEASE2}/Mutect2_somatic_pon_blood_Tonly_gatk4150_normalized_20200325/*_norm.vcf.gz $PWD/mesomics/release2/t-only
ln -s ${RELEASE2}/Mutect2_somatic_pon_blood_Tonly_gatk4150_normalized_20200325/*_norm.vcf.gz.tbi $PWD/mesomics/release2/t-only
#we remove the MESO_094 calls
rm -f  $PWD/mesomics/release2/t-only/MESO_094_T_filtered_PASS_norm.vcf.gz
rm -f  $PWD/mesomics/release2/t-only/MESO_094_T_filtered_PASS_norm.vcf.gz.tbi
#better calling for MESO_094
MESO_094=/data/gcs/mesomics/files/WGS/variant_calling/somatic_release2_MESO_094_28102020/intermediate_files/Mutect2-nf_results_MESOMICS_26102020_pon_blood_Tonly_gatk4150_MESO_094_normalized
ln -s ${MESO_094}/MESO_094_filtered_PASS_norm.vcf.gz $PWD/mesomics/release2/t-only/MESO_094_T_filtered_PASS_norm.vcf.gz
ln -s ${MESO_094}/MESO_094_filtered_PASS_norm.vcf.gz $PWD/mesomics/release2/t-only/MESO_094_T_filtered_PASS_norm.vcf.gz.tbi
```	

