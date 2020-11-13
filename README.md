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
 
 	
The script code/bash/get_mesomics_rel2.sh create the previous file structure.

## Training and optimizing a random forest classifier



### Random forest missing values

Typically, random forest methods/packages encourage two ways of handling missing values: a) drop data points with missing values (not recommended); b) fill in missing values with the median (for numerical values) or mode (for categorical values). For MPOS missing values we will use the median. [post](https://medium.com/airbnb-engineering/overcoming-missing-values-in-a-random-forest-classifier-7b1fc1fc03ba#:~:text=Typically%2C%20random%20forest%20methods%2Fpackages,mode%20(for%20categorical%20values))
 