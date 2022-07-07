# RF-mut-f
A random forest implementation to filter germline mutations from tumor-only samples


## Container
To ensure reproducibility, we created a docker/singularity container:

```
#fetch and create the singularity container
export TMPDIR=/tmp
singularity pull docker://adigenova/rf-mut-f:v2.0
```


## Directory structure

```
├── container						# container 
│   └── rf-mut-f_v2.0.sif	
├── databases						# external databases
│   ├── cosmic					
│   ├── genome
│   ├── gnomad
│   └── misc						# BED for centromers...
├── mesomics					   # mesomic data release
│   └── release2
└── README.md					   # main README.md	
```


## Mesomics 

The release2 directory contains three subdirectories:


1. **matched**
 
	Links to VCF files of matched tumors.
	
2. **matched-t-only**
 
   Links to VCF files of matched samples called as t-only.
   	  
3. **t-only**
  
 Links to VCF files of tumor-only samples.
 
 	
The script code/bash/get_mesomics_rel2.sh create the previous file structure.

## Step 1: Preparing the training data
This step creates a matrix of features for training the Random Forest Model to discriminate germline from somatic variants. The current list of features is the following:


| Feature           | type   | Description                                                            |
|:--------------------:|:----------:|-------------------------------------------------------------------------|
| COSMIC\_CENSUS\_GENE | factor   | var is in cosmic gene                                                   |
| COSMIC             | factor   | var is annotated in cosmic                                              |
| GNOMAD             | integer  | var is annotated in GNOMAD                                              |
| SNVS               | factor   | AT TC GA (Signatures)                                                   |
| BCSQ               | factor   | Variant impact                                                          |
| CENTROMER          | factor   | var is located in centromeric regions                                   |
| MPOS               | integer  | median distance from the end of read                                        |
| DP                 | integer  | Approximate read depth                                                  |
| GERMQ              | integer  | Phred-scaled quality that alt alleles are not germline variants       |
| SEQQ               | integer  | Phred-scaled quality that alt alleles are not sequencing errors       |
| STRANDQ            | integer  | Phred-scaled quality of strand bias artifact                            |
| TLOD               | integer  | Log 10 likelihood ratio score of variant existing versus not existing |
| NS                 | integer  | Number of samples with the var (Sample Freq)                            |
| AF                 | float    | Allele fractions of alternate alleles in the tumor                      |
| ADR                | integer  | Allelic depths for the ref allele                                       |
| ADA                | integer  | Allelic depths for the alt allele                                       |
| OR1                | integer  | Count of reads in F1R2 pair orientation supporting ref allele          |
| OR2                | integer  | Count of reads in F2R1 pair orientation supporting ref allele           |
| OA1                | integer  | Count of reads in F1R2 pair orientation supporting alt allele           |
| OA2                | integer  | Count of reads in F2R1 pair orientation supporting alt allele           |
| SOMATIC            | factor   | var is somatic or not                                                   |


A total of 20 features from 3 groups, including external databases (COSMIC, GNOMAD), genomic impact and composition, Mutect2 features related to sequencing errors, and read depth, among others, are used to build and train the RF classifier. 



### Code to build the feature matrices

The script **jobs/create\_matrix\_rel2.sh**, is a SLURM script that use the singularity container to create the following files: 

1. **Somatics.snv.matrix.txt**
    Matrix with somatic SNVs  
2. **Germline.snv.matrix.txt**
    Matrix with Germline SNVs 
3. **Somatics.indel.matrix.txt**
    Matrix with somatic INDELs 
4. **Germline.indel.matrix.txt** 
 	Matrix with Germline INDELs 
 	
These files are used for training and optimizing the RF classifier. 	
This job will allocate 50 CPUs in a single machine. 

```
sbatch jobs/create_matrix_rel2.sh
```

The above bash script call the ***code/makefiles/create_matrix_training.mk*** makefile that executes all the required within the container to generate the above files.

### Code to fetch database data

Under de ***databases***  directory, there are makefiles for fetching the required external databases.

## Step 2: Training and optimization of a random forest classifier

The bash script ***run-optimization-rel2.sh***  call the Rscript ***code/Rscripts/RF-SNV-optimization-dist-rel2.R*** to perform a grid search of the random forest model using a total of 48 parameter combinations (code/Rscripts/grid_search_parameters.txt) as described in Di Genova et al. 


## Step 3: Building the best RF model for SNVs and INDELs

The bash script ***jobs/train-best-model-rel2.sh*** call internally the ***code/Rscripts/RF-BEST-SNV-REL2.R*** to build the random forest models using the optimal parameters resulting from step 2. In the case of mesothelioma, they were 8 1000 5 for mtry, ntree, and node size, respectively. The resulting random forest models for SNPs and Indels are employed to classify tumor-only point mutations into germline or somatic.


## Step 4: Preparing and classifying tumor-only samples


The bash script ***jobs/predict_somatics_rel2.sh*** call the ***code/makefiles/predict_t_only.mk*** makefile to annotate the tumor-only variants, build the 20 features matrices, and finally use the optimal random forest to classify each variant of the tumor-only whole-genome sequence set into germline or somatics. The Rscript ***code/Rscripts/RF-APPLY-MODEL-CUTOFF-VOTES-PROB.R*** load the optimal models and performs variants classification.


## Step 5: Filtering tumor-only point mutations using RF votes


The bash script ***filter_vcfs_annovar_rel2.sh*** call the ***code/makefiles/filter_vcfs_using_model.mk*** makefile to filter tumor-only variants based on RF votes using a cuttoff of 0.75 and 0.5 for non-coding and coding variants,respectively.


After this final step, a set of VCFs with high-confidence somatic point mutations is generated for each tumor-only sample.

## Notes


### Missing MPOS values:
MPOS is a relevant variable for the model, the third one after Allele frequency and GERMQ. Some values of MPOS are set as  "." due to a bug of [Mutect2](https://github.com/broadinstitute/gatk/issues/6342). These missing values were replaced by the median MPOS of somatics (median=38) or germline (median=50) variants. A total of 152 and 346 variants were replaced by the median of somatic or germline variants, respectively.


### Random forest missing values

Typically, random forest methods/packages encourage two ways of handling missing values: a) drop data points with missing values (not recommended); b) fill in missing values with the median (for numerical values) or mode (for categorical values). For MPOS missing values, we will use the median. [post](https://medium.com/airbnb-engineering/overcoming-missing-values-in-a-random-forest-classifier-7b1fc1fc03ba#:~:text=Typically%2C%20random%20forest%20methods%2Fpackages,mode%20(for%20categorical%20values))

