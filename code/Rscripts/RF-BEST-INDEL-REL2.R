#library(dplyr)
library(randomForest)
library(caret)
library(e1071)
library(caTools)


args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}
#get the comand line parameters
#data.frame()
param_mtry=as.numeric(args[1])
param_ntree=as.numeric(args[2])
param_nodesize=as.numeric(args[3])  
param_snv_somatics=args[4]
param_snv_germline=args[5]
param_indel_somatics=args[6]
param_indel_germline=args[7]
param_model_out=args[8]
param_wdir=args[9]

print(paste("Tunning parameters:",param_mtry,param_ntree,param_nodesize))
#test1, with meso data
setwd(param_wdir)
#we load the somatics variants 
#ds=read.table("sample/matched/SOMATICS_SNPS.txt",h=T)
ds=read.table(param_snv_somatics,h=T)
vs=as.integer(ds[ds$MPOS != ".",]$MPOS)
#we compute the median of MPOS values
m_mpos_s=median(vs)
#we replace the "." by the median
ds$MPOS=replace(ds$MPOS,ds$MPOS==".",m_mpos_s)
print(paste0(m_mpos_s))




#we transform the variables to factor/integers when needed
ds <- transform(
  ds,
  FILE=as.factor(FILE), #sample name
  CHROM=as.factor(CHROM), #Chromosome
  POS=as.integer(POS), # Position
  SNVS=as.factor(SNVS), # join(REFATL) AT TC GA (Signatures)
  BCSQ=as.factor(BCSQ), #Variant impact 
  CODING=as.factor(CODING),#variant is coding
  COSMIC_CENSUS_GENE=as.factor(COSMIC_CENSUS_GENE), #var is in cosmic gene
  COSMIC=as.factor(COSMIC), #var is annotated in cosmic
  GNOMAD=as.integer(GNOMAD), # var is annotated in GNOMAD
  CENTROMER=as.factor(CENTROMER), # var is located in centromeric regions
  MPOS=as.integer(MPOS), # median distance from end of read
  DP=as.integer(DP), # Approximate read depth
  GERMQ=as.integer(GERMQ), # Phred-scaled quality that alt alleles are not germline variants
  SEQQ=as.integer(SEQQ), # Phred-scaled quality that alt alleles are not sequencing errors
  STRANDQ=as.integer(STRANDQ), # Phred-scaled quality of strand bias artifact
  TLOD=as.integer(TLOD), # Log 10 likelihood ratio score of variant existing versus not existing
  AF=AF,
  ADR=as.integer(ADR),
  ADA=as.integer(ADA),
  OR1=as.integer(OR1),
  OR2=as.integer(OR2),
  OA1=as.integer(OA1),
  OA2=as.integer(OA2),
  NS=as.integer(NS), # number of samples with the var
  SOMATIC=as.factor(SOMATIC) #var is somatic or not
)
#sapply(ds, class)

#germlines variants
dg=read.table(param_snv_germline,h=T)

vs=as.integer(dg[dg$MPOS != ".",]$MPOS)
m_mpos_s=median(vs)
dg$MPOS=replace(dg$MPOS,dg$MPOS==".",m_mpos_s)
print(paste0(m_mpos_s))



dg <- transform(
  dg,
  FILE=as.factor(FILE), #sample name
  CHROM=as.factor(CHROM), #Chromosome
  POS=as.integer(POS), # Position
  SNVS=as.factor(SNVS), # join(REFATL) AT TC GA (Signatures)
  BCSQ=as.factor(BCSQ), #Variant impact 
  CODING=as.factor(CODING),#variant is coding
  COSMIC_CENSUS_GENE=as.factor(COSMIC_CENSUS_GENE), #var is in cosmic gene
  COSMIC=as.factor(COSMIC), #var is annotated in cosmic
  GNOMAD=as.integer(GNOMAD), # var is annotated in GNOMAD
  CENTROMER=as.factor(CENTROMER), # var is located in centromeric regions
  MPOS=as.integer(MPOS), # median distance from end of read
  DP=as.integer(DP), # Approximate read depth
  GERMQ=as.integer(GERMQ), # Phred-scaled quality that alt alleles are not germline variants
  SEQQ=as.integer(SEQQ), # Phred-scaled quality that alt alleles are not sequencing errors
  STRANDQ=as.integer(STRANDQ), # Phred-scaled quality of strand bias artifact
  TLOD=as.integer(TLOD), # Log 10 likelihood ratio score of variant existing versus not existing
  AF=AF,
  ADR=as.integer(ADR),
  ADA=as.integer(ADA),
  OR1=as.integer(OR1),
  OR2=as.integer(OR2),
  OA1=as.integer(OA1),
  OA2=as.integer(OA2),
  NS=as.integer(NS), # number of samples with the var
  SOMATIC=as.factor(SOMATIC) #var is somatic or not
)

#we load indels
dsi=read.table(param_indel_somatics,h=T)
#we replace the missing mpos values by the median
vs=as.integer(dsi[dsi$MPOS != ".",]$MPOS)
m_mpos_s=median(vs)
dsi$MPOS=replace(dsi$MPOS,dsi$MPOS==".",m_mpos_s)
print(paste0(m_mpos_s))

#we transform the variables to factor/integers when needed
dsi <- transform(
  dsi,
  FILE=as.factor(FILE), #sample name
  CHROM=as.factor(CHROM), #Chromosome
  POS=as.integer(POS), # Position
  SNVS=as.factor(SNVS), # join(REFATL) AT TC GA (Signatures)
  BCSQ=as.factor(BCSQ), #Variant impact 
  CODING=as.factor(CODING),#variant is coding
  COSMIC_CENSUS_GENE=as.factor(COSMIC_CENSUS_GENE), #var is in cosmic gene
  COSMIC=as.factor(COSMIC), #var is annotated in cosmic
  GNOMAD=as.integer(GNOMAD), # var is annotated in GNOMAD
  CENTROMER=as.factor(CENTROMER), # var is located in centromeric regions
  MPOS=as.integer(MPOS), # median distance from end of read
  DP=as.integer(DP), # Approximate read depth
  GERMQ=as.integer(GERMQ), # Phred-scaled quality that alt alleles are not germline variants
  SEQQ=as.integer(SEQQ), # Phred-scaled quality that alt alleles are not sequencing errors
  STRANDQ=as.integer(STRANDQ), # Phred-scaled quality of strand bias artifact
  TLOD=as.integer(TLOD), # Log 10 likelihood ratio score of variant existing versus not existing
  AF=AF,
  ADR=as.integer(ADR),
  ADA=as.integer(ADA),
  OR1=as.integer(OR1),
  OR2=as.integer(OR2),
  OA1=as.integer(OA1),
  OA2=as.integer(OA2),
  NS=as.integer(NS), # number of samples with the var
  SOMATIC=as.factor(SOMATIC) #var is somatic or not
)
#sapply(ds, class)

#germlines variants
dgi=read.table(param_indel_germline,h=T)

#we replace the missing mpos values by the median
vs=as.integer(dgi[dsi$MPOS != ".",]$MPOS)
m_mpos_s=median(vs)
dgi$MPOS=replace(dgi$MPOS,dgi$MPOS==".",m_mpos_s)
print(paste0(m_mpos_s))

dgi <- transform(
  dgi,
  FILE=as.factor(FILE), #sample name
  CHROM=as.factor(CHROM), #Chromosome
  POS=as.integer(POS), # Position
  SNVS=as.factor(SNVS), # join(REFATL) AT TC GA (Signatures)
  BCSQ=as.factor(BCSQ), #Variant impact 
  CODING=as.factor(CODING),#variant is coding
  COSMIC_CENSUS_GENE=as.factor(COSMIC_CENSUS_GENE), #var is in cosmic gene
  COSMIC=as.factor(COSMIC), #var is annotated in cosmic
  GNOMAD=as.integer(GNOMAD), # var is annotated in GNOMAD
  CENTROMER=as.factor(CENTROMER), # var is located in centromeric regions
  MPOS=as.integer(MPOS), # median distance from end of read
  DP=as.integer(DP), # Approximate read depth
  GERMQ=as.integer(GERMQ), # Phred-scaled quality that alt alleles are not germline variants
  SEQQ=as.integer(SEQQ), # Phred-scaled quality that alt alleles are not sequencing errors
  STRANDQ=as.integer(STRANDQ), # Phred-scaled quality of strand bias artifact
  TLOD=as.integer(TLOD), # Log 10 likelihood ratio score of variant existing versus not existing
  AF=AF,
  ADR=as.integer(ADR),
  ADA=as.integer(ADA),
  OR1=as.integer(OR1),
  OR2=as.integer(OR2),
  OA1=as.integer(OA1),
  OA2=as.integer(OA2),
  NS=as.integer(NS), # number of samples with the var
  SOMATIC=as.factor(SOMATIC) #var is somatic or not
)

#sapply(dg, class)
#summary(dg)
#we sample randomly a number of germline variants from dgs (10% for demo)
dgs=dg[sample(nrow(dg), 0.75*dim(ds)*1),]
dss=ds[sample(nrow(ds), 0.75*dim(ds)*1),]
dim(dgs)
dim(dss)

dgsi=dgi[sample(nrow(dgi), 1*dim(dsi)*1),]
dssi=dsi[sample(nrow(dsi), 1*dim(dsi)*1),]
dim(dgsi)
dim(dssi)

# now we can merge both variant sets
meso_data=rbind(dss,dgs,dgsi,dssi)
head(meso_data)
#we drop unused columnns (vars)
df <- subset(meso_data, select = -c(FILE, CHROM, POS,CODING,SNVS))
#we train the randomForest (def values)
rf <- randomForest(
  SOMATIC ~ .,
  data=df,
  mtry = param_mtry, 
  ntree=param_ntree, 
  nodesize=param_nodesize,
  keep.forest = TRUE
)
varImp(rf)
#code to optimize training
print(rf) #default tree
cat('Model building complete. Saving model...\n')
saveRDS(rf, file=param_model_out)
cat('Model saved.\n')
