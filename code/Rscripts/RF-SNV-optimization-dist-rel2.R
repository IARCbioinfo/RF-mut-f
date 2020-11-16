#library(dplyr)
library(randomForest)
library(caret)
library(e1071)
library(caTools)

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("We expect : mtry ntree nodesize somatics germline working_directory", call.=FALSE)
}
#get the comand line parameters
#data.frame()
param_mtry=as.numeric(args[1])
param_ntree=as.numeric(args[2])
param_nodesize=as.numeric(args[3])  
param_somatics=args[4]
param_germline=args[5]
param_wdir=args[6]

print(paste("Tunning parameters:",param_mtry,param_ntree,param_nodesize))
#test1, with meso data
setwd(param_wdir)
#we load the somatics variants 
ds=read.table(param_somatics,h=T)
#table(ds$MPOS == ".")
#we handle the missing values of MPOS and replace them by median
vs=as.integer(ds[ds$MPOS != ".",]$MPOS)
#we compute the median of MPOS values
m_mpos_s=median(vs)
#we replace the "." by the median
ds$MPOS=replace(ds$MPOS,ds$MPOS==".",m_mpos_s)
print(paste0(m_mpos_s))

#summary(ds)
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

#summary(ds)
#germlines variants
dg=read.table(param_germline,h=T)
vg=as.integer(dg[dg$MPOS != ".",]$MPOS)
#table(dg$MPOS == ".")
#summary(va)
m_mpos_g=median(vg)
#we replace the missing MPOS values by median of germline variants
dg$MPOS=replace(dg$MPOS,dg$MPOS==".",m_mpos_g)
print(paste0(m_mpos_g))
#summary(dg)

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
#quit()
#sapply(dg, class)
#summary(dg)
#we sample randomly a number of germline variants from dgs (10% for demo)
dgs=dg[sample(nrow(dg), 1*dim(ds)*1),]
dss=ds[sample(nrow(ds), 1*dim(ds)*1),]
dim(dgs)
dim(dss)

# now we can merge both variant sets
meso_data=rbind(dss,dgs)
head(meso_data)
#we drop unused columnns (vars)
df <- subset(meso_data, select = -c(FILE, CHROM, POS,CODING))
#we split the data for testing and evaluation
sample = sample.split(df$SOMATIC, SplitRatio = .75)
train = subset(df, sample == TRUE) # 75% training
test  = subset(df, sample == FALSE) # 25% for evaluation
                                    # fitting datasset to avoid overfitting.
dim(train)
dim(test)

#we train the randomForest (def values)
rf <- randomForest(
  SOMATIC ~ .,
  data=train,
  mtry = param_mtry, 
  ntree=param_ntree, 
  nodesize=param_nodesize  
)
#we make the predictions
pred = predict(rf, newdata=test)
#we evaluate the prediction
confusionMatrix(pred, test$SOMATIC)
#we plot the relevant of variables
#varImpPlot(rf)
#
varImp(rf)
#code to optimize training
print(rf) #default tree
