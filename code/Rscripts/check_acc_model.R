#we load the current model
library(randomForest)
library(caret)
library(e1071)
library(caTools)
#we load the env for training model 8 1000 5
load("rf-8_1000_5_snv_meso61.r1.rds.RData")
RF <- readRDS("rf-8_1000_5_snv_meso61.r1.rds")
#total vars
dim(dg)
dim(ds)
#vars used in training
dim(dss)
dim(dgs)
#split in codding and non_coding
dg_coding=dg[dg$CODING==1,]
dg_non_coding=dg[dg$CODING==0,]
dim(dg_coding)
dim(dg_non_coding)

ds_coding=ds[ds$CODING==1,]
ds_non_coding=ds[ds$CODING==0,]

dim(ds_coding)
dim(ds_non_coding)

#coding
s_coding_somatic=table(ds_coding$FILE)
s_coding_germline=table(dg_coding$FILE)
a=rbind(s_coding_somatic,s_coding_germline)
a=t(a)
a[,1]<-ifelse(a[, 1] == 0, 1, a[,1])
summary(a[,2]/a[,1])
#total proportion
sum(a[,2])/sum(a[,1])
#non coding
s_nc_somatic=table(ds_non_coding$FILE)
s_nc_germline=table(dg_non_coding$FILE)
a=rbind(s_nc_somatic,s_nc_germline)
a=t(a)
a[,1]<-ifelse(a[, 1] == 0, 1, a[,1])
summary(a[,2]/a[,1])
#total proportion
sum(a[,2])/sum(a[,1])
#3.402437


#somatics not used in training 
nts=!(rownames(ds) %in% rownames(dss))
notrain_somatics=ds[nts,]
#germline not used in training
ntg=!(rownames(dg) %in% rownames(dgs))
notrain_germline=dg[ntg,]

#training data for model rf 8 1000 5
summary(df$SOMATIC)
#  1      0 
#163193 163193 
#we pick the coding variants not used in training for test the model
notrain_somatic_coding=notrain_somatics[notrain_somatics$CODING == 1, ]
notrain_germline_coding=notrain_germline[notrain_germline$CODING == 1,]
#we create the test  dataset for coding
cdgs=notrain_germline_coding[sample(nrow(notrain_germline_coding), 1*dim(notrain_somatic_coding)*3.5),]
cdss=notrain_somatic_coding[sample(nrow(notrain_somatic_coding), 1*dim(notrain_somatic_coding)*1),]
dim(cdgs)
dim(cdss)
# now we can merge both variant sets
meso_data_coding=rbind(cdss,cdgs)
summary(meso_data_coding$SOMATIC)
summary(meso_data_coding$CODING)
#  1   0 
#489 489 
dfc <- subset(meso_data_coding, select = -c(FILE, CHROM, POS,CODING))
dim(dfc)
#[1] 978  21
#pred = predict(RF, newdata=dfc)
#confusionMatrix(pred, dfc$SOMATIC)

pred.prob = predict(RF,newdata=dfc, type="prob")

pred.c<-ifelse(pred.prob[,1]>0.5, 1, 0)
pred.c=as.character(pred.c)
pred.c <- factor(pred.c, levels = levels(dfc$SOMATIC))
confusionMatrix(pred.c, dfc$SOMATIC)

pred.c<-ifelse(pred.prob[,1]>0.6, 1, 0)
pred.c=as.character(pred.c)
pred.c <- factor(pred.c, levels = levels(dfc$SOMATIC))
confusionMatrix(pred.c, dfc$SOMATIC)

pred.c<-ifelse(pred.prob[,1]>0.7, 1, 0)
pred.c=as.character(pred.c)
pred.c <- factor(pred.c, levels = levels(dfc$SOMATIC))
confusionMatrix(pred.c, dfc$SOMATIC)

pred.c<-ifelse(pred.prob[,1]>0.75, 1, 0)
pred.c=as.character(pred.c)
pred.c <- factor(pred.c, levels = levels(dfc$SOMATIC))
confusionMatrix(pred.c, dfc$SOMATIC)


pred.c<-ifelse(pred.prob[,1]>0.8, 1, 0)
pred.c=as.character(pred.c)
pred.c <- factor(pred.c, levels = levels(dfc$SOMATIC))
confusionMatrix(pred.c, dfc$SOMATIC)

pred.c<-ifelse(pred.prob[,1]>0.9, 1, 0)
pred.c=as.character(pred.c)
pred.c <- factor(pred.c, levels = levels(dfc$SOMATIC))
confusionMatrix(pred.c, dfc$SOMATIC)

#NO MESO_61

notrain_somatic_coding_no61=notrain_somatic_coding[notrain_somatic_coding$FILE != "MESO_061_filtered_PASS_norm.somatics.vcf.bgz", ]
notrain_germline_coding_no61=notrain_germline_coding[notrain_germline_coding$FILE != "MESO_061_filtered_PASS_norm.somatics.vcf.bgz",]


#we create the test  dataset for coding for meso61g
cdgsn61=notrain_germline_coding_no61[sample(nrow(notrain_germline_coding_no61), 1*dim(notrain_somatic_coding_no61)*1),]
cdssn61=notrain_somatic_coding_no61[sample(nrow(notrain_somatic_coding_no61), 1*dim(notrain_somatic_coding_no61)*1),]
dim(cdgsn61)
dim(cdssn61)

meso_data_coding_n61=rbind(cdssn61,cdgsn61)
summary(meso_data_coding_n61$SOMATIC)
summary(meso_data_coding_n61$CODING)
#  1   0 
#489 489 
dfcn61 <- subset(meso_data_coding_n61, select = -c(FILE, CHROM, POS,CODING))
dim(dfcn61)

pred.prob = predict(RF,newdata=dfcn61, type="prob")

pred.c<-ifelse(pred.prob[,1]>0.5, 1, 0)
pred.c=as.character(pred.c)
pred.c <- factor(pred.c, levels = levels(dfcn61$SOMATIC))
confusionMatrix(pred.c, dfcn61$SOMATIC)

pred.c<-ifelse(pred.prob[,1]>0.6, 1, 0)
pred.c=as.character(pred.c)
pred.c <- factor(pred.c, levels = levels(dfcn61$SOMATIC))
confusionMatrix(pred.c, dfcn61$SOMATIC)

pred.c<-ifelse(pred.prob[,1]>0.7, 1, 0)
pred.c=as.character(pred.c)
pred.c <- factor(pred.c, levels = levels(dfcn61$SOMATIC))
confusionMatrix(pred.c, dfcn61$SOMATIC)

pred.c<-ifelse(pred.prob[,1]>0.75, 1, 0)
pred.c=as.character(pred.c)
pred.c <- factor(pred.c, levels = levels(dfcn61$SOMATIC))
confusionMatrix(pred.c, dfcn61$SOMATIC)

pred.c<-ifelse(pred.prob[,1]>0.8, 1, 0)
pred.c=as.character(pred.c)
pred.c <- factor(pred.c, levels = levels(dfcn61$SOMATIC))
confusionMatrix(pred.c, dfcn61$SOMATIC)

pred.c<-ifelse(pred.prob[,1]>0.9, 1, 0)
pred.c=as.character(pred.c)
pred.c <- factor(pred.c, levels = levels(dfcn61$SOMATIC))
confusionMatrix(pred.c, dfcn61$SOMATIC)


#non-coding data not used in the training step 
notrain_somatics_nc=notrain_somatics[notrain_somatics$CODING==0,]
notrain_germline_nc=notrain_germline[notrain_germline$CODING==0,]

#we create the test  dataset for non-coding
ndgs=notrain_germline_nc[sample(nrow(notrain_germline_nc), 1*dim(notrain_somatics_nc)*3.5),]
ndss=notrain_somatics_nc[sample(nrow(notrain_somatics_nc), 1*dim(notrain_somatics_nc)*1),]
dim(ndgs)
dim(ndss)
# now we can merge both variant sets
meso_data_nc=rbind(ndss,ndgs)
summary(meso_data_nc$SOMATIC)
#    1     0 
#40310 40310 
summary(meso_data_nc$CODING)
#    0     1 
#80620     0 

dfnc <- subset(meso_data_nc, select = -c(FILE, CHROM, POS,CODING))
dim(dfnc)
#pred = predict(RF, newdata=dfnc)
#confusionMatrix(pred, dfnc$SOMATIC)
#1     0 
#40310 40310 
pred.prob = predict(RF,newdata=dfnc, type="prob")
#we transform the probalistic filter to Somatic or Germline
pred.c<-ifelse(pred.prob[,1]>0.5, 1, 0)
pred.c=as.character(pred.c)
pred.c <- factor(pred.c, levels = levels(dfnc$SOMATIC))
confusionMatrix(pred.c, dfnc$SOMATIC)

pred.c<-ifelse(pred.prob[,1]>0.6, 1, 0)
pred.c=as.character(pred.c)
pred.c <- factor(pred.c, levels = levels(dfnc$SOMATIC))
confusionMatrix(pred.c, dfnc$SOMATIC)
    
pred.c<-ifelse(pred.prob[,1]>0.7, 1, 0)
pred.c=as.character(pred.c)
pred.c <- factor(pred.c, levels = levels(dfnc$SOMATIC))
confusionMatrix(pred.c, dfnc$SOMATIC)
    
pred.c<-ifelse(pred.prob[,1]>0.75, 1, 0)
pred.c=as.character(pred.c)
pred.c <- factor(pred.c, levels = levels(dfnc$SOMATIC))
confusionMatrix(pred.c, dfnc$SOMATIC)
    
pred.c<-ifelse(pred.prob[,1]>0.8, 1, 0)
pred.c=as.character(pred.c)
pred.c <- factor(pred.c, levels = levels(dfnc$SOMATIC))
confusionMatrix(pred.c, dfnc$SOMATIC)

pred.c<-ifelse(pred.prob[,1]>0.9, 1, 0)
pred.c=as.character(pred.c)
pred.c <- factor(pred.c, levels = levels(dfnc$SOMATIC))
confusionMatrix(pred.c, dfnc$SOMATIC)





 
