.DELETE_ON_ERROR:


ROOT_D=/data/scratch/digenovaa/Somatic-reference-free/SNV-INDELs/RF-mut-f
RF_SNV=${ROOT_D}/mesomics/release2/t-only/rf_snv_m61.r1_predictions.txt
RF_INDEL=${ROOT_D}/mesomics/release2/t-only/rf_indels_m61.r1_predictions.txt

#code that filter the VCFs or annovar files using the random forest
PHOME=${ROOT_D}/code/perl

#we merge indels and snvs predictions
rf_m61.r1_predictions.txt:
	cat ${RF_SNV} ${RF_INDEL} > rf_m61.r1_predictions.txt

filter_m61=$(patsubst %_filtered_PASS_norm.vcf.hg38_multianno.vcf.gz,%.rf_rel2_m61.vcf.gz,$(wildcard *multianno.vcf.gz))
%.rf_rel2_m61.vcf.gz:%_filtered_PASS_norm.vcf.hg38_multianno.vcf.gz rf_m61.r1_predictions.txt
	perl ${PHOME}/filter-vcf.pl -a rf_m61.r1_predictions.txt -b $(subst _filtered_PASS_norm.vcf.hg38_multianno.vcf.gz,,$<) -c $< > $(subst .gz,,$@)
	bgzip -i $(subst .gz,,$@)
filter_meso61:$(filter_m61)

#annovar files
filter_an_m61=$(patsubst %_filtered_PASS_norm.vcf.hg38_multianno.txt,%.rf_rel2_m61.hg38_multianno.txt,$(wildcard *vcf.hg38_multianno.txt))
%.rf_rel2_m61.hg38_multianno.txt:%_filtered_PASS_norm.vcf.hg38_multianno.txt rf_m61.r1_predictions.txt
	perl ${PHOME}/filter-multianno.pl -a rf_m61.r1_predictions.txt -b $(subst _filtered_PASS_norm.vcf.hg38_multianno.txt,,$<) -c $< > $@
filter_an_meso61:$(filter_an_m61)

all: rf_m61.r1_predictions.txt filter_meso61 filter_an_meso61

