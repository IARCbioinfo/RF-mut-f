.DELETE_ON_ERROR:


ROOT_D=/data/scratch/digenovaa/Somatic-reference-free/SNV-INDELs/RF-mut-f
RF_SNV=${ROOT_D}/mesomics/release2/t-only/rf_snv_uncollapsed_m61.r1_predictions.txt
RF_INDEL=${ROOT_D}/mesomics/release2/t-only/rf_indel_uncollapsed_m61.r1_predictions.txt

#code that filter the VCFs or annovar files using the random forest
PHOME=${ROOT_D}/code/perl

filter_m61=$(patsubst %_filtered_PASS_norm.vcf.hg38_multianno.vcf.gz,%.rf_rel2_m61.vcf.gz,$(wildcard *multianno.vcf.gz))
%.rf_rel2_m61.vcf.gz:%_filtered_PASS_norm.vcf.hg38_multianno.vcf.gz 
	perl ${PHOME}/filter-vcf.pl -a ${RF_SNV} -b ${RF_INDEL} -d $< -s $(subst _filtered_PASS_norm.vcf.hg38_multianno.vcf.gz,,$<) -c 0.75 > $(subst .vcf.gz,.hc.vcf,$@) 2> $(subst .vcf.gz,.hc.vcf,$@).err 
	perl ${PHOME}/filter-vcf.pl -a ${RF_SNV} -b ${RF_INDEL} -d $< -s $(subst _filtered_PASS_norm.vcf.hg38_multianno.vcf.gz,,$<) > $(subst .gz,,$@) 2> $(subst .gz,,$@).err
	bgzip -i $(subst .vcf.gz,.hc.vcf,$@)
	bgzip -i $(subst .gz,,$@)
filter_meso61:$(filter_m61)


#annovar files
filter_an_m61=$(patsubst %_filtered_PASS_norm.vcf.hg38_multianno.txt,%.rf_rel2_m61.hg38_multianno.txt,$(wildcard *vcf.hg38_multianno.txt))
%.rf_rel2_m61.hg38_multianno.txt:%_filtered_PASS_norm.vcf.hg38_multianno.txt
	perl ${PHOME}/filter-multianno.pl -a ${RF_SNV} -b ${RF_INDEL} -d $< -s $(subst _filtered_PASS_norm.vcf.hg38_multianno.txt,,$<) -c 0.75 > $(subst .txt,.hc.txt,$@) 2>$(subst .txt,.hc.txt,$@).err
	perl ${PHOME}/filter-multianno.pl -a ${RF_SNV} -b ${RF_INDEL} -d $< -s $(subst _filtered_PASS_norm.vcf.hg38_multianno.txt,,$<) > $@ 2>$@.err

filter_an_meso61:$(filter_an_m61)

all: filter_meso61 filter_an_meso61

