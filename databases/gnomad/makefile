.DELETE_ON_ERROR:

URL=https://storage.googleapis.com/gnomad-public/release/3.0/vcf/genomes

#wget -qO- https://storage.googleapis.com/gnomad-public/release/3.0/vcf/genomes/gnomad.genomes.r3.0.sites.vcf.bgz | bgzip -dc | awk 'BEGIN{OFS="\t"}{if($0 ~/#/){print $0}else{split($8,a,";"); $8=""; print $0"\t"a[1]";"a[2]";"a[3]";"a[5]}}' > gnomad.genomes.r3.0.sites.vcf

CHRS=1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y
FILES=$(addprefix gnomad.genomes.r3.0.sites.chr,$(addsuffix .vcf.bgz,$(CHRS)))
#VCF:=$(addprefix gnomad.genomes.r3.0.sites.chr,$(addsuffix .vcf.bgz,$(CHRS)))
##"gnomad.genomes.r3.0.sites.chr"$1".vcf.bgz"
#vcf_files=$(patsubst %.vcf.bgz,%.vcf,$(FILES))
vcf_files=$(FILES)
#we donwload the variants and keep only those that pass all filters
%.vcf.bgz:
	wget -qO- $(URL)/$@ | bgzip -dc | awk 'BEGIN{OFS="\t"}{if($$0 ~/#/){print $$0}else{if($$0 ~/PASS/){split($$8,a,";"); $$8=""; print $$0""a[1]";"a[2]";"a[3]";"a[4]";"a[5]}}}' | bgzip >  $@
download: $(vcf_files)

#vcf_pass=$(patsubst %.vcf.bgz,%.pass.vcf.bgz,$(FILES))

#%.pass.vcf.bgz:%.vcf.bgz
#	zcat $< | awk 'BEGIN{OFS="\t";}{if($$0 ~/#/){print $$0}else{if($$0 ~/PASS/){print $$0}}}' | bgzip > $@
#filter_vars:$(vcf_pass)

#we compute stats
BCFTOOLS=/home/digenovaa/binaries/bcftools-1.10.2/BCFT/bin/bcftools
TABIX=/home/digenovaa/miniconda3/bin/tabix
vcf_stats=$(patsubst %.vcf.bgz,%.stats,$(vcf_files))
%.stats:%.vcf.bgz
	${TABIX} $<
	${BCFTOOLS} stats  $< > $@
stats:$(vcf_stats)
#we merge all the chromosomes
gnomad.genome.vcf.bgz: stats
	${BCFTOOLS} concat -o $@ -O z -n ${FILES}
	${TABIX} $@
	 ${BCFTOOLS} stats $@ > $@.stats
merge:gnomad.genome.vcf.bgz
#all: download filter_vars stats
all: download  stats merge
