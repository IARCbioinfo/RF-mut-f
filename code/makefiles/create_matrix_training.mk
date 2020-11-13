.DELETE_ON_ERROR:

#cromosomes considerated

CHRS=chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM
#directory with code and databases
ROOT_D=/data/scratch/digenovaa/Somatic-reference-free/SNV-INDELs/RF-mut-f
#we compute stats
BCFTOOLS=bcftools
TABIX=tabix

#annotate gnomad
GENOMAD_DB=${ROOT_D}/databases/gnomad/gnomad.genome.vcf.bgz
gnomad_annot_files=$(patsubst %.vcf.gz,%.gnomad.vcf.bgz,$(wildcard *.vcf.gz))
%.gnomad.vcf.bgz:%.vcf.gz
	${BCFTOOLS} annotate -r ${CHRS} -O z -o $@  -c "GNOMAD_AC:=AC,GNOMAD_AN:=AN" -m -NO_GNOMAD -a ${GENOMAD_DB} $<
	tabix $@
gnomad_annot:$(gnomad_annot_files)

#annotate cosmic vars
COSMIC_DB=${ROOT_D}/databases/cosmic/CosmicSNV_INDELS.gnomad.vcf.bgz
cosmic_annot_files=$(patsubst %.gnomad.vcf.bgz,%.cosmic.vcf.bgz,$(gnomad_annot_files))
%.cosmic.vcf.bgz:%.gnomad.vcf.bgz
	${BCFTOOLS} annotate -r ${CHRS} -O z -o $@ -c ID -m +COSMIC -a ${COSMIC_DB} $<
	tabix $@
comic_annot:$(cosmic_annot_files)

#annotate cosmic genes and centromers
add_centro_cgenes=$(patsubst %.cosmic.vcf.bgz,%.centro_cgenes.vcf.bgz,$(cosmic_annot_files))
%.centro_cgenes.vcf.bgz:%.cosmic.vcf.bgz
	${BCFTOOLS} annotate -a ${ROOT_D}/databases/misc/centromers.bed.gz -h ${ROOT_D}/databases/misc/centromers.hdr -c CNAME --columns CHROM,FROM,TO,CNAME -m +CENTROMER $< | \
	${BCFTOOLS} annotate -a ${ROOT_D}/databases/misc/cosmic_gene_census.bed.gz -h ${ROOT_D}/databases/misc/cosmic_gene_census.hdr -c CGENE --columns CHROM,FROM,TO,CGENE -O z -o $@ -m +COSMIC_CENSUS_GENE -
	${TABIX} $@
centro_genes:$(add_centro_cgenes)

#annotate variants effects
RFASTA=${ROOT_D}/databases/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa
RGFF3=${ROOT_D}/databases/genome/Homo_sapiens.GRCh38.99.chr.gff3.gz

add_effect=$(patsubst %.centro_cgenes.vcf.bgz,%.csq.vcf.bgz,$(add_centro_cgenes))
%.csq.vcf.bgz: %.centro_cgenes.vcf.bgz
	${BCFTOOLS} annotate --rename-chrs ${ROOT_D}/databases/misc/chr-names2.txt $< | ${BCFTOOLS} csq -p a -O z -o $@ -f ${RFASTA} -g ${RGFF3} -
add_csq:$(add_effect)

#we mark the SOMATICS Vars and annot common vars
annot_somatics=$(patsubst %.csq.vcf.bgz,%.annot_somatics.vcf.bgz,$(add_effect))

#Annot somatics
%.annot_somatics.vcf.bgz:%.csq.vcf.bgz
	${BCFTOOLS} annotate -a ../matched/$(subst .tonly.csq.vcf.bgz,.csq.vcf.bgz,$<) -m +SOMATIC -O z -o $@ $< 
som:$(annot_somatics)

#split somatics from germline
split_sg=$(patsubst %.annot_somatics.vcf.bgz,%.somatics.vcf.bgz,$(annot_somatics))
%.somatics.vcf.bgz:%.annot_somatics.vcf.bgz
	${BCFTOOLS}  view -e "SOMATIC==1" -O z -o $(subst .annot_somatics.vcf.bgz,.germline.vcf.bgz,$<) $<
	${BCFTOOLS}  view -i "SOMATIC==1" -O z -o $(subst .annot_somatics.vcf.bgz,.somatics.vcf.bgz,$<) $<
split:$(split_sg)
# we built the variants matrix
matrix=$(patsubst %.somatics.vcf.bgz,%.somatics.snv.matrix,$(split_sg))
#CMF=/Users/adigenova/Projects/IARC/bioinfo/iarcbioinfo/TON-Calling/VCF/extract_feature_from_VCF.pl
CMF=${ROOT_D}/code/perl/VCF/extract_feature_from_VCF.pl
%.somatics.snv.matrix:%.somatics.vcf.bgz
	perl ${CMF} -a $< > $@.log
	perl ${CMF} -a $(subst .somatics.vcf.bgz,.germline.vcf.bgz,$<) > $@.germ.log
to_matrix:$(matrix)

# we merge the matrix in a single file of variants for Somatics and Germline
COLUMS="FILE CHROM POS SNVS BCSQ CODING COSMIC_CENSUS_GENE COSMIC GNOMAD CENTROMER MPOS DP GERMQ SEQQ STRANDQ TLOD AF ADR ADA OR1 OR2 OA1 OA2 NS SOMATIC"
Somatics.snv.matrix.txt:$(matrix)
	cat $(matrix) |  grep -v "FILE"| awk -v colums=${COLUMS} 'BEGIN{print colums}{a[$$2"-"$$3]++; b[$$2"-"$$3]=$$0}END{for(i in a){print b[i]" "a[i]" 1"}}' | sed 's/_filtered_PASS_norm.tonly.somatics.vcf.bgz//g' > $@
	cat $(subst .snv.matrix,.indel.matrix,$(matrix)) |  grep -v "FILE"| awk -v colums=${COLUMS} 'BEGIN{print colums}{a[$$2"-"$$3]++; b[$$2"-"$$3]=$$0}END{for(i in a){print b[i]" "a[i]" 1"}}' | sed 's/_filtered_PASS_norm.tonly.somatics.vcf.bgz//g' > $(subst .snv.matrix.txt,.indel.matrix.txt,$@)
	cat $(subst .somatics.snv.matrix,.germline.snv.matrix,$(matrix)) |  grep -v "FILE"| awk -v colums=${COLUMS} 'BEGIN{print colums}{a[$$2"-"$$3]++; b[$$2"-"$$3]=$$0}END{for(i in a){print b[i]" "a[i]" 0"}}' | sed 's/_filtered_PASS_norm.tonly.germline.vcf.bgz//g' > Germline.snv.matrix.txt
	cat $(subst .somatics.snv.matrix,.germline.indel.matrix,$(matrix)) |  grep -v "FILE"| awk -v colums=${COLUMS} 'BEGIN{print colums}{a[$$2"-"$$3]++; b[$$2"-"$$3]=$$0}END{for(i in a){print b[i]" "a[i]" 0"}}' | sed 's/_filtered_PASS_norm.tonly.germline.vcf.bgz//g' > Germline.indel.matrix.txt

#remove MESO_061 from data
Somatics.snv.matrix.WMESO_061.txt:Somatics.snv.matrix.txt
	grep -v "MESO_061" Somatics.snv.matrix.txt   > Somatics.snv.matrix.WMESO_061.txt
	grep -v "MESO_061" Somatics.indel.matrix.txt > Somatics.indel.matrix.WMESO_061.txt
	grep -v "MESO_061" Germline.snv.matrix.txt > Germline.snv.matrix.WMESO_061.txt 
	grep -v "MESO_061" Germline.indel.matrix.txt > Germline.indel.matrix.WMESO_061.txt
	#files with only MESO_061
	egrep "^FILE|MESO_061" Somatics.snv.matrix.txt   > Somatics.snv.matrix.MESO_061.txt
	egrep "^FILE|MESO_061" Somatics.indel.matrix.txt > Somatics.indel.matrix.MESO_061.txt
	egrep "^FILE|MESO_061" Germline.snv.matrix.txt   > Germline.snv.matrix.MESO_061.txt 
	egrep "^FILE|MESO_061" Germline.indel.matrix.txt > Germline.indel.matrix.MESO_061.txt

all: gnomad_annot comic_annot centro_genes add_csq som split to_matrix Somatics.snv.matrix.txt Somatics.snv.matrix.WMESO_061.txt
