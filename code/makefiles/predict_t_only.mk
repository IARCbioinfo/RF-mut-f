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
	${TABIX} $@
gnomad_annot:$(gnomad_annot_files)

#annotate cosmic vars
COSMIC_DB=${ROOT_D}/databases/cosmic/CosmicSNV_INDELS.gnomad.vcf.bgz
cosmic_annot_files=$(patsubst %.gnomad.vcf.bgz,%.cosmic.vcf.bgz,$(gnomad_annot_files))
%.cosmic.vcf.bgz:%.gnomad.vcf.bgz
	${BCFTOOLS} annotate -r ${CHRS} -O z -o $@ -c ID -m +COSMIC -a ${COSMIC_DB} $<
	${TABIX} $@
comic_annot:$(cosmic_annot_files)

#annotate cosmic genes and centromers
add_centro_cgenes=$(patsubst %.cosmic.vcf.bgz,%.centro_cgenes.vcf.bgz,$(cosmic_annot_files))
%.centro_cgenes.vcf.bgz:%.cosmic.vcf.bgz
	${BCFTOOLS} annotate -a ${ROOT_D}/databases/misc/centromers.bed.gz -h ${ROOT_D}/databases/misc/centromers.hdr -c CNAME --columns CHROM,FROM,TO,CNAME -m +CENTROMER $< | \
	${BCFTOOLS} annotate -a ${ROOT_D}/databases/misc/cosmic_gene_census.bed.gz -h ${ROOT_D}/databases/misc/cosmic_gene_census.hdr  -c CGENE --columns CHROM,FROM,TO,CGENE -O z -o $@ -m +COSMIC_CENSUS_GENE -
	${TABIX} $@
centro_genes:$(add_centro_cgenes)

#annotate variants effects
RFASTA=${ROOT_D}/databases/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa
RGFF3=${ROOT_D}/databases/genome/Homo_sapiens.GRCh38.99.chr.gff3.gz

add_effect=$(patsubst %.centro_cgenes.vcf.bgz,%.csq.vcf.bgz,$(add_centro_cgenes))
%.csq.vcf.bgz: %.centro_cgenes.vcf.bgz
	${BCFTOOLS} annotate --rename-chrs ${ROOT_D}/databases/misc/chr-names2.txt $< | ${BCFTOOLS} csq -p a -f ${RFASTA} -g ${RGFF3} - | \
	${BCFTOOLS} annotate --rename-chrs ${ROOT_D}/databases/misc/number2names.txt -O z -o $@ -
	${TABIX} $@
add_csq:$(add_effect)

#build the matrix 
matrix=$(patsubst %.csq.vcf.bgz,%.csq.snv.matrix,$(add_effect))
CMF=${ROOT_D}/code/perl/VCF/extract_feature_from_VCF.pl
%.csq.snv.matrix:%.csq.vcf.bgz
	perl ${CMF} -a $< > $@.log
to_matrix:$(matrix)

# we merge the matrix in a single file of variants for Somatics and Germline
Mutations.snv.matrix.txt:$(matrix)
	cat $(matrix) |  grep -v "FILE"| awk 'BEGIN{print "FILE CHROM POS SNVS BCSQ CODING COSMIC_CENSUS_GENE COSMIC GNOMAD CENTROMER MPOS DP GERMQ SEQQ STRANDQ TLOD AF ADR ADA OR1 OR2 OA1 OA2 NS SOMATIC"}{a[$$2"-"$$3]++; b[$$2"-"$$3]=$$0}END{for(i in a){print b[i]" "a[i]" 1"}}' | sed 's/_filtered_PASS_norm.csq.vcf.bgz//g' > $@
	cat $(subst .snv.matrix,.indel.matrix,$(matrix)) |  grep -v "FILE"| awk 'BEGIN{print "FILE CHROM POS SNVS BCSQ CODING COSMIC_CENSUS_GENE COSMIC GNOMAD CENTROMER MPOS DP GERMQ SEQQ STRANDQ TLOD AF ADR ADA OR1 OR2 OA1 OA2 NS SOMATIC"}{a[$$2"-"$$3]++; b[$$2"-"$$3]=$$0}END{for(i in a){print b[i]" "a[i]" 1"}}' | sed 's/_filtered_PASS_norm.csq.vcf.bgz//g' > $(subst .snv.matrix.txt,.indel.matrix.txt,$@)

#random forest models
RFM=${ROOT_D}/mesomics/release2/matched-t-only
#we make the predictions using the best RF model for SNVs and INDELs
	#SNV
rf_snv_m61_predictions.txt:Mutations.snv.matrix.txt
	Rscript ${ROOT_D}/code/Rscripts/RF-APPLY-MODEL.R -i Mutations.snv.matrix.txt -o ${PWD} -m ${RFM}/rf-12_1500_25_snv_meso61.rds -s rf_snv_m61
	#Rscript RF-APPLY-MODEL.R -i Mutations.snv.matrix.txt -o ${PWD} -m ${RFM}/rf-12_1500_25_snv_wmeso61.rds -s rf_snv_wm61
	#INDELS
rf_indels_m61_predictions.txt:Mutations.snv.matrix.txt
	Rscript ${ROOT_D}/code/Rscripts/RF-APPLY-MODEL.R -i Mutations.indel.matrix.txt -o ${PWD} -m ${RFM}/rf-12_1500_25_indel_meso61.rds -s rf_indels_m61
	#Rscript RF-APPLY-MODEL.R -i Mutations.indel.matrix.txt -o ${PWD} -m ${RFM}/rf-12_1500_25_indel_wmeso61.rds -s rf_indels_wm61

all: centro_genes add_csq Mutations.snv.matrix.txt rf_snv_m61_predictions.txt rf_indels_m61_predictions.txt rf_indels_m61_predictions.txt



