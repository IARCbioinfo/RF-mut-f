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
	perl ${CMF} -a $< > $@.log 2>$@.err
to_matrix:$(matrix)

#random forest models
RFM=${ROOT_D}/mesomics/release2/matched-t-only

#we call the random forest on uncollapsed data
Mutations.snv.uncollapsed.matrix.txt:$(matrix)
	cat $(matrix) |  grep -v "FILE" | awk '{a[$$2"-"$$3]++;}END{for(i in a){print i" "a[i]}}' > snv.ns
	#the NF== 25 remove one var from MESO_084 that constains 24 columms due to the segfault error, the var is "chrM 14154"
	cat $(matrix) |  grep -v "FILE" | awk 'BEGIN{print "FILE CHROM POS SNVS BCSQ CODING COSMIC_CENSUS_GENE COSMIC GNOMAD CENTROMER MPOS DP GERMQ SEQQ STRANDQ TLOD AF ADR ADA OR1 OR2 OA1 OA2 NS SOMATIC"; while(getline <"snv.ns"){a[$$1]=$$2}}{print $$0" "a[$$2"-"$$3]" 1"}' | sed 's/_filtered_PASS_norm.csq.vcf.bgz//g' | awk '{if(NF == 25){print $$0}}' > $@

#we call the random forest on uncollapsed data
Mutations.indel.uncollapsed.matrix.txt:$(matrix)
	cat $(subst .snv.matrix,.indel.matrix,$(matrix)) |  grep -v "FILE" | awk '{a[$$2"-"$$3]++;}END{for(i in a){print i" "a[i]}}' > indel.ns
	cat $(subst .snv.matrix,.indel.matrix,$(matrix)) |  grep -v "FILE" | awk 'BEGIN{print "FILE CHROM POS SNVS BCSQ CODING COSMIC_CENSUS_GENE COSMIC GNOMAD CENTROMER MPOS DP GERMQ SEQQ STRANDQ TLOD AF ADR ADA OR1 OR2 OA1 OA2 NS SOMATIC"; while(getline <"indel.ns"){a[$$1]=$$2}}{print $$0" "a[$$2"-"$$3]" 1"}' | sed 's/_filtered_PASS_norm.csq.vcf.bgz//g' > $@

#apply the random forest to the uncollapsed data
rf_snv_uncollapsed_m61.r1_predictions.txt:Mutations.snv.uncollapsed.matrix.txt
	Rscript ${ROOT_D}/code/Rscripts/RF-APPLY-MODEL-CUTOFF-VOTES-PROB.R -i Mutations.snv.uncollapsed.matrix.txt -o ${PWD} -m ${RFM}/rf-8_1000_5_snv_meso61.r1.rds -s rf_snv_uncollapsed_m61.r1
#INDELS
rf_indel_uncollapsed_m61.r1_predictions.txt:Mutations.indel.uncollapsed.matrix.txt
	Rscript ${ROOT_D}/code/Rscripts/RF-APPLY-MODEL-CUTOFF-VOTES-PROB.R -i Mutations.indel.uncollapsed.matrix.txt -o ${PWD} -m ${RFM}/rf-8_1000_5_indel_meso61.r1.rds -s rf_indel_uncollapsed_m61.r1

all: centro_genes add_csq Mutations.snv.uncollapsed.matrix.txt Mutations.indel.uncollapsed.matrix.txt rf_snv_uncollapsed_m61.r1_predictions.txt rf_indel_uncollapsed_m61.r1_predictions.txt

