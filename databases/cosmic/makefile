.DELETE_ON_ERROR:


#we create a bed file of the cosmic genes
cancer_gene_census_curated.bed:
	awk '{print $$2"\t"$$1}' cancer_gene_census_curated.txt | sed 's/:/\t/; s/-/\t/' | grep -v "^Sym" | sort -k1,1n -k2,2n  > $@

#we remove duplicate SNV and INDELs and check the alleles with the reference
BCFT=/home/digenovaa/binaries/bcftools-1.10.2/BCFT/bin/bcftools
REF=/data/gcs/mesomics/work/alignment_references/hs38DH/hs38DH.fa

#we remove duplicated variants and check the ref bases and coordinates
CosmicCodingMuts.uniq.vcf.bgz:
	 zcat CosmicCodingMuts.vcf.gz | \
	 awk 'BEGIN{OFS="\t"}{if($$0~/#/){print $$0}else{if($$1 == "MT"){$$1="chrM";}else{$$1="chr"$$1;} print $$0}}' |  \
	  ${BCFT} norm -f ${REF} -O z -o $@ -c s --rm-dup all -  
	tabix $@

#we remove duplicated variants and check the ref bases
CosmicNonCodingVariants.uniq.vcf.bgz:
	 zcat CosmicNonCodingVariants.vcf.gz | \
	 awk 'BEGIN{OFS="\t"}{if($$0~/#/){print $$0}else{if($$1 == "MT"){$$1="chrM";}else{$$1="chr"$$1;} print $$0}}' |  \
	  ${BCFT} norm -f ${REF} -O z -o $@ -c s --rm-dup all -  
	 tabix $@

#we merge the variants in a single file
CosmicSNV_INDELS.vcf.bgz:CosmicNonCodingVariants.uniq.vcf.bgz CosmicCodingMuts.uniq.vcf.bgz
	${BCFT} merge -O z -o $@ CosmicNonCodingVariants.uniq.vcf.bgz CosmicCodingMuts.uniq.vcf.bgz 
	tabix $@

#we compute SNV_INDELS stats
CosmicSNV_INDELS.vcf.stats:CosmicSNV_INDELS.vcf.bgz
	${BCFT} stats CosmicSNV_INDELS.vcf.bgz > $@
	${BCFT} stats CosmicNonCodingVariants.uniq.vcf.bgz > CosmicNonCodingVariants.uniq.vcf.stats
	${BCFT} stats CosmicCodingMuts.uniq.vcf.bgz > CosmicCodingMuts.uniq.vcf.stats


#We annotate the cosmic variants present in genomad
GENOMAD_DB=/data/scratch/digenovaa/Somatic-reference-free/gnomad/gnomad.genome.vcf.bgz
CosmicSNV_INDELS.gnomad.vcf.bgz:CosmicSNV_INDELS.vcf.bgz
	bcftools annotate -O z -o $@  -c "GNOMAD_AC:=AC,GNOMAD_AN:=AN" -m -NO_GNOMAD -a ${GENOMAD_DB} $<
	tabix $@

#query examples 
#bcftools  query -f'%CHROM %POS %ID %GNOMAD_AC %CNT\n' -i 'GNOMAD_AC > 10' CosmicSNV_INDELS.gnomad.vcf.bgz
#

### zcat CosmicBreakpointsExport.tsv.gz | awk -F"\t" '{print $13" "$16" "$17" "$18" "$19" "$20" "$21" "$22" "$12" "$15}' | awk '{if($0~/Interchromosomal/){print $2"\t"$3"\t"$4"\tCOST"$1"-F"; print $6"\t"$7"\t"$8"\tCOST"$1"-T"}else{print $2"\t"$3"\t"$8"\tCOST"$1"-I"; print $0}}' | less

all: cancer_gene_census_curated.bed CosmicCodingMuts.uniq.vcf.bgz CosmicNonCodingVariants.uniq.vcf.bgz CosmicSNV_INDELS.vcf.bgz CosmicSNV_INDELS.vcf.stats CosmicSNV_INDELS.gnomad.vcf.bgz

