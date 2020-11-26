#create mesomics rel2 dataset
mkdir -p mesomics/release2
mkdir -p mesomics/release2/matched
mkdir -p mesomics/release2/matched-t-only
mkdir -p mesomics/release2/t-only
mkdir -p mesomics/release2/t-only-vcfs
#matched data (46)
RELEASE2=/data/gcs/mesomics/files/WGS/variant_calling/somatic_release2_26032020/intermediate_files/normalized_calling
ln -s ${RELEASE2}/Mutect2_somatic_pon_blood_TN_gatk4150_normalized_20200324/*_norm.vcf.gz $PWD/mesomics/release2/matched	
ln -s ${RELEASE2}/Mutect2_somatic_pon_blood_TN_gatk4150_normalized_20200324/*_norm.vcf.gz.tbi $PWD/mesomics/release2/matched
#matched-t-only data (46)
ln -s ${RELEASE2}/Mutect2_somatic_pon_blood_TN_Tonly_mode_gatk4150_normalized_20200324/*_norm.vcf.gz $PWD/mesomics/release2/matched-t-only
ln -s ${RELEASE2}/Mutect2_somatic_pon_blood_TN_Tonly_mode_gatk4150_normalized_20200324/*_norm.vcf.gz.tbi $PWD/mesomics/release2/matched-t-only
#t-only data (73)
ln -s ${RELEASE2}/Mutect2_somatic_pon_blood_Tonly_gatk4150_normalized_20200325/*_norm.vcf.gz $PWD/mesomics/release2/t-only
ln -s ${RELEASE2}/Mutect2_somatic_pon_blood_Tonly_gatk4150_normalized_20200325/*_norm.vcf.gz.tbi $PWD/mesomics/release2/t-only
#we remove the MESO_094 calls
rm -f  $PWD/mesomics/release2/t-only/MESO_094_T_filtered_PASS_norm.vcf.gz
rm -f  $PWD/mesomics/release2/t-only/MESO_094_T_filtered_PASS_norm.vcf.gz.tbi
#better calling for MESO_094
MESO_094=/data/gcs/mesomics/files/WGS/variant_calling/somatic_release2_MESO_094_28102020/intermediate_files/Mutect2-nf_results_MESOMICS_26102020_pon_blood_Tonly_gatk4150_MESO_094_normalized
ln -s ${MESO_094}/MESO_094_filtered_PASS_norm.vcf.gz $PWD/mesomics/release2/t-only/MESO_094_T_filtered_PASS_norm.vcf.gz
ln -s ${MESO_094}/MESO_094_filtered_PASS_norm.vcf.gz.tbi $PWD/mesomics/release2/t-only/MESO_094_T_filtered_PASS_norm.vcf.gz.tbi
#annovar annotations
REL2ANNOVAR=/data/gcs/mesomics/files/WGS/variant_calling/somatic_release2_26032020/intermediate_files/Tonly_before_RF
ln -s ${REL2ANNOVAR}/*.vcf.gz mesomics/release2/t-only-vcfs/
ln -s ${REL2ANNOVAR}/*.vcf.gz.tbi mesomics/release2/t-only-vcfs/
ln -s ${REL2ANNOVAR}/*.txt mesomics/release2/t-only-vcfs/
#we replace MESO_094
MESO_094_ANNOT=/data/gcs/mesomics/files/WGS/variant_calling/somatic_release2_MESO_094_28102020/intermediate_files/Mutect2-nf_results_MESOMICS_26102020_pon_blood_Tonly_gatk4150_MESO_094_normalized_annotated
#we delete the symbolic links
rm -f mesomics/release2/t-only-vcfs/MESO_094_T_filtered_PASS_norm.vcf.hg38_multianno.*
#we replace the new meso calling
ln -s ${MESO_094_ANNOT}/MESO_094_filtered_PASS_norm.vcf.hg38_multianno.txt mesomics/release2/t-only-vcfs/MESO_094_T_filtered_PASS_norm.vcf.hg38_multianno.txt
ln -s ${MESO_094_ANNOT}/MESO_094_filtered_PASS_norm.vcf.hg38_multianno.vcf.gz mesomics/release2/t-only-vcfs/MESO_094_T_filtered_PASS_norm.vcf.hg38_multianno.vcf.gz
ln -s ${MESO_094_ANNOT}/MESO_094_filtered_PASS_norm.vcf.hg38_multianno.vcf.gz.tbi mesomics/release2/t-only-vcfs/MESO_094_T_filtered_PASS_norm.vcf.hg38_multianno.vcf.gz.tbi
