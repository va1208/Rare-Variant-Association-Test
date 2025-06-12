#### calculate allele counts for each marker in whole plink file ####
module load plink/2.00a2.3_x86_64
cd /home/n11142006/TIwi_data/Genotypic_data/Final_vcf_file/6_Tiwi_updated/Genotype_data/Rare_Variant_Analysis/02_Run2

#### create a GRM for Tiwi samples ####
#### use the common variant to create grm for further analysis ####
plink1.9 --bfile ~/TIwi_data/Genotypic_data/Final_vcf_file/6_Tiwi_updated/Genotype_data/Merged_analysis/1_MAF0.05/Tiwi_data_filtered --keep rare_var_miss_hwe_het_fil.fam --keep-allele-order --make-bed --out common_variant_grm
Rscript ~/vig_Tools/SAIGE/extdata/createSparseGRM.R \
     --plinkFile=common_variant_grm \
     --nThreads=4  \
     --outputPrefix=03_SAIGE/sparseGRM \
     --numRandomMarkerforSparseKin=2000 \
     --relatednessCutoff=0.125
#### produce the grm for common variant ####
