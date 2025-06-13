########## write a pbs script for all variables to run ############
#!/bin/bash
#PBS -N my_job
#PBS -l walltime=10:00:00
#PBS -l select=1:ncpus=100:ht=true
module purge

# Set up variables
job_name="rare"
walltime="23:00:00"
memory="80GB"

#### input variable #####
geno_input="/mnt/home/n11142006/TIwi_data/Genotypic_data/Final_vcf_file/6_Tiwi_updated/Genotype_data/Rare_Variant_Analysis/02_Run2/03_SAIGE"    ## gneotype date ####
spareGRM="${geno_input}/sparseGRM_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx"    ### common variant GRM matrix
sparseGRM_sample="${geno_input}/sparseGRM_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"   #### GRM matrix sub file ####
null_plink="${geno_input}/rare_var_nullModel"  ### null model for the saige run ####
phenoFile="/mnt/home/n11142006/TIwi_data/Genotypic_data/Final_vcf_file/6_Tiwi_updated/Phenotype_data/Rare_variant_phenotype/rare_variant_edited_std.txt" ###phenotype data ####
out_null="${geno_input}/02_updated_grouping_AllVar/rare_var_null_model"   #### null model output #####
plinkfile="/mnt/home/n11142006/TIwi_data/Genotypic_data/Final_vcf_file/6_Tiwi_updated/Genotype_data/Rare_Variant_Analysis/02_Run2/rare_var_miss_hwe_het_fil"  #### plink file for the input ####
groupingFile="${geno_input}/saige_annotation_file_updated.txt"  ### saig grouping variable #####
out_result="${geno_input}/02_updated_grouping_AllVar/01_Results/Assoc_results"  #### output to save ######

x="/mnt/home/n11142006/TIwi_data/Genotypic_data/Final_vcf_file/6_Tiwi_updated/Genotype_data/Rare_Variant_Analysis/02_Run2/03_SAIGE/00_Script/" #### set working 



cd /mnt/home/n11142006/TIwi_data/Genotypic_data/Final_vcf_file/6_Tiwi_updated/Gevi notype_data/Rare_Variant_Analysis/02_Run2/03_SAIGE/00_Script/
word="weight height waist SBP_2 DBP_2  urine_creatinine_1 urine_albumin_mg_dL creatinine_1 albumin_1 hba1c_1 uric_acid_1 urine_osmolality_1 eGFR ACR BMI"
for line in $word; 
do echo $line;
   # Set up output file names
   output_files=${x}/${line}_rare.out
   
   ## Write PBS script fo file ##
   cat << EOF > ${line}_rare.pbs
   
   #!/bin/bash
   #PBS -N ${job_name}_${line}
   #PBS -l walltime=${walltime} 
   #PBS -l select=1:ncpus=4:mem=${memory}:ht=True
   #PBS -o ${output_file}
   #PBS -M vignesh.arunachalam@hdr.qut.edu.au
   #PBS -m e
   #PBS -V
 module purge
 source ~/.bashrc
 conda activate saige
 
 #### fit the null model ####
Rscript ~/vig_Tools/SAIGE/extdata/step1_fitNULLGLMM.R \
    --sparseGRMFile=${spareGRM}   \ 
	--sparseGRMSampleIDFile=${spareGRM_sample} \
    --plinkFile=${null_plink} \
    --useSparseGRMtoFitNULL=TRUE    \
    --phenoFile=${phenoFile} \
    --phenoCol=${line} \
    --covarColList=sex,Age,PC_1,PC_2,PC_3 \
	--sexCol=sex \
	--FemaleCode=2 \
	--MaleCode=1\
    --qCovarColList=Age,PC_1,PC_2,PC_3 \
    --sampleIDColinphenoFile=IID \
    --traitType=quantitative  \
    --isCateVarianceRatio=TRUE	\
    --outputPrefix=${out_null}_${line} \
    --IsOverwriteVarianceRatioFile=TRUE
	
Rscript ~/vig_Tools/SAIGE/extdata/step2_SPAtests.R \
    --sparseGRMFile=${spareGRM}   \ 
	--sparseGRMSampleIDFile=${spareGRM_sample} \
	--bedFile=${plinkfile}.bed \
	--bimFile=${plinkfile}.bim \
	--famFile=${plinkfile}.fam \
	--SAIGEOutputFile=${out_result}_${line}.txt \
	--groupFile=${groupingFile} \
	--varianceRatioFile=${out_null}_${line}.varianceRatio.txt \
	--GMMATmodelFile=${out_null}_${line}.rda \
	--is_output_markerList_in_groupTest=TRUE \
	--LOCO=FALSE \
	--is_fastTest=TRUE \
	--annotation_in_groupTest=intronic,exonic,nearby,UTR,intergenic,ncRNA,splicing

EOF
  #Submit PBS job
  qsub ${line}_rare.pbs
done
