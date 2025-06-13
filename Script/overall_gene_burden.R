#### set working directory ###
setwd("/../Rare_Variant_Analysis/03_Run_maf005")

#### get plink file to VCF file #### only run for first time 
system("plink1.9 --bfile likely_pathogenic_gene_var.bim --recode vcf-iid --out likely_pathogenic_gene_var_vcf")

### install the vcfR package to load VCF file ###
if (!requireNamespace("vcfR", quietly = TRUE)) {
  install.packages("vcfR")
}

### load packages ###
library(vcfR)
library(moments)
library(data.table)
library(broom)
library(gt)


#### import vcf file and extract the genotype ####
vcf <- read.vcfR("likely_pathogenic_gene_var_vcf.vcf", verbose = FALSE)
genotype_data <- extract.gt(vcf)  # Extracts only genotype calls
genotype_df <- as.data.frame(genotype_data)

### convert the genotype into 0s and 1s ###
convert_genotype <- function(gt) { 
    ifelse(gt %in% c("0|0", "0/0"), 0, 
        ifelse(gt %in% c("0/1", "0|1", "1/0", "1|0"), 1,
           ifelse(gt %in% c("1/1", "1|1"), 1, NA)))
}

binary_genotype <- apply(genotype_data, c(1, 2), convert_genotype)
binary_df <- as.data.frame(binary_genotype)

### transpose the data ###
binary_df_ind <- as.data.frame(t(as.matrix(binary_df)))

### calculate the cumulative sum for the individuals ###
binary_df_ind$Score <- rowSums(binary_df_ind, na.rm = TRUE)

### z score standardization ###
binary_df_ind$std_score <- scale(binary_df_ind$Score)

### min max std score ###
binary_df_ind$minmax_score <- (binary_df_ind$Score - min(binary_df_ind$Score))/(max(binary_df_ind$Score) - min(binary_df_ind$Score))

### save the data using fwrite ####
binary_df_ind$FID <- rownames(binary_df_ind)
fwrite(binary_df_ind, "/mnt/hpccs01/home/n11142006/TIwi_data/Genotypic_data/Final_vcf_file/6_Tiwi_updated/Genotype_data/Rare_Variant_Analysis/03_Run_maf005/04_TotalBurden/Pathogenic_variant_total_burden.tsv.gz", 
                sep = "\t", row.names = T, quote = F)

#### import the phenotye data ### 
pheno <- fread("/work/NagarajLab_BG/7_Tiwi_GWAS/Phenotype_data/Tiwi_phenotype_n492_full_updated_2024.csv")

### get required columns ###
required_col <- c("FID", "Age", "sex", "weight", "height", "waist", "urine_creatinine_1", "urine_albumin_mg_dL", "creatinine_1", "albumin_1", "hba1c_1", "uric_acid_1", "eGFR", "ACR", "diabetes_status", "dialysis_status", "ACR_status", "eGFR_status", "ACR_CKD", "eGFR_CKD")
pheno_sub <- pheno[,..required_col]
str(pheno_sub)
dim(pheno_sub)

dim(binary_df_ind)

### combine the genotype and phemnotype data ####
gt_pheno <- merge(binary_df_ind, pheno_sub, by = "FID")
str(gt_pheno)
dim(gt_pheno)
gt_pheno$sex <- as.factor(gsub("M", 0, gsub("F", 1, gt_pheno$sex)))

fwrite(gt_pheno, "/mnt/hpccs01/home/n11142006/TIwi_data/Genotypic_data/Final_vcf_file/6_Tiwi_updated/Genotype_data/Rare_Variant_Analysis/03_Run_maf005/04_TotalBurden/Pathogenic_variant_total_burden_gt_pheno.tsv.gz", 
                sep = "\t", row.names = F, quote = F)

### import the file for future ###
gt_pheno <- fread("/mnt/hpccs01/home/n11142006/TIwi_data/Genotypic_data/Final_vcf_file/6_Tiwi_updated/Genotype_data/Rare_Variant_Analysis/03_Run_maf005/04_TotalBurden/Pathogenic_variant_total_burden_gt_pheno.tsv.gz")

### determine the outliters
q1 <- quantile(gt_pheno$Score, probs = 0.25, na.rm = FALSE)
q3 <- quantile(gt_pheno$Score, probs = 0.75, na.rm = FALSE)

iqr_value <- q3-q1

lowerb <- q1 - 1.5 * iqr_value
lowerb
upperb <- q3 + 1.5 * iqr_value
upperb


gt_pheno_cleaned <- gt_pheno[gt_pheno$Score >= lowerb & gt_pheno$Score <= upperb,]
gt_pheno_cleaned$BMI <- gt_pheno_cleaned$weight/(gt_pheno_cleaned$height/100)^2

## calculate the std score again ###
gt_pheno_cleaned$std_score <- scale(gt_pheno_cleaned$Score)

### min max std score ###
gt_pheno_cleaned$minmax_score <- (gt_pheno_cleaned$Score - min(gt_pheno_cleaned$Score))/(max(gt_pheno_cleaned$Score) - min(gt_pheno_cleaned$Score))

gt_pheno_cleaned$sex <- as.factor(gt_pheno_cleaned$sex)
gt_pheno_cleaned$diabetes_status <- as.factor(gt_pheno_cleaned$diabetes_status)

# #### remove the ACR outliers ####
# q1 <- quantile(gt_pheno$ACR, probs = 0.25, na.rm = FALSE)
# q3 <- quantile(gt_pheno$ACR, probs = 0.75, na.rm = FALSE)

# iqr_value <- q3-q1

# lowerb <- q1 - 1.5 * iqr_value
# lowerb
# upperb <- q3 + 1.5 * iqr_value
# upperb

# gt_pheno_cleaned <- gt_pheno_cleaned[gt_pheno_cleaned$ACR >= lowerb & gt_pheno_cleaned$ACR <= upperb,]
# dim(gt_pheno_cleaned)


#gt_pheno$sex <- as.factor(gt_pheno$sex)

### fit a simple regresion with score as outcome ###
gt_pheno_cleaned$logACR <- log10(gt_pheno_cleaned$ACR)

fit <- lm(std_score ~ ACR + Age + sex, gt_pheno_cleaned)
fit <- lm(std_score ~ eGFR , gt_pheno_cleaned)

fit <- lm(std_score ~ eGFR + logACR + hba1c_1 + sex + Age + uric_acid_1 + waist + height + weight, 
                    gt_pheno_cleaned)
stepwise_model <- step(fit, direction = "both")
summary(stepwise_model)

subfit <- lm(logACR ~ std_score + hba1c_1 + sex + Age + uric_acid_1 + waist + height + weight + diabetes_status, gt_pheno_cleaned)
subfit_step <- step(subfit, direction = "backward")
summary(subfit_step)


#final <- lm(logACR ~ std_score + hba1c_1 + sex + Age + uric_acid_1 + weight + diabetes_status, gt_pheno_cleaned)
#summary_final <- summary(final)
#saveRDS(final, "04_TotalBurden/pathogenic_variant_linear_regression.RDS")

final <- lm(logACR ~ std_score + hba1c_1 + sex + Age + uric_acid_1 + weight, gt_pheno_cleaned)
summary_final <- summary(final)
saveRDS(final, "04_TotalBurden/pathogenic_variant_linear_regression.RDS")

#### create a table ####
coef_table <- summary_final$coefficients
coef_table
conf_int <- confint(final)
regres_table <- as.data.frame(cbind(coef_table, conf_int))
regres_table
colnames(regres_table) <- c("Effect_size", "Std_Error", "statistic", "p_value", "Lower_CI", "Upper_CI")
regres_table$variable <- rownames(regres_table)
regres_table <- regres_table[,c("variable", "Effect_size", "Std_Error", "statistic", "p_value", "Lower_CI", "Upper_CI")]
regres_table$variable <- gsub("[(]", "", gsub("[)]", "", regres_table$variable))
write.table(regres_table, "04_TotalBurden/pathogenic_variant_linear_regression.csv", quote = F, row.names = F, sep = ",")
