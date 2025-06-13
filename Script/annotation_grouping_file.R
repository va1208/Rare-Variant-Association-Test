library(data.table)
library(dplyr)
setwd("/.../Annovar_annotation/")
 
data <- fread("Annotation_prediction_reqCols.tsv")  ### VEP annotation file #####
head(data)
dim(data)

bim <- fread("/.../Rare_Variant_Analysis/02_Run2/rare_var_miss_hwe_het_fil.bim") ## load the final rare variant SNPs excluding intergenic region
length(intersect(bim$V2, data$SNP))

data <- subset(data, SNP %in% bim$V2)
length(unique(data$SNP))

sub <- data[,c("SNP", "Func.refGene", "Gene.refGene")]  ### Three required columns to create a grouping file ###
sub <- sub %>% distinct() 
### rewrite the column names based on your functional annotation####
sub <- sub %>% mutate(func = recode(Func.refGene, "downstream" = "downstream",
									"upstream" = "upstream", "upstream;downstream" = "upstream/downstream",
									"exonic;splicing" = "exonic", "intergenic;intergenic" = "intergenic",
									"intronic;intronic" = "intronic", "ncRNA_exonic" = "ncRNA",
									"ncRNA_intronic" = "ncRNA", "ncRNA_exonic;splicing" = "ncRNA",
									"ncRNA_splicing" = "ncRNA", "UTR5;UTR3" = "UTR",
									"UTR3" = "UTR", "UTR5" = "UTR"))
### set working directory ###
setwd("/.../Rare_Variant_Analysis/02_Run2/03_SAIGE/")

### function to do that
convert_to_saige_annotation <- function(input_df) {

  # Create var and anno rows for each gene
  result <- input_df %>% 
    group_by(Gene.refGene) %>%
    summarise(
      var_row = paste(Gene.refGene, "var", paste(SNP, collapse = "\t"), sep = "\t"),
      anno_row = paste(Gene.refGene, "anno", paste(func, collapse = "\t"), sep = "\t")
    )
  
  # Combine var and anno rows
  saige_annotation <- c(rbind(result$var_row, result$anno_row))
  gc()
  saige_annotation <- unique(saige_annotation)
  gc()
  return(saige_annotation)
  #Write to file
  #riteLines(saige_annotation, con = "saige_annotation_file.txt")
}

out <- convert_to_saige_annotation(sub)

writeLines(out, "saige_annotation_file_updated.txt") ## final required grouping variable for the saige association analysis ####
