bim <- fread("/mnt/home/n11142006/TIwi_data/Genotypic_data/Final_vcf_file/6_Tiwi_updated/Genotype_data/Rare_Variant_Analysis/02_Run2/rare_var_miss_hwe_het_fil.bim")
length(intersect(bim$V2, data$SNP))

data <- subset(data, SNP %in% bim$V2)
length(unique(data$SNP))

sub <- data[,c("SNP", "Func.refGene", "Gene.refGene")]
sub <- sub %>% distinct()
sub <- sub %>% mutate(func = recode(Func.refGene, "downstream" = "nearby",
									"upstream" = "nearby", "upstream;downstream" = "nearby",
									"exonic;splicing" = "exonic", "intergenic;intergenic" = "intergenic",
									"intronic;intronic" = "intronic", "ncRNA_exonic" = "ncRNA",
									"ncRNA_intronic" = "ncRNA", "ncRNA_exonic;splicing" = "ncRNA",
									"ncRNA_splicing" = "ncRNA", "UTR5;UTR3" = "UTR",
									"UTR3" = "UTR", "UTR5" = "UTR"))
### set working directory ###
setwd("/mnt/home/n11142006/TIwi_data/Genotypic_data/Final_vcf_file/6_Tiwi_updated/Genotype_data/Rare_Variant_Analysis/02_Run2/03_SAIGE/")

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

writeLines(out, "saige_annotation_file_updated.txt")
