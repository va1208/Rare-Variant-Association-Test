# Rare-Variant-Association-Test
This repository includes the analysis plan for the rare variant association testing (RVAT) in one of the Australian Indigenous populations.

The overall pipeline for the rare variant association analysis is given below

![pipeline_figure](https://github.com/user-attachments/assets/c8ea6968-2385-40b8-a5be-bca5b14305d4)


**Create a conda environment for saige**

    conda env create -f saige.yml


**Run an annotation using the PathVar pipeline**

Please refer to this pipeline to obtain information about the annotation pipeline - https://github.com/Mohammed-Alfayyadh/Pathogenic-variant-calling-pipeline-PathVar/tree/main

**Create multiple grouping files**

    Rscript annotation_grouping_file.R

**Grouping Strategies**

In the analysis, We have used multiple grouping strategies and includes,
        1. likely pathogenic, pathogenic variant genes
        2. Exonic region - grouped by genes
        3. Whole genome region - grouped all variant in the gene


**Run SAIGE-GENE pipeline**

Use the script to run Saige pipeline for different grouping strategies

    Rscript saige_run.sh ### Make sure all the corresponding input files are there ###
    
**Overall gene burden association**

This analysis is to check the overall burden of pathogenic or likely pathogenic variants are in association with kidney biomarkers after adjusting for covariates. Please use the script to run this analysis **overall_gene_burden.R**.  You will obtain the regression model coefficient table for this analysis.
