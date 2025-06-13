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
