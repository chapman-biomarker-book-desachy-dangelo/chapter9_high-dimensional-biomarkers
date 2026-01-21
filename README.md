# Chapter 9
High-Dimensional Genetic Biomarkers and Polygenic Risk Scores: Advanced Methods for Disease and Drug Response Prediction

## Overview

In this chapter, we select the most widely used PRS methods and introduce their software. These software tools are readily accessible for download from either GitHub or CRAN. To install and use them, users need to have R, Python, and/or Linux environment installed on their systems. Typically, these software tools require one or more inputs from disease GWAS summary statistics, PGx GWAS summary statistics, and individual-level disease or PGx GWAS data (training, validation, and testing). Additionally, users may need to provide the prefix of the bim file corresponding to the target dataset, as well as a reference panel with SNP LD information (e.g., from the 1000 Genomes Project), which can be supplied in either PLINK binary format (bed, bim, fam) or HDF5 format.

## Folder Instruction

- **code**: Main Rmarkdown files for the application of PRS software
- **function**: R scripts with key functions
- **data**:

## Usage

### Step 1: PolyPred

Update effect size estimates in disease GWAS summary statistics by conducting PolyPred (Weissbrod et al., 2022) with fine-mapping. Details can be found in **PolyPred** folder (**run_polypred_sim_example.sh**).

### Step 2: mtPRS-FIMEL

Run mtPRS-FIMEL with refined effect size estimates. Details can be found in **doc** folder for mtPRSFIMEL user manual (**mtPRSFIMEL_0.1.0.pdf**), and **vignettes** folder for a demo illustrating how to use our softwares (**README.Rmd**).
