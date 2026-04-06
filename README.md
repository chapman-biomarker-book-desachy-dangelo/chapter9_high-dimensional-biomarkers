# Chapter 9: High-Dimensional Genetic Biomarkers and Polygenic Risk Scores: Advanced Methods for Disease and Drug Response Prediction

## Overview

In this chapter, we select the most widely used PRS methods and introduce their software. These software tools are readily accessible for download from either GitHub or CRAN. To install and use them, users need to have R, Python, and/or Linux environment installed on their systems. Typically, these software tools require one or more inputs from disease GWAS summary statistics, PGx GWAS summary statistics, and individual-level disease or PGx GWAS data (training, validation, and testing). Additionally, users may need to provide the prefix of the bim file corresponding to the target dataset, as well as a reference panel with SNP LD information (e.g., from the 1000 Genomes Project), which can be supplied in either PLINK binary format (bed, bim, fam) or HDF5 format.

## Folder Instruction

- **code**: Main Rmarkdown files for the application of PRS software
- **function**: R scripts with key functions
- **data**: Include `base` for data bundle simulations, `input` as simulated inputs for specific methods (i.e., PRScs and PRS-PGx-TL), and `ref_panel` as reference panels for specific methods (i.e., Lassosum and PRScs)
- **result**: Pre-specified directories to save results from specific methods (i.e., LDpred2, PRScs, and PRS-PGx-TL)
- **software**: Software for PRScs and PRS-PGx-TL

## PRS Methods Instruction

- **Lassosum**: R package based method; Need to specify the reference panel
- **PRS-CS**: Python script based method; Run Python script in Terminal for PRS-CS analysis; Need to specify the reference panel, the bim file, the path to the software, and the path to save results
- **LDpred2**: R package based method; Need to specify the path to save intermediate results
- **PRS-PGx-L/GL/SGL**: R package based method; Require individual-level PGx data
- **PRS-PGx-Bayes**: R package based method; Need individual-level reference genotype data
- **PRS-PGx-TL**: R scripts based pipeline; Run multiple R scripts in Terminal for PRS-PGx-TL analysis - the next R script will take results from the previous R script as inputs
