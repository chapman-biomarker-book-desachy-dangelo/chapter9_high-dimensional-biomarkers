# Chapter 9
High-Dimensional Genetic Biomarkers and Polygenic Risk Scores: Advanced Methods for Disease and Drug Response Prediction

## Overview
mtPRSFIMEL package (Zhai et al., 2025) implements a novel multi-trait polygenic risk score (mtPRS) method, mtPRS-FIMEL. This method tackles key challenges in existing multi-trait PRS approaches by integrating information from three levels:
- variant-level: causal variants identified through fine-mapping analysis;
- trait-level: data from multiple genetically correlated traits;
- method-level: an ensemble learning framework that optimally combines complementary PRS methods, including mtPRS-PCA (Zhai et al., 2023), mtPRS-ML, and stPRSs.

## Installation

mtPRSFIMEL R package requires R version >= 4.0.3.

```
library(devtools)
devtools::install_github("yaowuliu/ACAT") # Need to install ACAT first through the Github
devtools::install_github("zhaiso1/mtPRSFIMEL")
```

## Usage

### Step 1: PolyPred

Update effect size estimates in disease GWAS summary statistics by conducting PolyPred (Weissbrod et al., 2022) with fine-mapping. Details can be found in **PolyPred** folder (**run_polypred_sim_example.sh**).

### Step 2: mtPRS-FIMEL

Run mtPRS-FIMEL with refined effect size estimates. Details can be found in **doc** folder for mtPRSFIMEL user manual (**mtPRSFIMEL_0.1.0.pdf**), and **vignettes** folder for a demo illustrating how to use our softwares (**README.Rmd**).
