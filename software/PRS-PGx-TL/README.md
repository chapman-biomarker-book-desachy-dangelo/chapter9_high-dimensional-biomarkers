# PRS-PGx-TL
A transfer learning (TL) based method to leverage large-scale disease GWAS summary statistics and individual-level pharmacogenomics (PGx) data to predict drug response

## Tutorial
### To use PRS-PGx-TL-M1, M2, M3, or M4
#### Step 1: inner layer CV
```
Rscript s1.R [file for PRS weights from disease GWAS] [file for phenotype and covariates] [file for genotype] [file for SNP information] [initial values] [number of total snps with non-zero effects] [path for outputs in s1] \
covar="age,gender"
```
- **file for PRS weights from disease GWAS:** The following columns are required: `chr`, `SNP`, `BP`, `Eff`, `Ref`, `Beta`, which represent chromosome, SNP, base pair position, effect allele, reference allele, and beta coefficient, respectively. Column names are needed.
```
  chr        SNP     BP Eff Ref          Beta        
1  19  rs2312724 266034   C   T  0.0006846221 
2  19 rs11084928 277776   A   G -0.0005535947 
3  19 rs11883060 335055   C   A -0.0003439436
4  19 rs12982646 499978   A   G  0.0014327104
5  19  rs4919908 508626   T   C  0.0002608665
6  19  rs2288951 538487   A   G -0.0007215167
...
```
- **file for phenotype and covariates:** The following columns are required: `ID`, `Y`, covariates (including `Tr`), `prs_g`, which represent the ID for each individual, the drug response phenotype to be predicted, other covariates (must have the treatment variable `Tr`, the PRS calculated from the weights from disease GWAS. Column names are needed.
```
   ID          Y Tr       prs_g
1 id1  0.6430834  0  0.33243675
2 id2 -0.2784456  0  0.67831252
3 id3  3.4273028  1 -0.01487689
4 id4 -0.7722807  0  1.46311777
5 id5 -1.3991453  0  0.62429906
6 id6 -0.1259501  0  0.32419528
...
```
- **file for genotype:** The genotype matrix of individual x SNP. Row names are ID and column names are SNP.
```
     rs2312724 rs11084928 rs11883060 rs12982646 rs4919908 rs2288951 ... 
id1          0          2          0          1         1         1          
id2          0          0          0          1         1         0          
id3          1          0          1          1         1         2          
id4          0          0          0          0         0         0          
id5          0          1          0          0         1         0          
id6          0          0          0          0         0         1          
id7          0          0          1          0         0         1          
id8          0          1          0          1         0         2          
id9          1          1          1          1         2         1          
id10         0          1          0          1         2         1          
...
```
- **file for SNP information:** Allele info for the SNPs in genotype file. The following columns are required: `SNP`, `Ref`. Column names are needed.
```
   Chromosome Position        SNP Eff Ref                       ID
1          19   266034  rs2312724   T   C  19_266034_rs2312724_T_C
2          19   277776 rs11084928   G   A 19_277776_rs11084928_G_A
3          19   335055 rs11883060   A   C 19_335055_rs11883060_A_C
4          19   499978 rs12982646   G   A 19_499978_rs12982646_G_A
5          19   508626  rs4919908   C   T  19_508626_rs4919908_C_T
6          19   538487  rs2288951   G   A  19_538487_rs2288951_G_A
...
```
- **initial values:** `PRS` or `zero`, indicating the starting values of the predictive effects.
- **number of total snps with non-zero effects:** The total number of SNPs in PRS weights file that have non-zero effects.
- **path for outputs in s1:** For example, `/Mypath/s1_result`. If step 1 is performed on each chromosome separately, please save the outputs as `/Mypath/s1_result_1`,`/Mypath/s1_result_2`...`/Mypath/s1_result_22`.
- **Optional arguments:** `covar` specifies covariates to be adjusted in the model (other than `Tr`). For example, specify `covar="age,gender"` to include `age` and `gender` in the model (`Tr` will be included automatically). If not specified, only `Tr` will be included.


#### Step 2: parameter tuning
```
Rscript s2.R [file for phenotype and covariates] [initial values] [path for outputs in s1] [path for outputs in s2] [indicator of whether files for all 22 chromosomes are separated] \
covar="age,gender"
```
- **file for phenotype and covariates:** The following columns are required: `ID`, `Y`, covariates (including `Tr`), which represent the ID for each individual, the drug response phenotype to be predicted, other covariates (must have the treatment variable `Tr`. Column names are needed.
- **initial values:** `PRS` or `zero`, indicating the starting values of the predictive effects.
- **path for outputs in s1:** The output from step 1. For example, `/Mypath/s1_result`. If step 1 is performed on each chromosome separately and the outputs are saved as `/Mypath/s1_result_1`,`/Mypath/s1_result_2`...`/Mypath/s1_result_22`, just need to input `/Mypath/s1_result` here and the algorithm will add `_chr` automatically.
- **path for outputs in s2:** The output from step 2. For example, `/Mypath/s2_result`.
- **indicator of whether files for all 22 chromosomes are separated:** `yes` or `no`, indicating whether step 1 is performed on all chromosomes or on each chromosome separately.
- **Optional arguments:** `covar` specifies covariates to be adjusted in the model (other than `Tr`). 


#### Step 3: rerun the algorithm and get the output with the best parameters
```
Rscript s3.R [file for PRS weights from disease GWAS] [file for phenotype and covariates] [file for genotype] [file for SNP information] [initial values] [number of total snps with non-zero effects] [criterion] [path for outputs in s2] [path for outputs in s3] \
covar="age,gender"
```
- The inputs `file for PRS weights from disease GWAS`, `file for phenotype and covariates`, `file for genotype`, `file for SNP information`, `initial values`, `number of total snps with non-zero effects`, and the optional argument `covar` are the same as those in `s1.R`. 
- **criterion:** `best_R2` or `best_condR2`, indicating whether to maximize overall R2 or the conditional R2 of GxT (conditional on G) when selecting the best parameters. 
- **path for outputs in s2:** The output from step 2. For example, `/Mypath/s2_result`.
- **path for outputs in s3:** The output from step 3. For example, `/Mypath/s3_result`. If step 3 is performed on each chromosome separately, please save the outputs as `/Mypath/s3_result_1`,`/Mypath/s3_result_2`...`/Mypath/s3_result_22`.

#### Step 4: sum the PRS over all chromosomes for the testing samples (optional)
If step 1 and 3 are performed on each chromosome separately, choose to run the step 4 to sum the PRS over all 22 chromosomes for the testing samples:
```
Rscript s4.R [initial values] [criterion] [path for outputs in s3] [path for outputs in s4]
```
- **initial values:** `PRS` or `zero`, indicating the starting values of the predictive effects.
- **criterion:** `best_R2` or `best_condR2`, indicating whether to maximize overall R2 or the conditional R2 of GxT (conditional on G) when selecting the best parameters. 
- **path for outputs in s3:** The output from step 3. If the outputs are saved as `/Mypath/s3_result_1`,`/Mypath/s3_result_2`...`/Mypath/s3_result_22`, just need to input `/Mypath/s3_result` here and the algorithm will add `_chr` automatically.
- **path for outputs in s4:** The output from step 4. For example, `/Mypath/s4_result`.


The choice of M1-M4 corresponds to the following arguments:
| Model | initial values | criterion
| :--- | :--- | :--- |
| M1 | `zero` | `best_R2`
| M2 | `PRS` | `best_R2`
| M3 | `zero` | `best_condR2`
| M4 | `PRS` | `best_condR2`

### To use PRS-PGx-TL-M5 or M6
#### Step 1: inner layer CV
```
Rscript s1_fixG.R [file for PRS weights from disease GWAS] [file for phenotype and covariates] [file for genotype] [file for SNP information] [initial values] [number of total snps with non-zero effects] [path for outputs in s1] \
covar="age,gender"
```
- **file for PRS weights from disease GWAS:** The following columns are required: `chr`, `SNP`, `BP`, `Eff`, `Ref`, `Beta`, which represent chromosome, SNP, base pair position, effect allele, reference allele, and beta coefficient, respectively. Column names are needed.
- **file for phenotype and covariates:** The following columns are required: `ID`, `Y`, covariates (including `Tr`), `prs_g`, which represent the ID for each individual, the drug response phenotype to be predicted, other covariates (must have the treatment variable `Tr`, the PRS calculated from the weights from disease GWAS. Column names are needed.
- **file for genotype:** The genotype matrix of individual x SNP. Row names are ID and column names are SNP.
- **file for SNP information:** Allele info for the SNPs in genotype file. The following columns are required: `SNP`, `Ref`. Column names are needed.
- **initial values:** `PRS` or `zero`, indicating the starting values of the predictive effects.
- **number of total snps with non-zero effects:** The total number of SNPs in PRS weights file that have non-zero effects.
- **path for outputs in s1:** For example, `/Mypath/s1_fixG_result`. If step 1 is performed on each chromosome separately, please save the outputs as `/Mypath/s1_fixG_result_1`,`/Mypath/s1_fixG_result_2`...`/Mypath/s1_fixG_result_22`.
- **Optional arguments:** `covar` specifies covariates to be adjusted in the model (other than `Tr`). For example, specify `covar="age,gender"` to include `age` and `gender` in the model (`Tr` will be included automatically). If not specified, only `Tr` will be included.


#### Step 2: parameter tuning
```
Rscript s2_fixG.R [file for phenotype and covariates] [initial values] [path for outputs in s1] [path for outputs in s2] [indicator of whether files for all 22 chromosomes are separated] \
covar="age,gender"
```
- **file for phenotype and covariates:** The following columns are required: `ID`, `Y`, covariates (including `Tr`), , `prs_g`, which represent the ID for each individual, the drug response phenotype to be predicted, other covariates (must have the treatment variable `Tr`, the PRS calculated from the weights from disease GWAS. Column names are needed.
- **initial values:** `PRS` or `zero`, indicating the starting values of the predictive effects.
- **path for outputs in s1:** The output from step 1. For example, `/Mypath/s1_fixG_result`. If step 1 is performed on each chromosome separately and the outputs are saved as `/Mypath/s1_fixG_result_1`,`/Mypath/s1_fixG_result_2`...`/Mypath/s1_fixG_result_22`, just need to input `/Mypath/s1_fixG_result` here and the algorithm will add `_chr` automatically.
- **path for outputs in s2:** The output from step 2. For example, `/Mypath/s2_fixG_result`.
- **indicator of whether files for all 22 chromosomes are separated:** `yes` or `no`, indicating whether step 1 is performed on all chromosomes or on each chromosome separately.
- **Optional arguments:** `covar` specifies covariates to be adjusted in the model (other than `Tr`). 


#### Step 3: rerun the algorithm and get the output with the best parameters
```
Rscript s3.R [file for PRS weights from disease GWAS] [file for phenotype and covariates] [file for genotype] [file for SNP information] [initial values] [number of total snps with non-zero effects] [path for outputs in s2] [path for outputs in s3] \
covar="age,gender"
```
- The inputs `file for PRS weights from disease GWAS`, `file for phenotype and covariates`, `file for genotype`, `file for SNP information`, `initial values`, `number of total snps with non-zero effects`, and the optional argument `covar` are the same as those in `s1_fixG.R`.  
- **path for outputs in s2:** The output from step 2. For example, `/Mypath/s2_fixG_result`.
- **path for outputs in s3:** The output from step 3. For example, `/Mypath/s3_fixG_result`. If step 3 is performed on each chromosome separately, please save the outputs as `/Mypath/s3_fixG_result_1`,`/Mypath/s3_fixG_result_2`...`/Mypath/s3_fixG_result_22`.

#### Step 4: sum the PRS over all chromosomes for the testing samples (optional)
If step 1 and 3 are performed on each chromosome separately, choose to run the step 4 to sum the PRS over all 22 chromosomes for the testing samples:
```
Rscript s4_fixG.R [initial values] [path for outputs in s3] [path for outputs in s4]
```
- **initial values:** `PRS` or `zero`, indicating the starting values of the predictive effects.
- **path for outputs in s3:** The output from step 3. If the outputs are saved as `/Mypath/s3_fixG_result_1`,`/Mypath/s3_fixG_result_2`...`/Mypath/s3_fixG_result_22`, just need to input `/Mypath/s3_fixG_result` here and the algorithm will add `_chr` automatically.
- **path for outputs in s4:** The output from step 4. For example, `/Mypath/s4_fixG_result`.

  
The choice of M5-M6 corresponds to the following arguments:
| Model | initial values 
| :--- | :--- |
| M5 | `zero` |
| M6 | `PRS` |




