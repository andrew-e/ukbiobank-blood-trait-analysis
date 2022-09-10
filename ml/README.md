# Steps for running lassosum

#### 1. `prepare_snp_range.sh` 

This script iterates through the original score file takes the surrounding 2000 SNPs for each entry, and extracts and recodes it for use in R. This is then used for both Bayesian and Linear regressions.

#### 2. `bayesian_regression_for_snp.R`

Takes data from step 1, and runs a bayesian regression using the package rstanarm.  Uses summary statistic data to calculate priors, and principal component data to add as covaites for the regresion

The output is a new file with one entry per SNP, which can be used to create a score file for PLINK

Called per ethnicity, marker and SNP, via `bayesian_regression.sh`

#### 3. `linear_regression_for_snp.R`

Takes data from step 1, and runs a linear regression.  Uses principal component data to add as covaites for the regresion

The output is a new file with one entry per SNP, which can be used to create a score file for PLINK

Called per ethnicity, marker and SNP, via `bayesian_regression.sh`

#### 4. `data_prep_for_plink.R`

Takes the data from steps 2 and 3, and consolidates the information to a format that can be ingested by PLINK --score parameter, so a new PGS can be calculated on the updated information 

#### 5. `plink_script.sh and sum_chromosome_scores.sh`

Recalculates the PGS using the new score file per chromosome, then sums up all the scores per chromosome so it can be used to find correlations to the existing GAM residuals, calculated in `biobank` directory
