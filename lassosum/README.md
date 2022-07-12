# Steps for running lassosum

#### 1. `lassosum_pipeline_only.R` 
This is the `lassosum.pipeline` function, and is run once per ethnicity, marker, and chromosome.  The data needed for this is:
* summary statistics (of marker), including position, major and minor allele, correlation (calculated from ss)
* ref, test, and keep subsets of the population.  These changed a bit and are a bit hardcoded currently for white ethnicity sizes.

This was called via `lassosum_pipeline.sh`

#### 2. `lassosum_merge_chromosomes.R`

Simple script that took the 22 results from lassosum.pipeline and merged them together

Called per ethnicity and marker, via `lassosum_merge.sh`

#### 3. `lassosum_validation_pheno.R`

This prepares the `pheno` (FID, IID, and phenotype which are the residuals as calculated previously from the GAM) and `covar` which is all 40 columns from the principal components data, and does not include the phenotype, as it's included in the `pheno`.

Once it runs the (first pass) validation, it takes a subst of the results based on the best performing `s` and `lambda`, as returned by the validate step.

Called via `lassosum_validation_pheno.sh`

#### 4. `lassosum_second_calibration.R`

This does very similar steps to 3.  Prepares the `pheno` and `covar` data frames, and runs validate on the subset calculated by 3.

Once that is complete, the results can be investigated.  There is a bunch of commented out code at the bottom that looked at the results.

Called via `lassosum_second_calibration.sh`

