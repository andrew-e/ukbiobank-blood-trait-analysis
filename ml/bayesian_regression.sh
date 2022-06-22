#PBS -l select=1:ncpus=1:mem=256gb
#PBS -l walltime=4:0:0

echo "Running job for $MARKER and $ETHNICITY"
SNP_DIR=CHR_1_SNP_rs174574

module load anaconda3/personal
Rscript /rds/general/user/are20/home/ukbiobank-blood-trait-analysis/ml/bayesian_regression_for_snp.R $MARKER $ETHNICITY
