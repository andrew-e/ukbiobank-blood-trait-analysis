#PBS -l select=1:ncpus=4:mem=12gb
#PBS -l walltime=2:0:0
#PBS -J 0-620

module load anaconda3/personal
source activate r4

cd /rds/general/user/are20/ephemeral/ml2000/$MARKER/
SNPS=($(ls))
SNP=${SNPS[$PBS_ARRAY_INDEX]}
SNP=${SNP:4}
cd -

echo "Running job for $MARKER and $ETHNICITY and $SNP"

Rscript /rds/general/user/are20/home/ukbiobank-blood-trait-analysis/ml/linear_regression_for_snp.R $MARKER $ETHNICITY $SNP

