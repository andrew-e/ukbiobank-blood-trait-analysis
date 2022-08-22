#PBS -l select=1:ncpus=1:mem=128gb
#PBS -l walltime=2:0:0

echo "Running job for $MARKER and $ETHNICITY"

module load anaconda3/personal
Rscript /rds/general/user/are20/home/ukbiobank-blood-trait-analysis/lassosum/lassosum_merge_chromosomes.R $MARKER $ETHNICITY
