#PBS -l select=1:ncpus=1:mem=256gb
#PBS -l walltime=4:0:0

echo "Running job for $MARKER and $ETHNICITY"

module load anaconda3/personal
Rscript /rds/general/user/are20/home/ukbiobank-blood-trait-analysis/lassosum_merge_chromosomes.R $MARKER $ETHNICITY
