#PBS -l select=1:ncpus=8:mem=384gb
#PBS -l walltime=6:0:0
#PBS -J 1-22

echo "Running job for $MARKER and $ETHNICITY"

module load anaconda3/personal
Rscript /rds/general/user/are20/home/ukbiobank-blood-trait-analysis/lassosum/lassosum_pipeline.R $PBS_ARRAY_INDEX $MARKER $ETHNICITY
