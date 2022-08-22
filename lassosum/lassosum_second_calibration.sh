#PBS -l select=1:ncpus=32:mem=1024gb
#PBS -l walltime=12:0:0

echo "Running job for $MARKER and $ETHNICITY"

module load anaconda3/personal
Rscript /rds/general/user/are20/home/ukbiobank-blood-trait-analysis/lassosum/lassosum_second_calibration.R $MARKER $ETHNICITY
