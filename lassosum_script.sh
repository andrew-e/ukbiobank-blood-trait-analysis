#PBS -l select=1:ncpus=16:mem=256gb
#PBS -l walltime=12:0:0

echo "Running job for $MARKER"

module load anaconda3/personal
for i in {1..22}; do
  Rscript /rds/general/user/are20/home/ukbiobank-blood-trait-analysis/lassosum_calculations.R $i $MARKER
done
