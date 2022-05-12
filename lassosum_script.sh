#PBS -l select=1:ncpus=10:mem=256gb
#PBS -l walltime=24:0:0

module load anaconda3/personal
Rscript /rds/general/user/are20/home/ukbiobank-blood-trait-analysis/lassosum_calculations.R

