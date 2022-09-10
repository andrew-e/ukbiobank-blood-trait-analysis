#PBS -l select=1:ncpus=2:mem=20gb
#PBS -l walltime=2:0:0

module load anaconda3/personal

Rscript /rds/general/user/are20/home/ukbiobank-blood-trait-analysis/ml/data_prep_for_plink.R

