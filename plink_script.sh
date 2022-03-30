## The environment variable $CHROMOSOME gives each
## subjob's index within the array
module load plink

UKB_DIR=/rds/general/project/uk-biobank-2018/live/reference/sdata_12032018
RESULTS_DIR=/rds/general/user/are20/home/ukbiobank-blood-trait-analysis/results/plink/$PHENOTYPE

mkdir -p $RESULTS_DIR 
for CHROMOSOME in $(seq 1 22)
do
  echo "Chromosome $CHROMOSOME"
  plink \
    --memory 96000 \
    --bed /rds/general/project/uk-biobank-2018/live/reference/sdata_12032018/ukb_imp_chr${CHROMOSOME}.bed \
    --bim /rds/general/project/uk-biobank-2018/live/reference/sdata_12032018/ukb_imp_chr${CHROMOSOME}.bim \
    --fam /rds/general/project/chadeau_ukbb_folder/live/data/project_data/Genetic_data_extraction/plink_files/ukb_imp.fam \
    --exclude /rds/general/user/are20/home/plink/duplicates/UKB_dupl_vars_to_exclude_chr${CHROMOSOME}.txt \
    --set-missing-var-ids @:#_\$1_\$2 \
    --score /rds/general/user/are20/home/plink/scores/reduced_${PHENOTYPE}_my.score_condind_common sum double-dosage \
    --out $RESULTS_DIR/chromosome_${CHROMOSOME} 

done
