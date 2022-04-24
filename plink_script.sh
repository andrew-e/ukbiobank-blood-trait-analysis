#PBS -lselect=1:ncpus=16:mem=10gb
#PBS -lwalltime=24:0:0
#PBS -J 1-22  


## The environment variable $CHROMOSOME gives each
## subjob's index within the array
module load plink

#for PHENOTYPE in baso baso_p eo eo_p hct hgb hlr hlr_p irf lymph lymph_p mchc mch mcv mono mono_p mpv mrv mscv neut neut_p pdw plt rbc rdw_cv ret ret_p wbc
for PHENOTYPE in baso eo hct hgb lymph mchc mch mcv mono mpv neut pct rbc rdw_cv wbc
do
  RESULTS_DIR=/rds/general/user/are20/home/ukbiobank-blood-trait-analysis/results/plink/$PHENOTYPE

  mkdir -p $RESULTS_DIR 
  for CHROMOSOME in $(seq 1 22)
  do
    echo "Phenotype: $PHENOTYPE, chromosome: $CHROMOSOME"

    if [ -f "$RESULTS_DIR/chromosome_$CHROMOSOME.profile" ]; then
      echo "Already calculated, skipping"
      continue
    fi

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
done
