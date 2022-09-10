#PBS -lselect=1:ncpus=4:mem=2gb
#PBS -lwalltime=6:0:0

module load plink

CHROMOSOME=1

#for MARKER in baso eo hct hgb lymph mchc mch mcv mono mpv neut pct rbc rdw_cv wbc
RESULTS_DIR=/rds/general/user/are20/home/ukbiobank-blood-trait-analysis/results/ml2000/plink/$MARKER/$ETHNICITY
mkdir -p $RESULTS_DIR

echo "Phenotype: $MARKER, chromosome: $CHROMOSOME"

#if [ -f "$RESULTS_DIR/chromosome_$CHROMOSOME.profile" ]; then
#  echo "Already calculated, skipping"
#  exit 0
#fi

plink \
  --memory 96000 \
  --bed /rds/general/project/uk-biobank-2018/live/reference/sdata_12032018/ukb_imp_chr${CHROMOSOME}.bed \
  --bim /rds/general/project/uk-biobank-2018/live/reference/sdata_12032018/ukb_imp_chr${CHROMOSOME}.bim \
  --fam /rds/general/user/are20/home/ukbiobank-blood-trait-analysis/data/working_fam_v6.fam \
  --exclude /rds/general/user/are20/home/plink/duplicates/UKB_dupl_vars_to_exclude_chr${CHROMOSOME}.txt \
  --set-missing-var-ids @:#_\$1_\$2 \
  --score /rds/general/user/are20/home/ukbiobank-blood-trait-analysis/results/ml2000/scores/${MARKER}_${ETHNICITY}.score_cleaned sum double-dosage \
  --out $RESULTS_DIR/chromosome_${CHROMOSOME} 

