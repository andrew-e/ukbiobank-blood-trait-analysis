#PBS -l select=1:ncpus=32:mem=512gb
#PBS -l walltime=6:0:0
#PBS -J 1-22

module load plink
#MARKERS=(baso eo hct hgb lymph mchc mch mcv mono mpv neut pct rbc rdw_cv wbc)
#MARKER="${MARKERS[PBS_ARRAY_INDEX]}"
#MARKER=pct

CHROMOSOME=$PBS_ARRAY_INDEX

BASE_DIR=/rds/general/user/are20/home/ukbiobank-blood-trait-analysis
DATA_DIR=$BASE_DIR/data/ml/snps

SNPS_PER_CHR=$DATA_DIR/${MARKER}/all_surrounding_snps_for_chr_${CHROMOSOME}.txt
OUTPUT_DIR=/rds/general/user/are20/ephemeral/ml/${MARKER}/extracted_${CHROMOSOME}

plink \
  --recode 12 AD \
  --bed /rds/general/project/uk-biobank-2018/live/reference/sdata_12032018/ukb_imp_chr${CHROMOSOME}.bed \
  --bim /rds/general/project/uk-biobank-2018/live/reference/sdata_12032018/ukb_imp_chr${CHROMOSOME}.bim \
  --fam $BASE_DIR/data/working_fam_v6.fam \
  --extract $SNPS_PER_CHR \
  --out $OUTPUT_DIR

