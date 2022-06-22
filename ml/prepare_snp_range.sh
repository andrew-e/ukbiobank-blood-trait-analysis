#PBS -l select=1:ncpus=1:mem=128gb
#PBS -l walltime=6:0:0
#PBS -J 0-14

module load plink

BASE_DIR=/rds/general/user/are20/home/ukbiobank-blood-trait-analysis
PHENOTYPES=(baso eo hct hgb lymph mchc mch mcv mono mpv neut pct rbc rdw_cv wbc)
PHENOTYPE="${PHENOTYPES[PBS_ARRAY_INDEX]}"
#PHENOTYPE=pct

echo "Starting $PHENOTYPE..."
EPH_DIR=/rds/general/user/are20/ephemeral/ml/${PHENOTYPE}

SCORE_FILE=/rds/general/user/are20/home/plink/scores/reduced_${PHENOTYPE}_my.score_condind_common
SCORE_SNP_LIST=$(wc -l < $SCORE_FILE)
echo "Calculating $SCORE_SNP_LIST SNPs"

for score_line in $(seq 1 $SCORE_SNP_LIST)
do
  SNP=$(sed "${score_line}q;d" $SCORE_FILE | awk '{print $1}')
  echo "Finding surrounding SNPs for $SNP..."

  if compgen -G "${EPH_DIR}/*${SNP}*/plink.raw" > /dev/null; then
    echo "$SNP is already calculated, skipping."
    continue
  fi

  for CHROMOSOME in $(seq 1 22)
  do
    BIM_FILE=/rds/general/project/uk-biobank-2018/live/reference/sdata_12032018/ukb_imp_chr${CHROMOSOME}.bim
    PADDED_SNP=" ${SNP} "
    SNPS=$(LC_ALL=C fgrep -A 500 -B 500 $PADDED_SNP $BIM_FILE)

    if [[ ! -z $SNPS ]]; then
      SNP_DIR=$EPH_DIR/SNP_${SNP}
      SURROUNDING_SNPS=$SNP_DIR/surrounding_snps.txt

      mkdir $SNP_DIR
      echo "$SNPS" > $SURROUNDING_SNPS

      echo "Preparing raw file for $SNP"
      plink \
        --recode 12 A \
        --bed /rds/general/project/uk-biobank-2018/live/reference/sdata_12032018/ukb_imp_chr${CHROMOSOME}.bed \
        --bim /rds/general/project/uk-biobank-2018/live/reference/sdata_12032018/ukb_imp_chr${CHROMOSOME}.bim \
        --fam $BASE_DIR/data/working_fam_v6.fam \
        --extract $SURROUNDING_SNPS \
        --out ${SNP_DIR}/plink
      break
    fi
  done
  echo "Completed prep for $SNP"
done


