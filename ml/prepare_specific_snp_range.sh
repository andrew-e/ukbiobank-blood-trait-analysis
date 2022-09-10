#PBS -l select=1:ncpus=1:mem=20gb
#PBS -l walltime=1:0:0

module load plink

BASE_DIR=/rds/general/user/are20/home/ukbiobank-blood-trait-analysis
PHENOTYPES=(baso eo hct hgb lymph mchc mch mcv mono mpv neut pct rbc rdw_cv wbc)
PHENOTYPE=rdw_cv
SNP=rs268
CHROMOSOME=8

EPH_DIR=/rds/general/user/are20/ephemeral/ml2000/${PHENOTYPE}

SCORE_FILE=/rds/general/user/are20/home/plink/scores/reduced_${PHENOTYPE}_my.score_condind_common
SCORE_SNP_LIST=$(wc -l < $SCORE_FILE)
echo "Calculating $SCORE_SNP_LIST SNPs"

echo "Finding surrounding SNPs for $SNP..."

if compgen -G "${EPH_DIR}/*${SNP}*/plink.raw" > /dev/null; then
  echo "$SNP is already calculated, skipping."
  continue
fi

BIM_FILE=/rds/general/project/uk-biobank-2018/live/reference/sdata_12032018/ukb_imp_chr${CHROMOSOME}.bim

SNPS=$(LC_ALL=C fgrep -A 1000 -B 1000 "${SNP}$(printf '\t')" $BIM_FILE)

if [[ ! -z $SNPS ]]; then
  SNP_DIR=$EPH_DIR/SNP_${SNP}
  SURROUNDING_SNPS=$SNP_DIR/surrounding_snps.txt

  mkdir -p $SNP_DIR
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

echo "Completed prep for $SNP"
