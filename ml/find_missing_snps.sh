#PBS -l select=1:ncpus=1:mem=128gb
#PBS -l walltime=6:0:0

module load plink

BASE_DIR=/rds/general/user/are20/home/ukbiobank-blood-trait-analysis
PHENOTYPES=(baso eo hct hgb lymph mchc mch mcv mono mpv neut pct rbc rdw_cv wbc)

for PHENOTYPE in $PHENOTYPES
do
  echo "For $PHENOTYPE"
  SCORE_FILE=/rds/general/user/are20/home/plink/scores/reduced_${PHENOTYPE}_my.score_condind_common
  SCORE_SNP_LIST=$(wc -l < $SCORE_FILE)
  for score_line in $(seq 1 $SCORE_SNP_LIST)
  do
    SNP=$(sed "${score_line}q;d" $SCORE_FILE | awk '{print $1}')
    BIM_FILES=/rds/general/project/uk-biobank-2018/live/reference/sdata_12032018/ukb_imp_chr*.bim
    
    SNPS=$(LC_ALL=C fgrep "${SNP}$(printf '\t')" $BIM_FILES)
    if [[ -z $SNPS ]]; then
      MISSING_SNPS+="$SNP\n"
    fi
  done

  echo "Missing SNPs for $PHENOTYPE:"
  echo $MISSING_SNPS
done

