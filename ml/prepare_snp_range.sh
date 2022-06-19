#PBS -l select=1:ncpus=1:mem=56gb
#PBS -l walltime=6:0:0

module load plink

BASE_DIR=/rds/general/user/are20/home/ukbiobank-blood-trait-analysis
DATA_DIR=$BASE_DIR/data/ml/snps

#for PHENOTYPE in baso eo hct hgb lymph mchc mch mcv mono mpv neut pct rbc rdw_cv wbc
#for PHENOTYPE in mcv mono mpv neut rbc rdw_cv wbc
for PHENOTYPE in pct 
do
  echo "Starting $PHENOTYPE..."

  SCORE_FILE=/rds/general/user/are20/home/plink/scores/reduced_${PHENOTYPE}_my.score_condind_common

  SCORE_SNP_LIST=$(wc -l < $SCORE_FILE)
  for score_line in $(seq 1 $SCORE_SNP_LIST)
  do
    SNP=$(sed "${score_line}q;d" $SCORE_FILE | awk '{print $1}')
    echo "Calculating $SNP..."
    if compgen -G "${DATA_DIR}/${PHENOTYPE}/*${SNP}*.txt" > /dev/null; then
      echo "$SNP is already calculated, skipping."
      continue
    fi

    for CHROMOSOME in $(seq 1 22)
    do
      SNPS_FILE=$DATA_DIR/${PHENOTYPE}/surrounding_snp_for_chr_${CHROMOSOME}_snp_${SNP}.txt
      BIM_FILE=/rds/general/project/uk-biobank-2018/live/reference/sdata_12032018/ukb_imp_chr${CHROMOSOME}.bim

      SNPS=$(LC_ALL=C fgrep -A 500 -B 500 $SNP $BIM_FILE)
      if [[ ! -z $SNPS ]]; then
        echo "$SNPS" > $SNPS_FILE
        echo "Completed file for SNP $SNP"
        break
      fi
    done
  done
  
  echo "Now, shoving them all together in one file..."
  for CHROMOSOME in $(seq 1 22)
  do
    SNPS_PER_CHR=$DATA_DIR/${PHENOTYPE}/all_surrounding_snps_for_chr_${CHROMOSOME}.txt
    cat $DATA_DIR/${PHENOTYPE}/surrounding_snp_for_chr_${CHROMOSOME}_snp_* | awk '!x[$0]++' > $SNPS_PER_CHR 
    rm $DATA_DIR/${PHENOTYPE}/surrounding_snp_for_chr_${CHROMOSOME}_snp_*    

    #plink \
    #  --recode 12 AD \
    #  --bed /rds/general/project/uk-biobank-2018/live/reference/sdata_12032018/ukb_imp_chr${CHROMOSOME}.bed \
    #  --bim /rds/general/project/uk-biobank-2018/live/reference/sdata_12032018/ukb_imp_chr${CHROMOSOME}.bim \
    #  --fam $BASE_DIR/data/working_fam_v6.fam \
    #  --extract $SNPS_PER_CHR \
    #  --out $DATA_DIR/${PHENOTYPE}/extracted_snps_for_chr_${CHROMOSOME}
  done

done

