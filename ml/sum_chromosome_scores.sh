#!/bin/bash
MARKER=$1
ETHNICITY=$2

BASE_DIR=/rds/general/user/are20/home/ukbiobank-blood-trait-analysis
RESULTS_DIR=$BASE_DIR/results/ml2000/plink/lm_train_${MARKER}/$ETHNICITY

MERGED_CHROMOSOMES=$RESULTS_DIR/all_chromosome_scores.txt
RESULT_FILE=$RESULTS_DIR/sum_chromosome_scores.txt

CHROMOSOME_FILES=$(ls $RESULTS_DIR/*.profile | wc -l)
if [[ $CHROMOSOME_FILES != 22 ]]; then
  echo "ERROR: $MARKER $ETHNICITY Some chromosome files missing, only $CHROMOSOME_FILES/22 present.  Exiting."
  exit 1
fi

if [ -f "$RESULT_FILE" ]; then
  echo "$MARKER $ETHNICITY Already calculated, skipping"
  exit 1
fi

cat $RESULTS_DIR/*.profile > $MERGED_CHROMOSOMES

echo "eid_19266 score snp_count" > $RESULT_FILE
awk '{b[$1]+=$4; a[$1]+=$6}END{for(i in a) print i,a[i],b[i]}' $MERGED_CHROMOSOMES >> $RESULT_FILE

rm $MERGED_CHROMOSOMES
