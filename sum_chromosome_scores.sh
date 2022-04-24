#!/bin/bash
PHENOTYPE=$1

MERGED_CHROMOSOMES=results/plink/$PHENOTYPE/all_chromosome_scores.txt
RESULT_FILE=results/plink/$PHENOTYPE/sum_chromosome_scores.txt

CHROMOSOME_FILES=$(ls results/plink/$PHENOTYPE/*.profile | wc -l)
if [[ $CHROMOSOME_FILES != 22 ]]; then
  echo "ERROR: $PHENOTYPE Some chromosome files missing"
  exit 1
fi
if [ -f "$RESULT_FILE" ]; then
  echo "$PHENOTYPE Already calculated, skipping"
  exit 1
fi

cat results/plink/$PHENOTYPE/*.profile > $MERGED_CHROMOSOMES

echo "eid_19266 score snp_count" > $RESULT_FILE
#TODO: add the CNT sum column here too, then divide the score by the CNT sum total number
awk '{b[$1]+=$4; a[$1]+=$6}END{for(i in a) print i,a[i],b[i]}' $MERGED_CHROMOSOMES >> $RESULT_FILE

rm $MERGED_CHROMOSOMES
