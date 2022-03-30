#!/bin/bash
PHENOTYPE=$1

MERGED_CHROMOSOMES=results/plink/$PHENOTYPE/all_chromosome_scores.txt
RESULT_FILE=results/plink/$PHENOTYPE/sum_chromosome_scores.txt

cat results/plink/$PHENOTYPE/*.profile > $MERGED_CHROMOSOMES

echo "eid_19266 score" > $RESULT_FILE
awk '{a[$1]+=$6}END{for(i in a) print i,a[i]}' $MERGED_CHROMOSOMES >> $RESULT_FILE

rm $MERGED_CHROMOSOMES
