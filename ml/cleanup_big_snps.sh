module load plink

PHENOTYPES="baso eo hct hgb lymph mchc mch mcv mono mpv neut pct rbc rdw_cv wbc"

for MARKER in $PHENOTYPES
do
  echo $MARKER
  cd /rds/general/user/are20/ephemeral/ml/$MARKER/
  SNPS=$(ls)

  for SNP in $SNPS
  do
    NUM_SNPS=$(wc -l $SNP/surrounding_snps.txt | awk -F' ' '{ print $1 }')
    if [ "$NUM_SNPS" -gt "1100" ]; then
      echo "should delete $SNP cause $NUM_SNPS";
      rm -r $SNP
    fi
  done
done

