library(data.table)
library(dbplyr)
library(tidyverse)

setwd("/rds/general/user/are20/home/ukbiobank-blood-trait-analysis")

EPHEMERAL_BASE <- "/rds/general/user/are20/ephemeral"

all_bims <- fread(sprintf("%s/bfile_symlinks/linked_%s_%s_%s.bim", EPHEMERAL_BASE, ethnicity, marker, 1))
for (chromosome in 2:22) {
  linked_file <- sprintf("%s/bfile_symlinks/linked_%s_%s_%s.bim", EPHEMERAL_BASE, ethnicity, marker, chromosome)
  bim <- fread(linked_file)
  all_bims <- rbind(all_bims, bim)
}

ethnicities <- c("asian", "black", "chinese", "mixed")
blood_markers <- c("baso", "eo",  "hct", "hgb", "lymph", "mch", "mchc", "mcv", "mono", "mpv", "neut", "pct", "rbc", "rdw_cv", "wbc")

for (marker in blood_markers) {
  for (ethnicity in ethnicities) {
    print(paste(marker, ethnicity))
    if (file.exists(sprintf("results/ml2000/scores/lm_train_%s_%s.score_cleaned", marker, ethnicity))) {
      print("already calculated, move on.")
      next
    }
    
    SCORE_FILENAME <- sprintf("results/ml2000/scores/lm_train_%s_%s.score", marker, ethnicity)
    scores <- fread(SCORE_FILENAME)
    scores$abs_mean_intercept <- abs(scores$mean_intercept)
    
    print(paste("Number of duplicated SNPs", sum(duplicated(scores$best_performing_snp))))
    
    scores <- unique(setDT(scores)[order(best_performing_snp, -abs_mean_intercept)], by = "best_performing_snp")
    scores_joined_with_bim <- merge(scores, all_bims, by.x="best_performing_snp", by.y="V2", sort = F)
    # Adding column based on other column:
    scores_joined_with_bim <- scores_joined_with_bim %>% mutate(allele = if_else(mean_intercept >= 0, V6, V5))
    
    cleaned_scores <- data.frame(scores_joined_with_bim$best_performing_snp, scores_joined_with_bim$allele, scores_joined_with_bim$abs_mean_intercept)
    
    write.table(cleaned_scores, sprintf("results/ml2000/scores/lm_train_%s_%s.score_cleaned", marker, ethnicity), quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
}

