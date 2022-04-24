base_dir = "/rds/general/user/are20/home/ukbiobank-blood-trait-analysis/results"

ethnicities = c("black", "asian", "chinese", "white", "mixed")
colours <- c("tomato1", "rosybrown3", "slateblue", "wheat", "skyblue1")
blood_markers = c("baso", "eo",  "hct", "hgb", "lymph", "mch", "mcv", "mono", "mpv", "neut", "pct", "rbc", "rdw_cv", "wbc") #"mchc"

for (marker in blood_markers) {
  markers_for_ethnicity <- sapply(ethnicities, function(ethnicity) {
    plink_results <- fread(sprintf("%s/plink/%s/sum_chromosome_scores.txt", base_dir, marker), sep=' ', header=TRUE, stringsAsFactors=FALSE)
    residuals <- fread(sprintf("%s/residuals/%s_%s_residuals.csv", base_dir, ethnicity, marker), sep=',', header=TRUE, stringsAsFactors=FALSE)
    plink_results$eid_19266 <- as.numeric(plink_results$eid_19266)
    plink_results$snp_count <- as.numeric(plink_results$snp_count)
    plink_results <- transform(plink_results, weighted_score = score / snp_count)
    merged_results <- merge(residuals, plink_results, by="eid_19266", all.x = TRUE)
    
    #weighted_correlation_result <- cor(merged_results$weighted_score, merged_results$residuals, use = "pairwise.complete.obs") #isn't much better.
    unweighted_correlation_result <- cor(merged_results$score, merged_results$residuals, use = "pairwise.complete.obs") #return this with ethnicity/marker combo
    
    signif(unweighted_correlation_result, digits = 3)
  })
  print(markers_for_ethnicity)
  barplot(markers_for_ethnicity, col = colours, main=paste("Correlation of", marker, "By Ethnicity"))
}

#MCHC errors: "Error: No valid entries in --score file."

#TODO: standardise the PRS score around CNT first of all, then around 0, and SD as well.
#also standardise the residuals as well in the same way.

#cor(merged_results$score, merged_results$residuals, use = "pairwise.complete.obs")
#Maybe run a linear regression on this, not just a correlation.
