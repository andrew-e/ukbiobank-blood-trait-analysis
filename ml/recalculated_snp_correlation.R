setwd("/rds/general/user/are20/home/ukbiobank-blood-trait-analysis")

ethnicities = c("black", "asian", "chinese", "mixed")
colours <- c("tomato1", "rosybrown3", "slateblue", "wheat", "skyblue1")
blood_markers = c("baso", "eo",  "hct", "hgb", "lymph", "mch", "mchc", "mcv", "mono", "mpv", "neut", "pct", "rbc", "rdw_cv", "wbc")

for (marker in blood_markers) {
  original_scores <- fread(sprintf("/rds/general/user/are20/home/plink/scores/reduced_%s_my.score_condind_common", marker))
  max_original_score <- max(original_scores$V3)
  
  markers_for_ethnicity <- sapply(ethnicities, function(ethnicity) {
    SCORE_FILENAME <- sprintf("results/ml2000/scores/%s_%s.score", marker, ethnicity)
    scores <- fread(SCORE_FILENAME)
    
    full_rows <- scores[complete.cases(scores),]
    dedup_rows <- unique(setDT(full_rows)[order(best_performing_snp)], by = "best_performing_snp")
    
    removed_rows <- (nrow(full_rows) - nrow(dedup_rows)) / nrow(full_rows)
    
    
    mean_regressions_skipped <- mean(scores$number_of_stans_skipped) / 2000
    sd_regressions_skipped <- sd(scores$number_of_stans_skipped)
    #effect_correlations <- cor(scores$mean_intercept, scores$original_effect_size, use = "pairwise.complete.obs")
    max_effect_diff <- max(abs(scores$mean_intercept)) - max_original_score
    
    plink_results <- fread(sprintf("results/ml2000/plink/lm_train_%s/%s/sum_chromosome_scores.txt", marker, ethnicity), sep=' ', header=TRUE, stringsAsFactors=FALSE)
    residuals <- fread(sprintf("results/residuals/train_%s_%s_residuals_scaled.csv", ethnicity, marker), sep=',', header=TRUE, stringsAsFactors=FALSE)
    
    plink_results$eid_19266 <- as.numeric(plink_results$eid_19266)
    plink_results$snp_count <- as.numeric(plink_results$snp_count)
    plink_results <- plink_results[complete.cases(plink_results),]
    merged_results <- merge(residuals, plink_results, by="eid_19266", all.x = TRUE)
    unweighted_correlation_result <- cor(merged_results$score, merged_results$scaled_residuals, use = "pairwise.complete.obs")
    
    return (c(cor = unweighted_correlation_result, mean_skipped = mean_regressions_skipped, sd_skipped = sd_regressions_skipped, max_effect_diff = max_effect_diff, na_rows = na_rows))
    return(removed_rows)
  })
  print(marker)
  
  write.csv(markers_for_ethnicity, sprintf("results/ml2000/final/lm_train_%s.csv", marker))
}





lm_skipped <- sapply(blood_markers, function(marker) {
  lm_effect_correlations <- fread(sprintf("results/ml2000/final/lm_%s.csv", marker))[2]
  lm_effect_list <- c(lm_effect_correlations$black, lm_effect_correlations$asian, lm_effect_correlations$chinese, lm_effect_correlations$mixed)
})
bayesian_cors <- sapply(blood_markers, function(marker) {
  bayesian_effect_correlations <- fread(sprintf("results/ml2000/final/%s.csv", marker))[2]
  bayesian_effect_list <- c(bayesian_effect_correlations$black, bayesian_effect_correlations$asian, bayesian_effect_correlations$chinese, bayesian_effect_correlations$mixed)
})

effect_cor_compare <- data.frame(matrix(NA, nrow = 4, ncol = 2))
effect_cor_compare[2] <- c(0.140190348, 0.078349813, 0.046022335, 0.077055078)
effect_cor_compare[1] <- c(0.21957299, 0.02700253, 0.24793637, 0.14293789)
rownames(effect_cor_compare) <- c("Asian", "Black", "Chinese", "Mixed")
colnames(effect_cor_compare) <- c("Bayesian", "Linear")
barplot(t(effect_cor_compare), beside=TRUE, ylab="Correlation Coefficient", col=c("rosybrown3", "slateblue"), main=sprintf("Correlation Between Original Effect Versus\n New Effect Across All Phenotypes"), legend = TRUE,  ylim = c(0,0.3))


colours <- c("skyblue1", "rosybrown3", "wheat", "slateblue")
blood_markers = c("baso", "eo",  "hct", "hgb", "lymph", "mch", "mchc", "mcv", "mono", "mpv", "neut", "pct", "rbc", "rdw_cv", "wbc")
for (ethnicity in ethnicities) {
  original_white_correlation <- fread("results/residuals/all_correlations.csv")$white
  original_correlations <- fread("results/residuals/all_correlations.csv")[[ethnicity]]
  bayesian_correlations <- fread(sprintf("results/ml2000/final/bayesian_all.csv"))[[ethnicity]]
  lm_correlations <- fread(sprintf("results/ml2000/final/lm_all_adjusted.csv"))[[ethnicity]]
  
  mean_correlations <- c("Baseline White" = mean(original_white_correlation), "Baseline Adj." = mean(original_correlations), "Bayesian" = mean(bayesian_correlations), "Linear" = mean(lm_correlations))
  
  max_value <- max(original_correlations, bayesian_correlations, lm_correlations)
  
  cor_compare <- data.frame(Linear = lm_correlations, Bayesian = bayesian_correlations, "Baseline Adj." = original_correlations, "Baseline White" = original_white_correlation, check.names = FALSE)
  rownames(cor_compare) <- blood_markers
  cor_compare <- cor_compare[order(nrow(cor_compare):1),]
  
  barplot(t(cor_compare), beside=TRUE, las = 2, horiz = TRUE, space = c(0, 2), legend = TRUE,  xlim = c(0,0.6),
          xlab="Correlation Coefficient",
          col=colours,
          main=sprintf("Comparing Correlations for %s Ethnicity per Phenotype", str_to_title(ethnicity))
  )
}

