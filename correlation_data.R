plink_results_file = "/rds/general/user/are20/home/ukbiobank-blood-trait-analysis/results/plink/pct/sum_chromosome_scores.txt"
white_pct_gam_file = "/rds/general/user/are20/home/ukbiobank-blood-trait-analysis/results/residuals/pct/white_platelet_count_residuals.csv"

white_pct_gam <- fread(white_pct_gam_file, sep=',', header=TRUE, stringsAsFactors=FALSE)
plink_results <- fread(plink_results_file, sep=' ', header=TRUE, stringsAsFactors=FALSE)
plink_results$eid_19266 <- as.numeric(plink_results$eid_19266)
merged_results <- merge(white_pct_gam, plink_results, by="eid_19266", all.x = TRUE)
