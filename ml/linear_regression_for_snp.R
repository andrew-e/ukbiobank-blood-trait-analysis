library(data.table)
library(rstanarm)
library(loo)
library(ggplot2)
options(mc.cores = parallel::detectCores())

setwd("/rds/general/user/are20/home/ukbiobank-blood-trait-analysis")
args = commandArgs(trailingOnly=TRUE)
marker <- args[1]
ethnicity <- args[2]
snp <- args[3]

EPHEMERAL_DIR <- paste0("/rds/general/user/are20/ephemeral/ml2000/", marker, "/SNP_", snp)
SCORE_FILENAME <- sprintf("results/ml2000/scores/lm_train_%s_%s.score", marker, ethnicity)

if (!file.exists(SCORE_FILENAME)) {
  file.create(SCORE_FILENAME)
  score_file_headers <- paste("best_performing_snp", "unadjusted_allele", "mean_intercept", "original_snp", "number_of_stans_skipped")
  write(score_file_headers, file = SCORE_FILENAME, append = TRUE)
}

#skip if already calculated...
score_file <- fread(SCORE_FILENAME)
if (any(score_file$original_snp == snp)) {
  print(paste(snp, "already calculated and stored in score file, skipping."))
  q()
}

principcal_components <- fread("data/plink/UKB_app_69328_genetic_PCs.txt")

bim_subset <- fread(paste0(EPHEMERAL_DIR, "/surrounding_snps.txt"))

plink_data <- fread(paste0(EPHEMERAL_DIR, "/plink.raw"))
plink_data = as.data.frame.matrix(plink_data)
plink_data <- plink_data[, !duplicated(colnames(plink_data))]
plink_data <- as.data.frame(plink_data)

residuals <- fread(sprintf("results/residuals/train_%s_%s_residuals_scaled.csv", ethnicity, marker), sep=',', header=TRUE, stringsAsFactors=FALSE)

calculateResults <- function(lm) {
  if (is.null(lm)) {
    return(c(NA, NA, NA))
  }
  
  r2 <- summary(lm)$r.squared
  p_value <- summary(lm)$coefficients[,4][2]
  effect_estimate <- summary(lm)$coefficients[,1][2]
  
  names(effect_estimate) <- NULL
  names(p_value) <- NULL
  
  return(c(effect_estimate, p_value, r2))
}

calculateLM <- function(surrounding_snp) {
  snp_no_trail <- gsub('.{2}$', '', surrounding_snp)
  
  plink_subset <- subset(plink_data_ethnicity, select = c("id", "eid_19266", "scaled_residuals", surrounding_snp))
  plink_subset <- merge(plink_subset, principcal_components, by.x="id", by.y="eid")
  
  #also add if there is less than 1% of individuals with non-0 (or differing snps)
  if(mean(plink_subset[[surrounding_snp]], na.rm=TRUE) < 0.02) {
    print(paste("Not enough variation in", snp_no_trail, "mean is", mean(plink_subset[[surrounding_snp]], na.rm=TRUE) ,"skipping."))
    return()
  }
  
  gc()
  pc_names <- colnames(plink_subset)[-4:-1]
  glm_formula = paste0("scaled_residuals ~ ", paste0("`", surrounding_snp, "`"), " + ", paste(pc_names, collapse=" + "))
  

  print(paste("Calculating for ", snp_no_trail))
  model <- lm(as.formula(glm_formula), data = plink_subset)

  return(model)
}

snp_list <- colnames(plink_data)[-6:-1]
plink_data_ethnicity <- merge(residuals, plink_data, by.x="eid_19266", by.y="FID")

stan_glm_pcs <- lapply(snp_list, calculateLM)
stan_results <- lapply(stan_glm_pcs, calculateResults)

stan_results <- data.frame( effect= sapply( stan_results, "[", 1), p_value = sapply( stan_results, "[", 2), 
                            r2 =sapply( stan_results, "[", 3) )

snp_list_no_trail <- lapply(snp_list, function(snp) {
  gsub('.{2}$', '', snp)
})

stan_results$snp <- unlist(snp_list_no_trail)
complete_stan_results <- stan_results[complete.cases(stan_results),]

good_p_values <- complete_stan_results[complete_stan_results$p_value < 0.05,]
lowest_p_performing_snp <- stan_results[which.min(stan_results$p_value),]

best_performing_snp <- good_p_values[which.max(good_p_values$r2),]
number_of_stans_skipped <- sum(is.na(stan_results$r2), na.rm = TRUE)

print("Best performing SNP with r2 (incl pcs): ")
print(best_performing_snp)

bim_entry <- bim_subset[bim_subset$V2 == best_performing_snp$snp]

score_file_entry <- paste(ifelse(is.na(best_performing_snp$snp), "NA", best_performing_snp$snp), ifelse(is.na(bim_entry$V6), "NA", bim_entry$V6), ifelse(is.na(best_performing_snp$effect), "NA", best_performing_snp$effect), snp, number_of_stans_skipped)
write(score_file_entry, file = SCORE_FILENAME, append = TRUE)
