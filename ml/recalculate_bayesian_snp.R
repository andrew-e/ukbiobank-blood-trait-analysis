library(data.table)
library(rstanarm)
library(loo)
library(ggplot2)
options(mc.cores = parallel::detectCores())
#You can read up on bayesian linear regression and stan_glm here:https://cran.r-project.org/web/packages/rstanarm/rstanarm.pdf#page.108

setwd("/rds/general/user/are20/home/ukbiobank-blood-trait-analysis")
args = commandArgs(trailingOnly=TRUE)
marker <- args[1]
ethnicity <- args[2]
snp <- args[3]

EPHEMERAL_DIR <- paste0("/rds/general/user/are20/ephemeral/ml2000/", marker, "/SNP_", snp)
SCORE_FILENAME <- sprintf("results/ml2000/scores/%s_%s.recalc_score", marker, ethnicity)


if (!file.exists(SCORE_FILENAME)) {
  file.create(SCORE_FILENAME)
  score_file_headers <- paste("best_performing_snp", "unadjusted_allele", "mean_intercept", "original_snp", "original_effect_size", "number_of_stans_skipped")
  write(score_file_headers, file = SCORE_FILENAME, append = TRUE)
}

#skip if already calculated...
score_file <- fread(SCORE_FILENAME)
if (any(score_file$original_snp == snp)) {
  print(paste(snp, "already calculated and stored in score file, skipping."))
  q()
}

#Get summary statistics, to be used for priors in the bayesian linear regression.
summary_stat_columns <- c("ID", "CHR", "REF", "ALT", "EFFECT", "SE")
summary_statistic <- fread(sprintf("data/plink/summary_statistics/%s.assoc", marker), select = summary_stat_columns)
principcal_components <- fread("data/plink/UKB_app_69328_genetic_PCs.txt")

bim_subset <- fread(paste0(EPHEMERAL_DIR, "/surrounding_snps.txt"))

plink_data <- fread(paste0(EPHEMERAL_DIR, "/plink.raw"))
plink_data = as.data.frame.matrix(plink_data)
plink_data <- plink_data[, !duplicated(colnames(plink_data))]
plink_data <- as.data.frame(plink_data)

residuals <- fread(sprintf("results/residuals/%s_%s_residuals_scaled.csv", ethnicity, marker), sep=',', header=TRUE, stringsAsFactors=FALSE)

calculateResults <- function(stan_glm) {
  if (is.null(stan_glm)) {
    return(c(NA, NA, NA))
  }
  
  mean_intercept <- summary(stan_glm)[, "mean"][2]
  names(mean_intercept) <- NULL
  
  posterior_interval <- posterior_interval(stan_glm)[2,]
  posterior_interval_range <- posterior_interval[2] - posterior_interval[1]
  names(posterior_interval_range) <- NULL
  
  ss_res <- var(residuals(stan_glm))
  ss_total <- var(fitted(stan_glm)) + var(residuals(stan_glm))
  r2 <- 1 - (ss_res / ss_total)
  return(c(mean_intercept, r2, posterior_interval_range))
}

calculateBayesianGLM <- function(surrounding_snp, includePCs) {
  snp_no_trail <- gsub('.{2}$', '', surrounding_snp)
  
  plink_subset <- subset(plink_data_ethnicity, select = c("id", "eid_19266", "scaled_residuals", surrounding_snp))
  plink_subset <- merge(plink_subset, principcal_components, by.x="id", by.y="eid")
  
  #also add if there is less than 1% of individuals with non-0 (or differing snps)
  if(mean(plink_subset[[surrounding_snp]], na.rm=TRUE) < 0.02) {
    print(paste("Not enough variation in", snp_no_trail, "mean is", mean(plink_subset[[surrounding_snp]], na.rm=TRUE) ,"skipping."))
    return()
  }
  
  gc()
  ss_snp <- summary_statistic[summary_statistic$ID == snp_no_trail]
  pc_names <- colnames(plink_subset)[-4:-1]
  glm_formula = paste0("scaled_residuals ~ ", paste0("`", surrounding_snp, "`"), " + ", paste(pc_names, collapse=" + "))
  
  scale <- ss_snp$SE * sqrt(410000) #410000 was the sample size of the white pop
  location <- ss_snp$EFFECT
  print(paste("Calculating for ", snp_no_trail, location, scale))
  
  stan_model <- NULL
  if (nrow(ss_snp) == 1) {
    stan_model <- stan_glm(as.formula(glm_formula), data = plink_subset, prior = normal(location = location, scale = scale, autoscale = FALSE), refresh=0, chains=4, cores=12)
  }
  else {
    stan_model <- stan_glm(as.formula(glm_formula), data = plink_subset, refresh=0, chains=4, cores=12)
  }

  if (any(summary(stan_model)[, "Rhat"] > 1.1)) {
    print("There is an Rhat > 1.1, therfore didn't converge: ignoring.")
    return()
  }

  return(stan_model)
}

snp_list <- colnames(plink_data)[-6:-1]
plink_data_ethnicity <- merge(residuals, plink_data, by.x="eid_19266", by.y="FID")


stan_glm_pcs <- readRDS(paste0(EPHEMERAL_DIR, "/", ethnicity, "_stan_results.rds"))
stan_results <- lapply(stan_glm_pcs, calculateResults)

stan_results <- data.frame( mean_intercept= sapply( stan_results, "[", 1), 
                            r2 =sapply( stan_results, "[", 2),
                            posterior_interval_range =sapply( stan_results, "[", 3)
                          )
snp_list_no_trail <- lapply(snp_list, function(snp) {
  gsub('.{2}$', '', snp)
})

stan_results$snp <- unlist(snp_list_no_trail)
stan_results <- merge(stan_results, summary_statistic[ , c("ID", "EFFECT")], by.x="snp", by.y="ID", all.x=TRUE, sort=FALSE)

write.csv(stan_results, paste0(EPHEMERAL_DIR, "/", ethnicity, "_stan_results.csv"), row.names = FALSE)

best_performing_snp <- stan_results[which.max(stan_results$r2),]
number_of_stans_run <- sum(is.na(stan_results$r2), na.rm = TRUE)

print("Best performing SNP with r2 (incl pcs): ")
print(best_performing_snp)

bim_entry <- bim_subset[bim_subset$V2 == best_performing_snp$snp]

score_file_entry <- paste(best_performing_snp$snp, bim_entry$V6, best_performing_snp$mean_intercept, snp, best_performing_snp$EFFECT, number_of_stans_run)
write(score_file_entry, file = SCORE_FILENAME, append = TRUE)
