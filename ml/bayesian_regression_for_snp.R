library(data.table)
library(rstanarm)
options(mc.cores = parallel::detectCores())
#You can read up on bayesian linear regression and stan_glm here:https://cran.r-project.org/web/packages/rstanarm/rstanarm.pdf#page.108

setwd("/rds/general/user/are20/home/ukbiobank-blood-trait-analysis")
args = commandArgs(trailingOnly=TRUE)
marker <- args[1]
ethnicity <- args[2]
snp <- args[3]

EPHEMERAL_DIR <- paste0("/rds/general/user/are20/ephemeral/ml/", marker, "/SNP_", snp)
SCORE_FILENAME <- sprintf("results/ml/scores/%s_%s.score", marker, ethnicity)

#skip if already calculated...
if (!file.exists(SCORE_FILENAME)) {
  file.create(SCORE_FILENAME)
}
score_file <- fread(SCORE_FILENAME)
if (any(score_file$V4 == snp)) {
  print(paste(snp, "already calculated and stored in score file, skipping."))
  q()
}

residuals <- fread(sprintf("results/residuals/%s_%s_residuals.csv", ethnicity, marker), sep=',', header=TRUE, stringsAsFactors=FALSE)

calculateR2 <- function(stan_glm) {
  if (is.null(stan_glm)) {
    return(NA)
  }
  
  ss_res <- var(residuals(stan_glm))
  ss_total <- var(fitted(stan_glm)) + var(residuals(stan_glm))
  r2 <- 1 - (ss_res / ss_total)
  return(r2)
}

calculateBayesianGLM <- function(surrounding_snp, includePCs) {
  snp_no_trail <- gsub('.{2}$', '', surrounding_snp)
  ss_snp <- summary_statistic[summary_statistic$ID == snp_no_trail]
  
  plink_subset <- subset(plink_data_ethnicity, select = c("id", "eid_19266", "residuals", surrounding_snp))
  scale <- ss_snp$SE * sqrt(410000) #410000 was the sample size of the white pop
  location <- ss_snp$EFFECT
  
  glm_formula <- paste0("residuals ~ ", paste0("`", surrounding_snp, "`"))
  if (includePCs == T) {
    plink_subset <- merge(plink_subset, principcal_components, by.x="id", by.y="eid")
    pc_names <- colnames(plink_subset)[-4:-1]
    glm_formula = paste0("residuals ~ ", paste0("`", surrounding_snp, "`"), " + ", paste(pc_names, collapse=" + "))
  }
  
  if(sd(plink_subset[[surrounding_snp]], na.rm=TRUE) == 0) {
    print(paste("No variation in", snp_no_trail, "skipping."))
    return()
  }
  print(paste("Calculating for ", snp_no_trail, location, scale))
  
  stan_model <- NULL
  if (nrow(ss_snp) == 1) {
    stan_model <- stan_glm(as.formula(glm_formula), data = plink_subset, prior = normal(location = location, scale = scale, autoscale = FALSE), refresh = 0)
  }
  else {
    stan_model <- stan_glm(as.formula(glm_formula), data = plink_subset, refresh = 0)
  }
  return(stan_model)
}


plink_data <- fread(paste0(EPHEMERAL_DIR, "/plink.raw"))
plink_data = as.data.frame.matrix(plink_data)
plink_data <- plink_data[, !duplicated(colnames(plink_data))]
plink_data <- as.data.frame(plink_data)

snp_list <- colnames(plink_data)[-6:-1]
plink_data_ethnicity <- merge(residuals, plink_data, by.x="eid_19266", by.y="FID")

#Get summary statistics, to be used for priors in the bayesian linear regression.
summary_stat_columns <- c("ID", "CHR", "EFFECT", "SE")
summary_statistic <- fread(sprintf("data/plink/summary_statistics/%s.assoc", marker), select = summary_stat_columns)

#get PCs and plink raw output
principcal_components <- fread("data/plink/UKB_app_69328_genetic_PCs.txt")

stan_glm_pcs <- lapply(snp_list, calculateBayesianGLM, includePCs=T)

r2_pcs <- lapply(stan_glm_pcs, calculateR2)
stan_comparison <- data.frame(snp = unlist(snp_list), r2_pcs = unlist(r2_pcs))
best_performing_snp <- stan_comparison[which.max(stan_comparison$r2_pcs),]

print("Best performing SNP with r2 (incl pcs): ")
print(best_performing_snp)

best_snp_no_trail <- gsub('.{2}$', '', best_performing_snp$snp)
ss_snp <- summary_statistic[summary_statistic$ID == best_snp_no_trail]


if (nrow(ss_snp) == 1) {
  bim <- fread(sprintf("/rds/general/apps/RDS_COMPATIBILITY_LINKFARM/groupvol/med-bio/uk-biobank-2017/release_12032018/converted_data/imp_bgen1.1_plink/ukb_imp_chr%s_v3.bim", ss_snp$CHR))
  minor_allele <- bim[bim$V2 == ss_snp$ID]$V5
  effect <- ss_snp$EFFECT
} else {
  minor_allele <- NA
  effect <- NA
}

score_file_entry <- paste(best_snp_no_trail, minor_allele, effect, snp)
write(score_file_entry, file = SCORE_FILENAME, append = TRUE)
