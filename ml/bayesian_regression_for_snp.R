library(data.table)
require(dplyr)
library(rstanarm)
options(mc.cores = parallel::detectCores())
#You can read up on bayesian linear regression and stan_glm here:https://cran.r-project.org/web/packages/rstanarm/rstanarm.pdf#page.108

setwd("/rds/general/user/are20/home/ukbiobank-blood-trait-analysis")
args = commandArgs(trailingOnly=TRUE)
marker <- args[1]
ethnicity <- args[2]
snp <- args[3]

EPHEMERAL_DIR <- paste0("/rds/general/user/are20/ephemeral/ml/", marker, "/SNP_", snp)
residuals <- fread(sprintf("results/residuals/%s_%s_residuals.csv", ethnicity, marker), sep=',', header=TRUE, stringsAsFactors=FALSE)

#Get summary statistics, to be used for priors in the bayesian linear regression.
summary_stat_columns <- c("ID", "CHR", "EFFECT", "SE")
summary_statistic <- fread(sprintf("data/plink/summary_statistics/%s.assoc", marker), select = summary_stat_columns)
ss_snp <- summary_statistic[summary_statistic$ID == snp]

#get PCs and plink raw output
principcal_components <- fread("data/plink/UKB_app_69328_genetic_PCs.txt")

plink_data <- fread(paste0(EPHEMERAL_DIR, "/plink.raw"))
plink_data <- plink_data %>% select(unique(colnames(.))) #removing duplicate columns, since multiple rows of the same SNP are found in the .bim files
plink_data_ethnicity <- merge(residuals, plink_data, by.x="eid_19266", by.y="FID")

snp_list <- colnames(plink_data)[-6:-1]
for (surrounding_snp in snp_list) {
  plink_subset <- subset(plink_data_ethnicity, select = c("id", "eid_19266", "residuals", surrounding_snp))
  snp_with_pcs <- merge(plink_subset, principcal_components, by.x="id", by.y="eid")
  names(snp_with_pcs)[4] <- "snp_dosage"

  #I think I have the relevant data..
  #The coded SNP for all the relevant ethnic background, annotated with the residuals of the blood marker

  sample_size <- nrow(snp_with_pcs)
  scale <- ss_snp$SE * sqrt(447409)
  dosage_breakdown <- table(snp_with_pcs$snp_dosage) #maybe this needs to be the snp we start with, not calculated for each snp...
  location <- as.numeric((ss_snp$EFFECT * dosage_breakdown[2] / sample_size) + (ss_snp$EFFECT * 2 * dosage_breakdown[3] / sample_size))
  
  #run stan model with and without PCs I guess?
  stan_model_with_prior <- stan_glm(residuals ~ snp_dosage, data = snp_with_pcs, prior = normal(location = location, scale = scale))
  
  pc_names <- colnames(snp_with_pcs)[-4:-1]
  formula_with_pcs = paste0("residuals ~ snp_dosage + ", paste(pc_names, collapse=" + "))
  stan_model_with_prior_and_pcs <- stan_glm(as.formula(formula_with_pcs), data = snp_with_pcs, prior = normal(location = location, scale = scale))
  
  prior_summary(stan_model_with_prior)
  prior_summary(stan_model_with_prior_and_pcs)
  
  summary(stan_model_with_prior)
  summary(stan_model_with_prior_and_pcs)
  
  #do something with the results and save them to a table?!?
}