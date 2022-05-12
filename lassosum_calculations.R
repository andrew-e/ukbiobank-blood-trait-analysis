library(devtools)
install_github("tshmak/lassosum")
library(lassosum)
library(data.table)
library(jsonlite)
library(parallel)

setwd("/rds/general/user/are20/home/")
#linking old and new biobank ids...
plink_file_prefix = "/rds/general/apps/RDS_COMPATIBILITY_LINKFARM/groupvol/med-bio/uk-biobank-2017/release_12032018/converted_data/imp_bgen1.1_plink/ukb_imp_chr"
biobank_id_map = "/rds/general/user/are20/home/ukbiobank-blood-trait-analysis/data/ids_and_ethnicity.csv"
patient_id_map <- fread(biobank_id_map, sep=',', header=TRUE,stringsAsFactors=FALSE)

#prepare PC
principal_components <- fread("/rds/general/user/are20/home/plink/UKB_app_69328_genetic_PCs.txt")
principal_components <- merge(principal_components, patient_id_map, by.x="eid", by.y="eid_69328")
principal_components <- cbind(phenotype = 1, principal_components) #is this right?
principal_components$IID <- principal_components$eid_19266
principal_components$FID <- principal_components$eid_19266

ethnicities = c("black", "asian", "chinese", "white", "mixed")
#blood_markers = c("baso", "eo",  "hct", "hgb", "lymph", "mch", "mchc", "mcv", "mono", "mpv", "neut", "pct", "rbc", "rdw_cv", "wbc")
blood_markers = c("pct")

#for (marker in blood_markers) {
marker = "pct"
ethnicity = "black"

summary_statistics <- fread(sprintf("/rds/general/user/are20/home/plink/summary_statistics/%s.assoc", marker))
summary_statistics <- summary_statistics[!P == 0]

pc_ethnicity <- principal_components[ethnic_background == ethnicity]
cor <- p2cor(p = as.numeric(summary_statistics$P),
             n = nrow(pc_ethnicity),
             sign = summary_statistics$EFFECT
)

summary_statistics$COR <- cor
summary_statistics <- summary_statistics[!is.na(COR)]

ld_blocks <- "AFR.hg19" # Change depending on ethnicity variable

lassosum_results <- list()
last_result <- list()
for (chromosome in 21:22) {
  if (exists("cl")) {
    rm(cl)
  }
  cl <- makeCluster(8)
  
  plink_chromosome_prefix = paste0(plink_file_prefix, chromosome)
  file.remove("linked.bim")
  file.symlink(paste0(plink_chromosome_prefix, "_v3.bim"), "linked.bim")
  file.remove("linked.bed")
  file.symlink(paste0(plink_chromosome_prefix, "_v3.bed"), "linked.bed")
  
  set.seed(1)
  current_lassosum_result <- lassosum.pipeline(
    cor = summary_statistics$COR,
    chr = summary_statistics$CHR,
    pos = summary_statistics$BP,
    A1 = summary_statistics$REF,
    A2 = summary_statistics$ALT,
    ref.bfile = "linked",
    sample = 5000,
    test.bfile = "linked",
    LDblocks = ld_blocks,
    cluster=cl
  )

  if (length(lassosum_results) == 0) {
    lassosum_results <- current_lassosum_result
  } else {
    lassosum_results <- merge(lassosum_results, current_lassosum_result)
    last_result <- current_lassosum_result
  }
}


#}


phenotypes <- fread(sprintf("/rds/general/user/are20/home/ukbiobank-blood-trait-analysis/results/residuals/%s_%s_residuals.csv", ethnicity, marker)) #is this the residuals?
target_pheno <- data.frame(FID = phenotypes$eid_19266, IID = phenotypes$eid_19266, residuals = phenotypes$residuals)
target.res <- validate(lassosum_results, pheno = as.data.frame(target_pheno))
r2 <- max(target.res$validation.table$value)^2
write(r2, file = sprintf("/rds/general/user/are20/home/ukbiobank-blood-trait-analysis/results/lassosum/%s_%s_r2.txt", ethnicity, marker))

print("just for my own sanity...")
target.res <- validate(last_result, pheno = as.data.frame(target_pheno))
single_r2 <- max(target.res$validation.table$value)^2
print(single_r2)

