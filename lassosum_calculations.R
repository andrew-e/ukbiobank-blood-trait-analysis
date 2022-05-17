library(devtools)
library(lassosum)
library(data.table)
library(jsonlite)
library(parallel)

#ethnicities = c("black", "asian", "chinese", "white", "mixed")
#blood_markers = c("baso", "eo",  "hct", "hgb", "lymph", "mch", "mchc", "mcv", "mono", "mpv", "neut", "pct", "rbc", "rdw_cv", "wbc")

args = commandArgs(trailingOnly=TRUE)
chromosome <- args[1]
marker <- args[2]

ethnicity = "black"
ld_blocks <- "AFR.hg19" # Change depending on ethnicity variable

print(sprintf("calculating score for r2_%s_%s_%s.txt...", ethnicity, marker, chromosome))
if (file.exists(sprintf("/rds/general/user/are20/home/ukbiobank-blood-trait-analysis/results/lassosum/r2_%s_%s_%s.txt", ethnicity, marker, chromosome))) {
  print(sprintf("r2_%s_%s_%s.txt already calculated, skipping", ethnicity, marker, chromosome))
  q()
}

start_time <- Sys.time()

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


summary_statistics <- fread(sprintf("/rds/general/user/are20/home/plink/summary_statistics/%s.assoc", marker))
summary_statistics <- summary_statistics[!P == 0]

pc_ethnicity <- principal_components[ethnic_background == ethnicity]
cor <- p2cor(p = as.numeric(summary_statistics$P),
             n = nrow(pc_ethnicity),
             sign = summary_statistics$EFFECT
)

summary_statistics$COR <- cor
summary_statistics <- summary_statistics[!is.na(COR)]

cl <- makeCluster(8)

plink_chromosome_prefix = paste0(plink_file_prefix, chromosome)
linked_file <- sprintf("linked_%s_%s", ethnicity, marker)
file.remove(paste0(linked_file, ".bim"))
file.symlink(paste0(plink_chromosome_prefix, "_v3.bim"), paste0(linked_file, ".bim"))
file.remove(paste0(linked_file, ".bed"))
file.symlink(paste0(plink_chromosome_prefix, "_v3.bed"), paste0(linked_file, ".bed"))

set.seed(1)
current_lassosum_result <- lassosum.pipeline(
  cor = summary_statistics$COR,
  chr = summary_statistics$CHR,
  pos = summary_statistics$BP,
  A1 = summary_statistics$REF,
  A2 = summary_statistics$ALT,
  ref.bfile = linked_file,
  sample = 5000,
  test.bfile = linked_file,
  LDblocks = ld_blocks,
  cluster=cl
)

print("now calculating the validation")

phenotypes <- fread(sprintf("/rds/general/user/are20/home/ukbiobank-blood-trait-analysis/results/residuals/%s_%s_residuals.csv", ethnicity, marker))
merged_phenotypes <- merge(patient_id_map, phenotypes, by.x="eid_19266", by.y="eid_19266", all.x=TRUE)
target_pheno <- data.frame(FID = merged_phenotypes$eid_19266, IID = merged_phenotypes$eid_19266, residuals = merged_phenotypes$residuals)
target.res <- lassosum::validate(current_lassosum_result, pheno = as.data.frame(target_pheno))
r2 <- max(target.res$validation.table$value)^2
write(r2, file = sprintf("/rds/general/user/are20/home/ukbiobank-blood-trait-analysis/results/lassosum/r2_%s_%s_%s.txt", ethnicity, marker, chromosome))

print(Sys.time() - start_time)
