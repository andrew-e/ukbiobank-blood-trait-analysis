library(lassosum)
library(data.table)
library(parallel)
setwd("/rds/general/user/are20/home/ukbiobank-blood-trait-analysis")

args = commandArgs(trailingOnly=TRUE)
marker <- args[1]
ethnicity <- args[2]
threads <- as.numeric(args[3])

print(sprintf("calculating lassosum validate for %s %s...", ethnicity, marker))
start_time <- Sys.time()

#patient_id_map <- fread("data/fam_filtered_patient_id_map.csv", sep=',', header=TRUE, stringsAsFactors=FALSE)
original_linked_ids <- fread("data/link_69328_19266_genotyped.txt")

#prepare PC
covar <- data.frame(FID = original_linked_ids$eid_19266, IID = original_linked_ids$eid_19266, eid_69328 = original_linked_ids$eid_69328)
principal_components <- fread("data/plink/UKB_app_69328_genetic_PCs.txt") #PCs are from the old data set, merge by secondary id
covar <- merge(covar, principal_components, by.x="eid_69328", by.y="eid", all.x = TRUE, sort = FALSE)
#pc_without_residuals <- cbind(covar, phenotype = 1)
#pc_without_residuals <- within(pc_without_residuals, phenotype[ethnic_background != 'white'] <- NA)
pc_without_residuals <- subset(covar, select = -c(eid_69328) )


#prepare phenotype
residuals <- fread(sprintf("results/residuals/%s_%s_residuals.csv", ethnicity, marker))
merged_phenotypes <- merge(original_linked_ids, residuals, by.x="eid_19266", by.y="eid_19266", all.x=TRUE, sort=FALSE)
#merged_phenotypes <- within(merged_phenotypes, residuals[ethnic_background != ethnicity] <- NA) #this should be true by default?
target_pheno <- data.frame(FID = merged_phenotypes$eid_19266, IID = merged_phenotypes$eid_19266, residuals = merged_phenotypes$residuals)

merged_pipeline <- readRDS(sprintf("results/lassosum/pipeline/merged_%s_%s.rds", ethnicity, marker))

EPHEMERAL_BASE <- "/rds/general/user/are20/ephemeral"
blood_markers = c("baso", "eo",  "hct", "hgb", "lymph", "mch", "mchc", "mcv", "mono", "mpv", "neut", "pct", "rbc", "rdw_cv", "wbc")
ethnicities = c("black", "asian", "chinese", "white", "mixed")

#for (marker in blood_markers) {
#  for (ethnicity in ethnicities) {
#    for (chromosome in 1:22) {
#      plink_chromosome_prefix = paste0("/rds/general/apps/RDS_COMPATIBILITY_LINKFARM/groupvol/med-bio/uk-biobank-2017/release_12032018/converted_data/imp_bgen1.1_plink/ukb_imp_chr", chromosome)
#      linked_file <- sprintf("%s/bfile_symlinks/linked_%s_%s_%s", EPHEMERAL_BASE, ethnicity, marker, chromosome)
#      file.symlink(paste0(plink_chromosome_prefix, "_v3.bim"), paste0(linked_file, ".bim"))
#      file.symlink(paste0(plink_chromosome_prefix, "_v3.bed"), paste0(linked_file, ".bed"))
#      file.symlink("/rds/general/user/are20/home/ukbiobank-blood-trait-analysis/data/working_fam_v6.fam", paste0(linked_file, ".fam"))
#    }
#  }
#}

#for (chromosome in 1:22) {
#  plink_chromosome_prefix = paste0("/rds/general/apps/RDS_COMPATIBILITY_LINKFARM/groupvol/med-bio/uk-biobank-2017/release_12032018/converted_data/imp_bgen1.1_plink/ukb_imp_chr", chromosome)
#  linked_file <- sprintf("%s/bfile_symlinks/linked_%s_%s_%s", EPHEMERAL_BASE, ethnicity, marker, chromosome)
#  file.symlink(paste0(plink_chromosome_prefix, "_v3.bim"), paste0(linked_file, ".bim"))
#  file.symlink(paste0(plink_chromosome_prefix, "_v3.bed"), paste0(linked_file, ".bed"))
#  file.symlink("/rds/general/user/are20/home/ukbiobank-blood-trait-analysis/data/working_fam_v6.fam", paste0(linked_file, ".fam"))
#}

cl <- makeCluster(threads)
lassosum_validation <- validate(merged_pipeline,
                                 pheno = as.data.frame(target_pheno),
                                 covar = as.data.frame(pc_without_residuals),
                                 cluster=cl
)

print("finished lassosum::validate...")
saveRDS(lassosum_validation, sprintf("results/lassosum/validation/first_pass_covar_pheno_%s_%s.rds", ethnicity, marker))

subset_lassosum_validation <- subset(merged_pipeline, s=lassosum_validation$best.s, lambda=lassosum_validation$best.lambda)
saveRDS(subset_lassosum_validation, sprintf("results/lassosum/validation/subset_covar_pheno_%s_%s.rds", ethnicity, marker))

print(paste("completed validation in", Sys.time() - start_time))

#for (chromosome in 1:22) {
#  file.remove(paste0(linked_file, ".bim"))
#  file.remove(paste0(linked_file, ".bed"))
#  file.remove(paste0(linked_file, ".fam"))
#}
