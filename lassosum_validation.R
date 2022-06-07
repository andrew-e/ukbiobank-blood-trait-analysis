library(lassosum)
library(data.table)
library(parallel)

start_time <- Sys.time()
#ethnicities = c("black", "asian", "chinese", "white", "mixed")
#blood_markers = c("baso", "eo",  "hct", "hgb", "lymph", "mch", "mchc", "mcv", "mono", "mpv", "neut", "pct", "rbc", "rdw_cv", "wbc")
args = commandArgs(trailingOnly=TRUE)
marker <- args[1]
ethnicity <- args[2]

print(sprintf("calculating lassosum validate for %s %s...", ethnicity, marker))

#setting bfiles up if they are not there already
for (chromosome in 1:22) {
  plink_file_prefix = "/rds/general/apps/RDS_COMPATIBILITY_LINKFARM/groupvol/med-bio/uk-biobank-2017/release_12032018/converted_data/imp_bgen1.1_plink/ukb_imp_chr"
  plink_chromosome_prefix = paste0(plink_file_prefix, chromosome)
  linked_file <- sprintf("linked_%s_%s_%s", ethnicity, marker, chromosome)
  
  file.remove(paste0(linked_file, ".bim"))
  file.symlink(paste0(plink_chromosome_prefix, "_v3.bim"), paste0(linked_file, ".bim"))
  file.remove(paste0(linked_file, ".bed"))
  file.symlink(paste0(plink_chromosome_prefix, "_v3.bed"), paste0(linked_file, ".bed"))
  file.remove(paste0(linked_file, ".fam"))
  file.symlink("/rds/general/user/are20/home/ukbiobank-blood-trait-analysis/data/working_fam_v6.fam", paste0(linked_file, ".fam"))
}

biobank_id_map = "/rds/general/user/are20/home/ukbiobank-blood-trait-analysis/data/fam_filtered_patient_id_map.csv"
patient_id_map <- fread(biobank_id_map, sep=',', header=TRUE,stringsAsFactors=FALSE)

#prepare PC
covar <- data.frame(FID = patient_id_map$id, IID = patient_id_map$id, ethnic_background = patient_id_map$ethnic_background, eid_69328 = patient_id_map$eid_69328)

principal_components <- fread("/rds/general/user/are20/home/plink/UKB_app_69328_genetic_PCs_altered_colnames.txt")
covar <- merge(covar, principal_components, by.x="eid_69328", by.y="eid", all.x = TRUE, sort = FALSE)
covar <- subset(covar, select = -c(ethnic_background, eid_69328) )

phenotypes <- fread(sprintf("/rds/general/user/are20/home/ukbiobank-blood-trait-analysis/results/residuals/%s_%s_residuals.csv", ethnicity, marker))
merged_phenotypes <- merge(patient_id_map, phenotypes, by.x="id", by.y="eid_19266", all.x=TRUE, sort=FALSE)
target_pheno <- data.frame(FID = merged_phenotypes$id, IID = merged_phenotypes$id, residuals = merged_phenotypes$residuals)

pc_with_residuals <- merge(covar, phenotypes, by.x="FID", by.y="eid_19266", all.x=TRUE, sort=FALSE)
pc_with_residuals <- subset(pc_with_residuals, select = -c(id) )

current_pipeline <- readRDS(sprintf("/rds/general/user/are20/home/ukbiobank-blood-trait-analysis/results/lassosum/pipeline/merged_%s_%s_v2.rds", ethnicity, marker))

cl <- makeCluster(8)
target.res <- lassosum::validate(current_pipeline,
                                 covar = as.data.frame(pc_with_residuals),
                                 cluster=cl
)

print("finished lassosum::validate...")
saveRDS(target.res, sprintf("/rds/general/user/are20/home/ukbiobank-blood-trait-analysis/results/lassosum/validation/covar_only_%s_%s.rds", ethnicity, marker))

r2 <- max(target.res$validation.table$value)^2
write(r2, file = sprintf("/rds/general/user/are20/home/ukbiobank-blood-trait-analysis/results/lassosum/r2_covar_%s_%s.txt", ethnicity, marker))

print(paste("completed validation in", Sys.time() - start_time))
