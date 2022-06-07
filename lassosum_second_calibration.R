library(lassosum)
library(data.table)
library(parallel)

#ethnicities = c("black", "asian", "chinese", "white", "mixed")
#blood_markers = c("baso", "eo",  "hct", "hgb", "lymph", "mch", "mchc", "mcv", "mono", "mpv", "neut", "pct", "rbc", "rdw_cv", "wbc")
args = commandArgs(trailingOnly=TRUE)
marker <- args[1]
ethnicity <- args[2]

print(sprintf("second validation and calibration of %s %s...", ethnicity, marker))
start_time <- Sys.time()

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
pc_without_residuals <- cbind(covar, phenotype = 1)

#prepare phenotype
phenotypes <- fread(sprintf("/rds/general/user/are20/home/ukbiobank-blood-trait-analysis/results/residuals/%s_%s_residuals.csv", ethnicity, marker))
merged_phenotypes <- merge(patient_id_map, phenotypes, by.x="id", by.y="eid_19266", all.x=TRUE, sort=FALSE)
target_pheno <- data.frame(FID = merged_phenotypes$id, IID = merged_phenotypes$id, residuals = merged_phenotypes$residuals)

pc_with_residuals <- merge(covar, phenotypes, by.x="FID", by.y="eid_19266", all.x=TRUE, sort=FALSE)
pc_with_residuals <- subset(pc_with_residuals, select = -c(id) )


current_pipeline <- readRDS(sprintf("/rds/general/user/are20/home/ukbiobank-blood-trait-analysis/results/lassosum/pipeline/subset_covar_pheno_%s_%s.rds", ethnicity, marker))

cl <- makeCluster(8)
second_validation <- validate(current_pipeline,
                                 pheno = as.data.frame(target_pheno),
                                 covar = as.data.frame(pc_without_residuals),
                                 cluster=cl
)


saveRDS(second_validation, sprintf("/rds/general/user/are20/home/ukbiobank-blood-trait-analysis/results/lassosum/validation/second_pass_covar_pheno_%s_%s.rds", ethnicity, marker))

#now, run an lm on the output of the second validation and the pc with residuals attached.
pgs = second_validation$results.table$best.pgs
names(pgs) = second_validation$results.table$IID

rownames(pc_with_residuals)=pc_with_residuals[,2]
pc_with_residuals=pc_with_residuals[,-c(1,2)]

pc_with_residuals=pc_with_residuals[names(pgs),]

formula=paste0("pgs~", paste0(paste0("pc_with_residuals$", colnames(pc_with_residuals)), collapse="+"))
print("Formula to adjust on technical covariates:")
print(formula)

adjusted_pgs = residuals(lm(as.formula(formula)))
print(head(adjusted_pgs))

print("Saving adjusted UK Biobank predictions...")
saveRDS(adjusted_pgs, sprintf("/rds/general/user/are20/home/ukbiobank-blood-trait-analysis/results/lassosum/validation/lm_adjusted_pgs_%s_%s.rds", ethnicity, marker))

print(paste("completed validation in", Sys.time() - start_time))
