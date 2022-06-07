library(lassosum)
library(data.table)
library(parallel)

#ethnicities = c("black", "asian", "chinese", "white", "mixed")
#blood_markers = c("baso", "eo",  "hct", "hgb", "lymph", "mch", "mchc", "mcv", "mono", "mpv", "neut", "pct", "rbc", "rdw_cv", "wbc")

args = commandArgs(trailingOnly=TRUE)
chromosome <- args[1]
marker <- args[2]
ethnicity <- args[3]

ld_blocks <- NULL
if (ethnicity == "white") {
  ld_blocks <- "EUR.hg19"
} else if (ethnicity == "black") {
  ld_blocks <- "AFR.hg19"
} else if (ethnicity == "chinese" || ethnicity == "asian") {
  ld_blocks <- "ASN.hg19"
} else {
  stop(paste("ethnicity", ethnicity, "is not supported"))
}

print(sprintf("running pipeline for chr_%s_%s_%s", chromosome, ethnicity, marker))

#if (file.exists(sprintf("/rds/general/user/are20/home/ukbiobank-blood-trait-analysis/results/lassosum/pipeline/new_chr_%s_%s_%s_v6.rds", chromosome, ethnicity, marker))) {
#  print(sprintf("new_chr_%s_%s_%s.rds already calculated, skipping", chromosome, ethnicity, marker))
#  q()
#}

start_time <- Sys.time()

setwd("/rds/general/user/are20/home/")
#linking old and new biobank ids...
plink_file_prefix = "/rds/general/apps/RDS_COMPATIBILITY_LINKFARM/groupvol/med-bio/uk-biobank-2017/release_12032018/converted_data/imp_bgen1.1_plink/ukb_imp_chr"
biobank_id_map = "/rds/general/user/are20/home/ukbiobank-blood-trait-analysis/data/fam_filtered_patient_id_map.csv"
patient_id_map <- fread(biobank_id_map, sep=',', header=TRUE,stringsAsFactors=FALSE)

ethnicity_list <- patient_id_map[patient_id_map$ethnic_background == "white"]

not_ethnicity_list <- patient_id_map[patient_id_map$ethnic_background != "white"]
remove_ref <- data.frame(FID=not_ethnicity_list$id, IID = not_ethnicity_list$id)

ethnicity_keep_ref <- head(ethnicity_list, n = 6000)
ethnicity_keep_ref <- data.frame(FID=ethnicity_keep_ref$id, IID = ethnicity_keep_ref$id)

ethnicity_keep_test <- tail(ethnicity_list, 6000)
ethnicity_keep_test <- data.frame(FID=ethnicity_keep_test$id, IID = ethnicity_keep_test$id)

#prepare PC
principal_components <- fread("/rds/general/user/are20/home/plink/UKB_app_69328_genetic_PCs_altered_colnames.txt")
principal_components <- merge(patient_id_map, principal_components, by.x="eid_69328", by.y="eid", all.x = TRUE, sort = FALSE)

summary_statistics <- fread(sprintf("/rds/general/user/are20/home/plink/summary_statistics/%s.assoc", marker))
summary_statistics <- summary_statistics[!P == 0]

cor <- p2cor(p = as.numeric(summary_statistics$P),
             n = nrow(principal_components),
             sign = summary_statistics$EFFECT
)

summary_statistics$COR <- cor
summary_statistics <- summary_statistics[!is.na(COR)]

plink_chromosome_prefix = paste0(plink_file_prefix, chromosome)

linked_file <- sprintf("linked_%s_%s_%s", ethnicity, marker, chromosome)

file.remove(paste0(linked_file, ".bim"))
file.symlink(paste0(plink_chromosome_prefix, "_v3.bim"), paste0(linked_file, ".bim"))
file.remove(paste0(linked_file, ".bed"))
file.symlink(paste0(plink_chromosome_prefix, "_v3.bed"), paste0(linked_file, ".bed"))
file.remove(paste0(linked_file, ".fam"))
file.symlink("/rds/general/user/are20/home/ukbiobank-blood-trait-analysis/data/working_fam_v6.fam", paste0(linked_file, ".fam"))

cl <- makeCluster(32)
current_lassosum_result <- lassosum.pipeline(
  cor = summary_statistics$COR,
  chr = summary_statistics$CHR,
  pos = summary_statistics$BP,
  A1 = summary_statistics$REF,
  A2 = summary_statistics$ALT,
  ref.bfile = linked_file,
  keep.ref = ethnicity_keep_ref,
  keep.test = ethnicity_keep_test,
  remove.ref = remove_ref,
  test.bfile = linked_file,
  LDblocks = ld_blocks,
  exclude.ambiguous=FALSE,
  max.ref.bfile.n=500000,
  cluster=cl
)

saveRDS(current_lassosum_result, sprintf("/rds/general/user/are20/home/ukbiobank-blood-trait-analysis/results/lassosum/pipeline/pipeline_chr_%s_%s_%s.rds", ethnicity, marker, chromosome))
print(paste("completed validation in", Sys.time() - start_time))