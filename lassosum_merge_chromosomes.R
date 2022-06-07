library(lassosum)
library(data.table)

#ethnicities = c("black", "asian", "chinese", "white", "mixed")
#blood_markers = c("baso", "eo",  "hct", "hgb", "lymph", "mch", "mchc", "mcv", "mono", "mpv", "neut", "pct", "rbc", "rdw_cv", "wbc")
args = commandArgs(trailingOnly=TRUE)
marker <- args[1]
ethnicity <- args[2]

print(sprintf("merging lassosum pipelines for %s %s...", ethnicity, marker))

lassosum_pipeline <- readRDS(sprintf("/rds/general/user/are20/home/ukbiobank-blood-trait-analysis/results/lassosum/pipeline/pipeline_chr_%s_%s_%s.rds", ethnicity, marker, chromosome))
for (chromosome in 2:22) {
  current_pipeline <- readRDS(sprintf("/rds/general/user/are20/home/ukbiobank-blood-trait-analysis/results/lassosum/pipeline/pipeline_chr_%s_%s_%s_v6.rds", ethnicity, marker, chromosome))

  plink_file_prefix = "/rds/general/apps/RDS_COMPATIBILITY_LINKFARM/groupvol/med-bio/uk-biobank-2017/release_12032018/converted_data/imp_bgen1.1_plink/ukb_imp_chr"
  plink_chromosome_prefix = paste0(plink_file_prefix, chromosome)
  linked_file <- sprintf("linked_%s_%s_%s", ethnicity, marker, chromosome)

  file.remove(paste0(linked_file, ".bim"))
  file.symlink(paste0(plink_chromosome_prefix, "_v3.bim"), paste0(linked_file, ".bim"))
  file.remove(paste0(linked_file, ".bed"))
  file.symlink(paste0(plink_chromosome_prefix, "_v3.bed"), paste0(linked_file, ".bed"))
  file.remove(paste0(linked_file, ".fam"))
  file.symlink("/rds/general/user/are20/home/ukbiobank-blood-trait-analysis/data/working_fam_v6.fam", paste0(linked_file, ".fam"))

  lassosum_pipeline <- merge(lassosum_pipeline, current_pipeline)
}

saveRDS(lassosum_pipeline, sprintf("/rds/general/user/are20/home/ukbiobank-blood-trait-analysis/results/lassosum/pipeline/merged_%s_%s.rds", ethnicity, marker))
print("finished lassosum:::merge")

