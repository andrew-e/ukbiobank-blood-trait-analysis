library(lassosum)
library(data.table)
library(parallel)
setwd("/rds/general/user/are20/home/ukbiobank-blood-trait-analysis")
start_time <- Sys.time()

EPHEMERAL_BASE <- "/rds/general/user/are20/ephemeral"

args = commandArgs(trailingOnly=TRUE)
chromosome <- args[1]
marker <- args[2]
ethnicity <- args[3]

ld_blocks <- NULL
if (ethnicity == "white" || ethnicity == "mixed") {
  ld_blocks <- "EUR.hg19"
} else if (ethnicity == "black") {
  ld_blocks <- "AFR.hg19"
} else if (ethnicity == "asian") {
  ld_blocks <- "ASN.hg19"
} else {
  stop(paste("ethnicity", ethnicity, "is not supported"))
}

print(sprintf("running pipeline for %s_%s_%s", ethnicity, marker, chromosome))

if (file.exists(sprintf("results/lassosum/pipeline/pipeline_chr_%s_%s_%s.rds", ethnicity, marker, chromosome))) {
  print("pipeline already exists, exiting.")
  q()
}

if (file.exists(sprintf("results/lassosum/pipeline/merged_%s_%s.rds", ethnicity, marker))) {
  print("pipelines already merged, exiting.")
  q()
}

#linking old and new biobank ids...
patient_id_map <- fread("data/fam_filtered_patient_id_map.csv", sep=',', header=TRUE,stringsAsFactors=FALSE)

ethnicity_list <- patient_id_map[patient_id_map$ethnic_background == ethnicity]
not_ethnicity_list <- patient_id_map[patient_id_map$ethnic_background != ethnicity]

remove_ref <- data.frame(FID=not_ethnicity_list$eid_19266, IID = not_ethnicity_list$eid_19266)

test_size <- floor(nrow(ethnicity_list) / 10)

ethnicity_keep_ref <- head(ethnicity_list, n = test_size)
ethnicity_keep_ref <- data.frame(FID=ethnicity_keep_ref$eid_19266, IID = ethnicity_keep_ref$eid_19266)

ethnicity_keep_test <- tail(ethnicity_list, n = test_size)
ethnicity_keep_test <- data.frame(FID=ethnicity_keep_test$eid_19266, IID = ethnicity_keep_test$eid_19266)

#prepare PC
principal_components <- fread("data/plink/UKB_app_69328_genetic_PCs.txt")
principal_components <- merge(patient_id_map, principal_components, by.x="eid_69328", by.y="eid", all.x = TRUE, sort = FALSE)

summary_statistics <- fread(sprintf("data/plink/summary_statistics/%s.assoc", marker))
summary_statistics <- summary_statistics[!P == 0]

cor <- p2cor(p = as.numeric(summary_statistics$P),
             n = nrow(principal_components),
             sign = summary_statistics$EFFECT
)

summary_statistics$COR <- cor
summary_statistics <- summary_statistics[!is.na(COR)]


linked_file <- sprintf("%s/bfile_symlinks/linked_%s_%s_%s", EPHEMERAL_BASE, ethnicity, marker, chromosome)

cl <- makeCluster(32)
current_lassosum_result <- lassosum.pipeline(
  cor = summary_statistics$COR,
  chr = summary_statistics$CHR,
  pos = summary_statistics$BP,
  A2 = summary_statistics$REF,
  A1 = summary_statistics$ALT,
  ref.bfile = linked_file,
  test.bfile = linked_file,
  keep.ref = ethnicity_keep_ref,
  keep.test = ethnicity_keep_test,
  remove.ref = remove_ref,
  LDblocks = ld_blocks,
  exclude.ambiguous=FALSE,
  max.ref.bfile.n=500000,
  cluster=cl
)

saveRDS(current_lassosum_result, sprintf("results/lassosum/pipeline/pipeline_chr_%s_%s_%s.rds", ethnicity, marker, chromosome))
print(paste("completed validation in", Sys.time() - start_time))
