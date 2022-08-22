library(lassosum)
library(data.table)
setwd("/rds/general/user/are20/home/ukbiobank-blood-trait-analysis")

args = commandArgs(trailingOnly=TRUE)
marker <- args[1]
ethnicity <- args[2]

print(sprintf("merging lassosum pipelines for %s %s...", ethnicity, marker))

lassosum_pipelines <- list()
for (chromosome in 1:22) {
  current_pipeline <- readRDS(sprintf("results/lassosum/pipeline/pipeline_chr_%s_%s_%s.rds", ethnicity, marker, chromosome))
  lassosum_pipelines <- append(lassosum_pipelines, list(current_pipeline))
}


merged_pipeline <- do.call(lassosum:::merge.lassosum.pipeline, lassosum_pipelines)

saveRDS(merged_pipeline, sprintf("results/lassosum/pipeline/merged_%s_%s.rds", ethnicity, marker))
print("finished lassosum:::merge")

for (chromosome in 1:22) {
  file.remove(sprintf("results/lassosum/pipeline/pipeline_chr_%s_%s_%s.rds", ethnicity, marker, chromosome))
}
