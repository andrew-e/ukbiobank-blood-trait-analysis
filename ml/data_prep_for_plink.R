library(data.table)

setwd("/rds/general/user/are20/home/ukbiobank-blood-trait-analysis")
args = commandArgs(trailingOnly=TRUE)
marker <- args[1]
ethnicity <- args[2]

SCORE_FILENAME <- sprintf("results/ml/scores/%s_%s.score", marker, ethnicity)
scores <- fread(SCORE_FILENAME)
scores$V4 <- NULL
scores <- unique(scores)
scores$V3 <- abs(scores$V3) #Just copying how the original score files were
scores <- scores[complete.cases(scores),]

write.table(scores, sprintf("results/ml/scores/%s_%s.score_cleaned", marker, ethnicity), quote = FALSE, row.names = FALSE, col.names = FALSE)

########## RANDOM SHIT  #########


#all_bims <- fread(sprintf("/rds/general/apps/RDS_COMPATIBILITY_LINKFARM/groupvol/med-bio/uk-biobank-2017/release_12032018/converted_data/imp_bgen1.1_plink/ukb_imp_chr1_v3.bim"))
#for (chr in 2:22) {
#  bim <- fread(sprintf("/rds/general/apps/RDS_COMPATIBILITY_LINKFARM/groupvol/med-bio/uk-biobank-2017/release_12032018/converted_data/imp_bgen1.1_plink/ukb_imp_chr%s_v3.bim", chr))
#  all_bims <- rbind(all_bims, bim)
#}



#chr12_bim <- fread("/rds/general/apps/RDS_COMPATIBILITY_LINKFARM/groupvol/med-bio/uk-biobank-2017/release_12032018/converted_data/imp_bgen1.1_plink/ukb_imp_chr12_v3.bim")


#head(chr12_bim)

#score_bim <- merge(scores, chr12_bim, by.x="V1", by.y="V2")

#scores_with_na <- fread("/rds/general/user/are20/home/ukbiobank-blood-trait-analysis/results/ml/scores/pct_black_with_NAs.score")
#scores_with_na_and_bim<- merge(scores_with_na, all_bims, by.x="V2", by.y="V2")

#original_score_file <- fread(sprintf("/rds/general/user/are20/home/plink/scores/reduced_%s_my.score_condind_common", marker))



