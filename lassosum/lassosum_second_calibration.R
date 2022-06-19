library(lassosum)
library(data.table)
library(parallel)
setwd("/rds/general/user/are20/home/ukbiobank-blood-trait-analysis")

args = commandArgs(trailingOnly=TRUE)
marker <- args[1]
ethnicity <- args[2]

print(sprintf("second validation and calibration of %s %s...", ethnicity, marker))
start_time <- Sys.time()

#patient_id_map <- fread("data/fam_filtered_patient_id_map.csv", sep=',', header=TRUE,stringsAsFactors=FALSE)
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

current_pipeline <- readRDS(sprintf("/rds/general/user/are20/home/ukbiobank-blood-trait-analysis/results/lassosum/validation/subset_covar_pheno_%s_%s.rds", ethnicity, marker))

cl <- makeCluster(32)
second_validation <- validate(current_pipeline,
                                 pheno = as.data.frame(target_pheno),
                                 covar = as.data.frame(pc_without_residuals),
                                 cluster=cl
)


saveRDS(second_validation, sprintf("/rds/general/user/are20/home/ukbiobank-blood-trait-analysis/results/lassosum/validation/second_pass_covar_pheno_%s_%s.rds", ethnicity, marker))



#TODO: here, you can take a look at a few things.  There is the result$validation.table
#   including the adjusted r^2 here: r2 <- max(result$validation.table$value)^2
#There is also the result$results.table, which has a pgs per FID, maybe this could be used to compare with the residuals via pearson correlation, in the same way we did with plink.
#There is also result$pgs, which is different than the data in the results.table ?
#> result$best.validation.result ?
# The other question is comparing like for like.  Should I be summing up the PGS per chromosome like before?  Could there be a difference? 


second_validation <- readRDS(sprintf("/rds/general/user/are20/home/ukbiobank-blood-trait-analysis/results/lassosum/validation/second_pass_covar_pheno_%s_%s.rds", ethnicity, marker))


#now, run an lm on the output of the second validation and the pc with residuals attached.
#pc_with_residuals <- merge(covar, residuals, by.x="FID", by.y="eid_19266", all.x=TRUE, sort=FALSE)
#pc_with_residuals <- subset(pc_with_residuals, select = -c(id) )

#pgs = second_validation$results.table$best.pgs
#names(pgs) = second_validation$results.table$IID

#compare_residuals_and_pgs <- merge(merged_phenotypes, second_validation$results.table, by.x="eid_19266", by.y="FID", sort=FALSE)
#complete_cases <- compare_residuals_and_pgs[complete.cases(compare_residuals_and_pgs[,c("best.pgs","residuals")]),]

#correlation_result <- cor(complete_cases$best.pgs, complete_cases$residuals, use = "pairwise.complete.obs")

#rownames(pc_with_residuals)=pc_with_residuals[,2]
#pc_with_residuals=pc_with_residuals[,-c(1,2)]

#pc_with_residuals=pc_with_residuals[names(pgs),]

#formula=paste0("pgs~", paste0(paste0("pc_with_residuals$", colnames(pc_with_residuals)), collapse="+"))
#print("Formula to adjust on technical covariates:")
#print(formula)

#lm_run <- lm(as.formula(formula))

#print("Saving adjusted UK Biobank predictions...")
#saveRDS(lm_run, sprintf("/rds/general/user/are20/home/ukbiobank-blood-trait-analysis/results/lassosum/validation/lm_adjusted_pgs_%s_%s.rds", ethnicity, marker))

#print(paste("completed validation in", Sys.time() - start_time))
