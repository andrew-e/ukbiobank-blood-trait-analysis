library(lassosum)
library(data.table)

setwd(system.file("data", package="lassosum")) # Directory where data and LD region files are stored

#linking old and new biobank ids...
biobank_id_map = "/rds/general/user/are20/home/ukbiobank-blood-trait-analysis/data/link_69328_19266_genotyped.txt"
patient_id_map <- fread(biobank_id_map, sep=' ', header=TRUE,stringsAsFactors=FALSE)


#### Get this data just for ethnicity map ###
biobank_csv = "/rds/general/project/chadeau_ukbb_folder/live/data/project_data/UKB_69328/basket_47946_14Sep21/ukb47946.csv"
covars_metadata <- fromJSON("ukbiobank-blood-trait-analysis/blood_type_covars.json", flatten=TRUE)
blood_trait_covar_ids <- c()
for (id in names(covars_metadata)) {
  blood_trait_covar_ids <- append(blood_trait_covar_ids, id)
}
blood_trait_covars <- fread(biobank_csv, sep=',', header=TRUE,stringsAsFactors=FALSE, select = blood_trait_covar_ids) #nrows = 100

covar_header_names = c()
for(column in colnames(blood_trait_covars)) {
  column_values <- covars_metadata[[column]]$values
  covar_header_names <<- append(covar_header_names, covars_metadata[[column]]$name)
  
  if(!is.null(column_values)) {
    #switch the ids with the names of the variables
    mapply(function(x, i) {
      blood_trait_covars[[column]][blood_trait_covars[[column]] == i] <<- x
    }, column_values, names(column_values)
    )
    
    #after swapping vars with names, turns column into a factor
    blood_trait_covars[[column]] <- as.factor(blood_trait_covars[[column]])
  }
  else {
    blood_trait_covars[[column]] <- as.numeric(blood_trait_covars[[column]]) #??
  }
}
#### End of this ethnicity map bullshit ###


### Read summary statistics file ###
ss <- fread("plink/summary_statistics/pct.assoc")
head(ss)

### Specify the PLINK file stub of the reference panel ###
ref.bfile <- "refpanel"

### Specify the PLINK file stub of the test data ###
test.bfile <- "testsample"

### Read LD region file ###
LDblocks <- "EUR.hg19" # This will use LD regions as defined in Berisa and Pickrell (2015) for the European population and the hg19 genome.
# Other alternatives available. Type ?lassosum.pipeline for more details. 

pcs <- fread("plink/UKB_app_69328_genetic_PCs.txt")
merged_pcs <- merge(pcs, patient_id_map, by.x="eid", by.y="eid_69328")

target.pheno <- fread("EUR.height")[,c("FID", "IID", "Height")]

#maybe don't have to do anything with covariages?
#covariate <- fread("EUR.cov")

