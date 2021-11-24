library(data.table)
library(jsonlite)
library(mltools)

phenotype_adjustments_from_json <- function() {
  biobank_csv = "/rds/general/project/chadeau_ukbb_folder/live/data/project_data/UKB_69328/basket_47946_14Sep21/ukb47946.csv"
  covars_metadata <- fromJSON("blood_type_covars.json", flatten=TRUE)
  
  blood_trait_covar_ids <- c()
  for (id in names(covars_metadata)) {
    blood_trait_covar_ids <- append(blood_trait_covar_ids, id)
  }

  blood_trait_covars <- fread(biobank_csv, sep=',', header=TRUE,stringsAsFactors=FALSE, select = blood_trait_covar_ids)

  header_names = c()
  for(column in colnames(blood_trait_covars)) {
    column_values <- covars_metadata[[column]]$values
    header_names <<- append(header_names, covars_test[[column]]$name)
    
    if(!is.null(column_values)) {
      #switch the ids with the names of the variables
      mapply(function(x, i) {
        blood_trait_covars[[column]][blood_trait_covars[[column]] == i] <<- x
      }, column_values, names(column_values)
      )

      #after swapping vars with names, turns column into a factor
      blood_trait_covars[[column]] <- as.factor(blood_trait_covars[[column]])
    }
  }

  colnames(blood_trait_covars) <- header_names

  #use the mltools.one_hot encoding function
  one_hot_blood_trait_covars <- one_hot(blood_trait_covars)
  head(one_hot_blood_trait_covars)
  
  #do scaling here?
  #Imputation?  Missing data: "For both datasets, where data-points were missing for a covariate, we imputed them by the mean covariate value and included a dummy variable to allow the mean of the index value for individuals with missing data to differ from the mean index value for individuals with non-missing data."
  
  #return(one_hot_blood_trait_covars)
}

