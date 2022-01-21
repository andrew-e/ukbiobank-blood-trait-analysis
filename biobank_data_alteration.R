library(data.table)
library(jsonlite)
library(mltools)
library(impute)
library(mgcv)

phenotype_adjustments_from_json <- function() {
  biobank_csv = "/rds/general/project/chadeau_ukbb_folder/live/data/project_data/UKB_69328/basket_47946_14Sep21/ukb47946.csv"
  covars_metadata <- fromJSON("ukbiobank-blood-trait-analysis/blood_type_covars.json", flatten=TRUE)
  
  blood_trait_covar_ids <- c()
  for (id in names(covars_metadata)) {
    blood_trait_covar_ids <- append(blood_trait_covar_ids, id)
  }

  blood_trait_covars <- fread(biobank_csv, sep=',', header=TRUE,stringsAsFactors=FALSE, select = blood_trait_covar_ids)

  header_names = c()
  for(column in colnames(blood_trait_covars)) {
    column_values <- covars_metadata[[column]]$values
    header_names <<- append(header_names, covars_metadata[[column]]$name)
    
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
  
  blood_trait_covars$pack_years_of_smoking[is.na(blood_trait_covars$pack_years_of_smoking)] <- 0
  blood_trait_covars$time_since_last_menstual_period[is.na(blood_trait_covars$time_since_last_menstual_period)] <- 0
  blood_trait_covars$alcoholic_drinks_yesterday[is.na(blood_trait_covars$alcoholic_drinks_yesterday)] <- 0
  blood_trait_covars$had_menopause[is.na(blood_trait_covars$had_menopause)] <- "no"
  
  blood_trait_covars <- impute.knn(blood_trait_covars, k = 3)
  

  
  
  #TODO: dont need to one hot encode
  #use the mltools.one_hot encoding function
  #one_hot_blood_trait_covars <- one_hot(blood_trait_covars)
  #head(one_hot_blood_trait_covars)
  
  #return(one_hot_blood_trait_covars)
  
  #ex. hemoglobin: so it changes what shows up in the blood sample.  So, they are known to influence blood measurements, so we have to take account for those differences.
  #look at the Figure 5A, and grab the data to put into this table
  #then you can 
  
  #save(list=c("blood_trait_covars"), file="blood_trait_covars.Rdata")
  return(blood_trait_covars)
}




doesthiswork <- phenotype_adjustments_from_json()
head(blood_trait_covars)

blood_markers_from_json <- function() {
  biobank_csv = "/rds/general/project/chadeau_ukbb_folder/live/data/project_data/UKB_69328/basket_47946_14Sep21/ukb47946.csv"
  blood_markers_metadata <- fromJSON("ukbiobank-blood-trait-analysis/blood_markers.json", flatten=TRUE)
  
  blood_marker_ids <- c()
  for (id in names(blood_markers_metadata)) {
    blood_marker_ids <- append(blood_marker_ids, id)
  }
  
  blood_markers <- fread(biobank_csv, sep=',', header=TRUE,stringsAsFactors=FALSE, select = blood_marker_ids)
  
  header_names = c()
  for(column in colnames(blood_markers)) {
    column_values <- blood_markers_metadata[[column]]$values
    header_names <<- append(header_names, blood_markers_metadata[[column]]$name)
    
    if(!is.null(column_values)) {
      #switch the ids with the names of the variables
      mapply(function(x, i) {
        blood_markers[[column]][blood_markers[[column]] == i] <<- x
      }, column_values, names(column_values)
      )
      
      #after swapping vars with names, turns column into a factor
      blood_markers[[column]] <- as.factor(blood_markers[[column]])
    }
  }
  colnames(blood_markers) <- header_names
  return(blood_markers)
  
}

try_this_out <- blood_markers_from_json()




run_gam_for_phenotype <- function(phenotype, covars) {
  covars$marker <- phenotype
  print(head(covars))
  
  #obvisouly don't use id, see what the adjust_scale_trait is all about in the other script
  gam_results = mgcv::gam(marker~s(age_at_recruitment, k=30, bs="ps") + 
                            s(ethnic_background) +
                            s(bmi, k=30, bs="tp") +
                            s(time_since_last_menstual_period, k=30, bs="ps")  + 
                            s(alcoholic_drinks_yesterday, k=30, bs="ps") + 
                            s(pack_years_of_smoking, k=30, bs="ps")
                          , optimizer=c("outer", "newton"), data=covars)
  
  return(gam_results)
}

checking_result = run_gam_for_phenotype(blood_markers$neutrophil_count, blood_trait_covars)
head(blood_markers$neutrophil_count)
