library(data.table)
library(jsonlite)
library(mltools)
library(mgcv)
library(dplyr)
library(bnstruct)
library(ggplot2)

biobank_csv = "/rds/general/project/chadeau_ukbb_folder/live/data/project_data/UKB_69328/basket_47946_14Sep21/ukb47946.csv"
biobank_id_map = "/rds/general/user/are20/home/ukbiobank-blood-trait-analysis/data/link_69328_19266_genotyped.txt"
covars_metadata <- fromJSON("ukbiobank-blood-trait-analysis/blood_type_covars.json", flatten=TRUE)
blood_markers_metadata <- fromJSON("ukbiobank-blood-trait-analysis/blood_markers.json", flatten=TRUE)
results_directory <- "ukbiobank-blood-trait-analysis/results/residuals"

######
# STEP 1: Preparing covars
#####

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

colnames(blood_trait_covars) <- covar_header_names
patient_id_map <- fread(biobank_id_map, sep=' ', header=TRUE,stringsAsFactors=FALSE)
blood_trait_covars <- merge(blood_trait_covars, patient_id_map, by.x="id", by.y="eid_69328")

blood_trait_covars$pack_years_of_smoking[is.na(blood_trait_covars$pack_years_of_smoking)] <- 0
blood_trait_covars$time_since_last_menstual_period[is.na(blood_trait_covars$time_since_last_menstual_period)] <- 0
blood_trait_covars$alcoholic_drinks_yesterday[is.na(blood_trait_covars$alcoholic_drinks_yesterday)] <- 0
blood_trait_covars$had_menopause[is.na(blood_trait_covars$had_menopause)] <- "no"

#TODO: do we need to do this?  NOT WORKING YETTTTT
#imputed <- impute.knn(as.matrix(blood_trait_covars), k = 3)
#imputed = knn.impute(as.matrix(blood_trait_covars), k=3, cat.var=c(2,3,5,9,11), to.impute = 1:nrow(blood_trait_covars), using = 2:nrow(blood_trait_covars))


#TODO: dont need to one hot encode
#use the mltools.one_hot encoding function
#one_hot_blood_trait_covars <- one_hot(blood_trait_covars)
#head(one_hot_blood_trait_covars)

#return(one_hot_blood_trait_covars)

#ex. hemoglobin: so it changes what shows up in the blood sample.  So, they are known to influence blood measurements, so we have to take account for those differences.
#look at the Figure 5A, and grab the data to put into this table
#then you can

#####
# STEP 2: prepare blood marker data,
#####

blood_marker_ids <- c()
for (id in names(blood_markers_metadata)) {
  blood_marker_ids <- append(blood_marker_ids, id)
}

blood_markers <- fread(biobank_csv, sep=',', header=TRUE,stringsAsFactors=FALSE, select = blood_marker_ids)

marker_header_names = c()
for(column in colnames(blood_markers)) {
  column_values <- blood_markers_metadata[[column]]$values
  marker_header_names <<- append(marker_header_names, blood_markers_metadata[[column]]$name)

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
colnames(blood_markers) <- marker_header_names
patient_id_map <- fread(biobank_id_map, sep=' ', header=TRUE,stringsAsFactors=FALSE)
blood_markers <- merge(blood_markers, patient_id_map, by.x="id", by.y="eid_69328")

#####
# STEP 4: run a gam on each marker compared to covars
# done in an extremely tedious manner because I don't understand R functions
####

run_gam <- function(df, name) {
  lower_bounds <- sd(df$marker)*-6
  upper_bounds <-  sd(df$marker)*6
  
  residual_col_names <- c("id", "eid_19266", "residuals")

  df_white <- subset(df, ethnic_background == "white")
  df_white_gam = mgcv::gam(marker ~ s(age_at_recruitment, by=as.factor(had_menopause), k=30, bs="ps") +
                                                              s(bmi, by=as.factor(had_menopause), k=30, bs="tp") +
                                                              s(time_since_last_menstual_period, k=30, bs="ps")  +
                                                              s(alcoholic_drinks_yesterday, k=30, bs="ps") +
                                                              s(pack_years_of_smoking, k=30, bs="ps")
                                                            , optimizer=c("outer", "newton"), data=df_white)

  to_save <- data.frame(df_white$id, df_white$eid_19266, df_white_gam$residuals)
  colnames(to_save) <- residual_col_names
  to_save <- to_save[to_save$residuals < upper_bounds & to_save$residuals > lower_bounds,]
  write.csv(to_save, sprintf("%s/white_%s_residuals.csv", results_directory, name), row.names = FALSE)


  df_asian <- subset(df, ethnic_background == "asian")
  df_asian_gam = mgcv::gam(marker ~ s(age_at_recruitment, by=as.factor(had_menopause), k=30, bs="ps") +
                             s(bmi, by=as.factor(had_menopause), k=30, bs="tp") +
                             s(time_since_last_menstual_period, k=30, bs="ps")  +
                             s(alcoholic_drinks_yesterday, k=30, bs="ps") +
                             s(pack_years_of_smoking, k=30, bs="ps")
                           , optimizer=c("outer", "newton"), data=df_asian)
  

  to_save <- data.frame( df_asian$id, df_asian$eid_19266, df_asian_gam$residuals)
  colnames(to_save) <- residual_col_names
  to_save <- to_save[to_save$residuals < upper_bounds & to_save$residuals > lower_bounds,]
  write.csv(to_save, sprintf("%s/asian_%s_residuals.csv", results_directory, name), row.names = FALSE)


  df_black <- subset(df, ethnic_background == "black")
  df_black_gam = mgcv::gam(marker ~ s(age_at_recruitment, by=as.factor(had_menopause), k=30, bs="ps") +
                             s(bmi, by=as.factor(had_menopause), k=30, bs="tp") +
                             s(time_since_last_menstual_period, k=30, bs="ps")  +
                             s(alcoholic_drinks_yesterday, k=30, bs="ps") +
                             s(pack_years_of_smoking, k=30, bs="ps")
                           , optimizer=c("outer", "newton"), data=df_black)
  
  to_save <- data.frame(df_black$id, df_black$eid_19266, df_black_gam$residuals)
  colnames(to_save) <- residual_col_names
  to_save <- to_save[to_save$residuals < upper_bounds & to_save$residuals > lower_bounds,]
  write.csv(to_save, sprintf("%s/black_%s_residuals.csv", results_directory, name), row.names = FALSE)


  df_chinese <- subset(df, ethnic_background == "chinese")
  df_chinese_gam = mgcv::gam(marker ~ s(age_at_recruitment, by=as.factor(had_menopause), k=30, bs="ps") +
                             s(bmi, by=as.factor(had_menopause), k=30, bs="tp") +
                             s(time_since_last_menstual_period, k=30, bs="ps")  +
                             s(alcoholic_drinks_yesterday, k=30, bs="ps") +
                             s(pack_years_of_smoking, k=30, bs="ps")
                           , optimizer=c("outer", "newton"), data=df_chinese)
  
  to_save <- data.frame( df_chinese$id, df_chinese$eid_19266, df_chinese_gam$residuals)
  colnames(to_save) <- residual_col_names
  to_save <- to_save[to_save$residuals < upper_bounds & to_save$residuals > lower_bounds,]
  write.csv(to_save, sprintf("%s/chinese_%s_residuals.csv", results_directory, name), row.names = FALSE)


  df_mixed <- subset(df, ethnic_background == "mixed")
  df_mixed_gam = mgcv::gam(marker ~ s(age_at_recruitment, by=as.factor(had_menopause), k=30, bs="ps") +
                               s(bmi, by=as.factor(had_menopause), k=30, bs="tp") +
                               s(time_since_last_menstual_period, k=30, bs="ps")  +
                               s(alcoholic_drinks_yesterday, k=30, bs="ps") +
                               s(pack_years_of_smoking, k=30, bs="ps")
                             , optimizer=c("outer", "newton"), data=df_mixed)
  
  to_save <- data.frame( df_mixed$id, df_mixed$eid_19266, df_mixed_gam$residuals)
  colnames(to_save) <- residual_col_names
  to_save <- to_save[to_save$residuals < upper_bounds & to_save$residuals > lower_bounds,]
  write.csv(to_save, sprintf("%s/mixed_%s_residuals.csv", results_directory, name), row.names = FALSE)

  plot(density(df_white_gam$residuals), col="wheat", main=paste("Adjusted Density plot of", name, "By Ethnicity"),
       xlab = name,
       ylab = "Density",
       lwd=2.0,
       #ylim=c(0,3.5)
  )
  lines(density(df_asian_gam$residuals), col="rosybrown3", lwd=2)
  lines(density(df_chinese_gam$residuals), col="slateblue", lwd=2)
  lines(density(df_black_gam$residuals), col="tomato1", lwd=2)
  lines(density(df_mixed_gam$residuals), col="skyblue1", lwd=2)

  legend(x="topright", legend=c("Black", "Asian", "Chinese", "White", "Mixed"),
                          col=c("tomato1", "rosybrown3", "slateblue", "wheat", "skyblue1"), lwd = 2)
}

descriptive_stats <- function(df, name) {
  lower_bounds <- mean(df$marker) - sd(df$marker)*6
  upper_bounds <- mean(df$marker) + sd(df$marker)*6

  white_marker <- subset(df, ethnic_background == "white")$marker
  white_marker <- white_marker[white_marker < upper_bounds & white_marker > lower_bounds]

  black_marker <- subset(df, ethnic_background == "black")$marker
  black_marker <- black_marker[black_marker < upper_bounds & black_marker > lower_bounds]

  asian_marker <- subset(df, ethnic_background == "asian")$marker
  asian_marker <- asian_marker[asian_marker < upper_bounds & asian_marker > lower_bounds]

  chinese_marker <- subset(df, ethnic_background == "chinese")$marker
  chinese_marker <- chinese_marker[chinese_marker < upper_bounds & chinese_marker > lower_bounds]

  mixed_marker <- subset(df, ethnic_background == "mixed")$marker
  mixed_marker <- mixed_marker[mixed_marker < upper_bounds & mixed_marker > lower_bounds]


  plot(density(white_marker), col="wheat", main=paste("Descriptive plot of", name, "By Ethnicity"),
       xlab = name,
       ylab = "Density",
       lwd=2.0,
       #ylim=c(0,0.3)
  )
  lines(density(asian_marker), col="rosybrown3", lwd=2)
  lines(density(black_marker), col="tomato1", lwd=2)
  lines(density(chinese_marker), col="slateblue", lwd=2)
  lines(density(mixed_marker), col="skyblue1", lwd=2)

  legend(x="topright", legend=c("Black", "Asian", "Chinese", "White", "Mixed"),
         col=c("tomato1", "rosybrown3", "slateblue", "wheat", "skyblue1"), lwd = 2)
}

#TODO: cut out the outliers after the gam is run, maybe at a 6 SD either side of the mean

#TODO: compare counts between ethnicities on raw values (for example, neutrophil counts are lower with black ethnicities)

#TODO: start writing up the basic descripctive info about the data set (age, sex, etc)
#TODO: plot all different blood type covars and overlay all different ethniticies for each covar.  This will help understand the baseline of the data
# we are working with.

#points.density(second data set)


#THIS WORKS NOW YAY
blood_trait_covars$marker <- blood_markers$red_blood_cell_distribution_width
red_blood_cell_distribution_width <- subset(blood_trait_covars, marker > 2) #change this
red_blood_cell_distribution_width <- red_blood_cell_distribution_width[complete.cases(red_blood_cell_distribution_width),]
descriptive_stats(red_blood_cell_distribution_width, "Red Blood Cell Distribution Width")
run_gam(red_blood_cell_distribution_width, "red_blood_cell_distribution_width")

blood_trait_covars$marker <- blood_markers$red_blood_cell_count
red_blood_cell_count <- subset(blood_trait_covars, marker > 1)
red_blood_cell_count <- red_blood_cell_count[complete.cases(red_blood_cell_count),]
descriptive_stats(red_blood_cell_count, "Red Blood Cell Count")
run_gam(red_blood_cell_count, "red_blood_cell_count")


blood_trait_covars$marker <- blood_markers$haemoglobin_concentration
haemoglobin_concentration <- subset(blood_trait_covars, marker < 22) #change this
haemoglobin_concentration <- haemoglobin_concentration[complete.cases(haemoglobin_concentration),]
descriptive_stats(haemoglobin_concentration, "Haemoglobin Concentration")
run_gam(haemoglobin_concentration, "haemoglobin_concentration")


blood_trait_covars$marker <- blood_markers$haematocrit_percentage
haematocrit_percentage <- subset(blood_trait_covars, marker > 10) #change this
haematocrit_percentage <- haematocrit_percentage[complete.cases(haematocrit_percentage),]
descriptive_stats(haematocrit_percentage, "Haematocrit Percentage")
run_gam(haematocrit_percentage, "haematocrit_percentage")


blood_trait_covars$marker <- blood_markers$mean_corpuscular_volume
mean_corpuscular_volume <- subset(blood_trait_covars, marker < 140) #change this
mean_corpuscular_volume <- mean_corpuscular_volume[complete.cases(mean_corpuscular_volume),]
descriptive_stats(mean_corpuscular_volume, "Mean Corpuscular Volume")
run_gam(mean_corpuscular_volume, "mean_corpuscular_volume")


blood_trait_covars$marker <- blood_markers$mean_corpuscular_haemoglobin
mean_corpuscular_haemoglobin <- subset(blood_trait_covars, marker > 7 & marker < 90) #change this
mean_corpuscular_haemoglobin <- mean_corpuscular_haemoglobin[complete.cases(mean_corpuscular_haemoglobin),]
descriptive_stats(mean_corpuscular_haemoglobin, "Mean Corpuscular Haemoglobin")
run_gam(mean_corpuscular_haemoglobin, "mean_corpuscular_haemoglobin")

blood_trait_covars$marker <- blood_markers$mean_corpuscular_haemoglobin_concentration
mean_corpuscular_haemoglobin_concentration <- subset(blood_trait_covars, marker <= 80) #change this
mean_corpuscular_haemoglobin_concentration <- mean_corpuscular_haemoglobin_concentration[complete.cases(mean_corpuscular_haemoglobin_concentration),]
descriptive_stats(mean_corpuscular_haemoglobin_concentration, "Mean Corpuscular Haemoglobin Concentration")
run_gam(mean_corpuscular_haemoglobin_concentration, "mean_corpuscular_haemoglobin_concentration")


blood_trait_covars$marker <- blood_markers$mean_platelet_volume
mean_platelet_volume <- blood_trait_covars
mean_platelet_volume <- mean_platelet_volume[complete.cases(mean_platelet_volume),]
descriptive_stats(mean_platelet_volume, "Mean Platelet Volume")
run_gam(mean_platelet_volume, "mean_platelet_volume")


#blood_trait_covars$marker <- blood_markers$platelet_count
#platelet_count <- subset(blood_trait_covars, marker <= 1000) #change this
#platelet_count <- platelet_count[complete.cases(platelet_count),]
#descriptive_stats(platelet_count, "Platelet Count")
#run_gam(platelet_count, "platelet_count")



blood_trait_covars$marker <- blood_markers$basophil_count
basophil_count <- subset(blood_trait_covars, marker < 2.5) #change this
basophil_count <- basophil_count[complete.cases(basophil_count),]
descriptive_stats(basophil_count, "Basophil Count")
run_gam(basophil_count, "basophil_count")

blood_trait_covars$marker <- blood_markers$neutrophil_count
neutrophil_count <- subset(blood_trait_covars, marker < 40)
neutrophil_count <- neutrophil_count[complete.cases(neutrophil_count),]
descriptive_stats(neutrophil_count, "Neutrophil Count")
run_gam(neutrophil_count, "neutrophil_count")


blood_trait_covars$marker <- blood_markers$monocyte_count
monocyte_count <- subset(blood_trait_covars, marker < 34) #change this
monocyte_count <- monocyte_count[complete.cases(monocyte_count),]
descriptive_stats(monocyte_count, "Monocyte Count")
run_gam(monocyte_count, "monocyte_count")


blood_trait_covars$marker <- blood_markers$eosinophill_count
eosinophill_count <- subset(blood_trait_covars, marker < 8) #change this
eosinophill_count <- eosinophill_count[complete.cases(eosinophill_count),]
descriptive_stats(eosinophill_count, "Eosinophill Count")
run_gam(eosinophill_count, "eosinophill_count")

blood_trait_covars$marker <- blood_markers$white_blood_cell_count
white_blood_cell_count <- subset(blood_trait_covars, marker < 200) #change this
white_blood_cell_count <- white_blood_cell_count[complete.cases(white_blood_cell_count),]
descriptive_stats(white_blood_cell_count, "White Blood Cell Count")
run_gam(white_blood_cell_count, "white_blood_cell_count")

blood_trait_covars$marker <- blood_markers$lymphocyte_count
lymphocyte_count <- subset(blood_trait_covars, marker < 140) #change this
lymphocyte_count <- lymphocyte_count[complete.cases(lymphocyte_count),]
descriptive_stats(lymphocyte_count, "Lymphocyte Count")
run_gam(lymphocyte_count, "lymphocyte_count")
