#####

library(populationmodels)
library(rwdcohort)
library(cobalt)
library(spatstat)
library(readr)

#######################################################################
# RWD data.
#######################################################################

rwd_data <- create_flatiron_data(treatment = "Crizotinib",
                                 data_cut_off_date = as.Date("2017-02-09"),
                                 censoring_strategy = "administrative_cutoff",
                                 SoC_date = as.Date("2011-08-26"))



summary(as.factor(rwd_data$smoking_status))
summary(as.factor(rwd_data$histology))
View(rwd_data)
#summary(as.factor(rwd_data))

ALEX_table1 <-
  tibble::tribble(
    ~characteristic, ~mean, ~sd, ~n,
    "age", 53.8, 13.5, 151,
    "sex_male", 0.424, NA, 64,
    "sex_female", 0.576, NA, 87,
    "race_asian", 0.457, NA,  69,
    "race_non_asian", 0.543, NA, 82,
    "ecog_0_1", 0.934, NA, 141,
    "ecog_2", 0.066, NA, 10,
    "smoking_status_ever", 0.351, NA, 53,
    "smoking_status_never", 0.649, NA, 98,
    "brain_mets", 0.384, NA, 58
  )

ALEX_variables_to_match =  c(
  "age",
  "sex_male",
  "ecog_0_1",
  "smoking_status_ever",   
  "brain_mets",
  "race_asian"
)

# Cohort with MAIC weights.

rwd_data_weighted <- create_weighted_cohort(
  cohort_data = rwd_data,
  reference_table1 = ALEX_table1,
  match_variables = ALEX_variables_to_match,
)[["cohort_data"]]  %>%
  rename(weight = patient_weight) %>%
  mutate(sex_male = 1*(sex == "Male"),
         sex_female = 1*(sex == "Female"),
         race_asian = 1*(race == "Asian"),
         race_non_asian = 1*(race != "Asian"),
         smoking_ever = 1*(smoking_status_ever == 1),
         smoking_never = 1*(smoking_status_ever == 0),
         ecog_0_or_1 = 1*ecog_0_1,
         ecog_2 = 1*!ecog_0_1,
         disease_local_advance = 1*(!is_stage_iv),
         disease_metastatic = 1*is_stage_iv,
         brain_mets_yes = 1*brain_mets,
         brain_mets_no = 1*(!brain_mets))

table1 <- tibble(variable = c("age_median", "age_range",
                    "sex_male", "sex_female",
                    "race_asian", "race_non_asian",
                    "smoking_ever", "smoking_never",
                    "ecog_0_or_1", "ecog_2",
                    "disease_local", "disease_meta",
                    "brain_mets_yes", "brain_mets_no")) %>%
  mutate(type = c("cts", "cts_range",
                  rep("binary",12))) %>%
  mutate(flatiron_column = c("age_at_lot_1_startdate", 
                             "age_at_lot_1_startdate",
                             "sex_male",
                             "sex_female",
                             "race_asian",
                             "race_non_asian",
                             "smoking_ever", 
                             "smoking_never",
                             "ecog_0_or_1", 
                             "ecog_2",
                             "disease_local_advance",
                             "disease_metastatic",
                             "brain_mets_yes",
                             "brain_mets_no")) %>%
  mutate(unweight_value =  0, 
        weight_value = 0)

for(i in 1:nrow(table1)){
  if(table1$type[i] == "binary"){

    variable_temp <- table1$flatiron_column[i]
    n_temp <- sum(rwd_data_weighted[variable_temp]) 
    perc_temp <- 100*sum(rwd_data_weighted[variable_temp])/nrow(rwd_data_weighted) 
    perc_temp <- round(perc_temp, 0)
    
    table1$unweight_value[i] <- paste0(n_temp, " (", perc_temp, "%)")
    
    temp_data <-  rwd_data_weighted[c(variable_temp,"weight")] 
    names(temp_data) <- c("variable", "weight")
    n_weight_temp <- sum(temp_data$variable*temp_data$weight) 
    n_weight_temp <- round(n_weight_temp, 1)
    perc_weight_temp <- 100*sum(temp_data$variable*temp_data$weight)/sum(temp_data$weight)
    perc_weight_temp <- round(perc_weight_temp, 0)
    
    table1$weight_value[i] <- paste0(n_weight_temp, " (", perc_weight_temp, "%)")
    
    
  } else if(table1$type[i] == "cts"){
      
    variable_temp <- table1$flatiron_column[i]
    
    temp_data <-  rwd_data_weighted[c(variable_temp,"weight")] 
    names(temp_data) <- c("variable", "weight")

    median_temp <- median(temp_data$variable)
    median_weight_temp <- weighted.median(temp_data$variable, temp_data$weight)

    table1$unweight_value[i] <- median_temp
    table1$weight_value[i] <- median_weight_temp
    
    }else if(table1$type[i] == "cts_range"){
      i <- 2
      variable_temp <- table1$flatiron_column[i]
      
      temp_data <-  rwd_data_weighted[c(variable_temp,"weight")] 
      names(temp_data) <- c("variable", "weight")
      
      range_temp <- paste0(min(temp_data$variable), "-",
                           max(temp_data$variable))
 
      table1$unweight_value[i] <- range_temp
      table1$weight_value[i] <- range_temp
 
    }
      
    
      
}


write_csv(table1, "Flatiron_data/Table1_rwe.csv")
