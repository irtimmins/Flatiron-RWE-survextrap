##############
library(populationmodels)
library(rwdcohort)
library(cobalt)
library(spatstat)
library(readr)
library(ggplot2)

#######################################################################
# RWD data.
#######################################################################
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

rwd_data <- create_flatiron_data(treatment = "Crizotinib",
                                 data_cut_off_date = as.Date("2017-02-09"),
                                 censoring_strategy = "administrative_cutoff",
                                 SoC_date = as.Date("2011-08-26"))


rwd_data_weighted <- create_weighted_cohort(
  cohort_data = rwd_data,
  reference_table1 = ALEX_table1,
  match_variables = ALEX_variables_to_match,
)[["cohort_data"]]  %>%
  rename(weight = patient_weight) 


rwd_data_weighted %>%
  ggplot()+
  geom_histogram(aes(x = weight), colour = "gray50", binwidth = 0.01)+
#  geom_histogram(aes(x = weight, colour = race_asian, fill = race_asian), binwidth = 0.01)+
  theme_classic()+
  scale_y_continuous("Frequency")+
  scale_x_continuous("Weight")









