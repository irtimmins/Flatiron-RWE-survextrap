
######################################################
#  Load packages.
######################################################

library(populationmodels)
library(rwdcohort)
library(ggplot2)
library(survival)
library(survminer)
library(maic)
library(azci)
library(bshazard)
library(dplyr)
library(tidyr)
library(ggplot2)
library(survival)
library(survminer)
library(survextrap)
library(ggh4x)
library(condsurv)
library(cowplot)
library(flexsurv)
library(utile.visuals)


######################################################
#  Derive Flatiron cohort for Crizotinib
######################################################

# Source dependencies.
source("Functions/Build_Flatiron_cohort.R")
source("Functions/Derive_conditional_survival.R")
source("Functions/Hybrid_model.R")

######################################################
#  Derive Flatiron cohort for Crizotinib
######################################################

# 262 - 39 ALK confirmed negative, 206 ALK confirmed positive. 17 unknown.
# Filter down to 206 ALK confirmed positive.

rwd_data <- create_flatiron_data(treatment = "Crizotinib",
                                 data_cut_off_date = as.Date("2017-02-09"),
                                 censoring_strategy = "administrative_cutoff",
                                 SoC_date = as.Date("2011-08-26"))


summary(rwd_data$age_at_lot_1_startdate)


######################################################
#  Apply MAIC weighting to Flatiron cohort
######################################################

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
  rename(weight = patient_weight)

sum(rwd_data_weighted$weight)^2/sum(rwd_data_weighted$weight^2)
#head(rwd_data_weighted$weight)
# Cohort without MAIC weights.
  
rwd_data_unweighted <- rwd_data_weighted %>%
  mutate(weight = 1,
         effective_sample_size = n())

######################################################
#  Plot RWD KM curves alongside trial data.
######################################################

trial_data <- readRDS("Data/trial_IPD_OS_ALEX_Feb_2017_v4.rds")
historic_trial_data <- readRDS("Data/trial_IPD_OS_PROFILE_1014_Nov_2016.rds")

km_trial_data  <- survfit(Surv(time, status) ~ trt, data=trial_data)
km_trial_plot <- ggsurvplot(km_trial_data, data=trial_data)["data.survplot"] [[1]] %>%
  select(time, surv, lower, upper, trt) %>%
  mutate(dataset = "ALEX") %>%
  as_tibble() %>%
  add_row(time = 0, surv = 1, trt = "Alectinib", dataset = "ALEX") %>%
  add_row(time = 0, surv = 1, trt = "Crizotinib", dataset = "ALEX") %>%
  mutate(trt = as.factor(trt))

km_historic_trial_data  <- survfit(Surv(time, status) ~ 1, data=historic_trial_data)
km_historic_trial_plot <- ggsurvplot(km_historic_trial_data, data=historic_trial_data)["data.survplot"][[1]] %>%
  mutate(trt = "Crizotinib",
         dataset = "PROFILE-1014")  %>%
  select(time, surv, lower, upper, trt, dataset) %>%
  as_tibble()

km_rwd_data_weighted <- survfit(Surv(time, status) ~ 1, data=rwd_data_weighted, 
                               weights = rwd_data_weighted$weight)
km_rwd_plot <- ggsurvplot(km_rwd_data_weighted, data=rwd_data_weighted)["data.survplot"][[1]] %>%
  as_tibble() %>%
  mutate(trt = "Crizotinib",
         dataset = "Flatiron RWE MAIC Weighted") %>%
  select(time, surv, lower, upper, trt, dataset)

km_rwd_data_unweighted <- survfit(Surv(time, status) ~ 1, data=rwd_data_unweighted)
km_rwd_plot_unweighted <- ggsurvplot(km_rwd_data_unweighted, data=rwd_data_unweighted)["data.survplot"][[1]] %>%
  mutate(trt = "Crizotinib",
         dataset = "Flatiron RWE Unweighted")  %>%
  select(time, surv, lower, upper, trt, dataset)

km_all <- bind_rows(km_trial_plot,
                    km_historic_trial_plot,
                    km_rwd_plot,
                    km_rwd_plot_unweighted)

km_all


######################################################
#  Create aggregate counts for RWD and historic data.
######################################################

increment_width <- 0.5
left_truncate_ALEX <- 2.5

rwd_maic_weighted_aggregate <- create_aggregate_counts(
  data = rwd_data_weighted %>%
    select(time, status, weight) %>%
    mutate(trt = "Crizotinib"),
  increment = increment_width,
  left_truncate = left_truncate_ALEX) %>%
  mutate(data = "Flatiron RWE")

rwd_unweighted_aggregate <- create_aggregate_counts(
  data = rwd_data_unweighted %>%
    select(time, status, weight) %>%
    mutate(trt = "Crizotinib"),
  increment = increment_width,
  left_truncate = left_truncate_ALEX) %>%
  mutate(data = "Flatiron RWE")

historic_trial_aggregate <- create_aggregate_counts(
  data = historic_trial_data %>%
    mutate(weight = 1),
  increment = increment_width,
  left_truncate = left_truncate_ALEX) %>% 
  mutate(data = "PROFILE-1014 trial") %>%
  rename(dataset = data)  %>%
  mutate(dataset = factor(dataset, 
                          levels = c("ALEX trial", "PROFILE-1014 trial")))

external_data_maic_weighted <- historic_trial_aggregate  %>%
  bind_rows(rwd_maic_weighted_aggregate) %>%
  filter(trt == "Crizotinib") %>%
  mutate(dataset = factor(dataset, 
                          levels = c("ALEX trial", "PROFILE-1014 trial", "Flatiron RWE")))

external_data_unweighted <- historic_trial_aggregate  %>%
  bind_rows(rwd_unweighted_aggregate) %>%
  filter(trt == "Crizotinib") %>%
  mutate(dataset = factor(dataset, 
                          levels = c("ALEX trial", "PROFILE-1014 trial", "Flatiron RWE")))

trial_data <- trial_data %>%
  mutate(dataset = "ALEX trial",
         dataset = factor(dataset, levels = c("ALEX trial", "PROFILE-1014 trial", "Flatiron RWE"))) 

saveRDS(trial_data, "Data/trial_data.rds")
saveRDS(historic_trial_aggregate, "Data/historic_trial_aggregate.rds")
saveRDS(external_data_maic_weighted, "Data/external_data_maic_weighted.rds")
saveRDS(external_data_unweighted, "Data/external_data_unweighted.rds")

