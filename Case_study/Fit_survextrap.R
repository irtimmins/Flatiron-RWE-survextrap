
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
library(purrr)


######################################################
#  Memoise functions.
######################################################

store_wd <- "/scratch/klvq491/case_study_nice_ta536/"

survextrap_mem <- memoise(survextrap, cache = cachem::cache_disk(dir = "../.cache",
                                                                 max_age = 31557600, # keep for 1-year
                                                                 max_size = 2e10)) # 20Gb of space 

rmst_mem <- memoise(rmst, cache = cachem::cache_disk(dir = "../.cache",
                                                          max_age = 31557600,
                                                          max_size = 2e10))

irmst_mem <- memoise(irmst, cache = cachem::cache_disk(dir = "../.cache",
                                                            max_age = 31557600,
                                                            max_size = 2e10))



######################################################
#  Get data
######################################################

trial_data <- readRDS("Data/trial_data.rds")
historic_trial_aggregate <- readRDS("Data/historic_trial_aggregate_no_overlap.rds")
external_data_aggregate <- readRDS("Data/external_data_no_overlap.rds")



######################################################
#  Models to fit.
######################################################
#?survextrap()

base_scenarios <- expand_grid(
  model = c("PH", "NON-PH", "Separate_arms"),
  datasets = c("trial_only", "trial_and_historic", "trial_and_all")) %>%
  mutate(
    df = 3,
    hsd_rate = 5,
    hrsd_rate = if_else(model == "NON-PH", 1, NA)) %>%
  mutate(
    store_file = paste0(store_wd, "base_model_", row_number(), ".rds"),
    hazard_survival_file = paste0(store_wd, "base_model_", row_number(), "_hs.rds"),
    rmst_file = paste0(store_wd, "base_model_", row_number(), "_rmst.rds"))
  
knots_sensitivity <- expand_grid(
  model = "PH",
  datasets = c("trial_only", "trial_and_historic", "trial_and_all"),
  df = c(3,6,10)) %>%
  mutate(
    hsd_rate = 5,
    hrsd_rate = if_else(model == "NON-PH", 1, NA)) %>%
  mutate(
    store_file = paste0(store_wd, "knots_model_", row_number(), ".rds"),
    hazard_survival_file = paste0(store_wd, "knots_model_", row_number(), "_hs.rds"),
    rmst_file = paste0(store_wd, "knots_model_", row_number(), "_rmst.rds"))


######################################################
#  Get survival, hazard and rmst.
######################################################



pmap(base_scenarios %>% 
       select(-hazard_survival_file, -rmst_file),
     fit_model)

pmap(base_scenarios %>% 
       select(store_file, hazard_survival_file) %>%
       rename(model_file = store_file, store_file =  hazard_survival_file),
     get_survival_and_hazard_survextrap)

pmap(base_scenarios %>% 
       select(store_file, rmst_file) %>%
       rename(model_file = store_file, store_file =  rmst_file),
     get_rmst_survextrap)

pmap(knots_sensitivity %>% 
       select(-hazard_survival_file, -rmst_file),
     fit_model)


# test <- readRDS("/scratch/klvq491/case_study_nice_ta536/base_model_1.rds")
# test2 <- readRDS("/scratch/klvq491/case_study_nice_ta536/base_model_7.rds")
# class(test)
# class(test2)
# test2$control

#test <- readRDS("/scratch/klvq491/case_study_nice_ta536/base_model_1_hs.rds")
#test2 <- readRDS("/scratch/klvq491/case_study_nice_ta536/base_model_7_hs.rds")

test <- readRDS("/scratch/klvq491/case_study_nice_ta536/base_model_1_rmst.rds")
test2 <- readRDS("/scratch/klvq491/case_study_nice_ta536/base_model_7_rmst.rds")
#head(test)
#tail(test)
#head(test2)
View(test)
View(test2)

##################################################
# Fit survextrap models to two arms.
##################################################

fit_model <- function(model, datasets, df , hsd_rate, hrsd_rate, store_file){
  
#  print(df)
  
  if(datasets == "trial_only") external_data <- NULL
  if(datasets == "trial_and_historic") external_data <- historic_trial_aggregate
  if(datasets == "trial_and_all") external_data <- external_data_aggregate
  
  if(model == "PH"){

    results <- survextrap_mem(Surv(time,status) ~ trt, 
                              data = trial_data , 
                              external = external_data,
                              mspline = list("df"=df),
                              prior_hsd = p_gamma(2, hsd_rate),
                              fit_method = "opt")
    
    
  } else if(model == "NON-PH"){
    
    results <- survextrap_mem(Surv(time,status) ~ trt, 
                              data = trial_data , 
                              nonprop = T,
                              external = external_data,
                              mspline = list("df"=df),
                              prior_hsd = p_gamma(2, hsd_rate),
                              prior_hrsd = p_gamma(2, hrsd_rate),
                              fit_method = "opt")
    
  }else if(model == "Separate_arms"){
    
    results_control <- survextrap_mem(Surv(time,status) ~ 1, 
                                      data = trial_data %>% filter(trt == "Crizotinib") , 
                                      external = external_data,
                                      mspline = list("df"=df),
                                      prior_hsd = p_gamma(2, hsd_rate),
                                      fit_method = "opt")
    
    
    results_active <- survextrap_mem(Surv(time,status) ~ 1, 
                                      data = trial_data %>% filter(trt == "Alectinib") , 
                                      mspline = list("df"=df),
                                      prior_hsd = p_gamma(2, hsd_rate),
                                      fit_method = "opt")
    
    results <- list(control = results_control, active = results_active)
    class(results) <- "two_models"
    
    }
  
    saveRDS(results, store_file)
  
}


##################################################
# Survival and hazard estimates from models.
##################################################

get_survival_and_hazard_survextrap <- function(model_file, store_file){
    
  #print(model_file)
  
  model <- readRDS(model_file)

  new_data_estimate <- tibble(trt = c("Alectinib", "Crizotinib"))
                           
  time_vec <- c(seq(from = 0, to = 5, length.out = 1e2),
                seq(from = 5, to = 20, length.out = 1e2))
    

  if(class(model) == "survextrap"){ # PH or NON-PH scenarios.
    
    hazard_all <- hazard(model, t = time_vec, newdata = new_data_estimate) 
   # print(hazard_all)
    survival_all <-survival(model, t = time_vec, newdata = new_data_estimate) 
  #  print(survival_all)

    results <- bind_rows(hazard_all,
                         survival_all)
    
  } else if(class(model) == "two_models"){ # separate arms scenarios.
    
    control_model <- model$control
    active_model <- model$active
    
    hazard_control <- hazard(control_model, t = time_vec) %>%
      mutate(trt = "Crizotinib")
    
    hazard_active <- hazard(active_model, t = time_vec) %>%
      mutate(trt = "Alectinib")
    
    survival_control <- survival(control_model, t = time_vec) %>%
      mutate(trt = "Crizotinib")
    
    survival_active <- survival(active_model, t = time_vec) %>%
      mutate(trt = "Alectinib")
    
    results <- bind_rows(hazard_active,
                         hazard_control,
                         survival_active,
                         survival_control) %>%
    select(variable, trt, t, median, lower, upper)
  }

  saveRDS(results, store_file)

}


##################################################
# rmst and difference in rmst (irmst) from models.
##################################################


get_rmst_survextrap <- function(model_file, store_file){
  
  print(model_file)
  
  model <- readRDS(model_file)
  
  new_data_estimate <- tibble(trt = c("Alectinib", "Crizotinib"))
  
  time_vec <- c(5, 10, 15, 20)
  
  
  if(class(model) == "survextrap"){ # PH or NON-PH scenarios.
    
    rmst_both_arms <- rmst(model, 
                     t = time_vec, 
                     newdata = new_data_estimate) 

    irmst_active_vs_control <- irmst(model, 
                                     t = time_vec, 
                                     newdata =   tibble(trt = c("Crizotinib", "Alectinib"))) %>%
      mutate(trt = NA) %>%
      select(variable, trt, t, median, lower, upper)
    
    results <- rmst_both_arms %>%
      bind_rows(irmst_active_vs_control)
    
  } else if(class(model) == "two_models"){ # separate arms scenarios.
    
    control_model <- model$control
    
    active_model <- model$active
    
    rmst_control <- rmst(control_model, t = time_vec)  %>%
      mutate(trt = "Crizotinib")

    rmst_active <- rmst(active_model, t = time_vec) %>%
      mutate(trt = "Alectinib")
    
    control_sample <- rmst(control_model, t = time_vec, sample = T)
    
    active_sample <- rmst(active_model, t = time_vec, sample = T)

    irmst_samples <- active_sample - control_sample
    
    irmst_active_vs_control <- survextrap:::summarise_output(
      irmst_samples,
      t =  time_vec, 
      summ_fns = list("median" = median, 
                      ~quantile(.x, probs=c(0.025, 0.975))),
      newdata = NULL, 
      summ_name = "irmst", 
      sample = F
      ) %>%
      rename(lower= "2.5%", upper = "97.5%") %>%
      mutate(trt = NA) %>%
      select(variable, trt, t, median, lower, upper)
    
    results <- bind_rows(
      rmst_active,
      rmst_control
      ) %>%
      select(variable, trt, t, median, lower, upper) %>%
      bind_rows(irmst_active_vs_control) 
  }
  
  saveRDS(results, store_file)
  
}


