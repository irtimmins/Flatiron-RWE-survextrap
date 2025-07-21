
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
library(rslurm)


######################################################
#  Get data
######################################################

trial_data <- readRDS("Data/trial_data.rds")
historic_trial_aggregate <- readRDS("Data/historic_trial_aggregate_no_overlap.rds")
external_data_aggregate <- readRDS("Data/external_data_no_overlap.rds")



######################################################
#  Specify scenarios.
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
# Fit models, using pmap or slurm.
######################################################

# using pmap.

pmap(base_scenarios %>% 
       select(-hazard_survival_file, -rmst_file),
     fit_model)

pmap(knots_sensitivity %>% 
       select(-hazard_survival_file, -rmst_file),
     fit_model)

######################################################

# or using slurm:

user <- Sys.info()["user"]
check_status <- paste0("sacct -S ", as.character(Sys.Date()-20),
                       " -u ", user ,
                       " --format=JobID,Jobname,partition,state,elapsed,ncpus -X")

objects_attach <- c("trial_data",
                    "historic_trial_aggregate",
                    "external_data_aggregate",
                    "survextrap_mem",
                    "rmst_mem", 
                    "irmst_mem")

package_attach <- c("dplyr", "tidyr", "readr",
                   "survextrap", "rstan", "survival",
                    "posterior")

fit_base_slurm <- slurm_apply(
  fit_model, 
  base_scenarios %>% 
    select(-hazard_survival_file, -rmst_file), 
  jobname = "fit_base",
  nodes = 4, 
  cpus_per_node = 4, 
  submit = T,
  global_objects = objects_attach,
  pkgs = package_attach,
  slurm_options = list(time='01:00:00',
                       partition='core',
                       "mem-per-cpu"= '16G'))

fit_sensistivity_slurm <- slurm_apply(
  fit_model, 
  knots_sensitivity %>% 
    select(-hazard_survival_file, -rmst_file), 
  jobname = "fit_knots_sensitivity",
  nodes = 4, 
  cpus_per_node = 4, 
  submit = T,
  global_objects = objects_attach,
  pkgs = package_attach,
  slurm_options = list(time='01:00:00',
                       partition='core',
                       "mem-per-cpu"= '16G'))



# check results.
system(check_status)
# results_0 <- readRDS("_rslurm_fit_base/results_0.RDS")
# results_0 <- readRDS("_rslurm_fit_knots_sensitivity/results_0.RDS")
# results_0

######################################################
#  Get survival, hazard and rmst.
######################################################

pmap(base_scenarios %>% 
       select(store_file, hazard_survival_file) %>%
       rename(model_file = store_file, store_file =  hazard_survival_file),
     get_survival_and_hazard_survextrap)

pmap(base_scenarios %>% 
       select(store_file, rmst_file) %>%
       rename(model_file = store_file, store_file =  rmst_file),
     get_rmst_survextrap)

pmap(knots_sensitivity %>% 
       select(store_file, hazard_survival_file) %>%
       rename(model_file = store_file, store_file =  hazard_survival_file),
     get_survival_and_hazard_survextrap)

pmap(knots_sensitivity %>% 
       select(store_file, rmst_file) %>%
       rename(model_file = store_file, store_file =  rmst_file),
     get_rmst_survextrap)

# test <- readRDS("/scratch/klvq491/case_study_nice_ta536/base_model_1_rmst.rds")
# test2 <- readRDS("/scratch/klvq491/case_study_nice_ta536/base_model_7_rmst.rds")
# 
# View(test)
# View(test2)


