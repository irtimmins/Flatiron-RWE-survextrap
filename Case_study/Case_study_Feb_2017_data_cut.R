
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
  select(time, surv, trt) %>%
  mutate(dataset = "ALEX") %>%
  as_tibble() %>%
  add_row(time = 0, surv = 1, trt = "Alectinib", dataset = "ALEX") %>%
  add_row(time = 0, surv = 1, trt = "Crizotinib", dataset = "ALEX") %>%
  mutate(trt = as.factor(trt))

km_historic_trial_data  <- survfit(Surv(time, status) ~ 1, data=historic_trial_data)
km_historic_trial_plot <- ggsurvplot(km_historic_trial_data, data=historic_trial_data)["data.survplot"][[1]] %>%
  mutate(trt = "Crizotinib",
         dataset = "PROFILE-1014") %>%
  select(time, surv, trt, dataset) %>%
  as_tibble()

km_rwd_data_weighted <- survfit(Surv(time, status) ~ 1, data=rwd_data_weighted, 
                               weights = rwd_data_weighted$weight)
km_rwd_plot <- ggsurvplot(km_rwd_data_weighted, data=rwd_data_weighted)["data.survplot"][[1]] %>%
  mutate(trt = "Crizotinib",
         dataset = "Flatiron RWE MAIC Weighted") %>%
  select(time, surv, trt, dataset)

km_rwd_data_unweighted <- survfit(Surv(time, status) ~ 1, data=rwd_data_unweighted)
km_rwd_plot_unweighted <- ggsurvplot(km_rwd_data_unweighted, data=rwd_data_unweighted)["data.survplot"][[1]] %>%
  mutate(trt = "Crizotinib",
         dataset = "Flatiron RWE Unweighted") %>%
  select(time, surv, trt, dataset)

km_all <- bind_rows(km_trial_plot,
                    km_historic_trial_plot,
                    km_rwd_plot,
                    km_rwd_plot_unweighted)

plot0 <- km_all %>%
  mutate(dataset = factor(dataset, 
                          levels = c("ALEX", "PROFILE-1014", 
                                     "Flatiron RWE MAIC Weighted", "Flatiron RWE Unweighted"))) %>%
  filter(dataset == "ALEX") %>%
  ggplot(aes(x = time, y = surv, colour = trt))+
  theme_classic()+
  geom_step(linewidth = 0.8)+
  scale_colour_discrete("Treatment")+
  scale_y_continuous("Overall Survival", 
                     limits = c(0,1), labels = scales::percent)+
  scale_x_continuous("Time (years)") + ##, limits = c(0,6))
  ggtitle("ALEX trial, Overall survival")

plot0

plot1 <- km_all %>%
  mutate(dataset = factor(dataset, 
                          levels = c("ALEX", "PROFILE-1014", 
                                     "Flatiron RWE MAIC Weighted", "Flatiron RWE Unweighted"))) %>%
  ggplot(aes(x = time, y = surv, colour = dataset , linetype = trt))+
  theme_classic()+
  geom_step(linewidth = 0.8)+
  scale_colour_discrete("Dataset")+
  scale_linetype_discrete("Treatment")+
  scale_y_continuous("Overall Survival", 
                     limits = c(0,1), labels = scales::percent)+
  scale_x_continuous("Time (years)") + ##, limits = c(0,6))
  ggtitle("ALEX, PROFILE-1014 and Flatiron RWE")

plot1

plot2 <- km_all %>%
  mutate(dataset = factor(dataset, 
                          levels = c("ALEX", "PROFILE-1014", 
                                     "Flatiron RWE MAIC Weighted", "Flatiron RWE Unweighted"))) %>%
  filter(dataset != "PROFILE 1014") %>%
  filter(trt == "Crizotinib") %>%
  ggplot(aes(x = time, y = surv, colour = dataset , linetype = trt))+
  theme_classic()+
  geom_step(linewidth = 0.8)+
  scale_colour_discrete("Dataset")+
  scale_linetype_discrete("Treatment")+
  scale_y_continuous("Overall Survival", 
                     limits = c(0,1), labels = scales::percent)+
  scale_x_continuous("Time (years)")

plot2

plot3 <- km_all %>%
  mutate(dataset = factor(dataset, 
                          levels = c("ALEX", "PROFILE-1014", 
                                     "Flatiron RWE MAIC Weighted", "Flatiron RWE Unweighted"))) %>%
  filter(trt == "Crizotinib") %>%
  ggplot(aes(x = time, y = surv, colour = dataset))+
  theme_classic()+
  geom_step(linewidth = 0.8)+
  scale_colour_discrete("Dataset")+
  scale_y_continuous("Overall Survival", 
                     limits = c(0,1), labels = scales::percent)+
  scale_x_continuous("Time (years)")+
  ggtitle("Crizotinib, Overall Survival, Feb 2017 Datacut")

plot3

######################################################
# Figure 1.
######################################################


plot_1a <- km_all %>%
  mutate(dataset = factor(dataset, 
                          levels = c("ALEX", "PROFILE-1014", 
                                     "Flatiron RWE MAIC Weighted", "Flatiron RWE Unweighted"))) %>%
  filter(dataset == "ALEX") %>%
  ggplot(aes(x = time, y = surv, colour = trt))+
  theme_classic()+
  geom_step(linewidth = 0.8)+
  scale_colour_discrete("Treatment")+
  scale_y_continuous("Overall Survival", 
                     limits = c(0,1), labels = scales::percent)+
  scale_x_continuous("Time (years)") + ##, limits = c(0,6))
  ggtitle("ALEX trial, Overall survival")

plot_1b <- km_all %>%
  mutate(dataset = factor(dataset, 
                          levels = c("ALEX", 
                                     "PROFILE-1014", 
                                     "Flatiron RWE MAIC Weighted", 
                                     "Flatiron RWE Unweighted"))) %>%
  filter(trt == "Crizotinib") %>%
  ggplot(aes(x = time, y = surv, colour = dataset))+
  theme_classic()+
  geom_step(linewidth = 0.8)+
  scale_colour_discrete("Dataset")+
  scale_y_continuous("Overall Survival", 
                     limits = c(0,1), labels = scales::percent)+
  scale_x_continuous("Time (years)")+
  ggtitle("Crizotinib, Overall Survival, Feb 2017 Datacut")







######################################################
# Conditional survival
######################################################

km_control_trial_data <- survfit(Surv(time, status) ~ 1, 
                                 data=trial_data %>%
                                    filter(trt == "Crizotinib"))


end_of_trial <- round(max(trial_data$time[trial_data$status == 1]), digits = 1)
prob_times1 <- c(0:1)
prob_times2 <- c(0:1, end_of_trial)
plot(km_trial_data)

trial_cond <- gg_conditional_surv(
  basekm = km_control_trial_data, 
  at = prob_times1,
  main = "Conditional survival",
  xlab = "Days"
  ) +
  labs(color = "Conditional time")
trial_cond$data

km_historic_trial_data <- survfit(Surv(time, status) ~ 1, 
                                 data=historic_trial_data %>%
                                   filter(trt == "Crizotinib"))

historic_trial_cond <- gg_conditional_surv(
  basekm = km_historic_trial_data, 
  at = prob_times2,
  main = "Conditional survival",
  xlab = "Days"
) +
  labs(color = "Conditional time")

km_rwd_unweighted_data <- survfit(Surv(time, status) ~ 1, 
                                  data=rwd_data_unweighted)

rwd_unweighted_cond <- gg_conditional_surv(
  basekm = km_rwd_unweighted_data , 
  at = prob_times2,
  main = "Conditional survival",
  xlab = "Days"
) +
  labs(color = "Conditional time")
#
#rwd_trial_cond$data
km_rwd_weighted_data <- survfit(Surv(time, status) ~ 1, data=rwd_data_weighted, 
                                weights = rwd_data_weighted$weight)

rwd_weighted_cond <- gg_conditional_surv_weight(
  basekm = km_rwd_weighted_data , 
  at = prob_times2,
  main = "Conditional survival",
  xlab = "Days"
) +
  labs(color = "Conditional time")

bind_rows(trial_cond$data %>% mutate(condtime = as.character(condtime),
                                     dataset = "ALEX"),
          historic_trial_cond$data %>% mutate(condtime = as.character(condtime),
                                              dataset = "PROFILE-1014"),
          rwd_unweighted_cond$data %>% mutate(condtime = as.character(condtime),
                                  dataset = "Flatiron RWE Unweighted"),
          rwd_weighted_cond$data %>% mutate(condtime = as.character(condtime),
                                         dataset = "Flatiron RWE MAIC Weighted")) %>%
  mutate(dataset = factor(dataset, 
                          levels = c("ALEX", "PROFILE-1014", 
                                     "Flatiron RWE MAIC Weighted", "Flatiron RWE Unweighted"))) %>%
  filter(condtime %in% c(0, 1.7)) %>%
  mutate(condtime = factor(condtime, levels = c(0,1.7),
                           labels = c("Overall Survival from Time = 0",
                                      "Conditional Survival from Time = 1.7-years \n (Last event time in ALEX trial)"))) %>%
  ggplot(aes(x = timept, y = prob, colour = dataset))+
  theme_classic()+
  geom_step()+
  scale_colour_discrete("Dataset")+
  facet_wrap(~condtime, nrow = 2)+
  scale_y_continuous("Survival", 
                     limits = c(0,1), labels = scales::percent)+
  scale_x_continuous("Time (years)")


######################################################
# Analyse empirical hazards for trial and RWD
######################################################

trial_smooth_bs_control <- bshazard(Surv(time, status) ~ 1, 
                                    data = trial_data  %>% filter(trt == "Crizotinib"))

trial_smooth_control <- trial_smooth_bs_control[c("time", "hazard")] %>%
  as_tibble() %>%
  mutate(data = "ALEX", trt = "Crizotinib")

trial_smooth_bs_active <- bshazard(Surv(time, status) ~ 1, 
                                   data = trial_data  %>% filter(trt == "Alectinib"))

trial_smooth_active <- trial_smooth_bs_active[c("time", "hazard")] %>%
  as_tibble() %>%
  mutate(data = "ALEX", trt = "Alectinib")

historic_trial_bs <- bshazard(Surv(time, status) ~ 1, 
                              data = historic_trial_data)
historic_smooth <- historic_trial_bs[c("time", "hazard")] %>%
  as_tibble() %>%
  mutate(data = "PROFILE-1014", trt = "Crizotinib")

rwd_bs <- bshazard(Surv(time, status) ~ 1, 
                   data = rwd_data)

rwd_smooth <- rwd_bs[c("time", "hazard")] %>%
  as_tibble() %>%
  mutate(data = "Flatiron RWE", trt = "Crizotinib")

bind_rows(trial_smooth_control,
          trial_smooth_active,
          historic_smooth,
          rwd_smooth) %>%
  ggplot(aes(x = time, y = hazard, colour = data, linetype = trt))+
  theme_classic()+
  scale_colour_discrete("Dataset")+
  scale_linetype_discrete("Treatment")+
  geom_line()+
  geom_vline(xintercept = max(trial_data$time)) +
  ylim(c(0,0.5))+
  ylab("Hazard")+
  xlab("Time (years)")+
  ggtitle("Empirical hazards")

bind_rows(trial_smooth_control,
          trial_smooth_active,
          historic_smooth,
          rwd_smooth) %>%
  ggplot(aes(x = time, y = log(hazard), colour = data, linetype = trt))+
  theme_classic()+
  scale_colour_discrete("Dataset")+
  scale_linetype_discrete("Treatment")+
  geom_line()+
  geom_vline(xintercept = max(trial_data$time)) +
 # ylim(c(0,0.5))+
  ylab("Log Hazard")+
  xlab("Time (years)")+
  ggtitle("Empirical log-hazards")

######################################################
#  Create aggregate counts for RWD and historic data.
######################################################

increment_width <- 0.5
left_truncate_ALEX <- max(trial_data$time[trial_data$status == 1])

rwd_aggregate_no_overlap <- create_aggregate_counts(
  data = rwd_data_weighted %>%
    select(time, status, weight) %>%
    mutate(trt = "Crizotinib"),
  increment = increment_width,
  left_truncate = left_truncate_ALEX) %>%
  mutate(data = "Flatiron RWE")

rwd_aggregate_partial <- create_aggregate_counts(
  data = rwd_data_weighted %>%
    select(time, status, weight) %>%
    mutate(trt = "Crizotinib"),
  increment = increment_width,
  left_truncate = left_truncate_ALEX-1) %>%
  mutate(data = "Flatiron RWE")

rwd_aggregate_full <- create_aggregate_counts(
  data = rwd_data_weighted %>%
    select(time, status, weight) %>%
    mutate(trt = "Crizotinib"),
  increment = increment_width,
  left_truncate = 0) %>%
  mutate(data = "Flatiron RWE")

historic_trial_aggregate_no_overlap <- create_aggregate_counts(
  data = historic_trial_data %>%
    mutate(weight = 1),
  increment = increment_width,
  left_truncate = left_truncate_ALEX) %>% 
  mutate(data = "PROFILE-1014 trial") %>%
  rename(dataset = data)  %>%
  mutate(dataset = factor(dataset, 
                          levels = c("ALEX trial", "PROFILE-1014 trial")))

historic_trial_aggregate_partial <- create_aggregate_counts(
  data = historic_trial_data %>%
    mutate(weight = 1),
  increment = increment_width,
  left_truncate = left_truncate_ALEX-1) %>% 
  mutate(data = "PROFILE-1014 trial") %>%
  rename(dataset = data)  %>%
  mutate(dataset = factor(dataset, 
                         levels = c("ALEX trial", "PROFILE-1014 trial")))

historic_trial_aggregate_full <- create_aggregate_counts(
  data = historic_trial_data %>%
    mutate(weight = 1),
  increment = increment_width,
  left_truncate = 0) %>% 
  mutate(data = "PROFILE-1014 trial")   %>%
  rename(dataset = data) %>%
  mutate(dataset = factor(dataset, 
         levels = c("ALEX trial", "PROFILE-1014 trial")))

external_data_no_overlap <- historic_trial_aggregate_no_overlap  %>%
  bind_rows(rwd_aggregate_no_overlap) %>%
  filter(trt == "Crizotinib") %>%
  mutate(dataset = factor(dataset, 
                          levels = c("ALEX trial", "PROFILE-1014 trial", "Flatiron RWE")))

external_data_partial <- historic_trial_aggregate_partial  %>%
  bind_rows(rwd_aggregate_partial) %>%
  filter(trt == "Crizotinib") %>%
  mutate(dataset = factor(dataset, 
                          levels = c("ALEX trial", "PROFILE-1014 trial", "Flatiron RWE")))

external_data_full <- historic_trial_aggregate_full  %>%
  bind_rows(rwd_aggregate_full) %>%
  filter(trt == "Crizotinib") %>%
  mutate(dataset = factor(dataset, 
                          levels = c("ALEX trial", "PROFILE-1014 trial", "Flatiron RWE"))) 

trial_data <- trial_data %>%
  mutate(dataset = "ALEX trial",
         dataset = factor(dataset, levels = c("ALEX trial", "PROFILE-1014 trial", "Flatiron RWE"))) 

trial_data_for_historic <- trial_data %>%
  mutate(dataset = "ALEX trial",
         dataset = factor(dataset, levels = c("ALEX trial", "PROFILE-1014 trial"))) 

#############################################################
# Model trial data alone
#############################################################

model_trial_data <- survextrap(Surv(time,status) ~ trt, 
                       data = trial_data , 
                       backhaz = cetux_bh,
                    #   nonprop = T,
                       df = 3 #,
                     #  add_knots = c(6)#,
                       #fit_method = "opt"
                       )
#model_trial_data$mspline$knots
plot(model_trial_data, tmax = 20)
saveRDS(model_trial_data, "/scratch/klvq491/case_study_nice_ta536/model_trial_data_Feb_2017.rds")
model_trial_data <- readRDS("/scratch/klvq491/case_study_nice_ta536/model_trial_data_Feb_2017.rds")
#model_trial_data

#############################################################
# Model including Historic trial data.
#############################################################

model_historic_no_overlap <- survextrap(Surv(time,status) ~ trt,#+dataset, 
                                     data = trial_data_for_historic, 
                                     external = historic_trial_aggregate_no_overlap,
                                     backhaz = cetux_bh,
                                     #    nonprop = ~trt,
                                     df = 3 #,
                                     # add_knots = c(10),
                                     #   fit_method = "opt"
                                      )
#model_historic_partial$mspline$knots
#View(summary(model_historic_partial2))
#plot(model_historic_no_overlap , tmax = 20)
saveRDS(model_historic_no_overlap, "/scratch/klvq491/case_study_nice_ta536/model_historic_no_overlap_Feb_2017.rds")
model_historic_no_overlap <- readRDS("/scratch/klvq491/case_study_nice_ta536/model_historic_no_overlap_Feb_2017.rds")

# model_historic_partial <- survextrap(Surv(time,status) ~ trt,#+dataset, 
#                                      data = trial_data_for_historic, 
#                                      external = historic_trial_aggregate_partial,
#                                      backhaz = cetux_bh,
#                                  #    nonprop = ~trt,
#                                      df = 3 #,
#                                     # add_knots = c(10),
#                                    #   fit_method = "opt"
#                                     )
# #model_historic_partial$mspline$knots
# #View(summary(model_historic_partial2))
# #plot(model_historic_partial, tmax = 20)
# saveRDS(model_historic_partial, "/scratch/klvq491/case_study_nice_ta536/model_historic_partial_Feb_2017.rds")
# model_historic_partial <- readRDS("/scratch/klvq491/case_study_nice_ta536/model_historic_partial_Feb_2017.rds")

# model_historic_full <- survextrap(Surv(time,status) ~trt,#+dataset, 
#                                      data = trial_data_for_historic, 
#                                      external = historic_trial_aggregate_full,
#                                      backhaz = cetux_bh,
#                                      nonprop = ~trt,
#                                      df = 3 #,
#                                    #  add_knots = c(10),
#                                 #  fit_method = "opt"
#                                   )
# #model_historic_full$mspline$knots
# #View(summary(model_historic_full))
# #plot(model_historic_full, tmax = 20)
# saveRDS(model_historic_full, "/scratch/klvq491/case_study_nice_ta536/model_historic_full_Feb_2017.rds")
# model_historic_full <- readRDS("/scratch/klvq491/case_study_nice_ta536/model_historic_full_Feb_2017.rds")

#############################################################
# Model dataset covariate to include RWE.
#############################################################

model_external_no_overlap <- survextrap(Surv(time,status) ~ trt, #+dataset, 
                                     data = trial_data, 
                                     external = external_data_no_overlap,
                                     backhaz = cetux_bh,
                                     nonprop = T,
                                     df = 3
                                     #add_knots = c(10).
                                     # fit_method = "opt"
)
#model_external_partial$mspline$knots
plot(model_external_no_overlap, tmax = 20)
saveRDS(model_external_no_overlap, "/scratch/klvq491/case_study_nice_ta536/model_external_no_overlap_Feb_2017.rds")
model_external_no_overlap <- readRDS("/scratch/klvq491/case_study_nice_ta536/model_external_no_overlap_Feb_2017.rds")

# # No dataset covariate? Since need overlapping time period for this.
# model_external_partial <- survextrap(Surv(time,status) ~ trt, #+dataset, 
#                            data = trial_data, 
#                            external = external_data_partial,
#                            backhaz = cetux_bh,
#                            nonprop = ~trt,
#                            df = 3
#                            #add_knots = c(10).
#                           # fit_method = "opt"
#                           )
# #model_external_partial$mspline$knots
# saveRDS(model_external_partial, "/scratch/klvq491/case_study_nice_ta536/model_external_partial_Feb_2017.rds")
# model_external_partial <- readRDS("/scratch/klvq491/case_study_nice_ta536/model_external_partial_Feb_2017.rds")
# 
# model_external_full <- survextrap(Surv(time,status) ~ trt, #+dataset, 
#                                         data = trial_data, 
#                                         external = external_data_full,
#                                         backhaz = cetux_bh,
#                                         nonprop =~trt,
#                                         df = 3,
#                                         #add_knots = c(10),
#                                         #fit_method = "opt"
#                                   )
# 
# #model_external_full$mspline$knots
# #plot(model_external_full, tmax = 20)
# saveRDS(model_external_full, "/scratch/klvq491/case_study_nice_ta536/model_external_full_Feb_2017.rds")
# model_external_full <- readRDS("/scratch/klvq491/case_study_nice_ta536/model_external_full_Feb_2017.rds")

new_data_estimate <- expand_grid(trt = c("Alectinib", "Crizotinib"),
                                 dataset =  c("ALEX trial", "PROFILE-1014 trial", "Flatiron RWE")) %>%
#  mutate(dataset = factor(dataset, levels =  c("ALEX trial", "PROFILE-1014 trial", "Flatiron RWE"))) %>%
  filter(trt == "Crizotinib" | dataset == "ALEX trial") 

time_vec <- c(seq(from = 0, to = 5, length.out = 1e2),
              seq(from = 5, to = 20, length.out = 1e2))

haz_mod1 <- hazard(model_trial_data, t = time_vec, newdata = new_data_estimate)  %>%
  mutate(model = "Trial data only") 

haz_mod2 <-  hazard(model_historic_no_overlap, 
                    t = time_vec,
                    newdata = new_data_estimate  %>%
                      mutate(dataset = factor(dataset, 
                                              levels = c("ALEX trial", "PROFILE-1014 trial")))) %>% 
  mutate(model = "Trial + Historic Trial")


haz_mod3 <-  hazard(model_external_no_overlap, 
                    t = time_vec,
                    newdata = new_data_estimate) %>%
  mutate(model = "Trial + Historic Trial + RWE")

# haz_mod4 <-  hazard(model_historic_full, 
#                     t = time_vec,
#                     newdata = new_data_estimate  %>%
#                       mutate(dataset = factor(dataset, 
#                                               levels = c("ALEX trial", "PROFILE-1014 trial")))) %>% 
#   mutate(model = "Trial + Historic Trial,\n full borrowing")

# haz_mod5 <-  hazard(model_external_full, 
#                     t = time_vec,
#                     newdata = new_data_estimate) %>%
#   mutate(model = "Trial + Historic Trial + RWE,\n full borrowing")

#####################################################

surv_mod1 <- survival(model_trial_data, t = time_vec, newdata = new_data_estimate)  %>%
  mutate(model = "Trial data only") 

surv_mod2 <-  survival(model_historic_no_overlap, 
                       t = time_vec,
                       newdata = new_data_estimate  %>%
                         mutate(dataset = factor(dataset, 
                                                 levels = c("ALEX trial", "PROFILE-1014 trial")))) %>% 
  mutate(model = "Trial + Historic Trial")

surv_mod3 <-  survival(model_external_no_overlap, 
                    t = time_vec,
                    newdata = new_data_estimate) %>%
  mutate(model = "Trial + Historic Trial + RWE")

# surv_mod4 <-  survival(model_historic_full, 
#                        t = time_vec,
#                        newdata = new_data_estimate  %>%
#                          mutate(dataset = factor(dataset, 
#                                                  levels = c("ALEX trial", "PROFILE-1014 trial")))) %>% 
#   mutate(model = "Trial + Historic Trial,\n full borrowing")
# surv_mod5 <-  survival(model_external_full, 
#                     t = time_vec,
#                     newdata = new_data_estimate) %>%
#   mutate(model = "Trial + Historic Trial + RWE,\n full borrowing")

#########################################
# Crizotinib arm.
#########################################

haz_mod1 %>%
  bind_rows(haz_mod2) %>%
  bind_rows(haz_mod3) %>%
#  bind_rows(haz_mod4) %>%
#  bind_rows(haz_mod5) %>%
  mutate(model = factor(model, 
                        levels = c("Trial data only", 
                                   "Trial + Historic Trial",
                                   "Trial + Historic Trial + RWE" #,
                                 #  "Trial + Historic Trial,\n full borrowing",
                                #   "Trial + Historic Trial + RWE,\n full borrowing"
                                ))) %>%  
  filter(dataset == "ALEX trial") %>%
  filter(trt == "Crizotinib") %>%
  filter(t > 0) %>%
  ggplot(aes(x = t, y = median, ymin = lower, ymax = upper, colour = model, fill = model))+
  theme_classic()+
  geom_line()+
  geom_ribbon(alpha = 0.4, colour = NA) +
  geom_vline(xintercept = max(trial_data$time)-1, alpha = 0.4)+
  geom_vline(xintercept = max(trial_data$time), alpha = 0.7)+
  scale_colour_discrete("Model")+
  scale_fill_discrete("Model")+
  facet_manual(~model, design = c("#AA#
                                  BBCC"))+
  xlab("Time (years)")+
  ylab("Hazard")+
  ggtitle("Hazard, Crizotinib, Feb 2017 Data cut")


haz_mod1 %>%
  bind_rows(haz_mod2) %>%
  bind_rows(haz_mod3) %>%
  # bind_rows(haz_mod4) %>%
  # bind_rows(haz_mod5) %>%
  mutate(model = factor(model, 
                        levels = c("Trial data only", 
                                   "Trial + Historic Trial",
                                   "Trial + Historic Trial + RWE"#,
                                 #  "Trial + Historic Trial,\n full borrowing",
                                  # "Trial + Historic Trial + RWE,\n full borrowing"
                                 ))) %>%  
  filter(dataset == "ALEX trial") %>%
  filter(trt == "Alectinib") %>%
  filter(t > 0) %>%
  ggplot(aes(x = t, y = median, ymin = lower, ymax = upper, colour = model, fill = model))+
  theme_classic()+
  geom_line()+
  geom_ribbon(alpha = 0.4, colour = NA) +
  #geom_vline(xintercept = max(trial_data$time)-1, alpha = 0.4)+
  geom_vline(xintercept = end_of_trial, alpha = 0.6)+
  scale_colour_discrete("Model")+
  scale_fill_discrete("Model")+
  facet_manual(~model, design = c("#AA#
                                  BBCC"))+
  xlab("Time (years)")+
  ylab("Hazard")+
  ggtitle("Hazard, Alectinib, Feb 2017 Data cut")

surv_mod1 %>%
  bind_rows(surv_mod2) %>%
  bind_rows(surv_mod3) %>%
#  bind_rows(surv_mod4) %>%
#  bind_rows(surv_mod5) %>%
  mutate(model = factor(model, 
                        levels = c("Trial data only", 
                                   "Trial + Historic Trial",
                                   "Trial + Historic Trial + RWE"#,
                                   #  "Trial + Historic Trial,\n full borrowing",
                                   # "Trial + Historic Trial + RWE,\n full borrowing"
                        ))) %>%  
  filter(dataset == "ALEX trial") %>%
  filter(trt == "Crizotinib") %>%
  filter(t > 0) %>%
  ggplot(aes(x = t, y = median, ymin = lower, ymax = upper, colour = model, fill = model))+
  theme_classic()+
   geom_ribbon(alpha = 0.4, colour =  NA)+
  geom_line()+
#  geom_ribbon(alpha = 0.2)+
#  geom_vline(xintercept = max(trial_data$time)-1, alpha = 0.4)+
  geom_vline(xintercept = end_of_trial, alpha = 0.6)+
  scale_colour_discrete("Model")+
  scale_fill_discrete("Model")+
  scale_y_continuous(limits = c(0,1),labels = scales::percent)+
  ylab("Overall Survival")+
  xlab("Time (years)")+
  facet_manual(~model, design = c("#AA#
                                  BBCC"))+
  ggtitle("Overall Survival, Crizotinib, Feb 2017 Data cut")

surv_mod1 %>%
  bind_rows(surv_mod2) %>%
  bind_rows(surv_mod3) %>%
# bind_rows(surv_mod4) %>%
# bind_rows(surv_mod5) %>%
  mutate(model = factor(model, 
                        levels = c("Trial data only", 
                                   "Trial + Historic Trial",
                                   "Trial + Historic Trial + RWE"#,
                                   #  "Trial + Historic Trial,\n full borrowing",
                                   # "Trial + Historic Trial + RWE,\n full borrowing"
                        ))) %>%  
  filter(dataset == "ALEX trial") %>%
  filter(trt == "Alectinib") %>%
  filter(t > 0) %>%
  ggplot(aes(x = t, y = median, ymin = lower, ymax = upper, colour = model, fill = model))+
  theme_classic()+
  geom_line()+
  geom_ribbon(alpha = 0.4, colour = NA)+
 # geom_ribbon(alpha = 0.2)+
#  geom_vline(xintercept = max(trial_data$time)-1, alpha = 0.4)+
  geom_vline(xintercept = end_of_trial, alpha = 0.7)+
  scale_colour_discrete("Model")+
  scale_fill_discrete("Model")+
  scale_y_continuous(limits = c(0,1),labels = scales::percent)+
  ylab("Overall Survival")+
  xlab("Time (years)")+
  facet_manual(~model, design = c("#AA#
                                  BBCC"))+
  ggtitle("Overall Survival, Alectinib, Feb 2017 Data cut")


####################################################
# Plot survextrap models with Kaplan-Meier data.
####################################################

km_data_all <- km_rwd_plot_unweighted %>%
  mutate(dataset = "Flatiron RWE") %>%
  bind_rows(km_trial_plot %>% mutate(dataset = "ALEX trial"),
            km_historic_trial_plot %>% mutate(dataset = "PROFILE-1014 trial")) %>%
  as_tibble() %>%
  mutate(model = "Kaplan-Meier") # %>%
 # filter(dataset == "ALEX") %>%
#  mutate(dataset = "ALEX trial")

summary(surv_mod1$dataset)

surv_mod1 %>%
  bind_rows(surv_mod2) %>%
  #bind_rows(surv_mod3) %>%
  bind_rows(km_data_all) %>%
  select(t, median, trt, dataset, model) %>%
  rename(time = t, surv = median) %>%
  mutate(dataset = as.character(dataset)) %>%
  bind_rows(km_data_all)  %>%
  mutate(model = factor(model, 
                        levels = c("Kaplan-Meier",
                                   "Trial data only", 
                                   "Trial + Historic Trial",
                                   "Trial + Historic Trial + RWE"#,
                                   #  "Trial + Historic Trial,\n full borrowing",
                                   # "Trial + Historic Trial + RWE,\n full borrowing"
                        ))) %>%  
 # filter(dataset == "PROFILE-1014 trial") %>%
#  filter(trt == "Crizotinib") %>%
  filter(time >= 0) %>%
  ggplot(aes(x = time, y = surv, linetype = trt, colour = model, fill = model))+
  theme_classic()+
  geom_line(linewidth = 0.8, alpha = 0.7)+
  #geom_ribbon(alpha = 0.2, colour = "white")+
  geom_vline(xintercept = max(trial_data$time)-1, alpha = 0.4)+
  geom_vline(xintercept = max(trial_data$time), alpha = 0.7)+
  scale_colour_discrete("Model")+
  scale_fill_discrete("Model")+
  scale_linetype_discrete("Treatment")+
  ylab("Survival")+
  scale_y_continuous(limits = c(0,1),labels = scales::percent)+
  xlab("Time (years)")+
  ggtitle("All models, both arms")+
  facet_wrap(~dataset)


surv_mod1 %>%
  bind_rows(surv_mod2) %>%
  bind_rows(surv_mod3) %>%
  bind_rows(km_data_all) %>%
  select(t, median, trt, dataset, model) %>%
  rename(time = t, surv = median) %>%
  mutate(dataset = as.character(dataset)) %>%
  bind_rows(km_data_all)  %>%
  mutate(model = factor(model, 
                        levels = c("Kaplan-Meier",
                                   "Trial data only", 
                                   "Trial + Historic Trial",
                                   "Trial + Historic Trial + RWE"#,
                                   #  "Trial + Historic Trial,\n full borrowing",
                                   # "Trial + Historic Trial + RWE,\n full borrowing"
                        ))) %>%  
  filter(dataset == "ALEX trial") %>%
  #  filter(trt == "Crizotinib") %>%
  filter(time >= 0) %>%
  ggplot(aes(x = time, y = surv, linetype = trt, colour = model, fill = model))+
  theme_classic()+
  geom_line(linewidth = 0.8, alpha = 0.7)+
  #geom_ribbon(alpha = 0.2, colour = "white")+
  geom_vline(xintercept = max(trial_data$time)-1, alpha = 0.4)+
  geom_vline(xintercept = max(trial_data$time), alpha = 0.7)+
  scale_colour_discrete("Model")+
  scale_fill_discrete("Model")+
  scale_linetype_discrete("Treatment")+
  ylab("Survival")+
  scale_y_continuous(limits = c(0,1),labels = scales::percent)+
  xlab("Time (years)")+
  ggtitle("All models, both arms")#+
  #facet_wrap(~dataset)


####################################################
# Hybrid model.
####################################################

trial_data <- readRDS("../Flatiron-RWE-survextrap/Data/trial_IPD_OS_ALEX_Feb_2017_v4.rds")

control <- trial_data %>%
  filter(trt == "Crizotinib") %>%
  select(time, status)

active <- trial_data %>%
  filter(trt == "Alectinib") %>%
  select(time, status)

set.seed(1212)

hybrid_control <- fit_hybrid_model(data = control, cut_point = 1.5)
hybrid_control_survival <- hybrid_model_survival_bootstrap(hybrid_control, 
                                                           t =seq(from = 0, to = 20, length.out = 1e3),
                                                           N = 1e3)

hybrid_active <- fit_hybrid_model(data = active, cut_point = 1.5)
hybrid_active_survival <- hybrid_model_survival_bootstrap(hybrid_active , 
                                                          t =seq(from = 0, to = 20, length.out = 1e3),
                                                             N = 1e3)
sum(active$status[active$time > 1.5])
sum(control$status[control$time > 1.5])

hybrid_active_survival <- hybrid_active_survival %>%
  mutate(treatment = "Alectinib")
  
hybrid_control_survival <- hybrid_control_survival %>%
  mutate(treatment = "Crizotinib")

hybrid_active_survival %>%
  bind_rows(hybrid_control_survival) %>%
  ggplot(aes(x = time, y = value, colour = treatment, fill = treatment))+
  theme_classic()+
  geom_ribbon(aes(x = time, ymin = lower, ymax = upper),
              alpha= 0.3, colour = NA)+
  geom_step()+
  geom_vline(xintercept = 1.5)+
 # geom_vline(xintercept = 0)+
  scale_y_continuous(limits = c(0,1))


###################################################
# RMST results.
####################################################

hybrid_model_rmst_bootstrap(hybrid_control , 
                            t =seq(from = 0, to = 20, length.out = 20),
                            N = 1e3)


hybrid_model_rmst_bootstrap(hybrid_active , 
                                t =seq(from = 0, to = 20, length.out = 20),
                                N = 1e3)


rmst(model_trial_data, 
     newdata = new_data_estimate, 
     t = 20)

rmst(model_historic_no_overlap, 
     newdata = new_data_estimate,
     t = 20)

rmst(model_external_no_overlap, 
     newdata = new_data_estimate,
     t = 20)


# Difference in RMST.







