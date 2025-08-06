
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


######################################################
#  Derive Flatiron cohort for Crizotinib
######################################################

rwd_data_control <- create_flatiron_data(treatment = "Crizotinib",
                                 SoC_date = as.Date("2011-08-26"))

rwd_data_active <- create_flatiron_data(treatment = "Alectinib",
                                         SoC_date = as.Date("2017-11-06"))


######################################################
#  Apply MAIC weighting to Flatiron cohort
######################################################

ALEX_table1_control <-
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

ALEX_table1_active <-
  tibble::tribble(
    ~characteristic, ~mean, ~sd, ~n,
    "age", 56.3, 12.0, 152,
    "sex_male", 0.447, NA, 68,
    "sex_female", 0.553, NA, 84,
    "race_asian", 0.454, NA,  69,
    "race_non_asian", 0.546, NA, 83,
    "ecog_0_1", 0.934, NA, 142,
    "ecog_2", 0.066, NA, 10,
    "smoking_status_ever", 0.395, NA, 60,
    "smoking_status_never", 0.605, NA, 92,
    "brain_mets", 0.421, NA, 64
  )

ALEX_variables_to_match =  c(
  "age",
  "sex_male",
  "ecog_0_1",
  "smoking_status_ever",   
  "brain_mets"
)

# Cohort with MAIC weights.

rwd_data_control_weighted <- create_weighted_cohort(
  cohort_data = rwd_data_control,
  reference_table1 = ALEX_table1_control,
  match_variables = ALEX_variables_to_match,
)[["cohort_data"]]  %>%
  rename(weight = patient_weight)

rwd_data_active_weighted <- create_weighted_cohort(
  cohort_data = rwd_data_active,
  reference_table1 = ALEX_table1_active,
  match_variables = ALEX_variables_to_match,
)[["cohort_data"]]  %>%
  rename(weight = patient_weight)

# Cohort without MAIC weights.

rwd_data_control_unweighted <- rwd_data_control_weighted %>%
  mutate(weight = 1,
         effective_sample_size = n())

rwd_data_active_unweighted <- rwd_data_active_weighted %>%
  mutate(weight = 1,
         effective_sample_size = n())

######################################################
#  Plot RWD KM curves alongside trial data.
######################################################

# trial_data <- readRDS("../Case_study_flatiron/data/trial_IPD_OS_ALEX_Nov_2019.rds")
# trial_data <- trial_data %>%
#   as_tibble() %>%
#   mutate(time = time/12) %>%
#   select(time, status, trt)
# saveRDS(trial_data, "Data/trial_IPD_OS_ALEX_Nov_2019.rds")  
trial_data <- readRDS("Data/trial_IPD_OS_ALEX_Nov_2019.rds")  

# historic_trial_data <- readRDS("../Case_study_flatiron/data/trial_IPD_OS_PROFILE_1014_Nov_2016.rds")
# historic_trial_data <- historic_trial_data  %>%
#   mutate(time = time /12) %>%
#   as_tibble() %>%
#   select(time, status, trt)
# summary(historic_trial_data)
# saveRDS(historic_trial_data, "Data/trial_IPD_OS_PROFILE_1014_Nov_2016.rds")
historic_trial_data <- readRDS("Data/trial_IPD_OS_PROFILE_1014_Nov_2016.rds")
#historic_trial_data

km_trial_data  <- survfit(Surv(time, status) ~ trt, data=trial_data)
km_trial_plot <- ggsurvplot(km_trial_data, data=trial_data)["data.survplot"][[1]] %>%
  select(time, surv, trt) %>%
  mutate(dataset = "ALEX") %>%
  as_tibble()

km_historic_trial_data  <- survfit(Surv(time, status) ~ 1, data=historic_trial_data)
km_historic_trial_plot <- ggsurvplot(km_historic_trial_data, data=historic_trial_data)["data.survplot"][[1]] %>%
  mutate(trt = "Crizotinib",
         dataset = "PROFILE 1014") %>%
  select(time, surv, trt, dataset) %>%
  as_tibble()

km_rwd_control_weighted_data <- survfit(Surv(time, status) ~ 1, 
                                data=rwd_data_control_weighted, 
                                weights = rwd_data_control_weighted$weight)
km_rwd_control_weighted_plot <- ggsurvplot(km_rwd_control_weighted_data, 
                                  data=rwd_data_control_weighted)["data.survplot"][[1]] %>%
  mutate(trt = "Crizotinib",
         dataset = "Flatiron RWE MAIC Weighted") %>%
  select(time, surv, trt, dataset) %>%
  as_tibble()

km_rwd_data_active_weighted <- survfit(Surv(time, status) ~ 1, 
                                        data=rwd_data_active_weighted, 
                                        weights = rwd_data_active_weighted$weight)
km_rwd_active_weighted_plot <- ggsurvplot(km_rwd_data_active_weighted, 
                                           data=rwd_data_active_weighted)["data.survplot"][[1]] %>%
  mutate(trt = "Alectinib",
         dataset = "Flatiron RWE MAIC Weighted") %>%
  select(time, surv, trt, dataset) %>%
  as_tibble()

km_rwd_data_control_unweighted <- survfit(Surv(time, status) ~ 1, 
                                          data=rwd_data_control_unweighted)
km_rwd_control_unweighted_plot <- ggsurvplot(km_rwd_data_control_unweighted, 
                                             data=rwd_data_control_unweighted)["data.survplot"][[1]] %>%
  mutate(trt = "Crizotinib",
         dataset = "Flatiron RWE Unweighted") %>%
  select(time, surv, trt, dataset) %>%
  as_tibble()

km_rwd_data_active_unweighted <- survfit(Surv(time, status) ~ 1, 
                                          data=rwd_data_active_unweighted)
km_rwd_active_unweighted_plot <- ggsurvplot(km_rwd_data_active_unweighted, 
                                             data=rwd_data_active_unweighted)["data.survplot"][[1]] %>%
  mutate(trt = "Alectinib",
         dataset = "Flatiron RWE Unweighted") %>%
  select(time, surv, trt, dataset) %>%
  as_tibble()

km_all <- bind_rows(km_trial_plot,
                    km_historic_trial_plot,
                    km_rwd_control_weighted_plot,
                    km_rwd_active_weighted_plot,
                    km_rwd_control_unweighted_plot,
                    km_rwd_active_unweighted_plot)


plot1 <- km_all %>%
  mutate(dataset = factor(dataset, 
                          levels = c("ALEX", "PROFILE 1014", 
                                     "Flatiron RWE MAIC Weighted", "Flatiron RWE Unweighted"))) %>%
  ggplot(aes(x = time, y = surv, colour = dataset , linetype = trt))+
  theme_classic()+
  geom_step(linewidth = 0.8)+
  scale_colour_discrete("Dataset")+
  scale_linetype_discrete("Treatment")+
  scale_y_continuous("Overall Survival", 
                     limits = c(0,1), labels = scales::percent)+
  scale_x_continuous("Time (years)")##, limits = c(0,6))

plot1

plot2 <- km_all %>%
  mutate(dataset = factor(dataset, 
                          levels = c("ALEX", "PROFILE 1014", 
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
  scale_x_continuous("Time (years)", limits = c(0,14))

plot2
plot3 <- km_all %>%
  mutate(dataset = factor(dataset, 
                          levels = c("ALEX", "PROFILE 1014", 
                                     "Flatiron RWE MAIC Weighted", "Flatiron RWE Unweighted"))) %>%
  filter(trt == "Crizotinib") %>%
  ggplot(aes(x = time, y = surv, colour = dataset))+
  theme_classic()+
  geom_step(linewidth = 0.8)+
  scale_colour_discrete("Dataset")+
  scale_y_continuous("Overall Survival", 
                     limits = c(0,1), labels = scales::percent)+
  scale_x_continuous("Time (years)", limits = c(0,14))

plot3


######################################################
# Analyse empirical hazards for trial and RWD
######################################################

trial_control_bs <- bshazard(Surv(time, status) ~ 1, 
                                    data = trial_data %>%
                               filter(trt == "Crizotinib"))

trial_control_smooth <- trial_control_bs[c("time", "hazard")] %>%
  as_tibble() %>%
  mutate(data = "ALEX", trt = "Crizotinib")

trial_active_bs <- bshazard(Surv(time, status) ~ 1, 
                                   data = trial_data %>%
                              filter(trt == "Alectinib"))

trial_active_smooth <- trial_active_bs[c("time", "hazard")] %>%
  as_tibble() %>%
  mutate(data = "ALEX", trt = "Alectinib")

historic_trial_bs <- bshazard(Surv(time, status) ~ 1, 
                              data = historic_trial_data)
historic_smooth <- historic_trial_bs[c("time", "hazard")] %>%
  as_tibble() %>%
  mutate(data = "PROFILE-1014", trt = "Crizotinib")

rwd_bs_control <- bshazard(Surv(time, status) ~ 1, 
                           data = rwd_data_control)

rwd_control_smooth <- rwd_bs_control[c("time", "hazard")] %>%
  as_tibble() %>%
  mutate(data = "Flatiron RWE", trt = "Crizotinib")

rwd_bs_active <- bshazard(Surv(time, status) ~ 1, 
                          data = rwd_data_active)

rwd_active_smooth <- rwd_bs_active[c("time", "hazard")] %>%
  as_tibble() %>%
  mutate(data = "Flatiron RWE", trt = "Alectinib")
bind_rows(trial_smooth_control,
          trial_smooth_active,
          historic_smooth,
          rwd_control_smooth,
          rwd_active_smooth) %>%
  ggplot(aes(x = time, y = hazard, colour = data, linetype = trt))+
  theme_classic()+
  geom_line()+
  geom_vline(xintercept = max(trial_data$time)) +
  ylim(c(0,0.5))+
  ylab("Hazard")+
  xlab("Time (years)")+
  ggtitle("Empirical log-hazards")

bind_rows(trial_smooth_control,
          trial_smooth_active,
          historic_smooth,
          rwd_control_smooth,
          rwd_active_smooth) %>%
  ggplot(aes(x = time, y = log(hazard), colour = data, linetype = trt))+
  theme_classic()+
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
left_truncate_ALEX <- max(trial_data$time)

rwd_control_aggregate <- create_aggregate_counts(
  data = rwd_data_control_unweighted %>%
    select(time, status, weight) %>%
    mutate(trt = "Crizotinib"),
  increment = increment_width,
  left_truncate = left_truncate_ALEX ) %>%
  mutate(data = "Flatiron RWE")

external_data <- rwd_control_aggregate %>%
  filter(trt == "Crizotinib") %>%
  rename(dataset = data) %>%
  mutate(dataset = factor(dataset, 
                          levels = c("ALEX trial", "PROFILE-1014 trial", "Flatiron RWE")))

trial_data <- trial_data %>%
  mutate(dataset = "ALEX trial",
         dataset = factor(dataset, levels = c("ALEX trial", "PROFILE-1014 trial", "Flatiron RWE")))


#############################################################
# Model trial data alone
#############################################################


test_mod <- survextrap(Surv(time,status) ~ trt, 
                       data = trial_data, 
                       backhaz = cetux_bh,
                       nonprop = T,
                       df = 6,
                       add_knots = c(10))


#############################################################
# Model dataset covariate.
#############################################################

test_mod_ext <- survextrap(Surv(time,status) ~ trt+dataset, 
                           data = trial_data, 
                           external = external_data,
                           backhaz = cetux_bh,
                           nonprop = ~trt,
                           df = 6,
                           add_knots = c(10))

new_data_estimate <- expand_grid(trt = c("Alectinib", "Crizotinib"),
                                 dataset =  c("ALEX trial", "PROFILE-1014 trial", "Flatiron RWE")) %>%
  mutate(dataset = factor(dataset, levels =  c("ALEX trial", "PROFILE-1014 trial", "Flatiron RWE"))) %>%
  filter(trt == "Crizotinib" | dataset == "ALEX trial")

haz_mod1 <- hazard(test_mod, t = seq(from = 0, to = 20, length.out = 1e2), newdata = new_data_estimate)  %>%
  mutate(model = "Trial data only") 

haz_mod2 <-  hazard(test_mod_ext, 
                    t = seq(from = 0, to = 20, length.out = 1e2),
                    newdata = new_data_estimate) %>%
  mutate(model = "Trial + RWE")


haz_mod1 %>%
  bind_rows(haz_mod2) %>%
  filter(dataset == "ALEX trial") %>%
  filter(trt == "Crizotinib") %>%
  filter(t > 0) %>%
  ggplot(aes(x = t, y = median, ymin = lower, ymax = upper, colour = model, fill = model))+
  theme_classic()+
  geom_line()+
  geom_ribbon(alpha = 0.2, colour = "white")


surv_mod1 <- survival(test_mod, t = seq(from = 0, to = 20, length.out = 1e2), 
                      newdata = new_data_estimate) %>%
  mutate(model = "Trial data only") 

surv_mod2 <-  survival(test_mod_ext, 
                       t = seq(from = 0, to = 20, length.out = 1e2),
                       newdata = new_data_estimate) %>%
  mutate(model = "Trial + RWE")

#########################################
# Crizotinib arm.
#########################################

surv_mod1 %>%
  bind_rows(surv_mod2) %>%
  filter(dataset == "ALEX trial") %>%
  filter(trt == "Crizotinib") %>%
  filter(t > 0) %>%
  ggplot(aes(x = t, y = median, ymin = lower, ymax = upper, colour = model, fill = model))+
  theme_classic()+
  geom_line()+
  geom_ribbon(alpha = 0.2, colour = "white")+
  geom_vline(xintercept = max(trial_data$time), alpha = 0.7)+
  scale_colour_discrete("Model")+
  scale_fill_discrete("Model")+
  ylim(c(0,1))+
  ylab("Survival")+
  scale_y_continuous(labels = scales::percent)+
  xlab("Time (years)")+
  ggtitle("Crizotinib")


rmst(test_mod, 
     newdata = new_data_estimate, 
     t = 20)

rmst(test_mod_ext, 
     newdata = new_data_estimate,
     t = 20)


surv_mod1 %>%
  bind_rows(surv_mod2) %>%
  filter(dataset == "ALEX trial") %>%
  filter(trt == "Alectinib") %>%
  filter(t > 0) %>%
  ggplot(aes(x = t, y = median, ymin = lower, ymax = upper, colour = model, fill = model))+
  theme_classic()+
  geom_line()+
  geom_ribbon(alpha = 0.2, colour = "white")+
  geom_vline(xintercept = max(trial_data$time), alpha = 0.7)+
  scale_colour_discrete("Model")+
  scale_fill_discrete("Model")+
  ylim(c(0,1))+
  scale_y_continuous(labels = scales::percent)+
  ylab("Survival")+
  xlab("Time (years)")+
  ggtitle("Alectinib")















