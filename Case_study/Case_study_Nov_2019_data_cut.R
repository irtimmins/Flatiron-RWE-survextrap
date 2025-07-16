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

rwd_data <- create_flatiron_data(treatment = "Crizotinib",
                                 data_cut_off_date = as.Date("2019-11-29"),
                                 censoring_strategy = "administrative_cutoff",
                                 SoC_date = as.Date("2011-08-26"))


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

# Cohort without MAIC weights.
  
rwd_data_unweighted <- rwd_data_weighted %>%
  mutate(weight = 1,
         effective_sample_size = n())

#sum(rwd_data_weighted$weight)^2/sum(rwd_data_weighted$weight^2)
#rwd_data_weighted$effective_sample_size[1:10]

######################################################
#  Plot RWD KM curves alongside trial data.
######################################################

trial_data <- readRDS("Data/trial_IPD_OS_ALEX_Nov_2019.rds")

historic_trial_data <- readRDS("Data/trial_IPD_OS_PROFILE_1014_Nov_2016.rds")

km_trial_data  <- survfit(Surv(time, status) ~ trt, data=trial_data)
km_trial_plot <- ggsurvplot(km_trial_data, data=trial_data)["data.survplot"] [[1]] %>%
  select(time, surv, trt) %>%
  mutate(dataset = "ALEX") %>%
  as_tibble()

km_historic_trial_data  <- survfit(Surv(time, status) ~ 1, data=historic_trial_data)
km_historic_trial_plot <- ggsurvplot(km_historic_trial_data, data=historic_trial_data)["data.survplot"][[1]] %>%
  mutate(trt = "Crizotinib",
         dataset = "PROFILE 1014") %>%
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
                          levels = c("ALEX", "PROFILE 1014", 
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
                          levels = c("ALEX", "PROFILE 1014", 
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
  scale_x_continuous("Time (years)")

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
  scale_x_continuous("Time (years)")+
  ggtitle("Crizotinib, Overall Survival")

plot3


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
left_truncate_ALEX <- max(trial_data$time)

rwd_aggregate_partial <- create_aggregate_counts(
  data = rwd_data_weighted %>%
    select(time, status, weight) %>%
    mutate(trt = "Crizotinib"),
  increment = increment_width,
  left_truncate = left_truncate_ALEX-3) %>%
  mutate(data = "Flatiron RWE")

rwd_aggregate_full <- create_aggregate_counts(
  data = rwd_data_weighted %>%
    select(time, status, weight) %>%
    mutate(trt = "Crizotinib"),
  increment = increment_width,
  left_truncate = 0) %>%
  mutate(data = "Flatiron RWE")

historic_trial_aggregate_partial <- create_aggregate_counts(
  data = historic_trial_data %>%
    mutate(weight = 1),
  increment = increment_width,
  left_truncate = left_truncate_ALEX-3) %>% 
  mutate(data = "PROFILE-1014 trial")

historic_trial_aggregate_full <- create_aggregate_counts(
  data = historic_trial_data %>%
    mutate(weight = 1),
  increment = increment_width,
  left_truncate = 0) %>% 
  mutate(data = "PROFILE-1014 trial")

external_data_partial <- historic_trial_aggregate_partial  %>%
  bind_rows(rwd_aggregate_full) %>%
  filter(trt == "Crizotinib") %>%
  rename(dataset = data) %>%
  mutate(dataset = factor(dataset, 
                          levels = c("ALEX trial", "PROFILE-1014 trial", "Flatiron RWE")))

external_data_full <- historic_trial_aggregate_full  %>%
  bind_rows(rwd_aggregate_full) %>%
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

model_trial_data <- survextrap(Surv(time,status) ~ trt, 
                       data = trial_data, 
                       backhaz = cetux_bh,
                       prior_hsd = p_gamma(2, 5),
                       nonprop = T,
                       df = 10,
                       add_knots = c(10))
model_trial_data$mspline$knots
plot(model_trial_data, tmax = 10)

saveRDS(model_trial_data, "/scratch/klvq491/case_study_nice_ta536/model_trial_data_Nov_2019.rds")
model_trial_data <- readRDS("/scratch/klvq491/case_study_nice_ta536/model_trial_data_Nov_2019.rds")

#############################################################
# Models with dataset covariate.
#############################################################

# No dataset covariate? Since need overlapping time period for this.
model_external_partial <- survextrap(Surv(time,status) ~ trt+dataset, 
                           data = trial_data, 
                           external = external_data_partial,
                           prior_hsd = p_gamma(2, 5),
                           backhaz = cetux_bh,
                           nonprop = ~trt,
                           df = 6,
                           add_knots = c(10))

plot(model_external_partial)

#?survextrap()
saveRDS(model_external_partial, "/scratch/klvq491/case_study_nice_ta536/model_external_partial_borrowing_Nov_2019.rds")

model_external_full <- survextrap(Surv(time,status) ~ trt+dataset, 
                                  data = trial_data, 
                                  external = external_data_full,
                                  prior_hsd = p_gamma(2, 5),
                               #   prior_hrsd = p_gamma(2, 5),
                                  backhaz = cetux_bh,
                                  nonprop = ~trt,
                                  df = 6,
                                  add_knots = c(10))

saveRDS(model_external_full, "/scratch/klvq491/case_study_nice_ta536/model_external_full_borrowing_Nov_2019.rds")
model_external_full <- readRDS("/scratch/klvq491/case_study_nice_ta536/model_external_full_borrowing_Nov_2019.rds")

new_data_estimate <- expand_grid(trt = c("Alectinib", "Crizotinib"),
                                 dataset =  c("ALEX trial", "PROFILE-1014 trial", "Flatiron RWE")) %>%
  mutate(dataset = factor(dataset, levels =  c("ALEX trial", "PROFILE-1014 trial", "Flatiron RWE"))) %>%
  filter(trt == "Crizotinib" | dataset == "ALEX trial") 


time_vec <- c(seq(from = 0, to = 5, length.out = 1e2),
              seq(from = 5, to = 20, length.out = 1e2))

haz_mod1 <- hazard(model_trial_data, t = time_vec, newdata = new_data_estimate)  %>%
  mutate(model = "Trial data only") 

haz_mod2 <-  hazard(model_external_partial, 
                    t = time_vec,
                    newdata = new_data_estimate) %>%
  mutate(model = "Trial + RWE, partial borrowing")

haz_mod3 <-  hazard(model_external_full, 
                    t = time_vec,
                    newdata = new_data_estimate) %>%
  mutate(model = "Trial + RWE, full borrowing")


haz_mod1 %>%
  bind_rows(haz_mod2) %>%
  bind_rows(haz_mod3) %>%
  filter(dataset == "ALEX trial") %>%
  filter(trt == "Crizotinib") %>%
  filter(t > 0) %>%
  ggplot(aes(x = t, y = median, ymin = lower, ymax = upper, colour = model, fill = model))+
  theme_classic()+
  geom_line()+
  geom_ribbon(alpha = 0.2, colour = "white")


surv_mod1 <- survival(model_trial_data, t = time_vec, 
                      newdata = new_data_estimate) %>%
  mutate(model = "Trial data only") 

surv_mod2 <-  survival(model_external_partial, 
                       t = time_vec,
                       newdata = new_data_estimate) %>%
  mutate(model = "Trial + RWE, partial borrowing")

surv_mod3 <-  survival(model_external_full, 
                       t = time_vec,
                       newdata = new_data_estimate) %>%
  mutate(model = "Trial + RWE, full borrowing")

#########################################
# Crizotinib arm.
#########################################

surv_mod1 %>%
  bind_rows(surv_mod2) %>%
  bind_rows(surv_mod3) %>%
  mutate(model = factor(model, 
                        levels = c("Trial data only", "Trial + RWE, partial borrowing","Trial + RWE, full borrowing" ))) %>%
  filter(dataset == "ALEX trial") %>%
  filter(trt == "Crizotinib") %>%
  filter(t > 0) %>%
  ggplot(aes(x = t, y = median, ymin = lower, ymax = upper))+
  theme_classic()+
  geom_line()+
  geom_ribbon(alpha = 0.1)+
  geom_vline(xintercept = max(trial_data$time), alpha = 0.7)+
  # scale_colour_discrete("Model")+
  # scale_fill_discrete("Model")+
  ylab("Overall Survival")+
  scale_y_continuous(limits = c(0,1), labels = scales::percent)+
  xlab("Time (years)")+
  facet_wrap(~model, ncol = 1)+
  ggtitle("Overall Survival, Crizotinib, Nov 2019 Data cut")


haz_mod1 %>%
  bind_rows(haz_mod2) %>%
  bind_rows(haz_mod3) %>%
  mutate(model = factor(model, 
                        levels = c("Trial data only", "Trial + RWE, partial borrowing","Trial + RWE, full borrowing" ))) %>%
  filter(dataset == "ALEX trial") %>%
  filter(trt == "Crizotinib") %>%
  filter(t > 0) %>%
  ggplot(aes(x = t, y = median, ymin = lower, ymax = upper, colour = model, fill = model))+
  theme_classic()+
  geom_line()+
  geom_ribbon(alpha = 0.1) +
  geom_vline(xintercept = max(trial_data$time), alpha = 0.7)+
#  facet_wrap(~model, ncol = 1)+
  ylab("Hazard")+
  ggtitle("Hazard, Crizotinib, Nov 2019 Data cut")

surv_mod1 %>%
  bind_rows(surv_mod2) %>%
  bind_rows(surv_mod3) %>%
  mutate(model = factor(model, 
                        levels = c("Trial data only", "Trial + RWE, partial borrowing","Trial + RWE, full borrowing" ))) %>%
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
  scale_y_continuous(limits = c(0,1),labels = scales::percent)+
  ylab("Overall Survival")+
  xlab("Time (years)")+
#  facet_wrap(~model, ncol = 1)+
  ggtitle("Overall Survival, Criztonib, Nov 2019 Data cut")

surv_mod1 %>%
  bind_rows(surv_mod2) %>%
  bind_rows(surv_mod3) %>%
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
  scale_y_continuous(limits = c(0,1),labels = scales::percent)+
  ylab("Overall Survival")+
  xlab("Time (years)")+
 # facet_wrap(~model, ncol = 1)+
  ggtitle("Overall Survival, Alectinib, Nov 2019 Data cut")


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
summary(as.factor(km_data_all$dataset))
summary(surv_mod1$dataset)

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
                                     "Trial + RWE, partial borrowing",
                                     "Trial + RWE, full borrowing")))  %>%
 # filter(dataset == "PROFILE-1014 trial") %>%
#  filter(trt == "Crizotinib") %>%
  filter(time >= 0) %>%
  ggplot(aes(x = time, y = surv, linetype = trt, colour = model, fill = model))+
  theme_classic()+
  geom_line(linewidth = 0.8, alpha = 0.7)+
  #geom_ribbon(alpha = 0.2, colour = "white")+
  geom_vline(xintercept = max(trial_data$time), alpha = 0.7)+
  scale_colour_discrete("Model")+
  scale_fill_discrete("Model")+
  scale_linetype_discrete("Treatment")+
  ylab("Survival")+
  scale_y_continuous(limits = c(0,1),labels = scales::percent)+
  xlab("Time (years)")+
  ggtitle("All models, both arms")+
  facet_wrap(~dataset)




####################################################
# RMST results.
####################################################


rmst(model_trial_data, 
     newdata = new_data_estimate, 
     t = 20)

rmst(model_external_no_overlap, 
     newdata = new_data_estimate,
     t = 20)

rmst(model_external_with_overlap, 
     newdata = new_data_estimate,
     t = 20)


# Difference in RMST.







