
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


plot_2a <- km_all %>%
  mutate(dataset = factor(dataset, 
                          levels = c("ALEX", "PROFILE-1014", 
                                     "Flatiron RWE MAIC Weighted", "Flatiron RWE Unweighted"))) %>%
  filter(dataset == "ALEX") %>%
  ggplot(aes(x = time, y = surv, colour = trt))+
  theme_classic()+
  theme(legend.position.inside = c(0.66,0.93),
        legend.position = "inside",
        legend.text=element_text(size=9),
        legend.key.spacing.y = unit(-8, "pt"),
        legend.box.spacing = unit(-10, "pt"),
        legend.title = element_text(size=12,
                                    margin = margin(l = 0, r = 0, b = 0, t = 0))) +
  geom_step(linewidth = 0.8)+
  scale_colour_manual("Treatment", values = c("#4393c3","#d6604d"))+
  scale_y_continuous("Overall Survival", 
                     limits = c(0,1), labels = scales::percent)+
  scale_x_continuous("Time (years)",breaks = 0:5, limits = c(0,3))
  #ggtitle("ALEX trial, Overall survival")

plot_2a


plot_2b <- km_all %>%
  mutate(dataset = factor(dataset, 
                          levels = c("ALEX", 
                                     "PROFILE-1014", 
                                     "Flatiron RWE MAIC Weighted", 
                                     "Flatiron RWE Unweighted"),
                          labels = c("ALEX", 
                                     "PROFILE-1014", 
                                     "Flatiron RWE MAIC weighted", 
                                     "Flatiron RWE unweighted"))) %>%
  filter(trt == "Crizotinib") %>%
  ggplot(aes(x = time, y = surv, colour = dataset))+
  theme_classic()+
  theme(legend.position.inside = c(0.66,0.93),
        legend.position = "inside",
        legend.text=element_text(size=9),
        legend.spacing.y = unit(-6, "pt"),
        legend.spacing.x = unit(0, "pt"),
        legend.key.spacing.y = unit(-8, "pt"),
        legend.box.spacing = unit(-10, "pt"),
        legend.title = element_text(size=12,
                                    margin = margin(l = 0, r = 0, b = 0, t = 0))) +
  geom_step(linewidth = 0.8)+
  scale_colour_manual("Dataset", 
                    values = c("#d6604d", "#7CAE00", "#00BFC4", "#C77CFF" ))+
  scale_y_continuous("Overall Survival", 
                     limits = c(0,1), labels = scales::percent)+
  scale_x_continuous("Time (years)",  breaks = 0:5, limits = c(0,5.5))
  #ggtitle("Crizotinib, Overall Survival, Feb 2017 Datacut")
plot_2b

figure_2 <- plot_grid(plot_2a+
            theme(  plot.title = element_text(hjust = -0.1,
                                              size=14, face="bold"))+
            labs(title = "(a)"),
          plot_2b+
            theme(  plot.title = element_text(hjust = -0.1,
                                              size=14, face="bold"))+
            labs(title = "(b)"),
          align = "h",
          rel_widths  =c(0.5,0.5),
          ncol = 2)
figure_2 

tiff(file = "Figures/Figure_2.tiff",   
     width = 7.6, 
     height = 3.7,
     units = 'in',  
     res = 300, 
     compression = "lzw")
print(figure_2)
dev.off()


######################################################
# MAIC wieghts
######################################################

sup_figure1 <- rwd_data_weighted %>%
  ggplot()+
  geom_histogram(aes(x = weight), colour = "gray50", binwidth = 0.01)+
  #  geom_histogram(aes(x = weight, colour = race_asian, fill = race_asian), binwidth = 0.01)+
  theme_classic()+
  scale_y_continuous("Frequency")+
  scale_x_continuous("Weight")


tiff(file = "Figures/Sup_figure_1.tiff",   
     width = 5.8, 
     height = 2.4,
     units = 'in',  
     res = 300, 
     compression = "lzw")
print(sup_figure1)
dev.off()



######################################################
# Conditional survival
######################################################

km_control_trial_data <- survfit(Surv(time, status) ~ 1, 
                                 data=trial_data %>%
                                    filter(trt == "Crizotinib"))


end_of_trial <- round(max(trial_data$time[trial_data$status == 1]), digits = 1)
prob_times1 <- c(0, end_of_trial)
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
  at = prob_times1,
  main = "Conditional survival",
  xlab = "Days"
) +
  labs(color = "Conditional time")

km_rwd_unweighted_data <- survfit(Surv(time, status) ~ 1, 
                                  data=rwd_data_unweighted)

rwd_unweighted_cond <- gg_conditional_surv(
  basekm = km_rwd_unweighted_data , 
  at = prob_times1,
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
  at = prob_times1,
  main = "Conditional survival",
  xlab = "Days"
) +
  labs(color = "Conditional time")

sup_figure2 <- bind_rows(trial_cond$data %>% mutate(condtime = as.character(condtime),
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
  mutate(condtime = factor(condtime, levels = c(0,1.8),
                           labels = c("Overall Survival from time zero",
                                      "Conditional Survival from Time = 1.8-years \n (Last event time in ALEX trial)"))) %>%
  ggplot(aes(x = timept, y = prob, colour = dataset))+
  theme_classic()+
  geom_step()+
  scale_colour_discrete("Dataset")+
  facet_wrap(~condtime, nrow = 2)+
  scale_y_continuous("Survival", 
                     limits = c(0,1), labels = scales::percent)+
  scale_x_continuous("Time (years)")


tiff(file = "Figures/Sup_figure_2.tiff",   
     width = 6.7, 
     height = 5.2,
     units = 'in',  
     res = 300, 
     compression = "lzw")
print(sup_figure2)
dev.off()


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

sup_figure3 <- bind_rows(trial_smooth_control,
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
  xlab("Time (years)")#+
  #ggtitle("Empirical hazards")

tiff(file = "Figures/Sup_figure_3.tiff",   
     width = 5.8, 
     height = 3.4,
     units = 'in',  
     res = 300, 
     compression = "lzw")
print(sup_figure3)
dev.off()


######################################################
#  Create aggregate counts for RWD and historic data.
######################################################

increment_width <- 0.5
left_truncate_ALEX <- end_of_trial

rwd_aggregate_no_overlap <- create_aggregate_counts(
  data = rwd_data_weighted %>%
    select(time, status, weight) %>%
    mutate(trt = "Crizotinib"),
  increment = increment_width,
  left_truncate = left_truncate_ALEX) %>%
  mutate(data = "Flatiron RWE")

# rwd_aggregate_partial <- create_aggregate_counts(
#   data = rwd_data_weighted %>%
#     select(time, status, weight) %>%
#     mutate(trt = "Crizotinib"),
#   increment = increment_width,
#   left_truncate = left_truncate_ALEX-1) %>%
#   mutate(data = "Flatiron RWE")
# 
# rwd_aggregate_full <- create_aggregate_counts(
#   data = rwd_data_weighted %>%
#     select(time, status, weight) %>%
#     mutate(trt = "Crizotinib"),
#   increment = increment_width,
#   left_truncate = 0) %>%
#   mutate(data = "Flatiron RWE")

historic_trial_aggregate_no_overlap <- create_aggregate_counts(
  data = historic_trial_data %>%
    mutate(weight = 1),
  increment = increment_width,
  left_truncate = left_truncate_ALEX) %>% 
  mutate(data = "PROFILE-1014 trial") %>%
  rename(dataset = data)  %>%
  mutate(dataset = factor(dataset, 
                          levels = c("ALEX trial", "PROFILE-1014 trial")))

# historic_trial_aggregate_partial <- create_aggregate_counts(
#   data = historic_trial_data %>%
#     mutate(weight = 1),
#   increment = increment_width,
#   left_truncate = left_truncate_ALEX-1) %>% 
#   mutate(data = "PROFILE-1014 trial") %>%
#   rename(dataset = data)  %>%
#   mutate(dataset = factor(dataset, 
#                          levels = c("ALEX trial", "PROFILE-1014 trial")))

# historic_trial_aggregate_full <- create_aggregate_counts(
#   data = historic_trial_data %>%
#     mutate(weight = 1),
#   increment = increment_width,
#   left_truncate = 0) %>% 
#   mutate(data = "PROFILE-1014 trial")   %>%
#   rename(dataset = data) %>%
#   mutate(dataset = factor(dataset, 
#          levels = c("ALEX trial", "PROFILE-1014 trial")))

external_data_no_overlap <- historic_trial_aggregate_no_overlap  %>%
  bind_rows(rwd_aggregate_no_overlap) %>%
  filter(trt == "Crizotinib") %>%
  mutate(dataset = factor(dataset, 
                          levels = c("ALEX trial", "PROFILE-1014 trial", "Flatiron RWE")))

# external_data_partial <- historic_trial_aggregate_partial  %>%
#   bind_rows(rwd_aggregate_partial) %>%
#   filter(trt == "Crizotinib") %>%
#   mutate(dataset = factor(dataset, 
#                           levels = c("ALEX trial", "PROFILE-1014 trial", "Flatiron RWE")))
# 
# external_data_full <- historic_trial_aggregate_full  %>%
#   bind_rows(rwd_aggregate_full) %>%
#   filter(trt == "Crizotinib") %>%
#   mutate(dataset = factor(dataset, 
#                           levels = c("ALEX trial", "PROFILE-1014 trial", "Flatiron RWE"))) 

trial_data <- trial_data %>%
  mutate(dataset = "ALEX trial",
         dataset = factor(dataset, levels = c("ALEX trial", "PROFILE-1014 trial", "Flatiron RWE"))) 

# trial_data_for_historic <- trial_data %>%
#   mutate(dataset = "ALEX trial",
#          dataset = factor(dataset, levels = c("ALEX trial", "PROFILE-1014 trial"))) 

saveRDS(trial_data, "Data/trial_data.rds")
saveRDS(historic_trial_aggregate_no_overlap, "Data/historic_trial_aggregate_no_overlap.rds")
saveRDS(external_data_no_overlap, "Data/external_data_no_overlap.rds")


#############################################################
# Model trial data alone
#############################################################

model_trial_data <- survextrap(Surv(time,status) ~ trt, 
                       data = trial_data , 
                       # backhaz = cetux_bh,
                      #  nonprop = T,
                        df = 3 ,
                       prior_hsd = p_gamma(2, 5),
                     #  prior_hrsd = p_gamma(2, 20)
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
                                     # backhaz = cetux_bh,
                                    # nonprop = T,
                                     prior_hsd = p_gamma(2, 5),
                                   #  prior_hrsd = p_gamma(2, 20),
                                     df = 3 #,
                                     # add_knots = c(10),
                                     #   fit_method = "opt"
                                      )
#model_historic_partial$mspline$knots
#View(summary(model_historic_partial2))
plot(model_historic_no_overlap , tmax = 20)
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
                                     #backhaz = cetux_bh,
                                     # nonprop = T,
                                     df = 3,
                                     prior_hsd = p_gamma(2, 5),
                                     #prior_hrsd = p_gamma(2, 20)
                                  
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

# #########################################
# # Crizotinib arm.
# #########################################
# 
# haz_mod1 %>%
#   bind_rows(haz_mod2) %>%
#   bind_rows(haz_mod3) %>%
# #  bind_rows(haz_mod4) %>%
# #  bind_rows(haz_mod5) %>%
#   mutate(model = factor(model, 
#                         levels = c("Trial data only", 
#                                    "Trial + Historic Trial",
#                                    "Trial + Historic Trial + RWE" #,
#                                  #  "Trial + Historic Trial,\n full borrowing",
#                                 #   "Trial + Historic Trial + RWE,\n full borrowing"
#                                 ))) %>%  
#   filter(dataset == "ALEX trial") %>%
#   filter(trt == "Crizotinib") %>%
#   filter(t > 0) %>%
#   ggplot(aes(x = t, y = median, ymin = lower, ymax = upper, colour = model, fill = model))+
#   theme_classic()+
#   geom_line()+
#   geom_ribbon(alpha = 0.4, colour = NA) +
#   geom_vline(xintercept = max(trial_data$time)-1, alpha = 0.4)+
#   geom_vline(xintercept = max(trial_data$time), alpha = 0.7)+
#   scale_colour_discrete("Model")+
#   scale_fill_discrete("Model")+
#   facet_manual(~model, design = c("#AA#
#                                   BBCC"))+
#   xlab("Time (years)")+
#   ylab("Hazard")+
#   ggtitle("Hazard, Crizotinib, Feb 2017 Data cut")
# 
# 
# haz_mod1 %>%
#   bind_rows(haz_mod2) %>%
#   bind_rows(haz_mod3) %>%
#   # bind_rows(haz_mod4) %>%
#   # bind_rows(haz_mod5) %>%
#   mutate(model = factor(model, 
#                         levels = c("Trial data only", 
#                                    "Trial + Historic Trial",
#                                    "Trial + Historic Trial + RWE"#,
#                                  #  "Trial + Historic Trial,\n full borrowing",
#                                   # "Trial + Historic Trial + RWE,\n full borrowing"
#                                  ))) %>%  
#   filter(dataset == "ALEX trial") %>%
#   filter(trt == "Alectinib") %>%
#   filter(t > 0) %>%
#   ggplot(aes(x = t, y = median, ymin = lower, ymax = upper, colour = model, fill = model))+
#   theme_classic()+
#   geom_line()+
#   geom_ribbon(alpha = 0.4, colour = NA) +
#   #geom_vline(xintercept = max(trial_data$time)-1, alpha = 0.4)+
#   geom_vline(xintercept = end_of_trial, alpha = 0.6)+
#   scale_colour_discrete("Model")+
#   scale_fill_discrete("Model")+
#   facet_manual(~model, design = c("#AA#
#                                   BBCC"))+
#   xlab("Time (years)")+
#   ylab("Hazard")+
#   ggtitle("Hazard, Alectinib, Feb 2017 Data cut")
# 
# surv_mod1 %>%
#   bind_rows(surv_mod2) %>%
#   bind_rows(surv_mod3) %>%
# #  bind_rows(surv_mod4) %>%
# #  bind_rows(surv_mod5) %>%
#   mutate(model = factor(model, 
#                         levels = c("Trial data only", 
#                                    "Trial + Historic Trial",
#                                    "Trial + Historic Trial + RWE"#,
#                                    #  "Trial + Historic Trial,\n full borrowing",
#                                    # "Trial + Historic Trial + RWE,\n full borrowing"
#                         ))) %>%  
#   filter(dataset == "ALEX trial") %>%
#   filter(trt == "Crizotinib") %>%
#   filter(t > 0) %>%
#   ggplot(aes(x = t, y = median, ymin = lower, ymax = upper, colour = model, fill = model))+
#   theme_classic()+
#    geom_ribbon(alpha = 0.4, colour =  NA)+
#   geom_line()+
# #  geom_ribbon(alpha = 0.2)+
# #  geom_vline(xintercept = max(trial_data$time)-1, alpha = 0.4)+
#   geom_vline(xintercept = end_of_trial, alpha = 0.6)+
#   scale_colour_discrete("Model")+
#   scale_fill_discrete("Model")+
#   scale_y_continuous(limits = c(0,1),labels = scales::percent)+
#   ylab("Overall Survival")+
#   xlab("Time (years)")+
#   facet_manual(~model, design = c("#AA#
#                                   BBCC"))+
#   ggtitle("Overall Survival, Crizotinib, Feb 2017 Data cut")
# 
# surv_mod1 %>%
#   bind_rows(surv_mod2) %>%
#   bind_rows(surv_mod3) %>%
# # bind_rows(surv_mod4) %>%
# # bind_rows(surv_mod5) %>%
#   mutate(model = factor(model, 
#                         levels = c("Trial data only", 
#                                    "Trial + Historic Trial",
#                                    "Trial + Historic Trial + RWE"#,
#                                    #  "Trial + Historic Trial,\n full borrowing",
#                                    # "Trial + Historic Trial + RWE,\n full borrowing"
#                         ))) %>%  
#   filter(dataset == "ALEX trial") %>%
#   filter(trt == "Alectinib") %>%
#   filter(t > 0) %>%
#   ggplot(aes(x = t, y = median, ymin = lower, ymax = upper, colour = model, fill = model))+
#   theme_classic()+
#   geom_line()+
#   geom_ribbon(alpha = 0.4, colour = NA)+
#  # geom_ribbon(alpha = 0.2)+
# #  geom_vline(xintercept = max(trial_data$time)-1, alpha = 0.4)+
#   geom_vline(xintercept = end_of_trial, alpha = 0.7)+
#   scale_colour_discrete("Model")+
#   scale_fill_discrete("Model")+
#   scale_y_continuous(limits = c(0,1),labels = scales::percent)+
#   ylab("Overall Survival")+
#   xlab("Time (years)")+
#   facet_manual(~model, design = c("#AA#
#                                   BBCC"))+
#   ggtitle("Overall Survival, Alectinib, Feb 2017 Data cut")
# 

# 
# summary(surv_mod1$dataset)
# 
# surv_mod1 %>%
#   bind_rows(surv_mod2) %>%
#   #bind_rows(surv_mod3) %>%
#   bind_rows(km_data_all) %>%
#   select(t, median, trt, dataset, model) %>%
#   rename(time = t, surv = median) %>%
#   mutate(dataset = as.character(dataset)) %>%
#   bind_rows(km_data_all)  %>%
#   mutate(model = factor(model, 
#                         levels = c("Kaplan-Meier",
#                                    "Trial data only", 
#                                    "Trial + Historic Trial",
#                                    "Trial + Historic Trial + RWE"#,
#                                    #  "Trial + Historic Trial,\n full borrowing",
#                                    # "Trial + Historic Trial + RWE,\n full borrowing"
#                         ))) %>%  
#  # filter(dataset == "PROFILE-1014 trial") %>%
# #  filter(trt == "Crizotinib") %>%
#   filter(time >= 0) %>%
#   ggplot(aes(x = time, y = surv, linetype = trt, colour = model, fill = model))+
#   theme_classic()+
#   geom_line(linewidth = 0.8, alpha = 0.7)+
#   #geom_ribbon(alpha = 0.2, colour = "white")+
#   geom_vline(xintercept = max(trial_data$time)-1, alpha = 0.4)+
#   geom_vline(xintercept = max(trial_data$time), alpha = 0.7)+
#   scale_colour_discrete("Model")+
#   scale_fill_discrete("Model")+
#   scale_linetype_discrete("Treatment")+
#   ylab("Survival")+
#   scale_y_continuous(limits = c(0,1),labels = scales::percent)+
#   xlab("Time (years)")+
#   ggtitle("All models, both arms")+
#   facet_wrap(~dataset)
# 
# 
# surv_mod1 %>%
#   bind_rows(surv_mod2) %>%
#   bind_rows(surv_mod3) %>%
#   bind_rows(km_data_all) %>%
#   select(t, median, trt, dataset, model) %>%
#   rename(time = t, surv = median) %>%
#   mutate(dataset = as.character(dataset)) %>%
#   bind_rows(km_data_all)  %>%
#   mutate(model = factor(model, 
#                         levels = c("Kaplan-Meier",
#                                    "Trial data only", 
#                                    "Trial + Historic Trial",
#                                    "Trial + Historic Trial + RWE"#,
#                                    #  "Trial + Historic Trial,\n full borrowing",
#                                    # "Trial + Historic Trial + RWE,\n full borrowing"
#                         ))) %>%  
#   filter(dataset == "ALEX trial") %>%
#   #  filter(trt == "Crizotinib") %>%
#   filter(time >= 0) %>%
#   ggplot(aes(x = time, y = surv, linetype = trt, colour = model, fill = model))+
#   theme_classic()+
#   geom_line(linewidth = 0.8, alpha = 0.7)+
#   #geom_ribbon(alpha = 0.2, colour = "white")+
#   geom_vline(xintercept = max(trial_data$time)-1, alpha = 0.4)+
#   geom_vline(xintercept = max(trial_data$time), alpha = 0.7)+
#   scale_colour_discrete("Model")+
#   scale_fill_discrete("Model")+
#   scale_linetype_discrete("Treatment")+
#   ylab("Survival")+
#   scale_y_continuous(limits = c(0,1),labels = scales::percent)+
#   xlab("Time (years)")+
#   ggtitle("All models, both arms")#+
#   #facet_wrap(~dataset)
# 

####################################################
# Create survival and hazard figure for manuscript.
####################################################

# Prepare Kaplan-Meier data to add to survival plot.

km_data_all <- km_rwd_plot_unweighted %>%
  mutate(dataset = "Flatiron RWE") %>%
  bind_rows(km_trial_plot %>% mutate(dataset = "ALEX trial"),
            km_historic_trial_plot %>% mutate(dataset = "PROFILE-1014 trial")) %>%
  as_tibble() %>%
  mutate(model = "Kaplan-Meier") # %>%
# filter(dataset == "ALEX") %>%
#  mutate(dataset = "ALEX trial")


figure3_a <- surv_mod1 %>%
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
                        ),
                        labels = c("ALEX trial data only", 
                                   "ALEX + PROFILE-1014",
                                   "ALEX + PROFILE-1014 + Flatiron RWE"#,
                                   #  "Trial + Historic Trial,\n full borrowing",
                                   # "Trial + Historic Trial + RWE,\n full borrowing"
                        ))) %>%  
  filter(dataset == "ALEX trial") %>%
 # filter(trt == "Crizotinib") %>%
  filter(t > 0) %>%
  ggplot()+
  theme_classic()+
  theme(legend.position = "bottom",
        legend.title = element_text(margin = margin(l = unit(4.0, 'cm'), r = unit(12.0, 'cm'))),
        legend.key.spacing.x =  unit(0.3, 'cm'),
        legend.box.spacing = unit(0, "inch"),
        plot.margin = unit(c(0.1,0,0,0), "cm"),
        legend.spacing = unit(c(0,0,0,0), "cm"),
        legend.spacing.x = unit(0, "mm"),
        legend.spacing.y = unit(0, "mm"))+
  geom_ribbon(aes(x = t, y = median, ymin = lower, ymax = upper, colour = trt, fill = trt),
              alpha = 0.25, colour =  NA)+
  geom_line(aes(x = t, y = median, colour = trt))+
  geom_line(data = km_data_all %>% 
              filter(dataset == "ALEX trial") %>% 
              select(-c(dataset,model)),
            aes(x = time, y = surv, colour = trt) )+
  #geom_ribbon(alpha = 0.2)+
  #geom_vline(xintercept = max(trial_data$time)-1, alpha = 0.4)+
  #geom_vline(xintercept = end_of_trial, alpha = 0.6)+
  scale_colour_discrete("Treatment")+
  scale_fill_discrete("Treatment")+
  scale_y_continuous(limits = c(0,1),labels = scales::percent)+
  ylab("Overall Survival")+
  xlab("Time (years)")+
  facet_wrap(~model, ncol = 1)

figure3_a

figure3_b <- haz_mod1 %>%
  bind_rows(haz_mod2) %>%
  bind_rows(haz_mod3) %>%
  #  bind_rows(surv_mod4) %>%
  #  bind_rows(surv_mod5) %>%
  mutate(model = factor(model, 
                        levels = c("Trial data only", 
                                   "Trial + Historic Trial",
                                   "Trial + Historic Trial + RWE"#,
                                   #  "Trial + Historic Trial,\n full borrowing",
                                   # "Trial + Historic Trial + RWE,\n full borrowing"
                        ),
                        labels = c("ALEX trial data only", 
                                   "ALEX + PROFILE-1014",
                                   "ALEX + PROFILE-1014 + Flatiron RWE"#,
                                   #  "Trial + Historic Trial,\n full borrowing",
                                   # "Trial + Historic Trial + RWE,\n full borrowing"
                        ))) %>%  
  filter(dataset == "ALEX trial") %>%
  # filter(trt == "Crizotinib") %>%
  filter(t > 0) %>%
  ggplot()+
  theme_classic()+
  theme(legend.box = "horizontal",
        legend.box.spacing = unit(0, "inch"),
        legend.title = element_text(margin = margin(b = 0.1)),
        plot.margin = unit(c(0.1,0,0,0), "cm"),
        legend.spacing = unit(c(0,0,0,0), "cm"),
        legend.spacing.x = unit(0, "mm"),
        legend.spacing.y = unit(0, "mm"))+
  geom_ribbon(aes(x = t, y = median, ymin = lower, ymax = upper, colour = trt, fill = trt),
              alpha = 0.25, colour =  NA)+
  geom_line(aes(x = t, y = median, colour = trt))+
  #geom_ribbon(alpha = 0.2)+
  #geom_vline(xintercept = max(trial_data$time), alpha = 0.4)+
  #geom_vline(xintercept = end_of_trial, alpha = 0.6)+
  scale_colour_discrete("Treatment")+
  scale_fill_discrete("Treatment")+
  scale_y_continuous("Hazard", limits = c(0, 0.4))+
  xlab("Time (years)")+
  facet_wrap(~model, ncol = 1)

figure3_b
get_legend_35 <- function(plot, legend_number = 1) {
  # find all legend candidates
  legends <- get_plot_component(plot, "guide-box", return_all = TRUE)
  # find non-zero legends
  idx <- which(vapply(legends, \(x) !inherits(x, "zeroGrob"), TRUE))
  # return either the chosen or the first non-zero legend if it exists,
  # and otherwise the first element (which will be a zeroGrob) 
  if (length(idx) >= legend_number) {
    return(legends[[idx[legend_number]]])
  } else if (length(idx) >= 0) {
    return(legends[[idx[1]]])
  } else {
    return(legends[[1]])
  }
}
plot_legend <- get_legend_35(figure3_a)

plot_figure_3_no_legend <- plot_grid(
  figure3_a +
    theme(legend.position="none",
          plot.title = element_text(hjust = -0.1,
                                    vjust = -0.02,
                                   size=14, face="bold"))+
    labs(title = "(a)"),
  NULL,
  figure3_b+
    theme(legend.position="none",
          plot.title = element_text(hjust = -0.1,
                                    vjust = -0.02,
                                   size=14, face="bold"))+
    labs(title = "(b)"), 
  rel_widths = c(0.51,0.01, 0.5),
  ncol = 3)

plot_figure_3 <- plot_grid(plot_figure_3_no_legend,
                           NULL,
                           plot_legend, 
                           rel_heights = c(0.9,-0.01, 0.1),
                           ncol = 1)
          
plot_figure_3

tiff(file = "Figures/Figure_3.tiff",   
     width = 7.2, 
     height = 6.0,
     units = 'in',  
     res = 300, 
     compression = "lzw")
print(plot_figure_3)
dev.off()



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
hybrid_control_survival <- hybrid_model_survival(hybrid_control, 
                                                           t =seq(from = 0, to = 20, length.out = 1e3),
                                                           N = 1e3)

hybrid_active <- fit_hybrid_model(data = active, cut_point = 1.5)
hybrid_active_survival <- hybrid_model_survival(hybrid_active , 
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
  ggplot()+
  theme_classic()+
  geom_line(data = km_data_all %>% 
              filter(dataset == "ALEX trial") %>% 
              select(-c(dataset,model)) %>%
              mutate(treatment = trt),
            aes(x = time, y = surv, colour = treatment),
            alpha = 0.5)+
  geom_ribbon(aes(x = time, ymin = lower, ymax = upper, fill = treatment),
              alpha= 0.25, colour = NA)+
  geom_step(aes(x = time, y = value, colour = treatment))+
  scale_colour_discrete("Treatment")+
  scale_fill_discrete("Treatment")+
  geom_vline(xintercept = 1.5, alpha = 0.4)+
  scale_y_continuous(limits = c(0,1.01),labels = scales::percent)+
  ylab("Overall Survival")+
  xlab("Time (years)")

#max(hybrid_active_survival$upper)
#max(hybrid_control_survival$upper)







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







