
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
# Log-cumulative hazard
######################################################


plot_sup_3a <-km_all %>%
  mutate(dataset = factor(dataset, 
                          levels = c("ALEX", "PROFILE-1014", 
                                     "Flatiron RWE MAIC Weighted", "Flatiron RWE Unweighted"))) %>%
  filter(dataset == "ALEX") %>%
  ggplot(aes(x = time, y = log(-log(surv)), colour = trt))+
  theme_classic()+
  theme(legend.position.inside = c(0.6,0.3),
        legend.position = "inside",
        legend.text=element_text(size=9),
        legend.key.spacing.y = unit(-8, "pt"),
        legend.box.spacing = unit(-10, "pt"),
        legend.title = element_text(size=12,
                                    margin = margin(l = 0, r = 0, b = 0, t = 0))) +
  geom_step(linewidth = 0.8)+
  scale_colour_manual("Treatment", values = c("#4393c3","#d6604d"))+
  scale_y_continuous("Log Cumulative Hazard")+
  scale_x_continuous("Time (years)",breaks = 0:5, limits = c(0,3))


plot_sup_3b <- km_all %>%
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
  ggplot(aes(x = time, y = log(-log(surv)), colour = dataset))+
  theme_classic()+
  theme(legend.position.inside = c(0.6,0.3),
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
   scale_y_continuous("Log Cumulative Hazard")+
   scale_x_continuous("Time (years)",  breaks = 0:5, limits = c(0,5.5))

plot_sup_3b

figure_sup_3 <- plot_grid(plot_sup_3a+
                        theme(  plot.title = element_text(hjust = -0.1,
                                                          size=14, face="bold"))+
                        labs(title = "(a)"),
                        plot_sup_3b+
                        theme(  plot.title = element_text(hjust = -0.1,
                                                          size=14, face="bold"))+
                        labs(title = "(b)"),
                      align = "h",
                      rel_widths  =c(0.5,0.5),
                      ncol = 2)
figure_sup_3


tiff(file = "Figures/Sup_figure_3.tiff",   
     width = 7.6, 
     height = 3.7,
     units = 'in',  
     res = 300, 
     compression = "lzw")
print(figure_sup_3)
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


tiff(file = "Figures/Sup_Figure_1.tiff",   
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
#end_of_trial <- round(max(trial_data$time[trial_data$trt == "Crizotinib"]), digits = 1)
# end_of_trial <- 2
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
print(sup_figure2)

tiff(file = "Figures/Sup_Figure_2.tiff",   
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






