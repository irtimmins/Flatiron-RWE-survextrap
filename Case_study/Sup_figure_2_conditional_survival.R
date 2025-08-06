

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