model_trial_data


model_trial_data_test <-  survextrap(Surv(time,status) ~ trt, 
                                     data = trial_data , 
                                     # backhaz = cetux_bh,
                                     #  nonprop = T,
                                     df = 3 ,
                                     prior_hsd = p_gamma(2, 5),
                                     #  prior_hrsd = p_gamma(2, 20)
                                     #  add_knots = c(6)#,
                                     fit_method = "opt"
)
  
rmst(model_trial_data_test, t = 5)

irmst(model_trial_data_test, t = 5, newdata =   tibble(trt = c("Crizotinib", "Alectinib")))



control_trial_data <- survextrap(Surv(time,status) ~ 1, 
                                 data = trial_data %>% 
                                   filter(trt == "Crizotinib") , 
                                 # backhaz = cetux_bh,
                                 #  nonprop = T,
                                 df = 3 ,
                                 prior_hsd = p_gamma(2, 5),
                                 #  prior_hrsd = p_gamma(2, 20)
                                 #  add_knots = c(6)#,
                                 fit_method = "opt")


active_trial_data <- survextrap(Surv(time,status) ~ 1, 
                                 data = trial_data %>% 
                                  filter(trt == "Alectinib") , 
                                 # backhaz = cetux_bh,
                                 #  nonprop = T,
                                 df = 3 ,
                                 prior_hsd = p_gamma(2, 5),
                                 #  prior_hrsd = p_gamma(2, 20)
                                 #  add_knots = c(6)#,
                                 fit_method = "opt")

rmst(control_trial_data, t = c(5,10,15,20))

rmst(active_trial_data, t = c(5,10,15,20))

control_sample <- rmst(control_trial_data, t = c(5,10,15,20), sample = T)

active_sample <- rmst(active_trial_data, t = c(5,10,15,20), sample = T)

irmst_samples <- active_sample - control_sample

survextrap:::summarise_output(irmst_samples , t =  c(5,10,15,20), 
                              summ_fns = list("median" = median, 
                                              ~quantile(.x, probs=c(0.025, 0.975))),
                              newdata = NULL, 
                              summ_name = "irmst", sample = F) %>%
  rename(lower= "2.5%", upper = "97.5%")




# irmst_sep <- posterior::summarise_draws(irmst_sep, median,
#                                         ~quantile(.x, probs=c(0.025, 0.975))) #%>%
#   rename(t=variable, lower="2.5%", ir_upper="97.5%") %>% 
#   select(ir_med, ir_lower, ir_upper)
# 
# #?irmst()
