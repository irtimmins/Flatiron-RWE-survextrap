











######################


base_scenarios
View(base_scenarios)

test_file <- base_scenarios %>%
  filter(model == "NON-PH",
         df == 6,
         hrsd_rate == 1,
         hsd_rate == 1,
         datasets == "trial_and_all")  %>%
  pull(store_file)

test_model <- readRDS(test_file)

test_model$mspline


test_model <- readRDS("/projects/aa/klvq491/Flatiron_ansclc/models/base_model_18.rds")

test_model$mspline


mspline_mod1 <-  mspline_spec(Surv(time,status) ~ trt, 
                         data = trial_data,
                         df=5,
                         add_knots = c(2.0, 3, 5.0))

test_mod1 <- survextrap(Surv(time,status) ~ trt, 
                          data = trial_data , 
                          nonprop = T,
                          external = external_data_aggregate, #%>% filter(start >= 3),
                          mspline = test_model$mspline,
                      prior_hsd = p_gamma(2,10))#,
#     fit_method = "opt")

plot(test_mod1, tmax = 20, ci = T)

##################################

mspline_mod2 <-  mspline_spec(Surv(time,status) ~ trt, 
                              data = trial_data,
                              df=6)#,
                              #add_knots = c(5.0))

test_mod2 <- survextrap(Surv(time,status) ~ trt, 
                        data = trial_data , 
                        nonprop = T,
                        mspline = mspline_mod2,
                        prior_hsd = p_gamma(2,10))#,
#     fit_method = "opt")

plot(test_mod2, tmax = 20, ci = T)


test_exp <- flexsurvreg(Surv(time, status) ~ 1, data = trial_data,  dist = "exp")
#test_exp
#test_exp$
  
