library(memoise)
library(survextrap)

######################################################
#  Memoise survextrap functions.
######################################################

store_wd <- "/projects/aa/klvq491/Flatiron_ansclc/cache"

survextrap_mem <- memoise(survextrap, 
                          cache = cachem::cache_disk(
                            dir = store_wd,
                            max_age = 31557600, # keep for 1-year
                            max_size = 2e10)) # 20Gb of space 

##################################################
# Fit survextrap models to two arms.
##################################################

fit_model <- function(model, datasets, df , hsd_rate, hrsd_rate, add_knots, fit_method, store_file){
  

  # Specify external dataset.
  
  add_knots <- get(add_knots)
  
  
  if(datasets == "trial_only") external_data <- NULL
  if(datasets == "trial_and_historic") external_data <- historic_trial_aggregate
  if(datasets == "trial_and_all_maic") external_data <- external_data_maic_weighted
  if(datasets == "trial_and_all_unweighted") external_data <-   external_data_unweighted

  
  # Specify M-spline, and add extra knots for external data.
  
  if(datasets == "trial_only"){
    
    mspline <-  mspline_spec(Surv(time,status) ~ trt, 
                             data = trial_data, 
                             df=df)
    
  } else if (datasets %in% c("trial_and_historic", 
                             "trial_and_all_maic",
                             "trial_and_all_unweighted")) {

    mspline <-  mspline_spec(Surv(time,status) ~ trt, 
                             data = trial_data,
                             df=df,
                             add_knots = add_knots)
        
  }

  # Define priors.
  
  prior_hscale <- p_meansurv(median=5, upper=20, mspline=mspline)
  prior_hsd <- p_gamma(2, hsd_rate)
  prior_loghr <- p_hr(median=1, upper=10)
  
  if(model == "NON-PH"){
    prior_hrsd <- p_gamma(2, hrsd_rate)
  }
  
    # Fit models.
  
  
  if(model == "PH"){
    
    results <- survextrap_mem(Surv(time,status) ~ trt, 
                              data = trial_data , 
                              external = external_data,
                              mspline = mspline,
                              prior_hscale = prior_hscale,
                              prior_loghr = prior_loghr,
                              prior_hsd = prior_hsd,
                              fit_method = fit_method)
    
    
  } else if(model == "NON-PH"){
    
    results <- survextrap_mem(Surv(time,status) ~ trt, 
                              data = trial_data , 
                              nonprop = T,
                              external = external_data,
                              mspline = mspline,
                              prior_hscale = prior_hscale,
                              prior_loghr = prior_loghr,
                              prior_hsd = prior_hsd,
                              prior_hrsd = prior_hrsd,
                              fit_method = fit_method)
    
  }else if(model == "Separate_arms"){
    
    results_control <- survextrap_mem(Surv(time,status) ~ 1, 
                                      data = trial_data %>% filter(trt == "Crizotinib") , 
                                      external = external_data,
                                      mspline = mspline,
                                      prior_hscale = prior_hscale,
                                      prior_loghr = prior_loghr,
                                      prior_hsd = prior_hsd,
                                      fit_method = fit_method)
    
    
    results_active <- survextrap_mem(Surv(time,status) ~ 1, 
                                     data = trial_data %>% filter(trt == "Alectinib") , 
                                     mspline = list("df"=df),
                                     prior_hscale = prior_hscale,
                                     prior_loghr = prior_loghr,
                                     prior_hsd = prior_hsd,
                                     fit_method = fit_method)
    
    results <- list(control = results_control, active = results_active)
    class(results) <- "two_models"
    
  }
  
  saveRDS(results, store_file)
  
}


##################################################
# Survival and hazard estimates from models.
##################################################

get_survival_and_hazard_survextrap <- function(model_file, store_file){
  
  print(model_file)
  
  model <- readRDS(model_file)
  
  new_data_estimate <- tibble(trt = c("Alectinib", "Crizotinib"))
  
  time_vec <- c(seq(from = 0, to = 5, length.out = 1e2),
                seq(from = 5, to = 20, length.out = 1e2))
  
  
  if(class(model) == "survextrap"){ # PH or NON-PH scenarios.
    
    hazard_all <- hazard(model, t = time_vec, newdata = new_data_estimate) 
    # print(hazard_all)
    survival_all <- survival(model, t = time_vec, newdata = new_data_estimate) 
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


##################################################
# further helper functions.
##################################################


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

filter_model <- function(.data, model = "PH", df_value = 3, hrsd_rate_value = NULL ){
  
  if(model == "PH"){
    
    .data %>%
      filter(model == "PH") %>%
      filter(df == df_value)
    
  } else if(model == "NON-PH"){
    
    .data %>%
      filter(model == "NON-PH") %>%
      filter(hrsd_rate == hrsd_rate_value, df == df_value) 
    
  } else if(model == "Separate_arms"){
    
    .data %>%
      filter(model == "Separate_arms") %>%
      filter(df == df_value)
    
  }
}
