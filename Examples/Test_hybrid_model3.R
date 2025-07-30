library(tibble)
library(dplyr)
library(tidyr)


################################################
# Function to fit hybrid model, also known as
# 'Liverpool approach'.
################################################

# data has time and status columns.
fit_hybrid_model <- function(data, cut_point){
  
  # data <- control
  # cut_point <- 1
  #data
  #  Step 1: Fit Kaplan-Meier on full data
  km_fit <- survfit(Surv(time, status) ~ 1, data = data)
  
  # Extract KM survival and standard error at cut_point
  # Find time point closest to cut_point
  km_summary <- summary(km_fit, times = cut_point)
  # survival at cut point
  survival_at_cut_point <- km_summary$surv[1]
  # s.e. of survival at cut point
  survival_at_cut_point_se <- km_summary$std.err[1]  
  
  # Step 2: Fit exponential model beyond cut_point 
  
  # Restrict data to those at risk at cut_point
  # and subtract cut_point from all remaining times
  data_tail <- data %>%
    filter(time > cut_point) %>%
    mutate(time_tail = time - cut_point)
  
  # Only consider events and censored times after cut point months
  # Fit exponential distribution to the tail
  exp_fit <- flexsurvreg(Surv(time_tail, status) ~ 1, data = data_tail,  dist = "exp")
  
  # Extract the exponential hazard rate:
  lambda <- exp_fit$res[,"est"]
  lambda_se <- exp_fit$res[,"se"]
  
  # Return key values in list as output
  results <- list("cut_point" = cut_point,
                  "km_fit" = km_fit,
                  "lambda" = lambda,
                  "lambda_se" = lambda_se,
                  "survival_at_cut_point" = survival_at_cut_point,
                  "survival_at_cut_point_se" = survival_at_cut_point_se)
  
  return(results)
  
  
}


################################################
# Get survival estimates from this model.
################################################


hybrid_model_survival <- function(hybrid_model, t, N) {
  
  
  # Combine all results in tibble.
  # add standard errors and confidence intervals
  results <- tibble("variable" = "Survival",
                    "time" = t,
                    "value" = 0,
                    "se" = 0) %>%     
    mutate(lower = value-1.96*se,
           upper = value+1.96*se) 
  
  
  # Unpack elements from hybrid_model list into environment
  list2env(hybrid_model, envir = environment())
  
  store_results <- expand_grid(time_id = t, iteration_id = 1:N) %>%
    mutate(value = NA)
  
  
  for(i in 1:N){
    
    
    survival_boot <- get_survival(km_fit = km_fit, 
                                  cut_point = cut_point, 
                                  survival_at_cut_point = survival_at_cut_point, 
                                  survival_at_cut_point_se = survival_at_cut_point_se, 
                                  lambda = lambda,
                                  lambda_se = lambda_se,
                                  t = t,
                                  with_error = T) 
    
    results_boot <- tibble(time_id = t,
                           iteration_id = i,
                           value = survival_boot)
    
    store_results <- store_results %>%
      left_join(results_boot, by = c("time_id", "iteration_id")) %>%
      mutate(value = coalesce(value.y, value.x)) %>%
      select(-value.y, -value.x)
    
    if(i %% 100 == 0){
      print(paste0(i, "/", N))
    }
    
    
  }
  
  surv_results <- get_survival(km_fit = km_fit, 
                               cut_point = cut_point, 
                               survival_at_cut_point = survival_at_cut_point,
                               survival_at_cut_point_se = survival_at_cut_point_se,
                               lambda = lambda, 
                               lambda_se = lambda_se,
                               t = t,
                               with_error = F)
  
  # Get s.e. from bootstrap
  se_results <- store_results %>%
    group_by(time_id) %>%
    summarise(lower = quantile(value, 0.025),
              upper = quantile(value, 0.975),
              se = sd(value))
  
  results <- results %>%
    mutate(
      value = surv_results, 
      se = se_results$se) %>%
    mutate(lower = se_results$lower,
           upper = se_results$upper)
  
  return(results)
}




hybrid_model_rmst <- function(hybrid_model, t, N) {
  
  
  # Combine all results in tibble.
  # add standard errors and confidence intervals
  results <- tibble("variable" = "rmst",
                    "time" = t,
                    "value" = 0,
                    "se" = 0) %>%     
    mutate(lower = value-1.96*se,
           upper = value+1.96*se) 
  
  
  # Unpack elements from hybrid_model list into environment
  list2env(hybrid_model, envir = environment())
  
  store_results <- expand_grid(time_id = t, iteration_id = 1:N) %>%
    mutate(value = NA)
  
  
  for(i in 1:N){
    
    rmst_boot <- get_rmst(
      km_fit = km_fit, 
      cut_point = cut_point, 
      survival_at_cut_point = survival_at_cut_point, 
      survival_at_cut_point_se = survival_at_cut_point_se, 
      lambda = lambda,
      lambda_se = lambda_se,
      t = t,
      with_error = T) 
    
    results_boot <- tibble(
      time_id = t,
      iteration_id = i,
      value = rmst_boot)
    
    store_results <- store_results %>%
      left_join(results_boot, by = c("time_id", "iteration_id")) %>%
      mutate(value = coalesce(value.y, value.x)) %>%
      select(-value.y, -value.x)
    
    if(i %% 100 == 0){
      print(paste0(i, "/", N))
    }
    
    
  }
  
  rmst_results <- get_rmst(
    km_fit = km_fit, 
    cut_point = cut_point, 
    survival_at_cut_point = survival_at_cut_point,
    survival_at_cut_point_se = survival_at_cut_point_se,
    lambda = lambda, 
    lambda_se = lambda_se,
    t = t,
    with_error = F)
  
  # Get s.e. from bootstrap
  se_results <- store_results %>%
    group_by(time_id) %>%
    summarise(lower = quantile(value, 0.025),
              upper = quantile(value, 0.975),
              se = sd(value))
  
  results <- results %>%
    mutate(
      value = rmst_results, 
      se = se_results$se) %>%
    mutate(lower = se_results$lower,
           upper = se_results$upper)
  
  return(results)
}


###################################################
## dependencies
###################################################

get_survival <- function(km_fit, 
                         cut_point, 
                         survival_at_cut_point, 
                         survival_at_cut_point_se, 
                         lambda, 
                         lambda_se,  
                         t, 
                         with_error = F){

  # Extract KM survival and standard error at cut_point
  # Find time point closest to cut_point
  km_summary <- summary(km_fit, times = cut_point)

  # Prepare vectors to hold survival and standard errors (s.e.).
  Surv <- numeric(length(t))

  # Divide time points before and after cut point.
  times_before_cut_point <- which(t <= cut_point)
  times_after_cut_point <- which(t > cut_point)
  
  # For time points before cut point, use Kaplan-Meier estimates,

  if(length(times_before_cut_point) > 0) {
   
    km_summary <- summary(km_fit, times=t[times_before_cut_point])
    
    Surv[times_before_cut_point] <- km_summary$surv
    
    # add error from Kaplan-Meier estimates before cut point.
    if(with_error == T){
      Surv[times_before_cut_point] <- Surv[times_before_cut_point] + rnorm(1)*km_summary$std.err
    }
  }
  
  # For time points after cut point
  if(length(times_after_cut_point) > 0) {
    
    times_subtract_cut_point <- t[times_after_cut_point] - cut_point
    
    if(with_error == T){
      # parametric bootstrap
      survival_at_cut_point <- survival_at_cut_point+rnorm(1)*survival_at_cut_point_se
      
      # on log scale for lambda, with delta approximation
      lambda <- exp(log(lambda)+rnorm(1)*lambda_se/lambda)
      
      Surv[times_after_cut_point] <-  survival_at_cut_point*exp(-lambda*times_subtract_cut_point)
      
    } else{
      
      # or get the mean.
      Surv[times_after_cut_point] <-  survival_at_cut_point*exp(-lambda*times_subtract_cut_point)
    }
    
    
  }
  
  return(Surv)
}

get_rmst <- function(km_fit, 
                     cut_point, 
                     survival_at_cut_point, 
                     survival_at_cut_point_se, 
                     lambda, 
                     lambda_se,  
                     t, 
                     with_error = F){
  
  # Extract KM survival and standard error at cut_point
  # Find time point closest to cut_point
  km_summary <- summary(km_fit, times = cut_point)
  
  # Prepare vectors to hold survival and standard errors (s.e.).
  RMST <- numeric(length(t))
  
  # Divide time points before and after cut point.
  times_before_cut_point <- which(t <= cut_point)
  times_after_cut_point <- which(t > cut_point)
  
  # For time points before cut point, use Kaplan-Meier estimates,
  
  if(length(times_before_cut_point) > 0) {
    
    km_summary <- summary(km_fit, times=t[times_before_cut_point])

    RMST[times_before_cut_point] <- km_summary$table["rmean"] %>% unname()
    
    # add error from Kaplan-Meier estimates before cut point.
    if(with_error == T){
      RMST[times_before_cut_point] <- RMST[times_before_cut_point] + rnorm(1)*(km_summary$table["se(rmean)"] %>% unname())
    }
  }
  
  # For time points after cut point
  if(length(times_after_cut_point) > 0) {
    
    times_subtract_cut_point <- t[times_after_cut_point] - cut_point
    
    if(with_error == T){
      # add standard error.  
      
      rmst_cut_point <- summary(km_fit, rmean = cut_point)$table["rmean"] %>% unname()
      rmst_cut_point_se <- summary(km_fit, rmean = cut_point)$table["se(rmean)"] %>% unname()
      
      rmst_cut_point <- rmst_cut_point + rnorm(1)*rmst_cut_point_se
      
      # area under survival curve after cut point.
      survival_at_cut_point <- survival_at_cut_point+rnorm(1)*survival_at_cut_point_se
      
      # on log scale for lambda, with delta approximation
      lambda <- exp(log(lambda)+rnorm(1)*lambda_se/lambda)
      
      rmst_after_cut_point <- survival_at_cut_point*(1-exp(-lambda * times_subtract_cut_point))/lambda
      
      survival_at_cut_point <- survival_at_cut_point+rnorm(1)*survival_at_cut_point_se
      lambda <- lambda+rnorm(1)*lambda_se
      
      RMST[times_after_cut_point] <- rmst_cut_point+rmst_after_cut_point
      
    } else{
      
      # or get the mean from the model.
      
      rmst_cut_point <- summary(km_fit, rmean = cut_point)$table["rmean"] %>% unname()

      # area under survival curve after cut point.
      rmst_after_cut_point <- survival_at_cut_point*(1-exp(-lambda * times_subtract_cut_point))/lambda

      RMST[times_after_cut_point] <- rmst_cut_point+rmst_after_cut_point
    }
    
    
  }
  
  return(RMST)
  
}


