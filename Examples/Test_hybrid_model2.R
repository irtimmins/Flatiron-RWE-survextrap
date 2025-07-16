# Load required packages
library(survival)
library(dplyr)
library(flexsurv)
library(ggplot2)
library(survextrap)
library(flexsurv)


################################################
# Function to fit hybrid model, also known as
# 'Liverpool approach'.
################################################

# data has time and status columns.
fit_hybrid_model <- function(data, cut_point){

  #  Step 1: Fit Kaplan-Meier on full data
  km_fit <- survfit(Surv(time, status) ~ 1, data = data)
  
  # Extract KM survival and standard error at cut_point
  # Find time point closest to cut_point
  km_summary <- summary(km_fit, times = cut_point)
  # survival at cut point
  S_km_cut_point <- km_summary$surv[1]
  # s.e. of survival at cut point
  SE_S_km_cut_point <- km_summary$std.err[1]  
  
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
  lambda_hat <- exp_fit$res[,"est"]
  var_lambda_hat <- exp_fit$res[,"se"]^2
  
  # Return key values in list as output
  results <- list("cut_point" = cut_point,
                  "km_fit" = km_fit,
                  "lambda_hat" = lambda_hat,
                  "var_lambda_hat" = var_lambda_hat,
                  "S_km_cut_point" = S_km_cut_point,
                  "SE_S_km_cut_point" = SE_S_km_cut_point)
  
  return(results)
  
  
}

################################################
# Function for evaluating survival probabilities
# from hybrid model.
################################################

hybrid_model_survival <- function(hybrid_model, t) {

  # Unpack elements from hybrid_model list into environment
  list2env(hybrid_model, envir = environment())
  
  # Prepare vectors to hold survival and standard errors (s.e.).
  Surv <- Surv_se <- numeric(length(t))
  
  # Divide time points before and after cut point.
  times_before_cut_point <- which(t <= cut_point)
  times_after_cut_point <- which(t > cut_point)
  
  # For time points before cut point, use Kaplan-Meier estimates,
  # and their standard errors.
  if(length(times_before_cut_point) > 0) {
    km_summary <- summary(km_fit, times=t[times_before_cut_point])
    Surv[times_before_cut_point] <- km_summary$surv
    Surv_se[times_before_cut_point] <- km_summary$std.err
  }
  
  # For time points after cut point, use delta method.
  if(length(times_after_cut_point) > 0) {
    
    # Take away the cut off point from the time-points,
    # and derive conditional survival from cut off point.
    times_subtract_cut_point <- t[times_after_cut_point] - cut_point
    Surv[times_after_cut_point] <- S_km_cut_point*exp(-lambda_hat * times_subtract_cut_point)
    
    # Apply the delta method to get the standard errors.
    # There is uncertainty in the estimate of survival at cut-point (A),
    # and uncertainty in the lambda from exponential model (B).
    # We take partial derivatives with respect to both of these estimators.
    
    # Let A = estimate of survival at cut-point from Kaplan-Meier.
    # and B = estimate of lambda from Exponential model.
    
    # After the time cut point we have:
    # Survival = f(A,B;t) = A*exp(-B*(t - cut_point))
    # Using the delta method, with partial derivatives:
    # Var[f(A,B;t)] = (df/dA)^2*Var[A]+(df/dB)^2*Var[B]
  
    # Note that df/dA = exp(-B*(t - cut_point))
    # and df/dB = -(t - cut_point)*A*exp(-B*(t - cut_point))
    # hence, the partial derivates df/dA and df/dB are:
    partial_wrt_A <- exp(-lambda_hat * times_subtract_cut_point)
    partial_wrt_B <- -(times_subtract_cut_point) * S_km_cut_point * exp(-lambda_hat * times_subtract_cut_point)
    
    # for Var[A] and Var[B] we have:
    var_A <- SE_S_km_cut_point^2
    var_B <- var_lambda_hat
    
    # Evaluate delta method using  
    # Var[f(A,B;t)] = (df/dA)^2*Var[A]+(df/dB)^2*Var[B].
    Surv_variance_from_A <- partial_wrt_A^2 * var_A
    Surv_variance_from_B <- partial_wrt_B^2 * var_B
    
    # Square root to get standard errors.
    Surv_se[times_after_cut_point] <- sqrt(Surv_variance_from_A + Surv_variance_from_B)

  }
  
  # Combine all results in tibble.
  # add standard errors and confidence intervals
  results <- tibble("variable" = "Survival",
                    "time" = t,
                    "value" = Surv,
                    "se" = NA) %>%     
    mutate(lower = NA,
           upper = NA) 
  
  return(results)
}

####################################################

hybrid_model_survival_semi_bootstrap <- function(data, cut_point, t, N){
  
  #  Step 1: Fit Kaplan-Meier on full data
  km_fit <- survfit(Surv(time, status) ~ 1, data = data)
  
  # Extract KM survival and standard error at cut_point
  # Find time point closest to cut_point
  km_summary <- summary(km_fit, times = cut_point)
  
  # survival at cut point
  S_km_cut_point <- km_summary$surv[1]
  # s.e. of survival at cut point
  SE_S_km_cut_point <- km_summary$std.err[1]  
  
  # Step 2: Fit exponential model beyond cut_point 
  
  # Restrict data to those at risk at cut_point
  # and subtract cut_point from all remaining times
  data_tail <- data %>%
    filter(time > cut_point) %>%
    mutate(time_tail = time - cut_point)
  
  # Only consider events and censored times after cut point months
  # Fit exponential distribution to the tail
  exp_fit <- flexsurvreg(Surv(time_tail, status) ~ 1, data = data_tail,  dist = "exp")
  
  # Extract the exponential hazard rate and s.e.:
  lambda <- exp_fit$res[,"est"]
  se_lambda <- exp_fit$res[,"se"]
  
  # Prepare vectors to hold survival probabilities.
  Surv <- numeric(length(N))
  
  #print(RMST)
  for(i in 1:N){
    
    #  print(i)
    # Set zero time points to NA.
    if(t == 0){
      Surv[i] <- NA
    }
    
    # For time points before cut point, use Kaplan-Meier estimates,
    # and their standard errors.
    
    
    if(t  <= cut_point) {
      
      km_summary <- summary(km_fit, times=t)
      mean_survival <- km_summary$surv %>% unname()
      se_survival <- km_summary$std.err %>% unname()
      
      Surv[i] <- rnorm(1, mean = mean_survival, sd = se_survival)
      
    }
    
    # For time points after cut point
    if(t > cut_point) {
      
      survival_boot <- rnorm(1, mean =   S_km_cut_point,
                             sd =   SE_S_km_cut_point)
      lambda_boot <- rnorm(1, mean = lambda, sd = se_lambda)
      
      # rmst = area under survival curve up to cut point, plus area after cut point.
      Surv[i] <- survival_boot*exp(-lambda_boot*(t-cut_point))
      
    }
  }
  
  return(sd(log(Surv)))
  
}

###################################################


hybrid_model_survival_log_scale <- function(hybrid_model, t) {
  
  # Unpack elements from hybrid_model list into environment
  list2env(hybrid_model, envir = environment())
  
  # Prepare vectors to hold survival and standard errors (s.e.).
  Surv <- Surv_se <- log_Surv <- log_Surv_se <-  numeric(length(t))
  
  # Divide time points before and after cut point.
  times_before_cut_point <- which(t <= cut_point)
  times_after_cut_point <- which(t > cut_point)
  
  # For time points before cut point, use Kaplan-Meier estimates,
  # and their standard errors.
  if(length(times_before_cut_point) > 0) {
    km_summary <- summary(km_fit, times=t[times_before_cut_point])
    Surv[times_before_cut_point] <- km_summary$surv
    Surv_se[times_before_cut_point] <- km_summary$std.err
  }
  
  # For time points after cut point, use delta method.
  if(length(times_after_cut_point) > 0) {
    
    # Take away the cut off point from the time-points,
    # and derive conditional survival from cut off point.
    times_subtract_cut_point <- t[times_after_cut_point] - cut_point
    Surv[times_after_cut_point] <- S_km_cut_point*exp(-lambda_hat * times_subtract_cut_point)
    

    log_Surv_se[times_after_cut_point] <- sqrt(SE_S_km_cut_point^2/S_km_cut_point^2 + var_lambda_hat*times_subtract_cut_point^2)
    
    # Square root to get standard errors.
    Surv_se[times_after_cut_point] <- S_km_cut_point*log_Surv_se[times_after_cut_point]
    
  }
  
  # Combine all results in tibble.
  # add standard errors and confidence intervals
  results <- tibble("variable" = "Survival",
                    "time" = t,
                    "value" = Surv,
                    "se" = Surv_se) %>%     
    mutate(lower = value-1.96*se,
           upper = value+1.96*se) 
  
  return(results)
}



################################################
# Function for evaluating restricted-mean
# survival from model.
################################################

hybrid_model_rmst <- function(hybrid_model, t) {
  
  # Unpack elements from hybrid_model list into environment
  list2env(hybrid_model, envir = environment())
  
  # Prepare vectors to hold survival and standard errors (s.e.).
  RMST <- RMST_se <- numeric(length(t))
  
  # Divide time points before and after cut point.
  times_zero <- which(t == 0)
  times_before_cut_point <- which(t <= cut_point & t > 0)
  times_after_cut_point <- which(t > cut_point)
  
  # Set zero time points to NA.
  if(!is.null(times_zero)){
    RMST[times_zero] <- NA
    RMST_se[times_zero] <- NA
  }
  
  # For time points before cut point, use Kaplan-Meier estimates,
  # and their standard errors.
  
  if(length(times_before_cut_point) > 0) {

    RMST[times_before_cut_point] <- sapply(
      t[times_before_cut_point], 
      function(t) {
        summary(km_fit, rmean = t)$table["rmean"] %>% unname()
      }) 
    
    RMST_se[times_before_cut_point] <- sapply(
      t[times_before_cut_point], 
      function(t) {
        summary(km_fit, rmean = t)$table["se(rmean)"] %>% unname()
      }) 

  }
  
  
  # For time points after cut point, use delta method.
  if(length(times_after_cut_point) > 0) {
    
    rmst_cut_point <- summary(km_fit, rmean = cut_point)$table["rmean"] %>% unname()
    se_rmst_cut_point <- summary(km_fit, rmean = cut_point)$table["se(rmean)"] %>% unname()
    
    times_subtract_cut_point <- t[times_after_cut_point] - cut_point
    
    # area under survival curve after cut point.
    rmst_after_cut_point <-  S_km_cut_point*(1-exp(-lambda_hat * times_subtract_cut_point))/lambda_hat
    
    # rmst = area under survival curve up to cut point, plus area after cut point.
    RMST[times_after_cut_point] <- rmst_cut_point+rmst_after_cut_point
    
    # for Var[A] and Var[B] we have:

    var_A <- SE_S_km_cut_point^2
    var_B <- var_lambda_hat
    
    partial_wrt_A <- (1-exp(-lambda_hat * times_subtract_cut_point))/lambda_hat
    partial_wrt_B <- S_km_cut_point *
      ((1-lambda_hat*times_subtract_cut_point)*exp(-lambda_hat * times_subtract_cut_point)-1)/lambda_hat^2
    
    # print(se_rmst_cut_point^2)
    # print(var_A*partial_wrt_A^2)
    # print(var_B*partial_wrt_B^2)
    
    var_all <- se_rmst_cut_point^2+var_A*partial_wrt_A^2+var_B*partial_wrt_B^2
    
    RMST_se[times_after_cut_point] <- sqrt(var_all)
    
    
  }
  
  
  # Combine all results in tibble.
  # add standard errors and confidence intervals
  results <- tibble("variable" = "rmst",
                    "time" = t,
                    "value" = RMST,
                    "se" = RMST_se) %>%     
    mutate(lower = value-1.96*se,
           upper = value+1.96*se) 
  
  return(results)
  
}


hybrid_model_rmst_semi_bootstrap <- function(data, cut_point, t, N){
  
  #  Step 1: Fit Kaplan-Meier on full data
  km_fit <- survfit(Surv(time, status) ~ 1, data = data)
  
  # Extract KM survival and standard error at cut_point
  # Find time point closest to cut_point
  km_summary <- summary(km_fit, times = cut_point)
  
  # survival at cut point
  S_km_cut_point <- km_summary$surv[1]
  # s.e. of survival at cut point
  SE_S_km_cut_point <- km_summary$std.err[1]  
  
  # Step 2: Fit exponential model beyond cut_point 
  
  # Restrict data to those at risk at cut_point
  # and subtract cut_point from all remaining times
  data_tail <- data %>%
    filter(time > cut_point) %>%
    mutate(time_tail = time - cut_point)
  
  # Only consider events and censored times after cut point months
  # Fit exponential distribution to the tail
  exp_fit <- flexsurvreg(Surv(time_tail, status) ~ 1, data = data_tail,  dist = "exp")
  
  # Extract the exponential hazard rate and s.e.:
  lambda_hat <- exp_fit$res[,"est"]
  var_lambda_hat <- exp_fit$res[,"se"]^2
  
  # Prepare vectors to hold RMST.
  RMST <- numeric(length(N))
  
  #print(RMST)
  for(i in 1:N){
    
    #  print(i)
    # Set zero time points to NA.
    if(t == 0){
      RMST[i] <- NA
    }
    
    # For time points before cut point, use Kaplan-Meier estimates,
    # and their standard errors.
    if(t  <= cut_point) {
      
      mean_RMST <- summary(km_fit, rmean = t)$table["rmean"] %>% unname()
      se_RMST <- summary(km_fit, rmean = t)$table["se(rmean)"] %>% unname()
      
      RMST[i] <- rnorm(1, mean = mean_RMST, sd = se_RMST)
      
    }
    
    # For time points after cut point
    if(t > cut_point) {
      
      rmst_cut_point <- summary(km_fit, rmean = cut_point)$table["rmean"] %>% unname()
      se_rmst_cut_point <- summary(km_fit, rmean = cut_point)$table["se(rmean)"] %>% unname()
      
      rmst_cut_point_boot <- rnorm(1, mean = rmst_cut_point, sd = se_rmst_cut_point)
      
      times_subtract_cut_point <- t - cut_point
      
      S_km_cut_point_boot <- rnorm(1, mean = S_km_cut_point, sd = SE_S_km_cut_point)
      lambda_hat_boot <- rnorm(1, mean = lambda_hat, sd = sqrt(var_lambda_hat))
      
      # area under survival curve after cut point.
      rmst_after_cut_point <-  S_km_cut_point_boot*(1-exp(-lambda_hat_boot * times_subtract_cut_point))/lambda_hat_boot
      
      # rmst = area under survival curve up to cut point, plus area after cut point.
      RMST[i] <- rmst_cut_point_boot+rmst_after_cut_point
      
    }
  }
  
  return(RMST)

}


hybrid_model_rmst_bootstrap <- function(data, cut_point, t, N) {
  
  
  mean <- get_rmst(data, cut_point, t)
  rmst_sample <- numeric(length(N))
  
  for(i in 1:N){
  
    data_boot <- data[sample(1:nrow(data), replace = T),] 
    #print(sum(data_boot$status[data_boot$time > cut_point]))
    rmst_boot <- get_rmst(data_boot, cut_point, t)
    rmst_sample[i] <- rmst_boot
    if(i %% 100 == 0){
    print(paste0(i, "/", N))
    }
  }
  
  
  return(rmst_sample)
  
  
  
}


get_rmst <- function(data, cut_point, t){
  
  #  Step 1: Fit Kaplan-Meier on full data
  km_fit <- survfit(Surv(time, status) ~ 1, data = data)
  
  # Extract KM survival and standard error at cut_point
  # Find time point closest to cut_point
  km_summary <- summary(km_fit, times = cut_point)
  # survival at cut point
  S_km_cut_point <- km_summary$surv[1]
  
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
  lambda_hat <- exp_fit$res[,"est"]
  
  # Prepare vectors to hold survival and standard errors (s.e.).
  RMST <- numeric(length(t))
  
  # Divide time points before and after cut point.
  times_zero <- which(t == 0)
  times_before_cut_point <- which(t <= cut_point & t > 0)
  times_after_cut_point <- which(t > cut_point)
  
  # Set zero time points to NA.
  if(!is.null(times_zero)){
    RMST[times_zero] <- NA
  }
  
  # For time points before cut point, use Kaplan-Meier estimates,
  # and their standard errors.
  
  if(length(times_before_cut_point) > 0) {
    
    RMST[times_before_cut_point] <- sapply(
      t[times_before_cut_point], 
      function(t) {
        summary(km_fit, rmean = t)$table["rmean"] %>% unname()
      }) 
    
  }
  
  
  # For time points after cut point, use delta method.
  if(length(times_after_cut_point) > 0) {
    
    rmst_cut_point <- summary(km_fit, rmean = cut_point)$table["rmean"] %>% unname()
    
    times_subtract_cut_point <- t[times_after_cut_point] - cut_point
    
    # area under survival curve after cut point.
    rmst_after_cut_point <-  S_km_cut_point*(1-exp(-lambda_hat * times_subtract_cut_point))/lambda_hat
    
    # rmst = area under survival curve up to cut point, plus area after cut point.
    RMST[times_after_cut_point] <- rmst_cut_point+rmst_after_cut_point
    
  }
  
  
  # Combine all results in tibble.
  # add standard errors and confidence intervals
  return(RMST)
  
  
}


hybrid_model_survival_bootstrap <- function(data, cut_point, t, N) {
  
  
  mean <- get_survival(data, cut_point, t)
  survival_sample <- numeric(length(N))
  
  for(i in 1:N){
    
    data_boot <- data[sample(1:nrow(data), replace = T),] 
    #print(sum(data_boot$status[data_boot$time > cut_point]))
    survival_boot <- get_survival(data_boot, cut_point, t)
    survival_sample[i] <- survival_boot
    if(i %% 100 == 0){
      print(paste0(i, "/", N))
    }
  }
  
  
  return(survival_sample)
  
  
  
}


get_survival <- function(data, cut_point, t){
  
  #  Step 1: Fit Kaplan-Meier on full data
  km_fit <- survfit(Surv(time, status) ~ 1, data = data)
  
  # Extract KM survival and standard error at cut_point
  # Find time point closest to cut_point
  km_summary <- summary(km_fit, times = cut_point)
  # survival at cut point
  S_km_cut_point <- km_summary$surv[1]
  
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
  lambda_hat <- exp_fit$res[,"est"]
  
  # Prepare vectors to hold survival and standard errors (s.e.).
  Surv <- numeric(1)
  
  # Set zero time points to NA.
  if(t == 0){
    Surv <- NA
  }
  
  # For time points before cut point, use Kaplan-Meier estimates,
  # and their standard errors.
  
  if(t <= cut_point) {
    
    Surv <- summary(km_fit, rmean = t)$table["rmean"] %>% unname()
    
  }
  
  
  # For time points after cut point
  
  if(t > cut_point) {
    
    # survival 
    Surv <-  S_km_cut_point*(exp(-lambda_hat * (t - cut_point)))
    
  }
  return(Surv)
  
}
