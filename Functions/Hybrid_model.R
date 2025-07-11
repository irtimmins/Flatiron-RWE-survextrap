# Load required packages
library(survival)
library(dplyr)
library(flexsurv)

test_data <- readRDS("Data/trial_IPD_OS_ALEX_Feb_2017.rds")
test_data <- test_data %>%
  filter(trt == "Alectinib") %>%
  select(time, status)

# Cut off point at 18 months (in years).
cut_point <- 1.5

hybrid_model_test <- fit_hybrid_model(data = test_data, 
                                      cut_point = 1.5)

hybrid_model_test
test_survival <- hybrid_model_survival(hybrid_model_test, t = seq(from = 0, to = 10, length.out = 100))

test_survival %>%
  ggplot(aes(x = time, y = value))+
  theme_classic()+
  geom_ribbon(aes(x = time, ymin = lower, ymax = upper),
              alpha= 0.2)+
  geom_line()+
  geom_vline(xintercept = 1.5)+
  scale_y_continuous(limits = c(0,1))

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
    
    # for Var[A] and Var[B] we have:
    var_A <- SE_S_km_cut_point^2
    var_B <- var_lambda_hat
    
    # the partial derivates df/dA and df/dB are:
    partial_wrt_A <- exp(-lambda_hat * times_subtract_cut_point)
    partial_wrt_B <- -(times_subtract_cut_point) * S_km_cut_point * exp(-lambda_hat * times_subtract_cut_point)
    
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
  RMST[times_zero] <- NA
  RMST_se[time_zero] <- NA
  
  # For time points before cut point, use Kaplan-Meier estimates,
  # and their standard errors.
  
  if(length(times_before_cut_point) > 0) {
    
    # test_km <- survfit(Surv(time, status) ~ 1, data = test_data)
    # summary(test_km, rmean = 40)$table["rmean"] %>% unname()
    # 
    # 
    # sapply(c(0.3,0.8, 1.5), function(t) {
    #   summary(test_km, rmean = t)$table["rmean"] %>% unname()
    # }) 
    # 
    # sapply(c(0.3,0.8, 1.5), function(t) {
    #   summary(test_km, rmean = t)$table["se(rmean)"] %>% unname()
    # }) 
    
    # Use sapply to vectorise.
    
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
    
    ## .....
    
    RMST1_cut_point <- RMST1_se_cut_point <- numeric(length(times_after_cut_point))
    
    ## .....
    
    
    # for Var[A] and Var[B] we have:
    var_A <- SE_S_km_cut_point^2
    var_B <- var_lambda_hat
    
    partial_wrt_A <- (1-exp(-lambda_hat * times_subtract_cut_point))/lambda_hat
    partial_wrt_B <- S_km_cut_point *
      ((1-lambda_hat*times_subtract_cut_point)*exp(-lambda_hat * times_subtract_cut_point)-1)/lambda_hat^2
    
    
  }
  
  
}

