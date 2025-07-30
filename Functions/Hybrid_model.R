library(tibble)
library(dplyr)
library(tidyr)
library(survRM2)


################################################
# Function to fit hybrid model, also known as
# 'Liverpool approach'.
################################################

# data has time and status columns.
fit_hybrid_model <- function(data, cut_point){
  
  #  Fit Kaplan-Meier on full data
  km_fit <- survfit(Surv(time, status) ~ 1, data = data)
  
  # Then fit exponential distribution to use for tail.
  exp_fit <- flexsurvreg(Surv(time, status) ~ 1, data = data,  dist = "exp")
  smallest_survival <- min(data$time[data$status == 1])

  # Return key values in list as output
  results <- list("km_fit" = km_fit,
                  "exp_fit" = exp_fit,
                  "cut_point" = cut_point,
                  "smallest_survival" = smallest_survival)
  
  return(results)
  
  
}


################################################
# Get survival estimates from this model.
################################################

hybrid_model_survival <- function(hybrid_model, t, N) {
  
  # Unpack elements from hybrid_model list into environment
  list2env(hybrid_model, envir = environment())
  
  # Prepare vectors to hold survival and standard errors (s.e.).
  Surv <- Surv_SE <- Surv_Upper <- Surv_Lower <- numeric(length(t))
  
  # Divide time points before and after cut point.
  times_before_cut_point <- which(t <= cut_point)
  times_after_cut_point <- which(t > cut_point)
  
  km_summary <- summary(km_fit, times=t[times_before_cut_point])
  
  Surv[times_before_cut_point] <- km_summary$surv
  Surv_SE[times_before_cut_point] <- km_summary$std.err
  Surv_Lower[times_before_cut_point] <- km_summary$surv - 1.96*km_summary$std.err
  Surv_Upper[times_before_cut_point] <- km_summary$surv + 1.96*km_summary$std.err
  
  
  # Then use exponential tail.
  
  tail_results <- summary(exp_fit, 
          t = t[times_after_cut_point],
          se = T)[[1]] %>%
    as_tibble() %>%
    rename(value = est, lower = lcl, upper = ucl) %>%
    select(time, value, se, lower, upper)

    Surv[times_after_cut_point] <- tail_results$value
    Surv_SE[times_after_cut_point] <- tail_results$se
    Surv_Lower[times_after_cut_point] <- tail_results$lower
    Surv_Upper[times_after_cut_point] <- tail_results$upper
  
  results <- tibble(
      time = t,
      value = Surv, 
      se = Surv_SE,
      lower = Surv_Lower,
      upper = Surv_Upper)
  
  return(results)
}



hybrid_model_rmst <- function(hybrid_model, t, N) {
  
  # Unpack elements from hybrid_model list into environment
  list2env(hybrid_model, envir = environment())
  
  # Prepare vectors to hold survival and standard errors (s.e.).
  RMST <- RMST_SE <- RMST_Upper <- RMST_Lower <- numeric(length(t))
  
  # Divide time points before and after cut point.
  times_before_cut_point <- which(t <= cut_point)
  times_after_cut_point <- which(t > cut_point)

  if(length(times_before_cut_point) > 0){
    for(i in 1:length(times_before_cut_point)){
      if(t[times_before_cut_point[i]] > smallest_survival){
     # print(t[times_before_cut_point[i]])  
      km_summary <- summary(km_fit, rmean=t[times_before_cut_point[i]])
      #print(km_summary$table["rmean"] %>% unname())
      #print(km_summary$table["se(rmean)"] %>% unname())
      RMST[times_before_cut_point[i]] <- km_summary$table["rmean"] %>% unname()
      RMST_SE[times_before_cut_point[i]] <- km_summary$table["se(rmean)"] %>% unname()
      } else {
      RMST[times_before_cut_point[i]] <- t[times_before_cut_point[i]]
      RMST_SE[times_before_cut_point[i]] <- 0
      }
    }
  }  
  
  RMST_Lower[times_before_cut_point] <- RMST[times_before_cut_point] - RMST_SE[times_before_cut_point]
  RMST_Upper[times_before_cut_point] <- RMST[times_before_cut_point] + RMST_SE[times_before_cut_point]
  
  # After cut point use both Kaplan-Meier and exponential.
  
  
  # Get Kaplan-Meier estimate of RMST up to cut point.
  km_summary_cut_point <- summary(km_fit, rmean=cut_point)
  RMST_cut_point <- km_summary_cut_point$table["rmean"] %>% unname()
  RMST_se_cut_point <- km_summary_cut_point$table["se(rmean)"] %>% unname()
  # print(RMST_cut_point)
  # print(RMST_se_cut_point)
  
  # Then get closed form of conditional RMST from cut point.
  # and its standard error
  lambda <- exp_fit$res[,"est"]

  cut_point_RMST_from_exp <- summary(exp_fit, 
          t = cut_point,
          type = "rmst",
          se = T)[[1]] %>% 
    as_tibble() %>%
    rename(value = est, lower = lcl, upper = ucl)

  
  RMST_from_exp <- summary(exp_fit, 
                           t = t[times_after_cut_point],
                           type = "rmst",
                           se = T)[[1]] %>% 
    as_tibble() %>%
    rename(value = est, lower = lcl, upper = ucl)  %>%
    select(-c(lower, upper)) %>%
    mutate(value= value-cut_point_RMST_from_exp$value,
           se = sqrt(se^2-cut_point_RMST_from_exp$se^2))
  
  # Break variance into components.
  # Var_term1 <- (1/lambda^2)*(exp(-lambda*tc)*exp(-2*lambda*ta)-exp(-2*lambda*tc)*exp(-lambda*ta))
  # Var_term2 <- (1/lambda^2)*(-(1/3)*exp(-3*lambda*ta)+(1/3)*exp(-3*lambda*tc))
  # Var_term3 <- (1/lambda^2)*(-(3/2)*exp(-2*lambda*tc)*exp(-2*lambda*ta)-(1/4)*exp(-4*lambda*tc)-(1/4)*exp(-4*lambda*ta))
  # Var_term4 <- (1/lambda^2)*(exp(-lambda*tc)*exp(-3*lambda*ta)+exp(-3*lambda*tc)*exp(-lambda*ta))
  
  Var_RMST_at_cut_point <- RMST_se_cut_point^2
  # Var_RMST_after_cut_point <- Var_term1 + Var_term2 + Var_term3 + Var_term4

  RMST[times_after_cut_point] <- RMST_cut_point+RMST_from_exp$value
  RMST_SE[times_after_cut_point] <- sqrt(Var_RMST_at_cut_point+RMST_from_exp$se^2)
  
  RMST_Lower[times_after_cut_point] <- RMST[times_after_cut_point] - 1.96*RMST_SE[times_after_cut_point]
  RMST_Upper[times_after_cut_point] <- RMST[times_after_cut_point] + 1.96*RMST_SE[times_after_cut_point]
  
  results <- tibble(
    time = t,
    value = RMST, 
    se = RMST_SE,
    lower = RMST_Lower,
    upper = RMST_Upper)
  
  return(results)
}

