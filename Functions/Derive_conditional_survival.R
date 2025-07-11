
# Function to create aggregate counts for conditional survival estimates.

create_aggregate_counts <- function(
    data = NULL,
    increment = 1,
    left_truncate = 0){
  
  min_t <- max(0, left_truncate)
  max_t <- ceiling(max(data$time))
  
  data_survSplit <- survSplit(Surv(time,status) ~.,
                              data,
                              cut=seq(from = min_t, to = max_t, by = increment), 
                              episode ="timegroup")
  
  data_aggregate <- 
    data_survSplit %>%
    as_tibble() %>%
    rename(start = tstart) %>%
    group_by(start, trt) %>%
    summarise(n = sum(weight),
              r = sum(weight)-sum(weight*status),
              .groups = "keep") %>%
    mutate(hazard = -log(r/n)/increment) %>% # average hazard.
    ungroup() %>%
    group_by(trt)  %>%
    mutate(stop = lead(start, default = max_t)) %>%
    relocate(c(hazard, start, stop, trt), .after = r) %>%
    filter(start >= min_t) %>%
    ungroup()  %>%
    mutate(n = ceiling(n), r = round(r))
  
  return(data_aggregate)
  
}


gg_conditional_surv_weight <- function (basekm, at, main = NULL, xlab = "Years", ylab = "Survival probability", 
          lwd = 1) 
{
  if (class(basekm) != "survfit") {
    stop("Argument to basekm must be of class survfit")
  }
  if (max(at) > max(basekm$time)) {
    stop(paste("Argument to at specifies value(s) outside the range of observed times;", 
               "the maximum observed time is", round(max(basekm$time), 
                                                     2)))
  }
  nt <- length(at)
  fitkm <- list()
  fitkmdat <- list()
  for (i in 1:nt) {
    fitkm[[i]] <- survival::survfit(formula = stats::as.formula(basekm$call$formula), 
                                    data = eval(basekm$call$data), 
                                    weights = eval(basekm$call$weights),
                                    start.time = at[i])
    fitkmdat[[i]] <- tibble::tibble(timept = fitkm[[i]]$time, 
                                    prob = fitkm[[i]]$surv)
  }
  condsurvdat <- fitkmdat %>% purrr::map_df(`[`, .id = "which_at") %>% 
    dplyr::mutate(condtime = factor(which_at, levels = seq(1, 
                                                           nt), labels = at))
  ggplot2::ggplot(condsurvdat, ggplot2::aes(x = timept, y = prob, 
                                            color = condtime)) + ggplot2::geom_step(lwd = lwd) + 
    ggplot2::ylim(0, 1) + ggplot2::labs(x = xlab, y = ylab, 
                                        title = main) + ggplot2::labs(color = "x") + ggplot2::theme_bw()
}


