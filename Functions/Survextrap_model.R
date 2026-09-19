# Model fitting and estimate extraction for the survextrap analysis.
#
# Priors come from Choose_priors.R and Check_prior_packages.R as four
# numbers carried in the scenario table, so a scenario is reproducible
# from its row. Background mortality is supplied as a time/hazard data
# frame, which survextrap calls form (a), so predictions describe total
# survival. The smoothing model is the package default random walk, set
# explicitly because the prior summary functions default otherwise.
#
# Knots vary by data cut and by evidence set: with no external data the
# boundary sits at the end of the trial, since nothing informs a knot
# beyond it. They travel as the name of a global vector, a numeric
# vector not being something a table column can hold.
#
# RMST is computed at 10 and 30 years, undiscounted and discounted at
# 3.5 per cent, and the survival grid runs to 30, matching TA536.

library(memoise)
library(survextrap)

######################################################
#  Memoise survextrap.
######################################################
#
# The cache directory is relative to the project, so this works on a
# laptop as well as on the cluster. Set survextrap_cache before sourcing
# to override.

if (!exists("survextrap_cache"))
  survextrap_cache <- file.path(getwd(), "cache")
dir.create(survextrap_cache, showWarnings = FALSE, recursive = TRUE)

survextrap_mem <- memoise(survextrap,
                          cache = cachem::cache_disk(
                            dir = survextrap_cache,
                            max_age = 31557600,   # keep for a year
                            max_size = 2e10))     # 20Gb


##################################################
# Fit survextrap models to two arms.
##################################################
#
# Paths are absolute, built by the driver, so they resolve the same here
# or inside an rslurm working directory. external_file is NA for the
# trial only scenarios; tau_rate is used only by NON-PH.

fit_model <- function(model, datasets, df, add_knots,
                      hsd_rate, eta_med, eta_up, tau_rate,
                      fit_method, trial_file, external_file,
                      backhaz_file, store_file) {
  
  trial_data <- readRDS(trial_file)
  backhaz <- readRDS(backhaz_file)
  external_data <- if (is.na(external_file)) NULL else readRDS(external_file)
  
  add_knots <- get(add_knots)
  
  mspline <- mspline_spec(Surv(time, status) ~ trt, data = trial_data,
                          df = df, add_knots = add_knots)
  
  prior_hscale <- p_meansurv(median = eta_med, upper = eta_up,
                             mspline = mspline)
  prior_hsd <- p_gamma(2, hsd_rate)
  prior_loghr <- p_hr(median = 1, upper = 2)
  smooth_model <- "random_walk"
  
  if (model == "PH") {
    
    results <- survextrap_mem(Surv(time, status) ~ trt,
                              data = trial_data,
                              external = external_data,
                              mspline = mspline,
                              prior_hscale = prior_hscale,
                              prior_loghr = prior_loghr,
                              prior_hsd = prior_hsd,
                              backhaz = backhaz,
                              smooth_model = smooth_model,
                              fit_method = fit_method)
    
  } else if (model == "NON-PH") {
    
    results <- survextrap_mem(Surv(time, status) ~ trt,
                              data = trial_data,
                              nonprop = ~ trt,
                              external = external_data,
                              mspline = mspline,
                              prior_hscale = prior_hscale,
                              prior_loghr = prior_loghr,
                              prior_hsd = prior_hsd,
                              prior_hrsd = p_gamma(2, tau_rate),
                              backhaz = backhaz,
                              smooth_model = smooth_model,
                              fit_method = fit_method)
    
  } else if (model == "Separate_arms") {
    
    # the control arm gets the external data, since PROFILE-1014 and
    # Flatiron are both crizotinib; the active arm has none, which is
    # what makes this specification behave badly and is the point of
    # showing it. Both arms use the same spline, so they differ only in
    # the evidence.
    results_control <- survextrap_mem(Surv(time, status) ~ 1,
                                      data = trial_data %>%
                                        filter(trt == "Crizotinib"),
                                      external = external_data,
                                      mspline = mspline,
                                      prior_hscale = prior_hscale,
                                      prior_hsd = prior_hsd,
                                      backhaz = backhaz,
                                      smooth_model = smooth_model,
                                      fit_method = fit_method)
    
    results_active <- survextrap_mem(Surv(time, status) ~ 1,
                                     data = trial_data %>%
                                       filter(trt == "Alectinib"),
                                     mspline = mspline,
                                     prior_hscale = prior_hscale,
                                     prior_hsd = prior_hsd,
                                     backhaz = backhaz,
                                     smooth_model = smooth_model,
                                     fit_method = fit_method)
    
    results <- list(control = results_control, active = results_active)
    class(results) <- "two_models"
  }
  
  saveRDS(results, store_file)
  invisible(store_file)
}


##################################################
# Survival and hazard estimates from models.
##################################################

get_survival_and_hazard_survextrap <- function(model_file, store_file,
                                               t_max = 30) {
  
  model <- readRDS(model_file)
  
  new_data_estimate <- tibble(trt = c("Alectinib", "Crizotinib"))
  
  # dense over the trial period, coarser over the extrapolation
  time_vec <- c(seq(from = 0, to = 5, length.out = 1e2),
                seq(from = 5, to = t_max, length.out = 1.5e2))
  
  if (inherits(model, "survextrap")) {
    
    results <- bind_rows(
      hazard(model, t = time_vec, newdata = new_data_estimate),
      survival(model, t = time_vec, newdata = new_data_estimate))
    
  } else if (inherits(model, "two_models")) {
    
    results <- bind_rows(
      hazard(model$active, t = time_vec) %>% mutate(trt = "Alectinib"),
      hazard(model$control, t = time_vec) %>% mutate(trt = "Crizotinib"),
      survival(model$active, t = time_vec) %>% mutate(trt = "Alectinib"),
      survival(model$control, t = time_vec) %>% mutate(trt = "Crizotinib")) %>%
      select(variable, trt, t, median, lower, upper)
  }
  
  saveRDS(results, store_file)
  invisible(store_file)
}


##################################################
# rmst and difference in rmst, undiscounted and discounted.
##################################################

# Does the installed survextrap take a discount rate directly.
takes_disc <- function(f) "disc_rate" %in% names(formals(f))

# Summarise a matrix of samples the way survextrap does, so the
# fallback route returns the same shape as irmst().
summarise_samples <- function(samples, time_vec) {
  survextrap:::summarise_output(
    samples, t = time_vec,
    summ_fns = list("median" = median,
                    ~quantile(.x, probs = c(0.025, 0.975))),
    newdata = NULL, summ_name = "irmst", sample = FALSE) %>%
    rename(lower = "2.5%", upper = "97.5%") %>%
    mutate(trt = NA) %>%
    select(variable, trt, t, median, lower, upper)
}


get_rmst_survextrap <- function(model_file, store_file,
                                times = c(10, 30),
                                disc_rates = c(0, 0.035)) {
  
  model <- readRDS(model_file)
  
  new_data_estimate <- tibble(trt = c("Alectinib", "Crizotinib"))
  
  one_rate <- function(dr) {
    
    if (inherits(model, "survextrap")) {
      
      rmst_both <- rmst(model, t = times, newdata = new_data_estimate,
                        disc_rate = dr)
      
      if (takes_disc(irmst)) {
        irmst_diff <- irmst(model, t = times,
                            newdata = tibble(trt = c("Crizotinib", "Alectinib")),
                            disc_rate = dr) %>%
          mutate(trt = NA) %>%
          select(variable, trt, t, median, lower, upper)
      } else {
        # difference the discounted samples arm by arm
        s <- rmst(model, t = times, newdata = new_data_estimate,
                  disc_rate = dr, sample = TRUE)
        act <- s[, which(new_data_estimate$trt == "Alectinib"), , drop = TRUE]
        con <- s[, which(new_data_estimate$trt == "Crizotinib"), , drop = TRUE]
        irmst_diff <- summarise_samples(act - con, times)
      }
      
      out <- bind_rows(rmst_both %>%
                         select(variable, trt, t, median, lower, upper),
                       irmst_diff)
      
    } else if (inherits(model, "two_models")) {
      
      rmst_control <- rmst(model$control, t = times, disc_rate = dr) %>%
        mutate(trt = "Crizotinib")
      rmst_active <- rmst(model$active, t = times, disc_rate = dr) %>%
        mutate(trt = "Alectinib")
      
      irmst_diff <- summarise_samples(
        rmst(model$active, t = times, disc_rate = dr, sample = TRUE) -
          rmst(model$control, t = times, disc_rate = dr, sample = TRUE),
        times)
      
      out <- bind_rows(rmst_active, rmst_control) %>%
        select(variable, trt, t, median, lower, upper) %>%
        bind_rows(irmst_diff)
    }
    
    mutate(out, disc_rate = dr)
  }
  
  results <- bind_rows(lapply(disc_rates, one_rate))
  
  saveRDS(results, store_file)
  invisible(store_file)
}


##################################################
# Convergence diagnostics from models.
##################################################
#
# One small file per model, so this can be run over the grid with slurm
# rather than by opening every fitted model in one session. Returns the
# worst rhat, the smallest effective sample size and the number of
# divergent transitions. A separate arms model gives one row per arm.

get_diagnostics_survextrap <- function(model_file, store_file){
  
  model <- readRDS(model_file)
  
  if(inherits(model, "two_models")){
    fits <- list(control = model$control$stanfit,
                 active = model$active$stanfit)
  } else {
    fits <- list(model = model$stanfit)
  }
  
  results <- NULL
  
  for(arm in names(fits)){
    
    fit <- fits[[arm]]
    
    temp <- tryCatch({
      
      summ <- rstan::summary(fit)$summary
      sampler <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
      
      divergences <- 0
      for(chain in sampler){
        divergences <- divergences + sum(chain[, "divergent__"])
      }
      
      tibble(arm = arm,
             max_rhat = max(summ[, "Rhat"], na.rm = TRUE),
             min_ess = min(summ[, "n_eff"], na.rm = TRUE),
             divergences = divergences)
      
    },
    error = function(e)
      tibble(arm = arm, max_rhat = NA_real_, min_ess = NA_real_,
             divergences = NA_real_))
    
    results <- bind_rows(results, temp)
  }
  
  saveRDS(results, store_file)
  invisible(store_file)
}


##################################################
# further helper functions.
##################################################

get_legend_35 <- function(plot, legend_number = 1) {
  legends <- get_plot_component(plot, "guide-box", return_all = TRUE)
  idx <- which(vapply(legends, \(x) !inherits(x, "zeroGrob"), TRUE))
  if (length(idx) >= legend_number) {
    return(legends[[idx[legend_number]]])
  } else if (length(idx) >= 0) {
    return(legends[[idx[1]]])
  } else {
    return(legends[[1]])
  }
}

# the grid no longer sweeps priors, so a model is picked by its name and
# its scenario label rather than by prior values
filter_model <- function(.data, model_name = "PH",
                         scenario_value = "primary") {
  .data %>% filter(model == model_name, scenario == scenario_value)
}