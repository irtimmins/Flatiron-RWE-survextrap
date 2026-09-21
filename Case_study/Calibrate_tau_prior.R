# Calibrating the prior for tau, the flexibility of the hazard ratio in
# the non-proportional hazards model.
#
# What tau controls, and what it does not. In the non-proportional
# hazards model the spline coefficients for the treated arm are
# perturbed by delta ~ Normal(0, tau) in each region of time, so tau
# governs how far the hazard ratio departs from a constant. At tau = 0
# the model is exactly proportional hazards. Where the hazard ratio sits
# is a separate question, set by prior_loghr; tau is only about how much
# it moves. A judgement of the form "the hazard ratio might be as low as
# 0.5 and as high as 1.5" is therefore about variation, and the quantity
# that expresses it is the ratio between the two, 1.5 / 0.5 = 3.
#
# That ratio is what prior_hr_sd returns as hrr: the high over low value
# of the hazard ratio across time, at the 90th and 10th percentiles by
# default. The vignette says trial and error is needed to find a Gamma
# prior matching a given judgement, so this sweeps the rate and reports
# which value lands nearest the target rather than trying them one at a
# time.
#
# Everything else has to be supplied as it is in the real fit, since the
# simulation draws from the joint prior: change the mspline, the hazard
# scale prior or the sigma prior and the answer for tau changes with it.

library(dplyr)
library(survextrap)
library(survival)

set.seed(1234)

# The judgement being calibrated against. Change these two numbers and
# the target ratio follows.
hr_low <- 0.5
hr_high <- 1.5
target_hrr <- hr_high / hr_low

# Rates to try for the Gamma(2, rate) prior on tau. Larger rate means a
# smaller tau on average, so a stiffer hazard ratio. The range runs well
# past the values in the scenario grid because the vague end is vaguer
# than it looks: p_gamma(2,1) gives a median high over low hazard ratio
# of about 16, and an upper limit in the millions.
rates <- c(1, 2, 3, 5, 7, 10, 15, 20, 30, 40, 60, 80)

nsim <- 4000
cut_wanted <- "feb_2017"

cat("target: hazard ratio varying between", hr_low, "and", hr_high,
    "over time,\n        a high over low ratio of", target_hrr, "\n\n")


######################################################
#  Rebuild the spline and the other priors, as fitted.
######################################################
#
# These must match Fit_survextrap.R, or the calibration is for a model
# that is not the one being fitted.

trial_data <- readRDS("Data/trial_data.rds")

df_value <- 6
add_knots <- c(2, 3, 5)          # add_knots_feb
hsd_rate <- 1                    # the sigma prior in the base case

sp <- mspline_spec(Surv(time, status) ~ trt, data = trial_data,
                   df = df_value, add_knots = add_knots)

prior_hscale <- p_meansurv(median = 5, upper = 20, mspline = sp)
prior_hsd <- p_gamma(2, hsd_rate)
prior_loghr <- p_hr(median = 1, upper = 10)

cat("spline knots:", paste(round(sp$knots, 2), collapse = ", "), "\n")

# The boundary knot is the last element of knots, not a field of its own;
# bsmooth is a logical flag for the smoothness constraint, not a knot.
boundary_knot <- max(sp$knots)
cat("boundary knot:", round(boundary_knot, 2), "years\n")
cat("  beyond this the hazard is held constant, so with knots at",
    paste(add_knots, collapse = ", "), "the model\n",
    "  assumes a flat hazard from", round(boundary_knot, 1),
    "years to the end of the horizon. That is what adding\n",
    "  a knot at 10 years changes, and it is worth keeping in view while",
    "reading\n  the calibration below: tau governs how the hazard ratio",
    "moves within the\n  splined region, not beyond it.\n\n")


######################################################
#  Sweep the rate.
######################################################

# Both arms have to carry the full set of factor levels. Passing a bare
# string gives a factor with one level, and model.matrix cannot build
# contrasts from that: it is the error this hit first. The vignette's
# example uses a numeric 0/1 covariate and so does not run into it.
trt_levels <- levels(trial_data$trt)
if (is.null(trt_levels)) trt_levels <- sort(unique(as.character(trial_data$trt)))
cat("treatment levels, reference first:",
    paste(trt_levels, collapse = ", "), "\n\n")

as_trt <- function(x) factor(x, levels = trt_levels)

sweep <- bind_rows(lapply(rates, function(r) {
  
  s <- prior_hr_sd(mspline = sp,
                   prior_hsd = prior_hsd,
                   prior_hscale = prior_hscale,
                   prior_loghr = prior_loghr,
                   prior_hrsd = p_gamma(2, r),
                   # tau only exists in the non-proportional hazards
                   # model, and this argument says which covariates have
                   # one. It takes a formula, not TRUE: passing a logical
                   # gives "$ operator is invalid for atomic vectors",
                   # the code having tried to read terms off it. Left
                   # off altogether it defaults to proportional hazards,
                   # where prior_hrsd has nothing to act on and the sweep
                   # returns the same numbers at every rate.
                   nonprop = ~ trt,
                   formula = ~ trt,
                   newdata = data.frame(trt = as_trt("Alectinib")),
                   newdata0 = data.frame(trt = as_trt("Crizotinib")),
                   quantiles = c(0.025, 0.5, 0.975),
                   nsim = nsim)
  
  # prior_hr_sd returns quantities as columns (sd_hr, hrr) and quantiles
  # as rows, named "2.5%", "50%" and so on. hrr is the high over low
  # hazard ratio across time, which is the quantity being calibrated.
  if (r == rates[1]) {
    cat("prior_hr_sd returns columns:", paste(colnames(s), collapse = ", "),
        "\n                    rows:", paste(rownames(s), collapse = ", "),
        "\n\n")
  }
  
  if (!"hrr" %in% colnames(s)) {
    print(s)
    stop("no hrr column in prior_hr_sd's output")
  }
  
  q <- function(want) {
    hit <- match(want, rownames(s))
    if (is.na(hit)) {
      print(s)
      stop("no row called '", want, "' in prior_hr_sd's output; rows are ",
           paste(rownames(s), collapse = ", "))
    }
    s[hit, "hrr"]
  }
  
  tibble(rate = r,
         hrr_median = q("50%"),
         hrr_lower = q("2.5%"),
         hrr_upper = q("97.5%"))
}))

sweep <- sweep %>%
  mutate(distance = abs(hrr_median - target_hrr))

# The rate must actually change the answer. If the medians are flat
# across an eightyfold range then tau is not entering the simulation and
# the sweep is comparing one prior with itself, which looks like a
# result but is not one. Better to stop than to report the noisiest draw
# as a calibration.
spread <- diff(range(sweep$hrr_median)) / min(sweep$hrr_median)
if (spread < 0.5) {
  print(as.data.frame(sweep), row.names = FALSE)
  stop("the median high over low hazard ratio barely changes across ",
       "rates ", min(rates), " to ", max(rates), ", so prior_hrsd is not ",
       "reaching the simulation. Check nonprop is passed as a formula, ",
       "~ trt rather than TRUE, and that prior_hrsd matches the argument ",
       "name in formals(prior_hr_sd).")
}

cat("prior distribution of the high over low hazard ratio, by rate\n")

# the upper tail runs to very large values at the vague end, so it is
# formatted rather than printed at full width
fmt_hrr <- function(x)
  ifelse(x >= 1000, format(signif(x, 2), scientific = TRUE),
         sprintf("%.2f", x))

print(as.data.frame(sweep %>%
                      transmute(rate,
                                median = fmt_hrr(hrr_median),
                                interval = paste0("(", fmt_hrr(hrr_lower), ", ",
                                                  fmt_hrr(hrr_upper), ")"),
                                from_target = round(distance, 2))), row.names = FALSE)

best <- sweep %>% slice_min(distance, n = 1, with_ties = FALSE)

cat("\nclosest to a median ratio of", target_hrr, "is p_gamma(2,",
    best$rate, ")\n")
cat("  median ", fmt_hrr(best$hrr_median),
    ", 95% interval (", fmt_hrr(best$hrr_lower), ", ",
    fmt_hrr(best$hrr_upper), ")\n", sep = "")

cat("\nRead the median as the ratio this prior considers typical and the\n")
cat("interval as the range it allows. A prior centred exactly on the\n")
cat("target still permits a good deal more and a good deal less, which\n")
cat("is intended: it is a prior, not a constraint. If the upper end looks\n")
cat("too permissive, take the next rate up and accept a median below the\n")
cat("target.\n")


######################################################
#  What the chosen rate implies, alongside the ones in use.
######################################################

cat("\nthe rates currently in the scenario grid, for comparison\n")
print(as.data.frame(sweep %>%
                      filter(rate %in% c(1, 5, 10)) %>%
                      transmute(rate,
                                median = fmt_hrr(hrr_median),
                                interval = paste0("(", fmt_hrr(hrr_lower), ", ",
                                                  fmt_hrr(hrr_upper), ")"))),
      row.names = FALSE)

cat("\nIf the calibrated rate is not among these, it is worth adding to\n")
cat("hrsd_rate in Fit_survextrap.R rather than picking the nearest one\n")
cat("already there.\n")