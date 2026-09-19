# Calibrating the priors against what they imply for survival, with
# background mortality included.
#
# The earlier calibration used the package's summary functions directly,
# which simulate the M-spline alone. With an additive background hazard
# in the fitted model that is the excess hazard, and calibrating without
# the background is calibrating a different model. The fix is small: the
# package exports prior_sample_hazard, which returns the excess hazard
# for every draw at every time, and background mortality is added to
# that before anything is summarised. The hazard itself comes from the
# package's own code, not a reconstruction.
#
# The judgements are also different in kind. Rather than a ratio of high
# to low hazard, which is hard to hold an opinion about, they are stated
# as what fraction of this cohort could plausibly be alive at 5, 10, 20
# and 30 years, with limits on both sides: implausibly good scenarios,
# such as a third alive at 30 years, and implausibly bad ones, such as
# hazards several times what advanced cancer produces. Ten year survival
# is left able to be poor, since for this population it could be.
#
# eta and sigma are swept together, since what a given scale implies
# depends on how much the hazard is allowed to move. tau is then
# calibrated with both in force, as before.
#
# The trial's own knots and the extra ones at 2, 3, 5, 10 and 20 are
# taken as settled.

library(dplyr)
library(tidyr)
library(survextrap)
library(survival)

set.seed(1)


######################################################
#  The judgements.
######################################################
#
# For each landmark, the range that most of the prior should fall in:
# the 5th percentile of prior draws should be at least surv_low, and the
# 95th at most surv_high. A value of 0 for surv_low means no lower
# limit. haz_high is the most the total hazard should plausibly be at
# that time, again at the 95th percentile.
#
# These are starting points to argue with, not findings.

landmarks <- tribble(
  ~time, ~surv_low, ~surv_high, ~haz_high,
  5,      0.15,       0.85,      0.50,
  10,      0.03,       0.60,      0.40,
  20,      0.00,       0.30,      0.40,
  30,      0.00,       0.20,      0.40)

judgement_text <- c(
  "five year survival is plausibly between 15 and 85 per cent",
  "ten year survival could be poor, but not above 60 per cent",
  "fewer than 30 per cent alive at 20 years, fewer than 20 at 30",
  "the total hazard does not plausibly exceed 0.4 to 0.5 per year")

tau_target <- c(low = 1.5, high = 2.5)

nsim <- 1000
horizon <- 30


######################################################
#  Basis, background, and the fixed pieces.
######################################################

trial_data <- readRDS("Data/trial_data.rds")
trt_levels <- levels(trial_data$trt)
as_trt <- function(x) factor(x, levels = trt_levels)

sp <- mspline_spec(Surv(time, status) ~ trt, data = trial_data,
                   df = 6, add_knots = c(2, 3, 5, 10, 20))

if (!file.exists("Data/backhaz.rds"))
  stop("no Data/backhaz.rds; run Build_background_mortality.R first")
backhaz <- readRDS("Data/backhaz.rds")

cat("knots:", paste(round(sp$knots, 2), collapse = ", "), "\n")
cat("background hazard at 0, 10, 20, 30 years:",
    paste(sprintf("%.3f", approx(backhaz$time, backhaz$hazard,
                                 xout = c(0, 10, 20, 30), rule = 2,
                                 method = "constant")$y), collapse = ", "),
    "\n\n")

prior_loghr <- p_hr(median = 1, upper = 10)
criz <- data.frame(trt = as_trt("Crizotinib"))
alec <- data.frame(trt = as_trt("Alectinib"))

# the same smoothing structure survextrap fits with; prior_haz_sd's
# default is this, prior_sample's is not, so it is set explicitly
smooth_model <- "exchangeable"


######################################################
#  Drawing hazard curves with the background added.
######################################################
#
# prior_sample_hazard gives the excess hazard per draw per time. The
# background is added, survival integrated, and both returned per draw.

draw_total <- function(prior_hsd, prior_hscale, prior_hrsd = NULL,
                       nonprop = NULL, newdata = criz, newdata0 = NULL,
                       nsim = 1000) {
  
  pred <- prior_sample_hazard(
    knots = sp$knots, degree = sp$degree, bsmooth = sp$bsmooth,
    prior_hsd = prior_hsd, prior_hscale = prior_hscale,
    smooth_model = smooth_model, prior_loghr = prior_loghr,
    formula = ~ trt, nonprop = nonprop,
    newdata = newdata, newdata0 = newdata0, prior_hrsd = prior_hrsd,
    tmin = 0, tmax = horizon, nsim = nsim)
  
  # the columns are found rather than assumed; time is whichever numeric
  # column is not rep, haz or hr
  if (!all(c("rep", "haz") %in% names(pred))) {
    print(utils::head(pred))
    stop("prior_sample_hazard's output has no rep and haz columns; ",
         "it is printed above")
  }
  time_col <- setdiff(names(pred)[vapply(pred, is.numeric, logical(1))],
                      c("rep", "haz", "hr", "haz0"))
  if (!length(time_col)) {
    print(utils::head(pred))
    stop("cannot find a time column in prior_sample_hazard's output")
  }
  pred$time <- pred[[time_col[1]]]
  
  pred %>%
    mutate(backhaz = approx(backhaz$time, backhaz$hazard, xout = time,
                            rule = 2, method = "constant")$y,
           total = haz + backhaz) %>%
    group_by(rep) %>%
    arrange(time, .by_group = TRUE) %>%
    mutate(dt = c(time[1], diff(time)),
           surv = exp(-cumsum(total * dt))) %>%
    ungroup()
}

# value at a landmark, per draw, by nearest time
at_landmark <- function(d, tt, what) {
  d %>%
    group_by(rep) %>%
    slice_min(abs(time - tt), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    pull({{ what }})
}

# does this prior meet the landmark judgements
assess <- function(d) {
  bind_rows(lapply(seq_len(nrow(landmarks)), function(i) {
    tt <- landmarks$time[i]
    s <- at_landmark(d, tt, surv)
    h <- at_landmark(d, tt, total)
    tibble(time = tt,
           surv_5 = quantile(s, 0.05), surv_50 = median(s),
           surv_95 = quantile(s, 0.95),
           haz_50 = median(h), haz_95 = quantile(h, 0.95),
           ok_low = surv_5 >= landmarks$surv_low[i],
           ok_high = surv_95 <= landmarks$surv_high[i],
           ok_haz = haz_95 <= landmarks$haz_high[i])
  }))
}


######################################################
#  Step 1: eta and sigma together.
######################################################

cat("step 1, eta and sigma\n")
for (j in judgement_text) cat("  ", j, "\n", sep = "")
cat("\n")

grid <- expand_grid(
  sigma_mean = c(0.6, 0.8, 1.0),
  eta_median = c(2, 3, 4, 5, 6),
  eta_upper_ratio = c(1.5, 2, 3))

sigma_shape <- 10

sweep <- bind_rows(lapply(seq_len(nrow(grid)), function(i) {
  g <- grid[i, ]
  d <- draw_total(
    prior_hsd = p_gamma(sigma_shape, sigma_shape / g$sigma_mean),
    prior_hscale = p_meansurv(median = g$eta_median,
                              upper = g$eta_median * g$eta_upper_ratio,
                              mspline = sp),
    nsim = nsim)
  a <- assess(d)
  tibble(g,
         meets = all(a$ok_low & a$ok_high & a$ok_haz),
         n_fail = sum(!(a$ok_low & a$ok_high & a$ok_haz)),
         s10 = a$surv_50[a$time == 10], s10_lo = a$surv_5[a$time == 10],
         s10_hi = a$surv_95[a$time == 10],
         s20_hi = a$surv_95[a$time == 20],
         s30_hi = a$surv_95[a$time == 30],
         h10_hi = a$haz_95[a$time == 10],
         h20_hi = a$haz_95[a$time == 20])
}))

cat("prior survival and hazard at landmarks, with background added\n")
cat("s10 is median (5th to 95th) at ten years; the rest are 95th\n\n")
print(as.data.frame(sweep %>%
                      transmute(sigma = sprintf("Gamma(%g, %.1f)", sigma_shape,
                                                sigma_shape / sigma_mean),
                                eta = sprintf("meansurv(%g, %g)", eta_median,
                                              eta_median * eta_upper_ratio),
                                s10 = sprintf("%.0f%% (%.0f-%.0f)", 100 * s10, 100 * s10_lo,
                                              100 * s10_hi),
                                s20_95 = sprintf("%.0f%%", 100 * s20_hi),
                                s30_95 = sprintf("%.0f%%", 100 * s30_hi),
                                h10_95 = sprintf("%.2f", h10_hi),
                                h20_95 = sprintf("%.2f", h20_hi),
                                fails = n_fail,
                                meets = ifelse(meets, "yes", ""))), row.names = FALSE)

if (any(sweep$meets)) {
  # among those meeting every judgement, prefer the vaguest: the widest
  # eta and the largest sigma, since a prior should not be tighter than
  # the judgement requires
  pick <- sweep %>% filter(meets) %>%
    arrange(desc(eta_upper_ratio), desc(sigma_mean)) %>% slice(1)
  cat("\n  ", sum(sweep$meets), " combinations meet every judgement;",
      " taking the vaguest\n", sep = "")
} else {
  pick <- sweep %>% arrange(n_fail, desc(eta_upper_ratio)) %>% slice(1)
  cat("\n  nothing meets every judgement; taking the fewest failures.\n")
  cat("  Which landmark fails most often:\n")
  fails <- bind_rows(lapply(seq_len(nrow(grid)), function(i) {
    g <- grid[i, ]
    d <- draw_total(
      prior_hsd = p_gamma(sigma_shape, sigma_shape / g$sigma_mean),
      prior_hscale = p_meansurv(median = g$eta_median,
                                upper = g$eta_median * g$eta_upper_ratio,
                                mspline = sp), nsim = 300)
    assess(d) %>% mutate(i = i)
  }))
  print(as.data.frame(fails %>%
                        group_by(time) %>%
                        summarise(too_low = mean(!ok_low), too_high = mean(!ok_high),
                                  hazard_too_high = mean(!ok_haz), .groups = "drop") %>%
                        mutate(across(-time, ~sprintf("%.0f%%", 100 * .x)))),
        row.names = FALSE)
}

chosen_hsd <- p_gamma(sigma_shape, sigma_shape / pick$sigma_mean)
chosen_hscale <- p_meansurv(median = pick$eta_median,
                            upper = pick$eta_median * pick$eta_upper_ratio,
                            mspline = sp)

cat("\n  chosen: sigma Gamma(", sigma_shape, ", ",
    round(sigma_shape / pick$sigma_mean, 2), "), eta p_meansurv(median = ",
    pick$eta_median, ", upper = ", pick$eta_median * pick$eta_upper_ratio,
    ")\n", sep = "")

final <- draw_total(chosen_hsd, chosen_hscale, nsim = 2000)
cat("\n  what it implies, crizotinib, with background\n")
print(as.data.frame(assess(final) %>%
                      transmute(time,
                                survival = sprintf("%.0f%% (%.0f-%.0f)", 100 * surv_50,
                                                   100 * surv_5, 100 * surv_95),
                                hazard = sprintf("%.3f (95th %.3f)", haz_50, haz_95))),
      row.names = FALSE)


######################################################
#  Step 2: tau.
######################################################
#
# On the total hazard ratio, which is the one a clinician would
# recognise. tau acts on the excess; as background grows the total
# ratio is pulled towards 1 whatever tau does, so this is calibrated
# over the trial-relevant period rather than the whole horizon.

cat("\nstep 2, tau: the effect may attenuate towards nothing but is",
    "unlikely to reverse\n\n")

tau_out <- bind_rows(lapply(c(1, 2, 3, 5, 10), function(r) {
  d <- draw_total(chosen_hsd, chosen_hscale,
                  prior_hrsd = p_gamma(2, r), nonprop = ~ trt,
                  newdata = alec, newdata0 = criz, nsim = nsim)
  if (!"hr" %in% names(d)) stop("no hr column with newdata0 supplied")
  # total hazard ratio needs the comparator's total too; the package
  # returns the excess ratio as hr, so the total is rebuilt from it:
  # excess for the comparator is haz / hr
  d <- d %>% mutate(haz0 = haz / hr, total0 = haz0 + backhaz,
                    hr_total = total / total0) %>%
    filter(time <= 10)
  hrr <- d %>% group_by(rep) %>%
    summarise(r = quantile(hr_total, 0.9) / quantile(hr_total, 0.1),
              above1 = mean(hr_total > 1), .groups = "drop")
  tibble(rate = r, ratio = median(hrr$r), upper = quantile(hrr$r, 0.975),
         favours_comparator = mean(hrr$above1))
}))

print(as.data.frame(tau_out %>%
                      transmute(prior = sprintf("Gamma(2, %g)", rate),
                                hr_varies_by = sprintf("%.2f", ratio),
                                hr_range = sprintf("%.2f to %.2f", 1 / sqrt(ratio), sqrt(ratio)),
                                upper_975 = sprintf("%.0f", upper),
                                favours_comparator = sprintf("%.0f%%", 100 * favours_comparator),
                                meets = ifelse(ratio >= tau_target["low"] &
                                                 ratio <= tau_target["high"], "yes", ""))),
      row.names = FALSE)

tau_pick <- tau_out %>%
  filter(ratio >= tau_target["low"], ratio <= tau_target["high"]) %>%
  slice_min(upper, n = 1, with_ties = FALSE)
if (!nrow(tau_pick))
  tau_pick <- tau_out %>% mutate(d = abs(ratio - 2)) %>% slice_min(d, n = 1)

cat("\n  chosen: Gamma(2,", tau_pick$rate, ")\n")


######################################################
#  The set.
######################################################

cat("\n\nthe calibrated set, with background mortality\n")
cat("  knots  6 df on event-time quantiles, plus 2, 3, 5, 10, 20\n")
cat("  sigma  Gamma(", sigma_shape, ", ",
    round(sigma_shape / pick$sigma_mean, 2), ")\n", sep = "")
cat("  eta    p_meansurv(median = ", pick$eta_median, ", upper = ",
    pick$eta_median * pick$eta_upper_ratio, ")\n", sep = "")
cat("  tau    Gamma(2, ", tau_pick$rate, ")\n", sep = "")
cat("\nCalibrated by simulating excess hazards from the prior, adding the\n")
cat("cohort's background mortality, and checking implied survival at 5,\n")
cat("10, 20 and 30 years against judgements that rule out both\n")
cat("implausibly good and implausibly bad outcomes.\n")