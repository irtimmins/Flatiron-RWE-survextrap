# Choosing the priors for the survextrap models, and showing what they
# assume. Four steps: sigma, then the scale, then the hazard ratio's
# width, then tau. Sigma first because its summary is a ratio and so
# does not depend on the scale.
#
# Outputs Figures/prior_predictive.jpg and Data/calibrated_priors.rds

library(dplyr)
library(tidyr)
library(ggplot2)
library(survextrap)
library(survival)


#  Judgements and settings

# What was knowable at TA536 in 2018. Chemotherapy-era ALK positive
# advanced NSCLC had a median survival of 1.5 to 2 years (Shaw 2011);
# PROFILE 1014 gave four year survival on crizotinib of about 57 per
# cent (Solomon 2018); no cohort had been followed past about six years.
#
# The ceilings rule out implausibly good extrapolations. The two year
# floor is the one constraint on pessimism, and comes from the
# chemotherapy era: a targeted agent should not be assumed worse than
# what it replaced. Without it nothing stops the scale collapsing.
# 0 means no floor, 1 no ceiling; the 5th and 95th percentiles of prior
# draws are checked.
landmarks <- tribble(
  ~time, ~surv_low, ~surv_high,
  2,      0.25,       1.00,
  10,      0.00,       0.50,
  20,      0.00,       0.30,
  30,      0.00,       0.15)

judgement_text <- c(
  "at least a quarter alive at two years, the chemotherapy era floor",
  "no more than half alive at ten years",
  "fewer than 30 per cent alive at 20 years, fewer than 15 at 30",
  "otherwise the observed period is left to the data")

# Among priors that all meet a judgement, how much slack to take.
# Taking the loosest at every step compounds into intervals too wide to
# inform anything.
prefer <- "narrow"        # or "wide"

# sigma: how much the hazard moves over the splined period. The median
# is what a typical draw does, the 90th what the prior commonly permits.
sigma_ratio <- c(low = 2.5, high = 8)
sigma_p90_max <- 60

# tau: how far the hazard ratio moves over time. A plausible span of
# 0.5 to 1.25 is a ratio of about 2.5, read as what the prior permits
# rather than what a typical draw does, so the median is lower and the
# tail capped. tau is what widens intervals in non-PH models.
tau_target <- c(low = 1.4, high = 2.6)
tau_p975_max <- 12

min_pass_rate <- 0.8
n_seeds <- 5
nsim <- 800
horizon <- 30


#  Basis, background and fixed priors

trial_data <- readRDS("Data/trial_data.rds")
trt_levels <- levels(trial_data$trt)
as_trt <- function(x) factor(x, levels = trt_levels)

# Boundary at 5, so the excess hazard is constant from there: an
# exponential tail as in TA536. A later boundary gives intervals too
# wide to inform a decision, especially under non-proportional hazards.
# The cost is that the shape past five years is assumed, not estimated.
sp <- mspline_spec(Surv(time, status) ~ trt, data = trial_data,
                   df = 5, add_knots = c(2, 5))

backhaz <- readRDS("Data/backhaz.rds")
criz <- data.frame(trt = as_trt("Crizotinib"))
alec <- data.frame(trt = as_trt("Alectinib"))

# survextrap fits with random_walk; the prior functions default to
# exchangeable, so it is passed explicitly to every call.
smooth_model <- "random_walk"

# centred on one, so nothing here assumes either drug is better; the
# width covers the plausible span of 0.5 to 1.25
prior_loghr <- p_hr(median = 1, upper = 2)
loghr_upper <- 2


#  Drawing from the prior, with background added

draw_total <- function(prior_hsd, prior_hscale, prior_hrsd = NULL,
                       nonprop = NULL, newdata = criz, newdata0 = NULL,
                       loghr = prior_loghr, nsim = 800) {
  
  pred <- prior_sample_hazard(
    knots = sp$knots, degree = sp$degree, bsmooth = sp$bsmooth,
    prior_hsd = prior_hsd, prior_hscale = prior_hscale,
    smooth_model = smooth_model, prior_loghr = loghr,
    formula = ~ trt, nonprop = nonprop,
    newdata = newdata, newdata0 = newdata0, prior_hrsd = prior_hrsd,
    tmin = 0, tmax = horizon, nsim = nsim)
  
  time_col <- setdiff(names(pred)[vapply(pred, is.numeric, logical(1))],
                      c("rep", "haz", "hr", "haz0"))
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

assess <- function(d) {
  bind_rows(lapply(seq_len(nrow(landmarks)), function(i) {
    tt <- landmarks$time[i]
    at <- d %>% group_by(rep) %>%
      slice_min(abs(time - tt), n = 1, with_ties = FALSE) %>% ungroup()
    s5 <- quantile(at$surv, 0.05); s95 <- quantile(at$surv, 0.95)
    tibble(time = tt, surv_5 = s5, surv_50 = median(at$surv),
           surv_95 = s95, haz_50 = median(at$total),
           haz_95 = quantile(at$total, 0.95),
           too_high = s95 > landmarks$surv_high[i],
           ok = !(s5 < landmarks$surv_low[i] |
                    s95 > landmarks$surv_high[i]))
  }))
}


#  Step 1: sigma
#
# Gamma(2, rate), only the rate swept: the package default family and
# what Timmins et al vary. Measured over the splined period, since past
# the boundary the hazard is flat and every prior would look alike.

cat("step 1, sigma: hazard varies ", sigma_ratio["low"], " to ",
    sigma_ratio["high"], " fold to ", max(sp$knots),
    " years, 90th under ", sigma_p90_max, "\n\n", sep = "")

sigma_out <- bind_rows(lapply(c(0.5, 1, 1.5, 2, 3, 4, 5, 8, 12, 20),
                              function(rt) {
                                s <- prior_haz_sd(mspline = sp, prior_hsd = p_gamma(2, rt),
                                                  prior_hscale = p_meansurv(median = 3, upper = 6,
                                                                            mspline = sp),
                                                  prior_loghr = prior_loghr, smooth_model = smooth_model,
                                                  formula = ~ trt, newdata = criz,
                                                  tmin = 0, tmax = max(sp$knots),
                                                  quantiles = c(0.5, 0.9, 0.975), nsim = 2000)
                                tibble(rate = rt, ratio = s["50%", "hr"], p90 = s["90%", "hr"],
                                       p975 = s["97.5%", "hr"],
                                       sig_lo = qgamma(0.025, 2, rt), sig_hi = qgamma(0.975, 2, rt))
                              })) %>%
  mutate(meets = ratio >= sigma_ratio["low"] & ratio <= sigma_ratio["high"] &
           p90 <= sigma_p90_max)

# mutate leaves existing columns where they were, so the order is set
# explicitly rather than by the order written above
print(as.data.frame(sigma_out %>%
                      mutate(.keep = "none", prior = sprintf("Gamma(2, %g)", rate),
                             sigma = sprintf("%.2f-%.2f", sig_lo, sig_hi),
                             varies_by = sprintf("%.1f", ratio), p90 = sprintf("%.1f", p90),
                             meets = ifelse(meets, "yes", "")) %>%
                      select(prior, sigma, varies_by, p90, meets)), row.names = FALSE)

sigma_ok <- sigma_out %>% filter(meets)

if (nrow(sigma_ok)) {
  sigma_choice <- if (prefer == "narrow")
    sigma_ok %>% slice_max(rate, n = 1, with_ties = FALSE)
  else sigma_ok %>% slice_min(rate, n = 1, with_ties = FALSE)
} else {
  cat("\n  nothing meets both; taking the nearest on typical variation\n")
  sigma_choice <- sigma_out %>%
    mutate(d = abs(ratio - mean(sigma_ratio))) %>%
    slice_min(d, n = 1, with_ties = FALSE)
}

sigma_rate <- sigma_choice$rate
chosen_hsd <- p_gamma(2, sigma_rate)

cat("\n  chosen Gamma(2, ", sigma_rate, "): varies ",
    sprintf("%.1f", sigma_choice$ratio), " fold typically, ",
    sprintf("%.1f", sigma_choice$p90), " at the 90th\n\n", sep = "")


#  Step 2: the scale, given that sigma
#
# The scale has to meet the ceilings on its own, which is why the sweep
# reaches so low. Replicated across seeds so the choice is not a fluke.

grid <- expand_grid(eta_median = c(0.75, 1, 1.25, 1.5, 2, 2.5, 3),
                    eta_upper_ratio = c(1.5, 2, 2.5, 3))

cat("step 2, the scale:", nrow(grid), "priors across", n_seeds, "seeds\n")
for (j in judgement_text) cat("  ", j, "\n", sep = "")

reps <- bind_rows(lapply(seq_len(n_seeds), function(s) {
  set.seed(100 + s)
  bind_rows(lapply(seq_len(nrow(grid)), function(i) {
    g <- grid[i, ]
    a <- assess(draw_total(
      prior_hsd = chosen_hsd,
      prior_hscale = p_meansurv(median = g$eta_median,
                                upper = g$eta_median * g$eta_upper_ratio,
                                mspline = sp), nsim = nsim))
    tibble(g, passes = all(a$ok),
           high_10 = a$too_high[a$time == 10],
           high_20 = a$too_high[a$time == 20],
           high_30 = a$too_high[a$time == 30],
           s10 = a$surv_50[a$time == 10], s10_hi = a$surv_95[a$time == 10],
           s20_95 = a$surv_95[a$time == 20],
           s30_95 = a$surv_95[a$time == 30])
  }))
}))

summary_grid <- reps %>%
  group_by(eta_median, eta_upper_ratio) %>%
  summarise(pass_rate = mean(passes), across(s10:s30_95, mean),
            .groups = "drop") %>%
  arrange(desc(pass_rate), eta_upper_ratio, desc(eta_median))

# a block of passers is robust, a lone cell is a fluke
cat("\npass rate, scale median across and prior width down\n")
print(as.data.frame(summary_grid %>%
                      mutate(cell = sprintf("%3.0f", 100 * pass_rate)) %>%
                      select(eta_upper_ratio, eta_median, cell) %>%
                      pivot_wider(names_from = eta_median, values_from = cell,
                                  names_prefix = "median=") %>%
                      arrange(eta_upper_ratio)), row.names = FALSE)

cat("\nthe passers, and which landmark the rest fail\n")
print(as.data.frame(summary_grid %>%
                      filter(pass_rate > 0) %>%
                      mutate(.keep = "none",
                             eta = sprintf("ms(%g, %g)", eta_median,
                                           eta_median * eta_upper_ratio),
                             pass = sprintf("%.0f%%", 100 * pass_rate),
                             `10y` = sprintf("%.0f (to %.0f)", 100 * s10, 100 * s10_hi),
                             `20y 95th` = sprintf("%.0f", 100 * s20_95),
                             `30y 95th` = sprintf("%.0f", 100 * s30_95))), row.names = FALSE)
cat("  failing, 10y/20y/30y: ",
    paste(sprintf("%.0f%%", 100 * c(mean(reps$high_10), mean(reps$high_20),
                                    mean(reps$high_30))),
          collapse = " "), "\n", sep = "")

eligible <- summary_grid %>% filter(pass_rate >= min_pass_rate)

if (nrow(eligible)) {
  # narrowest prior and, within that, the least pessimistic scale
  pick <- if (prefer == "narrow")
    eligible %>% arrange(eta_upper_ratio, desc(eta_median)) %>% slice(1)
  else
    eligible %>% arrange(desc(eta_upper_ratio), desc(eta_median)) %>% slice(1)
} else {
  pick <- summary_grid %>% slice(1)
  cat("\n  none passes robustly, best ",
      sprintf("%.0f%%", 100 * max(summary_grid$pass_rate)),
      " of seeds; loosen the binding landmark\n", sep = "")
}

chosen_hscale <- p_meansurv(median = pick$eta_median,
                            upper = pick$eta_median * pick$eta_upper_ratio,
                            mspline = sp)

cat("\n  chosen p_meansurv(", pick$eta_median, ", ",
    pick$eta_median * pick$eta_upper_ratio, ")\n\n", sep = "")


#  Step 3: tau

cat("step 3, tau: ratio ", tau_target["low"], " to ", tau_target["high"],
    ", 97.5th under ", tau_p975_max, "\n\n", sep = "")

set.seed(99)
tau_out <- bind_rows(lapply(c(1, 2, 3, 5, 10), function(r) {
  hrr <- draw_total(chosen_hsd, chosen_hscale, prior_hrsd = p_gamma(2, r),
                    nonprop = ~ trt, newdata = alec, newdata0 = criz,
                    nsim = nsim) %>%
    mutate(hr_total = total / (haz / hr + backhaz)) %>%
    filter(time <= 10) %>%
    group_by(rep) %>%
    summarise(r = quantile(hr_total, 0.9) / quantile(hr_total, 0.1),
              .groups = "drop")
  tibble(rate = r, ratio = median(hrr$r), upper = quantile(hrr$r, 0.975))
})) %>%
  mutate(meets = ratio >= tau_target["low"] & ratio <= tau_target["high"] &
           upper <= tau_p975_max)

print(as.data.frame(tau_out %>%
                      mutate(.keep = "none",
                             prior = sprintf("Gamma(2, %g)", rate),
                             varies_by = sprintf("%.2f", ratio),
                             spans = sprintf("%.2f to %.2f", 1 / sqrt(ratio), sqrt(ratio)),
                             p975 = sprintf("%.0f", upper),
                             meets = ifelse(meets, "yes", "")) %>%
                      select(prior, varies_by, spans, p975, meets)), row.names = FALSE)

tau_ok <- tau_out %>% filter(meets)

if (nrow(tau_ok)) {
  tau_pick <- if (prefer == "narrow")
    tau_ok %>% arrange(desc(rate)) %>% slice(1)
  else tau_ok %>% arrange(rate) %>% slice(1)
} else {
  cat("\n  nothing meets both; taking the tightest inside the cap\n")
  under_cap <- tau_out %>% filter(upper <= tau_p975_max)
  tau_pick <- if (nrow(under_cap))
    under_cap %>% arrange(desc(ratio)) %>% slice(1)
  else tau_out %>% arrange(upper) %>% slice(1)
}

cat("\n  chosen Gamma(2, ", tau_pick$rate, "): ratio ",
    sprintf("%.2f", tau_pick$ratio), ", spans ",
    sprintf("%.2f to %.2f", 1 / sqrt(tau_pick$ratio),
            sqrt(tau_pick$ratio)), ", ", sprintf("%.0f", tau_pick$upper),
    " at the 97.5th\n\n", sep = "")


#  What the set assumes, and whether it fights the data

set.seed(7)
final <- draw_total(chosen_hsd, chosen_hscale, nsim = 2000)

cat("crizotinib under the chosen priors, background included\n")
print(as.data.frame(assess(final) %>%
                      mutate(.keep = "none", time,
                             survival = sprintf("%.0f%% (%.0f-%.0f)", 100 * surv_50,
                                                100 * surv_5, 100 * surv_95),
                             hazard = sprintf("%.3f (95th %.3f)", haz_50, haz_95),
                             within = ifelse(ok, "yes", "no"))), row.names = FALSE)

saveRDS(list(hscale = chosen_hscale, hsd = chosen_hsd,
             hrsd = p_gamma(2, tau_pick$rate), loghr = prior_loghr,
             mspline = sp, landmarks = landmarks, choice = pick,
             smooth_model = smooth_model),
        "Data/calibrated_priors.rds")

# Band inclusion alone is too weak: a wide band contains almost
# anything. The gap between the prior's median and the data is what
# says whether the prior is centred where the evidence is.
check_t <- 2
km_at <- summary(survfit(Surv(time, status) ~ trt, data = trial_data),
                 times = check_t, extend = TRUE)
obs <- setNames(km_at$surv, sub("^trt=", "", km_at$strata))

pb <- final %>%
  group_by(rep) %>%
  slice_min(abs(time - check_t), n = 1, with_ties = FALSE) %>%
  ungroup() %>% pull(surv) %>% quantile(c(0.05, 0.5, 0.95))

gap <- vapply(obs, function(o) o - pb[2], numeric(1))
max_gap <- 0.20

cat("\nat ", check_t, " years: prior ",
    sprintf("%.0f%% (%.0f to %.0f)", 100 * pb[2], 100 * pb[1], 100 * pb[3]),
    ", observed ",
    paste(sprintf("%s %.0f%%", names(obs), 100 * obs), collapse = ", "),
    "\n  gap to the prior median: ",
    paste(sprintf("%+.0f pts", 100 * gap), collapse = ", "), "\n", sep = "")

if (any(obs < pb[1] | obs > pb[3])) {
  cat("  outside the prior's band: the prior treats the data as",
      "unlikely\n")
} else if (any(abs(gap) > max_gap)) {
  cat("  inside the band but the median sits ",
      sprintf("%.0f", 100 * max(abs(gap))),
      " points off: the prior is centred where the data is not\n", sep = "")
} else {
  cat("  inside the band and within ", round(100 * max_gap),
      " points: not fighting the data\n", sep = "")
}


#  Figure

bands <- final %>%
  group_by(time) %>%
  summarise(lo = quantile(surv, 0.05), md = median(surv),
            hi = quantile(surv, 0.95), .groups = "drop")

limits <- landmarks %>%
  filter(surv_low > 0 | surv_high < 1) %>%
  mutate(.keep = "none", time, ymin = surv_low, ymax = surv_high)

km <- survfit(Surv(time, status) ~ trt, data = trial_data)
km_df <- tibble(time = km$time, surv = km$surv,
                trt = rep(trt_levels, km$strata))

fig <- ggplot(bands, aes(time)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#2C6FB3", alpha = 0.18) +
  geom_line(aes(y = md), colour = "#2C6FB3", linewidth = 0.9) +
  geom_linerange(data = limits, aes(x = time, ymin = ymin, ymax = ymax),
                 colour = "#B3452C", linewidth = 1.2, inherit.aes = FALSE) +
  geom_line(data = km_df, aes(time, surv, group = trt),
            colour = "grey25", linetype = "22", linewidth = 0.5) +
  scale_y_continuous("Overall survival", limits = c(0, 1),
                     labels = scales::percent) +
  scale_x_continuous("Time (years)", breaks = seq(0, horizon, 5)) +
  theme_classic()

print(fig)

dir.create("Figures", showWarnings = FALSE)
jpeg("Figures/prior_predictive.jpg", width = 6, height = 4, units = "in",
     res = 600, quality = 100)
print(fig)
dev.off()

cat("\nthe calibrated set\n")
cat("  knots  ", paste(round(sp$knots, 2), collapse = ", "),
    ", boundary ", max(sp$knots), "\n", sep = "")
cat("  sigma  Gamma(2, ", sigma_rate, ")\n", sep = "")
cat("  eta    p_meansurv(", pick$eta_median, ", ",
    pick$eta_median * pick$eta_upper_ratio, ")\n", sep = "")
cat("  loghr  p_hr(1, ", loghr_upper, ")\n", sep = "")
cat("  tau    Gamma(2, ", tau_pick$rate, ")\n", sep = "")
cat("\nwritten to Figures/prior_predictive.jpg and",
    "Data/calibrated_priors.rds\n")