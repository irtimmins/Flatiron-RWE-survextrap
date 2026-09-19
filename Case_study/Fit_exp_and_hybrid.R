# Exponential and piecewise models, as used in NICE TA536.
#
# Two purposes, and they are different quantities that have been
# conflated before.
#
# First, reproduce the appraisal. TA536 reports life years: survival
# integrated over a 30 year horizon in weekly cycles with a half cycle
# correction, discounted at 3.5 per cent. The company's original base
# case used the February 2017 cut with Kaplan-Meier to 18 months and an
# exponential tail, giving 5.14 alectinib and 4.32 crizotinib. After the
# appraisal consultation the committee preferred the December 2017 cut
# with Kaplan-Meier to 25 months, giving 5.46 and 4.42. Those are the
# numbers this script has to land on.
#
# Second, produce benchmarks for the paper: restricted mean survival at
# 5, 10, 20 and 30 years, undiscounted and discounted, on the same
# footing as the survextrap models, plus survival and hazard curves for
# the figures. RMST is not the same as a life year and will not match the
# appraisal figures, which is fine as long as the two are not put in one
# column.
#
# The exponential tail is rebased onto the Kaplan-Meier at the cut point,
# so survival is S_km(cut) times exp(-lambda (t - cut)). Carrying the
# exponential's own level forward instead is what TA536 does not do and
# is worth about 0.06 life years.
#
# Uncertainty comes from a bootstrap rather than from combining analytic
# standard errors, since there is no closed form for a discounted
# piecewise quantity and subtracting variances can go negative.
#
# Hazards. The exponential's hazard is exactly lambda, constant over the
# whole horizon, and its interval is the bootstrap distribution of
# lambda. The piecewise model has no continuous hazard before the cut:
# it uses the Kaplan-Meier there, whose hazard is a set of spikes at the
# event times, so something has to be chosen to represent it. A penalised
# spline from bshazard at a fixed penalty is used over the trial period,
# and exactly lambda after the cut. bshazard is the same estimator at the
# same penalty as the validation hazard in the figures, so the two are
# drawn on comparable terms; it should be read as the shape of the
# observed hazard rather than as a modelled quantity. The jump at the cut
# point is real and belongs to the model, not to the smoothing.

library(dplyr)
library(tidyr)
library(tibble)
library(survival)
library(bshazard)
library(ggplot2)

set.seed(1212)

# Smoothing for the hazard, matching the value in
# Figure3_survival_hazard.R, since the piecewise model's pre-cut hazard
# and the validation hazard appear on the same plot and should be
# smoothed alike. If one changes, change both.
#
# Note this sits inside the bootstrap: one spline fit per replicate per
# arm per cut, so roughly 4 * B fits in total. bshazard is slower than a
# kernel smoother, though fixing lambda avoids the REML search and is
# much faster than letting it choose. If the run becomes too slow, B is
# the thing to reduce, not the smoothing.
hazard_lambda <- 1000

cfg <- list(
  cuts = list(
    feb_2017 = list(file = "Data/trial_data.rds",
                    cut = 18 / 12,
                    label = "February 2017",
                    target = c(Alectinib = 5.14, Crizotinib = 4.32,
                               increment = 0.83)),
    dec_2017 = list(file = "Data/trial_data_dec_2017.rds",
                    cut = 25 / 12,
                    label = "December 2017",
                    target = c(Alectinib = 5.46, Crizotinib = 4.42,
                               increment = 1.05))),
  arms      = c("Alectinib", "Crizotinib"),
  tau       = c(5, 10, 20, 30),
  horizon   = 30,
  disc_rate = 0.035,
  cycle     = 7,          # days, the appraisal's cycle length
  B         = 100,        # bootstrap replicates
  grid      = seq(0, 30, length.out = 601))


# ---- the fitted pieces ----------------------------------------------------
#
# The exponential rate is the maximum likelihood estimate for the whole
# arm, deaths over person time, which is a closed form and so cheap
# enough to sit inside a bootstrap loop.

smooth_hazard <- function(d, at, max_time) {
  # bshazard has no max.time, so the window is imposed by censoring
  # everyone at the cut point and fitting that. Unlike a kernel smoother
  # it copes with the block of tied censoring times this creates, and the
  # same lambda applies whatever the window length, so the pre-cut hazard
  # and the validation curve drawn beside it are smoothed alike. A
  # bandwidth could not do that here: the pre-cut windows are only 1.5
  # and 2.08 years, far shorter than the validation data.
  d_pre <- d %>%
    mutate(status = ifelse(time > max_time, 0, status),
           time = pmin(time, max_time))
  fit <- try(bshazard(Surv(time, status) ~ 1, data = d_pre,
                      lambda = hazard_lambda, verbose = FALSE),
             silent = TRUE)
  if (inherits(fit, "try-error")) return(rep(NA_real_, length(at)))
  approx(fit$time, fit$hazard, xout = at, rule = 2)$y
}

surv_fns <- function(d, cut) {
  
  km <- survfit(Surv(time, status) ~ 1, data = d)
  st <- stepfun(km$time, c(1, km$surv))
  lambda <- sum(d$status) / sum(d$time)
  s_cut <- st(cut)
  
  list(
    lambda = lambda,
    
    exponential = function(t) exp(-lambda * t),
    
    hybrid = function(t)
      ifelse(t <= cut, st(t), s_cut * exp(-lambda * (t - cut))),
    
    # constant everywhere, exactly the fitted rate
    exponential_hazard = function(t) rep(lambda, length(t)),
    
    # smoothed observed hazard to the cut, then exactly lambda
    hybrid_hazard = function(t) {
      h <- rep(lambda, length(t))
      pre <- t <= cut
      if (any(pre))
        h[pre] <- smooth_hazard(d, at = t[pre], max_time = cut)
      h
    })
}


# ---- the estimands --------------------------------------------------------
#
# Life years use the appraisal's machinery: weekly cycles, a half cycle
# correction, and discounting applied in whole years. Restricted mean
# survival is a plain integral to a horizon, with discounting optional so
# the same function serves both columns of the results table.

life_years <- function(S, horizon = cfg$horizon, disc = cfg$disc_rate,
                       cycle_days = cfg$cycle, hcc = TRUE) {
  step <- cycle_days / 365.25
  t <- seq(0, horizon, by = step)
  w <- rep(1, length(t))
  if (hcc) { w[1] <- 0.5; w[length(w)] <- 0.5 }
  sum(S(t) * (1 + disc)^(-floor(t)) * w) * step
}

rmst_at <- function(S, tau, disc = 0) {
  g  <- seq(0, tau, length.out = 20001)
  dx <- g[2] - g[1]
  sum(head(S(g) * (1 + disc)^(-g), -1)) * dx
}


# ---- one arm, point estimate and bootstrap --------------------------------

one_fit <- function(f, m) {
  S <- f[[m]]
  H <- f[[paste0(m, "_hazard")]]
  list(surv = S(cfg$grid),
       hazard = H(cfg$grid),
       ly = life_years(S),
       rmst = sapply(cfg$tau, function(tt) rmst_at(S, tt, 0)),
       rmst_disc = sapply(cfg$tau, function(tt) rmst_at(S, tt, cfg$disc_rate)))
}

analyse_arm <- function(d, cut, B = cfg$B) {
  
  point <- surv_fns(d, cut)
  models <- c("exponential", "hybrid")
  
  est <- lapply(models, function(m) one_fit(point, m))
  names(est) <- models
  
  boot <- lapply(seq_len(B), function(b) {
    db <- d[sample(nrow(d), nrow(d), replace = TRUE), ]
    if (sum(db$status) < 2) return(NULL)
    fb <- surv_fns(db, cut)
    setNames(lapply(models, function(m) one_fit(fb, m)), models)
  })
  boot <- boot[!sapply(boot, is.null)]
  
  list(est = est, boot = boot)
}


# ---- run both cuts --------------------------------------------------------

collect <- function(cut_key) {
  
  cc <- cfg$cuts[[cut_key]]
  trial <- readRDS(cc$file) %>%
    mutate(trt = as.character(trt)) %>%
    select(trt, time, status)
  
  res <- lapply(cfg$arms, function(a)
    analyse_arm(filter(trial, trt == a), cc$cut))
  names(res) <- cfg$arms
  
  models <- c("exponential", "hybrid")
  
  # survival and hazard curves, with bootstrap percentile bands. One
  # object with a variable column, the same shape base_model_all.rds
  # uses, so the figure scripts can treat the two sources alike.
  curves <- bind_rows(lapply(cfg$arms, function(a)
    bind_rows(lapply(models, function(m)
      bind_rows(lapply(c("surv", "hazard"), function(key) {
        bs <- sapply(res[[a]]$boot, function(b) b[[m]][[key]])
        tibble(data_cut = cc$label, model = m, trt = a,
               variable = ifelse(key == "surv", "survival", "hazard"),
               t = cfg$grid,
               value = res[[a]]$est[[m]][[key]],
               lower = apply(bs, 1, quantile, 0.025, na.rm = TRUE),
               upper = apply(bs, 1, quantile, 0.975, na.rm = TRUE))
      }))))))
  
  # life years, the quantity the appraisal reports
  ly <- bind_rows(lapply(models, function(m) {
    v <- sapply(cfg$arms, function(a) res[[a]]$est[[m]]$ly)
    bs <- sapply(cfg$arms, function(a)
      sapply(res[[a]]$boot, function(b) b[[m]]$ly))
    inc <- bs[, "Alectinib"] - bs[, "Crizotinib"]
    tibble(data_cut = cc$label, model = m,
           quantity = c(cfg$arms, "increment"),
           value = c(v[["Alectinib"]], v[["Crizotinib"]],
                     v[["Alectinib"]] - v[["Crizotinib"]]),
           lower = c(quantile(bs[, "Alectinib"], 0.025),
                     quantile(bs[, "Crizotinib"], 0.025),
                     quantile(inc, 0.025)),
           upper = c(quantile(bs[, "Alectinib"], 0.975),
                     quantile(bs[, "Crizotinib"], 0.975),
                     quantile(inc, 0.975)),
           published = unname(cc$target[c("Alectinib", "Crizotinib",
                                          "increment")]))
  }))
  
  # restricted mean survival, undiscounted and discounted
  rmst <- bind_rows(lapply(models, function(m)
    bind_rows(lapply(c("rmst", "rmst_disc"), function(key) {
      dr <- if (key == "rmst") 0 else cfg$disc_rate
      arm_tab <- bind_rows(lapply(cfg$arms, function(a) {
        bs <- sapply(res[[a]]$boot, function(b) b[[m]][[key]])
        tibble(variable = "rmst", trt = a, time = cfg$tau,
               value = res[[a]]$est[[m]][[key]],
               lower = apply(bs, 1, quantile, 0.025),
               upper = apply(bs, 1, quantile, 0.975))
      }))
      bs_a <- sapply(res[["Alectinib"]]$boot, function(b) b[[m]][[key]])
      bs_c <- sapply(res[["Crizotinib"]]$boot, function(b) b[[m]][[key]])
      inc <- bs_a - bs_c
      diff_tab <- tibble(
        variable = "irmst", trt = NA_character_, time = cfg$tau,
        value = res[["Alectinib"]]$est[[m]][[key]] -
          res[["Crizotinib"]]$est[[m]][[key]],
        lower = apply(inc, 1, quantile, 0.025),
        upper = apply(inc, 1, quantile, 0.975))
      bind_rows(arm_tab, diff_tab) %>%
        mutate(data_cut = cc$label, model = m, disc_rate = dr)
    }))))
  
  list(curves = curves, ly = ly, rmst = rmst)
}

out <- lapply(names(cfg$cuts), collect)
names(out) <- names(cfg$cuts)

curves_all <- bind_rows(lapply(out, `[[`, "curves"))
ly_all     <- bind_rows(lapply(out, `[[`, "ly"))
rmst_all   <- bind_rows(lapply(out, `[[`, "rmst"))


# ---- does it reproduce the appraisal --------------------------------------
#
# The piecewise model is the comparison that matters: Kaplan-Meier to the
# cut then an exponential tail is what TA536 used. The pure exponential
# is a separate comparator the company also fitted, and has no published
# life years to check against.

cat("life years against NICE TA536\n")
cat("30 year horizon, weekly cycles, half cycle correction,",
    cfg$disc_rate * 100, "per cent discounting\n\n")
print(as.data.frame(ly_all %>%
                      filter(model == "hybrid") %>%
                      transmute(data_cut, quantity,
                                reconstructed = round(value, 2), published,
                                difference = round(value - published, 2),
                                ci = sprintf("(%.2f, %.2f)", lower, upper))), row.names = FALSE)

cat("\nexponential model, same conventions, no published value to check\n")
print(as.data.frame(ly_all %>%
                      filter(model == "exponential") %>%
                      transmute(data_cut, quantity, life_years = round(value, 2),
                                ci = sprintf("(%.2f, %.2f)", lower, upper))), row.names = FALSE)


# ---- the fitted hazards ----------------------------------------------------

cat("\nfitted exponential rates, per year\n")
print(as.data.frame(curves_all %>%
                      filter(variable == "hazard", model == "exponential") %>%
                      group_by(data_cut, trt) %>%
                      summarise(rate = round(first(value), 4),
                                ci = sprintf("(%.4f, %.4f)", first(lower), first(upper)),
                                .groups = "drop")), row.names = FALSE)

peak <- curves_all %>%
  filter(variable == "hazard", model == "hybrid") %>%
  summarise(peak = max(value, na.rm = TRUE)) %>% pull(peak)
cat(sprintf("\npiecewise hazard peaks at %.3f per year\n", peak))
if (peak > 0.45)
  cat("  this is above the 0.45 limit used on the hazard figures and\n",
      "  will be clipped there\n", sep = "")


# ---- benchmarks for the results table --------------------------------------

cat("\nrestricted mean survival, for comparison with the survextrap models\n")
print(as.data.frame(rmst_all %>%
                      filter(time %in% c(10, 30)) %>%
                      transmute(data_cut, model, disc_rate, variable,
                                trt = ifelse(is.na(trt), "difference", trt), time,
                                value = sprintf("%.2f (%.2f, %.2f)", value, lower, upper)) %>%
                      pivot_wider(names_from = c(variable, time), values_from = value)),
      row.names = FALSE)

dir.create("Data", showWarnings = FALSE)
saveRDS(curves_all, "Data/exp_hybrid_curves.rds")
saveRDS(ly_all,     "Data/exp_hybrid_life_years.rds")
saveRDS(rmst_all,   "Data/exp_hybrid_rmst.rds")

# survival only, under the name and column shape the earlier figure
# scripts read, so nothing that has not been updated yet breaks
saveRDS(curves_all %>%
          filter(variable == "survival") %>%
          select(data_cut, model, trt, time = t, value, lower, upper),
        "Data/exp_hybrid_survival.rds")

# the shapes the existing Table 2 script reads, February cut, undiscounted
saveRDS(rmst_all %>%
          filter(data_cut == "February 2017", model == "exponential",
                 disc_rate == 0) %>%
          select(variable, trt, time, value, lower, upper),
        "Data/exp_model_rmst.rds")
saveRDS(rmst_all %>%
          filter(data_cut == "February 2017", model == "hybrid",
                 disc_rate == 0) %>%
          select(variable, trt, time, value, lower, upper),
        "Data/hybrid_model_rmst.rds")


# ---- supplementary figure --------------------------------------------------
#
# Survival and hazard for the two company models, February cut, with the
# cut point marked. The jump in the piecewise hazard at that line is the
# model switching from the observed data to its exponential tail.

km_data_all <- readRDS("Data/km_data_all.rds")

model_labels <- c(exponential = "Exponential",
                  hybrid = "Kaplan-Meier to 18 months,\nexponential tail")

feb <- curves_all %>%
  filter(data_cut == "February 2017") %>%
  mutate(model = factor(model, levels = names(model_labels),
                        labels = model_labels))

sup_a <- feb %>%
  filter(variable == "survival") %>%
  ggplot() +
  theme_classic() +
  theme(legend.position = "bottom") +
  geom_line(data = km_data_all %>%
              filter(dataset == "ALEX trial") %>%
              select(-c(dataset, model)),
            aes(x = time, y = surv, colour = trt), alpha = 0.5) +
  geom_ribbon(aes(x = t, ymin = lower, ymax = upper, fill = trt),
              alpha = 0.25, colour = NA) +
  geom_vline(xintercept = cfg$cuts$feb_2017$cut, alpha = 0.4) +
  geom_line(aes(x = t, y = value, colour = trt)) +
  scale_colour_discrete("Treatment") +
  scale_fill_discrete("Treatment") +
  scale_y_continuous("Overall Survival", limits = c(0, 1.01),
                     labels = scales::percent) +
  xlab("Time (years)") +
  facet_wrap(~model)

sup_b <- feb %>%
  filter(variable == "hazard") %>%
  ggplot() +
  theme_classic() +
  theme(legend.position = "bottom") +
  geom_ribbon(aes(x = t, ymin = lower, ymax = upper, fill = trt),
              alpha = 0.25, colour = NA) +
  geom_vline(xintercept = cfg$cuts$feb_2017$cut, alpha = 0.4) +
  geom_line(aes(x = t, y = value, colour = trt)) +
  scale_colour_discrete("Treatment") +
  scale_fill_discrete("Treatment") +
  scale_y_continuous("Hazard") +
  coord_cartesian(ylim = c(0, 0.45)) +
  xlab("Time (years)") +
  facet_wrap(~model)

print(sup_a)
print(sup_b)

dir.create("Figures", showWarnings = FALSE)
jpeg(file = "Figures/Sup_figure_4_exp_hybrid.jpg",
     width = 7.2, height = 6.6, units = "in", res = 600, quality = 100)
print(cowplot::plot_grid(sup_a + labs(title = "(a)"),
                         sup_b + labs(title = "(b)"),
                         ncol = 1))
dev.off()