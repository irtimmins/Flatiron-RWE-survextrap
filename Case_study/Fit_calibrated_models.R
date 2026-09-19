# Fitting the calibrated model under three evidence sets, with
# background mortality. Same priors and basis throughout; only what the
# model sees differs. Where the three agree the data decides, where they
# differ the extrapolation is carried by whatever was added.
#
# Priors and the smoothing model come from Choose_priors.R rather than
# being retyped. Background is supplied as a time/hazard data frame,
# which survextrap calls form (a): it then knows the background at all
# times, so predictions describe total survival.
#
# Outputs Figures/fitted_models.jpg and Data/fitted_models.rds

library(dplyr)
library(tidyr)
library(ggplot2)
library(survextrap)
library(survival)
library(cowplot)

fit_method <- "mcmc"       # "opt" for a quick look, mcmc for the answer
horizon <- 30
data_dir <- "Data"


#  Priors, data and background

cal <- readRDS(file.path(data_dir, "calibrated_priors.rds"))
backhaz <- readRDS(file.path(data_dir, "backhaz.rds"))
trial_data <- readRDS(file.path(data_dir, "trial_data.rds"))
trt_levels <- levels(trial_data$trt)

if (!all(c("time", "hazard") %in% names(backhaz)) || backhaz$time[1] != 0)
  stop("backhaz needs time and hazard columns with time starting at 0")

# the calibration was done under this, so the fit has to match
smooth_model <- if (is.null(cal$smooth_model)) "random_walk" else
  cal$smooth_model

cat("priors from Choose_priors.R\n")
cat("  knots ", paste(round(cal$mspline$knots, 2), collapse = ", "),
    ", boundary ", max(cal$mspline$knots), "\n", sep = "")
cat("  sigma Gamma(", cal$hsd$shape, ", ", round(cal$hsd$rate, 2),
    "), tau Gamma(", cal$hrsd$shape, ", ", round(cal$hrsd$rate, 2), ")\n",
    sep = "")
cat("  smoothing ", smooth_model, ", method ", fit_method, "\n\n", sep = "")


#  The three evidence sets

evidence_order <- c("ALEX only", "ALEX + PROFILE-1014",
                    "ALEX + PROFILE-1014 + Flatiron")

evidence <- tibble(
  label = evidence_order,
  file = c(NA_character_,
           file.path(data_dir, "historic_trial_aggregate.rds"),
           file.path(data_dir, "external_data_maic_weighted.rds")))

missing <- evidence$file[!is.na(evidence$file) & !file.exists(evidence$file)]
if (length(missing)) stop("missing: ", paste(missing, collapse = ", "))

read_external <- function(file) {
  if (is.na(file)) return(NULL)
  ext <- readRDS(file)
  needed <- c("start", "stop", "n", "r", "trt")
  if (!all(needed %in% names(ext)))
    stop(file, " needs ", paste(needed, collapse = ", "))
  ext$trt <- factor(as.character(ext$trt), levels = trt_levels)
  ext
}

fit_one <- function(label, file) {
  external <- read_external(file)
  cat("fitting: ", label,
      if (!is.null(external)) sprintf("  (%d external rows)",
                                      nrow(external)) else "", "\n", sep = "")
  
  args <- list(formula = Surv(time, status) ~ trt, data = trial_data,
               mspline = cal$mspline, nonprop = ~ trt,
               prior_hscale = cal$hscale, prior_hsd = cal$hsd,
               prior_loghr = cal$loghr, prior_hrsd = cal$hrsd,
               backhaz = backhaz, smooth_model = smooth_model,
               fit_method = fit_method)
  if (!is.null(external)) args$external <- external
  do.call(survextrap, args)
}

fits <- lapply(seq_len(nrow(evidence)),
               function(i) fit_one(evidence$label[i], evidence$file[i]))
names(fits) <- evidence$label


#  Convergence, before any number is read off

if (fit_method == "mcmc") {
  cat("\nconvergence\n")
  for (l in names(fits)) {
    sf <- fits[[l]]$stanfit
    cat(sprintf("  %-32s max rhat %.3f  divergences %d\n", l,
                max(rstan::summary(sf)$summary[, "Rhat"], na.rm = TRUE),
                sum(rstan::get_divergent_iterations(sf))))
  }
  cat("  rhat above 1.01, or any divergences, means that fit's intervals",
      "should not be reported\n")
}


#  Survival, hazard and the hazard ratio

t_grid <- seq(0, horizon, length.out = 301)
newdata <- data.frame(trt = factor(trt_levels, levels = trt_levels))

curves <- bind_rows(lapply(names(fits), function(l) {
  f <- fits[[l]]
  bind_rows(
    survival(f, newdata = newdata, t = t_grid) %>%
      mutate(.keep = "none", trt, t, median, lower, upper,
             variable = "survival"),
    hazard(f, newdata = newdata, t = t_grid) %>%
      mutate(.keep = "none", trt, t, median, lower, upper,
             variable = "hazard")) %>%
    mutate(evidence = l)
})) %>%
  mutate(evidence = factor(evidence, levels = evidence_order),
         trt = factor(trt, levels = trt_levels))

# the hazard ratio is what non-PH models widen, so it is reported.
# hazard_ratio's column names have varied between versions, so they are
# matched rather than assumed.
hr_curves <- bind_rows(lapply(names(fits), function(l) {
  h <- hazard_ratio(fits[[l]], t = t_grid)
  tcol <- if ("t" %in% names(h)) "t" else
    setdiff(names(h)[vapply(h, is.numeric, logical(1))],
            c("median", "lower", "upper"))[1]
  tibble(t = h[[tcol]], median = h$median, lower = h$lower,
         upper = h$upper, evidence = l)
})) %>%
  mutate(evidence = factor(evidence, levels = evidence_order))

saveRDS(list(fits = fits, curves = curves, hr = hr_curves, priors = cal,
             backhaz = backhaz, smooth_model = smooth_model,
             fit_method = fit_method),
        file.path(data_dir, "fitted_models.rds"))

cat("\nsurvival at selected times\n")
print(as.data.frame(curves %>%
                      filter(variable == "survival", t %in% c(5, 10, 20, 30)) %>%
                      mutate(.keep = "none", evidence, trt, t,
                             survival = sprintf("%.3f (%.3f, %.3f)", median, lower, upper))),
      row.names = FALSE)

cat("\nhazard ratio, alectinib against crizotinib\n")
print(as.data.frame(hr_curves %>%
                      filter(t %in% c(2, 5, 10, 30)) %>%
                      mutate(.keep = "none", evidence, t,
                             hr = sprintf("%.2f (%.2f, %.2f)", median, lower, upper))),
      row.names = FALSE)


#  Figure

km <- survfit(Surv(time, status) ~ trt, data = trial_data)
km_df <- tibble(t = km$time, median = km$surv,
                trt = factor(rep(trt_levels, km$strata),
                             levels = trt_levels),
                variable = "survival")

arm_colours <- setNames(c("#F8766D", "#00BFC4"), trt_levels)
panel <- theme_classic() +
  theme(strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text = element_text(size = 8))

surv_plot <- ggplot(filter(curves, variable == "survival"),
                    aes(t, median, colour = trt, fill = trt)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 0.85) +
  geom_line(data = km_df, linewidth = 0.4, linetype = "22") +
  scale_colour_manual(NULL, values = arm_colours) +
  scale_fill_manual(NULL, values = arm_colours) +
  scale_y_continuous("Overall survival", limits = c(0, 1),
                     labels = scales::percent) +
  scale_x_continuous("Time (years)", breaks = seq(0, horizon, by = 10)) +
  facet_wrap(~evidence, nrow = 1) + panel +
  theme(legend.position = "bottom")

# grey step is the background alone, which the total cannot fall below
haz_plot <- ggplot(filter(curves, variable == "hazard"),
                   aes(t, median, colour = trt, fill = trt)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 0.85) +
  geom_step(data = tibble(t = backhaz$time, median = backhaz$hazard),
            aes(t, median), inherit.aes = FALSE, colour = "grey30",
            linewidth = 0.5, linetype = "dashed") +
  scale_colour_manual(NULL, values = arm_colours) +
  scale_fill_manual(NULL, values = arm_colours) +
  scale_y_continuous("Total hazard (per year)") +
  coord_cartesian(ylim = c(0, 0.6)) +
  scale_x_continuous("Time (years)", breaks = seq(0, horizon, by = 10)) +
  facet_wrap(~evidence, nrow = 1) + panel +
  theme(legend.position = "none")

hr_plot <- ggplot(hr_curves, aes(t, median)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey60",
              alpha = 0.25) +
  geom_line(linewidth = 0.85) +
  geom_hline(yintercept = 1, colour = "grey40", linetype = "dashed",
             linewidth = 0.4) +
  scale_y_continuous("Hazard ratio", trans = "log10") +
  scale_x_continuous("Time (years)", breaks = seq(0, horizon, by = 10)) +
  facet_wrap(~evidence, nrow = 1) + panel

fig <- plot_grid(surv_plot, haz_plot, hr_plot, ncol = 1,
                 rel_heights = c(1.15, 1, 1))
print(fig)

dir.create("Figures", showWarnings = FALSE)
jpeg("Figures/fitted_models.jpg", width = 9, height = 9, units = "in",
     res = 600, quality = 100)
print(fig)
dev.off()

cat("\nwritten to Figures/fitted_models.jpg and Data/fitted_models.rds\n")
if (fit_method == "opt")
  cat("fitted with opt, so the intervals are approximate; refit with",
      "mcmc before taking any number from this\n")