# Figure 3: survival and hazard, company models against survextrap under
# increasing evidence, for each ALEX data cut.
#
# Five rows, one per model or evidence set, each with survival on the
# left and hazard on the right, under a title spanning both:
#   1. Exponential (company)
#   2. Piecewise: Kaplan-Meier + exponential tail (company)
#   3. survextrap, ALEX trial only
#   4. survextrap, ALEX + PROFILE-1014
#   5. survextrap, ALEX + PROFILE-1014 + Flatiron RWE
#
# Built from ten self-contained plots rather than one faceted plot. A
# facet cannot carry a different y axis title per column, since axis
# titles are drawn once for the whole plot, and having "Overall survival"
# and "Hazard" as axis labels rather than strips is the point of this
# layout. Each cell therefore keeps its own axes, so the time axis
# repeats down the rows.
#
# Row 5 always uses the MAIC weighted Flatiron cohort; the title does
# not say so, since the figure caption is where that belongs. Flatiron
# without MAIC weighting is left out entirely.
#
# Rows 3 to 5 use one survextrap model type at a time, PH, Non-PH or
# separate arms, so the function is called once per model type per cut,
# giving six output files. Rows 1 and 2 are the same in all three, since
# the company models have no PH, Non-PH or separate arms versions.
#
# The priors are the calibrated set, one per scenario, matching Table 2.
# manuscript's stated base case (Non-PH, tau ~ Gamma(2,10)).
#
# Four series carry the legend: each arm's fitted curve, and each arm's
# Kaplan-Meier from the final overall survival analysis, 28 April 2025,
# which no curve here has seen. The distance between a fitted curve and
# its dashed counterpart past the end of the trial data is the point of
# the figure. The hazard cells use a penalised spline from bshazard at a
# fixed penalty on that same final analysis data.
#
# The legend merges colour and linetype into one key per series, which
# only happens when both scales carry the same title, limits and labels.
# Giving them different titles and then hiding the title text does not
# count and produces a second key of black dashed lines.
#
# The Kaplan-Meier of the cut being fitted to is drawn thin and solid in
# the arm's colour, with a fixed colour rather than a mapped one, so it
# stays out of the legend. It sits under the curve fitted to it.
#
# Output files
#   February: Figures/Sup_figure_5_ph.jpg
#             Figures/Figure_3_nonph.jpg
#             Figures/Sup_figure_6_sep_arms.jpg
#   December: Figures/Sup_figure_10_ph_dec.jpg
#             Figures/Sup_figure_11_nonph_dec.jpg
#             Figures/Sup_figure_12_sep_arms_dec.jpg
# The December file numbers are placeholders.

library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(survival)
library(bshazard)

source("Functions/Survextrap_model.R")

# The grid no longer sweeps priors: one calibrated set per scenario, so
# a model is picked by its name and its scenario label.
scenario_value <- "primary"
x_max <- 30
hazard_max <- 0.45

# Smoothing for the hazard estimates, shared with Fit_exp_and_hybrid.R.
# bshazard fits a penalised spline; lambda is the penalty, and a larger
# one gives a smoother curve. Fixing it rather than letting bshazard
# choose by REML keeps the two arms and the three data cuts smoothed
# alike, which a bandwidth in years could not do, since their follow-up
# lengths differ by a factor of four. If one script changes, change both.
hazard_lambda <- 1000

km_label <- "Kaplan-Meier, ALEX trial\nApr 2025 data cut"

series_levels <- c("Alectinib (model)",
                   "Crizotinib (model)",
                   paste0("Alectinib (", km_label, ")"),
                   paste0("Crizotinib (", km_label, ")"))

series_colours <- setNames(
  c("#F8766D", "#00BFC4", "firebrick", "darkcyan"), series_levels)

# "22" is a short dash, short gap; the default dashed pattern is a fixed
# length and reads as a handful of long strokes over thirty years
series_linetypes <- setNames(c("solid", "solid", "22", "22"), series_levels)

fitted_colours <- c(Alectinib = "#F8766D", Crizotinib = "#00BFC4")

model_series <- function(trt)
  factor(paste0(trt, " (model)"), levels = series_levels)
valid_series <- function(trt)
  factor(paste0(trt, " (", km_label, ")"), levels = series_levels)

panel_levels <- c("Exponential",
                  "Piecewise: Kaplan-Meier + exponential tail",
                  "survextrap, ALEX trial only",
                  "survextrap, ALEX + PROFILE-1014",
                  "survextrap, ALEX + PROFILE-1014 + Flatiron RWE")

survextrap_dataset_levels <- c("trial_only", "trial_and_historic",
                               "trial_and_all_maic")


######################################################
#  Read what the figure needs.
######################################################

base_results <- readRDS("Data/base_model_all.rds")
km_data_all <- readRDS("Data/km_data_all.rds")
alex_2025 <- readRDS("Data/trial_data_apr_2025.rds")
company_curves <- readRDS("Data/exp_hybrid_curves.rds")

if (max(company_curves$t) < 29)
  stop("exp_hybrid_curves.rds only reaches ",
       round(max(company_curves$t), 1), " years. Re-run ",
       "Fit_exp_and_hybrid.R first.")

validation_hazard <- bind_rows(lapply(unique(as.character(alex_2025$trt)),
                                      function(a) {
                                        d <- alex_2025 %>% filter(trt == a)
                                        # A penalised spline at a fixed penalty rather than a kernel with a
                                        # fixed window. A window in years cannot serve both arms here:
                                        # crizotinib's follow-up is about a third of alectinib's, so any
                                        # single bandwidth is either too wide for one or too narrow for the
                                        # other's sparse tail. The penalty degrades more gracefully where
                                        # events thin out.
                                        fit <- bshazard(Surv(time, status) ~ 1, data = d,
                                                        lambda = hazard_lambda, verbose = FALSE)
                                        tibble(trt = a, time = fit$time, value = fit$hazard,
                                               variable = "hazard")
                                      }))

validation_survival <- km_data_all %>%
  filter(dataset == "ALEX trial, final analysis") %>%
  transmute(trt, time, value = surv, variable = "survival")

if (nrow(validation_survival) == 0)
  stop("no Kaplan-Meier rows for the final analysis in km_data_all.rds")

validation_lines <- bind_rows(validation_survival, validation_hazard) %>%
  mutate(series = valid_series(trt))

cat("smoothed validation hazard estimated to",
    round(max(validation_hazard$time), 1), "years\n")


######################################################
#  The two company rows, for one data cut.
######################################################

company_rows <- function(cut_wanted) {
  cs <- company_curves %>%
    filter(data_cut == cut_wanted) %>%
    mutate(panel = ifelse(model == "exponential", panel_levels[1],
                          panel_levels[2]))
  if (nrow(cs) == 0)
    stop("no company rows for ", cut_wanted, " in exp_hybrid_curves.rds")
  cs %>% transmute(panel, trt, t, median = value, lower, upper, variable)
}


######################################################
#  One cell: survival or hazard, for one row.
######################################################

make_cell <- function(model_lines, trial_km, panel_wanted, variable_wanted,
                      y_label, y_scale, side) {
  
  d <- filter(model_lines, panel == panel_wanted, variable == variable_wanted)
  v <- filter(validation_lines, variable == variable_wanted)
  
  p <- ggplot() +
    geom_ribbon(data = d, aes(x = t, ymin = lower, ymax = upper, fill = trt),
                alpha = 0.20, colour = NA)
  
  # the cut being fitted to, only on the survival cells
  if (variable_wanted == "survival")
    for (a in names(fitted_colours))
      p <- p + geom_line(data = filter(trial_km, trt == a),
                         aes(x = time, y = value),
                         colour = fitted_colours[[a]], linewidth = 0.4)
  
  p +
    geom_line(data = v, aes(x = time, y = value, colour = series,
                            linetype = series), linewidth = 0.9) +
    geom_line(data = d, aes(x = t, y = median, colour = series,
                            linetype = series), linewidth = 0.85) +
    scale_colour_manual(NULL, values = series_colours, limits = series_levels) +
    scale_linetype_manual(NULL, values = series_linetypes,
                          limits = series_levels) +
    scale_fill_manual(values = fitted_colours, guide = "none") +
    # the default 5 per cent expansion is what left bare axis running on
    # past 30 years; a small additive pad keeps the last curve off the
    # frame without the overrun
    scale_x_continuous("Time (years)", breaks = seq(0, x_max, by = 5),
                       expand = c(0, 0)) +
    y_scale +
    coord_cartesian(xlim = c(-0.3, x_max + 0.4), expand = FALSE) +
    ylab(y_label) +
    theme_classic() +
    theme(legend.position = "none",
          axis.text = element_text(size = 7),
          axis.title = element_text(size = 8),
          # titles pulled in towards their axes
          axis.title.x = element_text(margin = margin(t = 0.8)),
          axis.title.y = element_text(margin = margin(r = 0.8)),
          # tight on the inner edge so the pair sits close to flush
          plot.margin = if (side == "left")
            unit(c(0.09, 0, 0.015, 0.04), "cm")
          else unit(c(0.09, 0.06, 0.015, 0), "cm"))
}


######################################################
#  A row title spanning both cells.
######################################################

row_title <- function(label) {
  ggplot() +
    annotate("text", x = 0.004, y = 0.62, label = label, hjust = 0,
             vjust = 0.5, size = 3.5, fontface = "bold") +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    theme_void() +
    theme(plot.margin = unit(c(0.04, 0.1, 0.04, 0.04), "cm"))
}


######################################################
#  One figure, for one cut and one survextrap model type.
######################################################

make_figure3 <- function(cut_wanted, km_dataset, model_type, figure_file) {
  
  # survival and hazard only: base_model_all.rds also holds rmst and
  # irmst rows, which have no place on a curve and which previously
  # arrived here as extra rows with a missing panel label
  survextrap_rows <- base_results %>%
    filter(cut_label == cut_wanted, variable %in% c("survival", "hazard"),
           t > 0, datasets %in% survextrap_dataset_levels) %>%
    filter_model(model_name = model_type,
                 scenario_value = scenario_value) %>%
    mutate(panel = as.character(factor(datasets,
                                       levels = survextrap_dataset_levels,
                                       labels = panel_levels[3:5]))) %>%
    transmute(panel, trt, t, median, lower, upper, variable)
  
  if (nrow(survextrap_rows) == 0)
    stop("no survextrap rows for ", cut_wanted, ", ", model_type,
         ". Check cut_label, scenario and datasets against the values ",
         "actually present in base_model_all.rds")
  
  model_lines <- bind_rows(company_rows(cut_wanted), survextrap_rows) %>%
    filter(variable %in% c("survival", "hazard")) %>%
    mutate(series = model_series(trt))
  
  trial_km <- km_data_all %>%
    filter(dataset == km_dataset) %>%
    transmute(trt, time, value = surv)
  
  if (nrow(trial_km) == 0)
    stop("no Kaplan-Meier rows for '", km_dataset, "' in km_data_all.rds. ",
         "Available: ", paste(unique(km_data_all$dataset), collapse = ", "))
  
  surv_scale <- scale_y_continuous(limits = c(0, 1), labels = scales::percent)
  haz_scale <- scale_y_continuous(limits = c(0, hazard_max),
                                  oob = scales::oob_keep)
  
  rows <- lapply(seq_along(panel_levels), function(i) {
    pl <- panel_levels[i]
    left <- make_cell(model_lines, trial_km, pl, "survival",
                      "Overall survival", surv_scale, "left")
    right <- make_cell(model_lines, trial_km, pl, "hazard",
                       "Hazard", haz_scale, "right")
    plot_grid(row_title(paste0("(", letters[i], ") ", pl)),
              plot_grid(left, right, nrow = 1, rel_widths = c(1, 1),
                        align = "h", axis = "tb"),
              ncol = 1, rel_heights = c(0.115, 0.885))
  })
  
  # one legend, built from a cell with all four series present and the
  # position set, then taken by name; a search across guide boxes can
  # return an empty one in ggplot2 3.5 and later
  legend_source <- make_cell(model_lines, trial_km, panel_levels[1],
                             "survival", "Overall survival", surv_scale,
                             "left") +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE),
           linetype = guide_legend(nrow = 2, byrow = TRUE)) +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 8),
          legend.key.width = unit(26, "pt"),
          # the two rows sit closer; the default key spacing leaves a
          # visible band between the model and Kaplan-Meier entries
          legend.key.spacing.y = unit(0, "pt"),
          legend.key.height = unit(12, "pt"),
          legend.background = element_rect(colour = "grey65", fill = NA,
                                           linewidth = 0.35),
          legend.margin = margin(1, 6, 4, 6),
          # the default box margin (about 5.5pt on every side) sits
          # outside legend.margin and is not touched by it; left alone it
          # adds space the explicit spacer rows below cannot see or trim
          legend.box.margin = margin(0, 0, 0, 0))
  
  plot_legend <- cowplot::get_plot_component(legend_source, "guide-box-bottom")
  if (inherits(plot_legend, "zeroGrob"))
    stop("no legend found; check colour and linetype are mapped")
  
  body <- plot_grid(plotlist = rows, ncol = 1)
  
  # spacers above and below the legend: above so it is not crowded by
  # the last "Time (years)", below so its border is off the figure edge
  whole_figure <- plot_grid(body, plot_legend, NULL, ncol = 1,
                            rel_heights = c(0.94, 0.063, 0.002))
  
  jpeg(file = figure_file, width = 7.6, height = 9.6, units = "in",
       res = 600, quality = 100)
  print(whole_figure)
  dev.off()
  
  cat("wrote", figure_file, "\n")
}


######################################################
#  Every model type, for both cuts.
######################################################

cuts <- tibble(
  cut_label = c("February 2017", "December 2017"),
  km_dataset = c("ALEX trial", "ALEX trial, December 2017"),
  files = list(
    c("Figures/Sup_figure_5_ph.jpg",
      "Figures/Figure_3_nonph.jpg",
      "Figures/Sup_figure_6_sep_arms.jpg"),
    c("Figures/Sup_figure_10_ph_dec.jpg",
      "Figures/Sup_figure_11_nonph_dec.jpg",
      "Figures/Sup_figure_12_sep_arms_dec.jpg")))

model_types <- c("PH", "NON-PH", "Separate_arms")

for (i in seq_len(nrow(cuts))) {
  for (j in seq_along(model_types)) {
    make_figure3(cut_wanted = cuts$cut_label[i],
                 km_dataset = cuts$km_dataset[i],
                 model_type = model_types[j],
                 figure_file = cuts$files[[i]][j])
  }
}