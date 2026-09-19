# Supplementary sensitivity figures: survival and hazard as one model
# specification is varied at a time.
#
# Three figures, one per parameter, all in the layout of Figure 3: rows
# of survival on the left and hazard on the right under a spanning
# title, with the trial's own Kaplan-Meier drawn thin underneath and the
# final analysis Kaplan-Meier dashed over the top as the validation.
# There are no manufacturer model rows here; every row is survextrap.
#
#   sigma   the prior on the baseline hazard's smoothing, Gamma(2, r)
#           for r in 1 and 3
#   tau     the prior on the hazard ratio's smoothing in the
#           non-proportional hazards model, Gamma(2, r) for r in 1, 5
#           and 10
#   df      the M-spline degrees of freedom, 6 and 10
#
# Everything not being varied is held at the base case: the
# non-proportional hazards model on ALEX plus PROFILE-1014 plus
# Flatiron, six degrees of freedom, sigma Gamma(2,1) and tau Gamma(2,10),
# February 2017 cut. So the base case appears as one row of every figure,
# which is what makes the rows comparable.
#
# On knot positions. These were asked for as a third sensitivity and are
# not here, because they are not in the fitted grid: Fit_survextrap.R
# sweeps df, hsd_rate and hrsd_rate but gives add_knots exactly one value
# per data cut, so base_model_all.rds holds no alternative placements to
# plot. The third figure is degrees of freedom alone. Varying positions
# as well means adding further add_knots vectors to the scenario grid and
# refitting, after which this script draws them with one more entry in
# the settings list below.
#
# Output files
#   Figures/Sup_figure_7_sigma.jpg
#   Figures/Sup_figure_8_tau.jpg
#   Figures/Sup_figure_9_df.jpg

library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(survival)
library(bshazard)

source("Functions/Survextrap_model.R")

x_max <- 30
hazard_max <- 0.45
hazard_lambda <- 1000          # as in Figure3_survival_hazard.R

cut_wanted <- "February 2017"
km_dataset <- "ALEX trial"

# everything held fixed unless it is the parameter being varied
base_spec <- list(model = "NON-PH", df = 6, hsd_rate = 1, hrsd_rate = 10,
                  datasets = "trial_and_all_maic")

km_label <- "Kaplan-Meier, ALEX trial\nApr 2025 data cut"

series_levels <- c("Alectinib (model)",
                   "Crizotinib (model)",
                   paste0("Alectinib (", km_label, ")"),
                   paste0("Crizotinib (", km_label, ")"))

series_colours <- setNames(
  c("#F8766D", "#00BFC4", "firebrick", "darkcyan"), series_levels)

series_linetypes <- setNames(c("solid", "solid", "22", "22"), series_levels)

fitted_colours <- c(Alectinib = "#F8766D", Crizotinib = "#00BFC4")

model_series <- function(trt)
  factor(paste0(trt, " (model)"), levels = series_levels)
valid_series <- function(trt)
  factor(paste0(trt, " (", km_label, ")"), levels = series_levels)


######################################################
#  What varies in each figure.
######################################################
#
# Greek letters are written as escapes so the file stays plain ASCII.

settings <- list(
  sigma = list(
    param = "hsd_rate",
    values = c(1, 3),
    label = function(v) paste0("\u03c3 ~ Gamma(2, ", v, ")"),
    file = "Figures/Sup_figure_7_sigma.jpg"),
  
  tau = list(
    param = "hrsd_rate",
    values = c(1, 5, 10),
    label = function(v) paste0("\u03c4 ~ Gamma(2, ", v, ")"),
    file = "Figures/Sup_figure_8_tau.jpg"),
  
  df = list(
    param = "df",
    values = c(6, 10),
    label = function(v) paste0("M-spline with ", v, " degrees of freedom"),
    file = "Figures/Sup_figure_9_df.jpg"))


######################################################
#  Read what the figures need.
######################################################

base_results <- readRDS("Data/base_model_all.rds")
km_data_all <- readRDS("Data/km_data_all.rds")
alex_2025 <- readRDS("Data/trial_data_apr_2025.rds")

validation_hazard <- bind_rows(lapply(unique(as.character(alex_2025$trt)),
                                      function(a) {
                                        d <- alex_2025 %>% filter(trt == a)
                                        fit <- bshazard(Surv(time, status) ~ 1, data = d,
                                                        lambda = hazard_lambda, verbose = FALSE)
                                        tibble(trt = a, time = fit$time, value = fit$hazard, variable = "hazard")
                                      }))

validation_survival <- km_data_all %>%
  filter(dataset == "ALEX trial, final analysis") %>%
  transmute(trt, time, value = surv, variable = "survival")

if (nrow(validation_survival) == 0)
  stop("no Kaplan-Meier rows for the final analysis in km_data_all.rds")

validation_lines <- bind_rows(validation_survival, validation_hazard) %>%
  mutate(series = valid_series(trt))

trial_km <- km_data_all %>%
  filter(dataset == km_dataset) %>%
  transmute(trt, time, value = surv)

if (nrow(trial_km) == 0)
  stop("no Kaplan-Meier rows for '", km_dataset, "' in km_data_all.rds")


######################################################
#  One cell: survival or hazard, for one row.
######################################################

make_cell <- function(model_lines, panel_wanted, variable_wanted,
                      y_label, y_scale, side) {
  
  d <- filter(model_lines, panel == panel_wanted, variable == variable_wanted)
  v <- filter(validation_lines, variable == variable_wanted)
  
  p <- ggplot() +
    geom_ribbon(data = d, aes(x = t, ymin = lower, ymax = upper, fill = trt),
                alpha = 0.20, colour = NA)
  
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
    scale_x_continuous("Time (years)", breaks = seq(0, x_max, by = 5),
                       expand = c(0, 0)) +
    y_scale +
    coord_cartesian(xlim = c(-0.3, x_max + 0.4), expand = FALSE) +
    ylab(y_label) +
    theme_classic() +
    theme(legend.position = "none",
          axis.text = element_text(size = 7),
          axis.title = element_text(size = 8),
          axis.title.x = element_text(margin = margin(t = 0.8)),
          axis.title.y = element_text(margin = margin(r = 0.8)),
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
#  One sensitivity figure.
######################################################

make_sensitivity <- function(setting) {
  
  param <- setting$param
  values <- setting$values
  panel_levels <- vapply(values, setting$label, character(1))
  
  # one block of rows per value of the parameter, with everything else
  # held at the base case
  model_lines <- bind_rows(lapply(seq_along(values), function(i) {
    spec <- base_spec
    spec[[param]] <- values[i]
    
    d <- base_results %>%
      filter(cut_label == cut_wanted,
             variable %in% c("survival", "hazard"), t > 0,
             model == spec$model,
             datasets == spec$datasets,
             df == spec$df,
             hsd_rate == spec$hsd_rate,
             hrsd_rate == spec$hrsd_rate | is.na(hrsd_rate))
    
    if (nrow(d) == 0)
      stop("no rows for ", param, " = ", values[i],
           ". Check base_model_all.rds contains model ", spec$model,
           ", datasets ", spec$datasets, ", df ", spec$df,
           ", hsd_rate ", spec$hsd_rate, ", hrsd_rate ", spec$hrsd_rate)
    
    d %>% transmute(panel = panel_levels[i], trt, t, median, lower, upper,
                    variable)
  })) %>%
    mutate(panel = factor(panel, levels = panel_levels),
           series = model_series(trt))
  
  surv_scale <- scale_y_continuous(limits = c(0, 1), labels = scales::percent)
  haz_scale <- scale_y_continuous(limits = c(0, hazard_max),
                                  oob = scales::oob_keep)
  
  rows <- lapply(seq_along(panel_levels), function(i) {
    pl <- panel_levels[i]
    left <- make_cell(model_lines, pl, "survival",
                      "Overall survival", surv_scale, "left")
    right <- make_cell(model_lines, pl, "hazard",
                       "Hazard", haz_scale, "right")
    plot_grid(row_title(paste0("(", letters[i], ") ", pl)),
              plot_grid(left, right, nrow = 1, rel_widths = c(1, 1),
                        align = "h", axis = "tb"),
              ncol = 1, rel_heights = c(0.115, 0.885))
  })
  
  legend_source <- make_cell(model_lines, panel_levels[1], "survival",
                             "Overall survival", surv_scale, "left") +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE),
           linetype = guide_legend(nrow = 2, byrow = TRUE)) +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 8),
          legend.key.width = unit(26, "pt"),
          legend.key.spacing.y = unit(0, "pt"),
          legend.key.height = unit(12, "pt"),
          legend.background = element_rect(colour = "grey65", fill = NA,
                                           linewidth = 0.35),
          legend.margin = margin(1, 6, 4, 6),
          legend.box.margin = margin(0, 0, 0, 0))
  
  plot_legend <- cowplot::get_plot_component(legend_source, "guide-box-bottom")
  if (inherits(plot_legend, "zeroGrob"))
    stop("no legend found; check colour and linetype are mapped")
  
  body <- plot_grid(plotlist = rows, ncol = 1)
  
  # the legend takes a fixed share of a figure whose height grows with
  # the number of rows, so its share shrinks as rows are added; Figure 3
  # has five rows in 9.6 inches, which is the same row height as this
  n <- length(panel_levels)
  fig_height <- 1.78 * n + 0.9
  legend_share <- 0.6 / fig_height
  
  whole_figure <- plot_grid(body, plot_legend, NULL, ncol = 1,
                            rel_heights = c(1 - legend_share - 0.002,
                                            legend_share, 0.002))
  
  jpeg(file = setting$file, width = 7.6, height = fig_height, units = "in",
       res = 600, quality = 100)
  print(whole_figure)
  dev.off()
  
  cat("wrote", setting$file, "with", n, "rows at",
      round(fig_height, 1), "inches\n")
}


######################################################
#  All three.
######################################################

cat("base case held fixed:",
    paste(names(base_spec), unlist(base_spec), sep = " = ", collapse = ", "),
    "\n")
cat("data cut:", cut_wanted, "\n\n")

for (s in settings) make_sensitivity(s)