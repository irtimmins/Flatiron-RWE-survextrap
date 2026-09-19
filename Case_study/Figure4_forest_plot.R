# Figure 4: forest plot of restricted mean survival at 10 and 30 years,
# company models against survextrap under increasing evidence, for each
# ALEX data cut.
#
# Company models (Exponential; Piecewise, Kaplan-Meier + exponential
# tail) sit alongside survextrap's three model types (PH, Non-PH,
# separate arms), each under three evidence sets: ALEX trial only, ALEX
# + PROFILE-1014, ALEX + PROFILE-1014 + Flatiron RWE. The Flatiron
# cohort is always the MAIC weighted one; that is stated here and in the
# caption rather than in the row label, matching Table 2 and Figure 3.
# Flatiron without MAIC weighting is left out, for the same reason as in
# Figure 3: it is a sensitivity comparison, not part of this figure's
# question.
#
# Company models have no evidence beyond the trial, so they are drawn at
# the "ALEX trial data only" tier alongside survextrap's own trial only
# fit: same evidence, different model, directly comparable. They do not
# appear at the two tiers with external data, since no such version of
# them exists.
#
# Two visual channels, read independently. Shape is the model:
# exponential is a cross, piecewise a star, and the three survextrap
# models are a circle, a square and a triangle. Fill is the evidence:
# empty for the trial alone, lightly filled once PROFILE-1014 is added,
# solid once Flatiron is added as well. So a solid triangle is the
# separate arms model with everything, and a hollow circle is PH on the
# trial alone. The cross and star have no interior, which is consistent
# rather than an omission: those two only ever appear at the empty tier.
#
# Four columns per horizon, matching Table 2: restricted mean survival
# for the control arm, for the active arm, the difference, and the
# difference again discounted at 3.5 per cent. Arm level figures stay
# undiscounted, as in the table, since discounting is what changes a
# decision on the difference rather than on the levels.
#
# The 10 year panels carry a red reference line at the restricted mean
# survival actually observed in the final analysis, 28 April 2025, with
# a light band one standard error either side. No model here has seen
# that data, so how close each estimate sits to that line is the
# validation. The discounted column gets the discounted validation
# figure, so the two are compared on the same basis. The 30 year panels
# have no such line: follow-up reaches 10.6 years, so there is nothing
# to check a 30 year extrapolation against, and the absence is the point
# rather than an omission.
#
# The priors are the calibrated set, one per scenario, matching Table 2.
# manuscript's stated base case (Non-PH, tau ~ Gamma(2,10)).
#
# Reference gridlines and axis limits are computed from the data in each
# panel rather than fixed in advance. The original figure had these
# hardcoded for a five and twenty year horizon with no company rows; both
# the horizons and the row set have changed enough that those numbers no
# longer describe the data.
#
# Within a horizon, the three panels (control, active, difference) share
# a y-axis ordering built once and reused, so a given row sits at the
# same height in all three. That depends on arranging by dataset then
# model, in that order, identically in every panel; if that arrangement
# ever differs between panels the rows will not line up and the figure
# will still render without warning, so this is worth checking by eye
# after any change here.
#
# Output files
#   February: Figures/Figure_4.jpg
#   December: Figures/Sup_figure_13_forest_dec.jpg
# The December figure number is a placeholder; slot it into the actual
# supplementary numbering once the rest of the supplement is finalised.

library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)

source("Functions/Survextrap_model.R")

# The grid no longer sweeps priors: one calibrated set per scenario.
scenario_value <- "primary"
time_points <- c(10, 30)

survextrap_dataset_levels <- c("trial_only", "trial_and_historic",
                               "trial_and_all_maic")

dataset_labels <- c("ALEX trial data only",
                    "ALEX + PROFILE-1014",
                    "ALEX + PROFILE-1014 +\nFlatiron RWE")

model_labels <- c("Exponential",
                  "Piecewise: Kaplan-Meier +\nexponential tail",
                  "survextrap, PH",
                  "survextrap, Non-PH",
                  "survextrap, Separate arms")

# Shape carries the model, fill carries the evidence. Only shapes 21 to
# 25 take a fill, so the three survextrap models get circle, square and
# triangle; the two company models get a cross and a star, which have no
# interior and need none, since they only ever appear at the trial only
# tier where the fill would be empty anyway.
model_shapes <- c(21, 22, 24)                       # PH, Non-PH, Separate arms
names(model_shapes) <- model_labels[3:5]
model_shapes <- c(setNames(c(4, 8), model_labels[1:2]), model_shapes)

# Empty, lightly filled, solid: the more evidence, the more filled in.
dataset_fills <- setNames(c("white", "#a6bddb", "#045a8d"), dataset_labels)

# The validation reference gets its own key. Mapping colour rather than
# setting it is what makes a legend appear at all, and colour is free
# here since the markers and intervals both use fixed colours.
validation_colour <- "#1b7837"
validation_alpha <- 0.12
validation_label <- paste0("Kaplan-Meier, ALEX trial,\n",
                           "April 2025 data cut,\n",
                           "RMST at 10-y (\u00b1 1 s.e.)")

# The key shows the band and the line together, as they appear on the
# panel. The band cannot be mapped to fill, which datasets already own,
# so the glyph draws it directly: a rectangle at the band's transparency
# with the line through its middle.
`%||%` <- function(a, b) if (is.null(a)) b else a

draw_key_validation <- function(data, params, size) {
  # both the band and the line run to 65 per cent of the key height,
  # centred: on the panel each spans the full height, but in a key this
  # small they read as overlong against the label
  key_frac <- 0.58
  grid::grobTree(
    grid::rectGrob(height = key_frac,
                   gp = grid::gpar(fill = scales::alpha(data$colour,
                                                        validation_alpha),
                                   col = NA)),
    grid::segmentsGrob(x0 = 0.5, x1 = 0.5,
                       y0 = (1 - key_frac) / 2, y1 = 1 - (1 - key_frac) / 2,
                       gp = grid::gpar(col = data$colour,
                                       lwd = (data$linewidth %||% 0.7) * .pt,
                                       alpha = data$alpha %||% 1)))
}


######################################################
#  Read and combine, one cut at a time.
######################################################

base_results <- readRDS("Data/base_model_all.rds")
exp_hybrid_rmst <- readRDS("Data/exp_hybrid_rmst.rds")

# Restricted mean survival from the final analysis, 28 April 2025, which
# no model here has seen. It exists at 5, 7 and 10 years only: follow-up
# reaches 10.6 years, so there is nothing to validate the 30 year panel
# against and none is drawn there.
#
# The band is one standard error either side of the estimate, not a 95
# per cent interval. It is a reference for the eye rather than a formal
# test, and at 1.96 standard errors it is wide enough to overlap almost
# every model estimate, which would tell the reader nothing.
validation_rmst <- readRDS("Data/validation_alex_apr_2025.rds")

get_validation <- function(time_point, trt_wanted, disc_wanted) {
  v <- validation_rmst %>%
    filter(t == time_point, disc_rate == disc_wanted)
  v <- if (is.null(trt_wanted)) filter(v, is.na(trt)) else filter(v, trt == trt_wanted)
  if (nrow(v) == 0) return(NULL)
  list(value = v$rmst[1],
       lower = v$rmst[1] - v$se[1],
       upper = v$rmst[1] + v$se[1])
}

get_rmst_data <- function(cut_wanted) {
  
  survextrap_rows <- base_results %>%
    filter(cut_label == cut_wanted, variable %in% c("rmst", "irmst"),
           t %in% time_points, scenario == scenario_value,
           datasets %in% survextrap_dataset_levels) %>%
    transmute(variable, trt, t, disc_rate, median, lower, upper,
              model = recode(model,
                             "PH" = "survextrap, PH",
                             "NON-PH" = "survextrap, Non-PH",
                             "Separate_arms" = "survextrap, Separate arms"),
              datasets = recode(datasets,
                                trial_only = "ALEX trial data only",
                                trial_and_historic = "ALEX + PROFILE-1014",
                                trial_and_all_maic = "ALEX + PROFILE-1014 +\nFlatiron RWE"))
  
  if (nrow(survextrap_rows) == 0)
    stop("no survextrap rmst rows for ", cut_wanted,
         ". Check cut_label, scenario and datasets against the values ",
         "actually present in base_model_all.rds")
  
  company_rows <- exp_hybrid_rmst %>%
    filter(data_cut == cut_wanted,
           variable %in% c("rmst", "irmst"), time %in% time_points) %>%
    transmute(variable, trt, t = time, disc_rate,
              median = value, lower, upper,
              model = ifelse(model == "exponential", "Exponential",
                             "Piecewise: Kaplan-Meier +\nexponential tail"),
              datasets = "ALEX trial data only")
  
  if (nrow(company_rows) == 0)
    stop("no company rmst rows for ", cut_wanted, " in exp_hybrid_rmst.rds")
  
  bind_rows(company_rows, survextrap_rows) %>%
    mutate(model = factor(model, levels = model_labels),
           datasets = factor(datasets, levels = dataset_labels))
}


# Axis limits shared across a pair of panels. Built from the estimates
# and the validation reference, padded and rounded to whole break points,
# rather than from the full intervals: one wide interval would otherwise
# set the scale for the whole row. Intervals past the edge get an
# arrowhead instead.
#
# To fix a range by hand instead, put it in axis_override, keyed by
# horizon and then "arms" or "diff".
axis_override <- list("10" = list(arms = c(3, 8),
                                  diff = c(-2.5, 3.8)),
                      "30" = list(arms = c(0, 20),
                                  diff = c(-11, 12)))

# Tick spacing, same keys as axis_override. Where a step is given the
# ticks are whole multiples of it across the range; where it is not,
# pretty() chooses. The 10 year difference panels tick every year, since
# a year of life is the unit the reader is weighing there and the range
# is narrow enough to show them all without crowding.
axis_step <- list("10" = list(diff = 1))

shared_breaks <- function(x_limits, time_point, key) {
  step <- axis_step[[as.character(time_point)]][[key]]
  brk <- if (!is.null(step)) {
    seq(ceiling(x_limits[1] / step) * step,
        floor(x_limits[2] / step) * step, by = step)
  } else {
    pretty(x_limits, n = 4)
  }
  brk[brk >= x_limits[1] & brk <= x_limits[2]]
}

shared_limits <- function(data, variable_wanted, time_point, disc_rates,
                          symmetric = FALSE, extra = NULL, pad_frac = 0.35) {
  
  key <- if (variable_wanted == "irmst") "diff" else "arms"
  fixed <- axis_override[[as.character(time_point)]][[key]]
  if (!is.null(fixed)) return(fixed)
  
  sub <- data %>%
    filter(variable == variable_wanted, t == time_point,
           disc_rate %in% disc_rates)
  
  centre <- range(c(sub$median, extra), na.rm = TRUE)
  pad <- max(diff(centre) * pad_frac, 0.5)
  lim <- centre + c(-pad, pad)
  
  if (symmetric) {
    m <- max(abs(lim))
    lim <- c(-m, m)
  }
  
  brk <- pretty(lim, n = 4)
  c(min(brk), max(brk))
}


######################################################
#  One forest panel: control, active or the difference, one horizon.
######################################################

one_forest <- function(data, variable_wanted, trt_wanted, time_point,
                       disc_wanted, x_label, x_limits, zero_line = FALSE,
                       validation = NULL) {
  
  d <- data %>%
    filter(variable == variable_wanted, t == time_point,
           disc_rate == disc_wanted) %>%
    { if (!is.null(trt_wanted)) filter(., trt == trt_wanted) else . } %>%
    arrange(datasets, model) %>%
    mutate(group_id = -row_number(),
           lo_clip = pmax(lower, x_limits[1]),
           hi_clip = pmin(upper, x_limits[2]),
           lo_out = lower < x_limits[1],
           hi_out = upper > x_limits[2])
  
  if (nrow(d) == 0)
    stop("no rows for ", variable_wanted, " at t = ", time_point,
         ", disc_rate = ", disc_wanted,
         if (!is.null(trt_wanted)) paste0(", trt = ", trt_wanted) else "")
  
  x_breaks <- shared_breaks(
    x_limits, time_point,
    if (variable_wanted == "irmst") "diff" else "arms")
  
  y_limits <- c(min(d$group_id) - 0.5, max(d$group_id) + 0.5)
  
  p <- ggplot(d) +
    theme_classic() +
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.line.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 8),
          axis.title.x = element_text(size = 10),
          legend.background = element_blank(),
          legend.box.background = element_rect(colour = "black"),
          # gap between the Model, Datasets and Validation blocks; the
          # ggplot default here is 11pt
          legend.spacing.y = unit(5.5, "pt"),
          legend.key.spacing.y = unit(3, "pt"),
          legend.text = element_text(size = 8)) +
    geom_vline(xintercept = x_breaks, colour = "gray70", alpha = 0.4)
  
  if (zero_line)
    p <- p + geom_vline(xintercept = 0, colour = "gray20",
                        linewidth = 0.6, alpha = 0.75)
  
  # The validation reference: a band one standard error either side and a
  # line at the estimate, drawn first so the model estimates sit on top
  # of it rather than behind it. colour is mapped, not set, so the line
  # earns a legend key of its own.
  if (!is.null(validation))
    p <- p +
    annotate("rect", xmin = validation$lower, xmax = validation$upper,
             ymin = -Inf, ymax = Inf, fill = validation_colour,
             alpha = validation_alpha) +
    geom_vline(data = data.frame(v = validation$value),
               aes(xintercept = v, colour = validation_label),
               linewidth = 0.7, linetype = "solid", alpha = 0.85,
               key_glyph = draw_key_validation)
  
  p +
    # intervals as two segments from the estimate outwards, so each side
    # can be capped independently; an arrowhead means the interval runs
    # past the edge of the panel rather than ending there
    geom_segment(data = filter(d, !lo_out),
                 aes(x = median, xend = lo_clip, y = group_id, yend = group_id),
                 colour = "gray35", linewidth = 0.5) +
    geom_segment(data = filter(d, !hi_out),
                 aes(x = median, xend = hi_clip, y = group_id, yend = group_id),
                 colour = "gray35", linewidth = 0.5) +
    geom_segment(data = filter(d, lo_out),
                 aes(x = median, xend = lo_clip, y = group_id, yend = group_id),
                 colour = "gray35", linewidth = 0.5,
                 arrow = arrow(length = unit(0.06, "in"), type = "open",
                               ends = "last")) +
    geom_segment(data = filter(d, hi_out),
                 aes(x = median, xend = hi_clip, y = group_id, yend = group_id),
                 colour = "gray35", linewidth = 0.5,
                 arrow = arrow(length = unit(0.06, "in"), type = "open",
                               ends = "last")) +
    geom_point(aes(x = median, y = group_id, shape = model, fill = datasets),
               colour = "black", stroke = 0.6, size = 2.4) +
    scale_x_continuous(x_label, breaks = x_breaks) +
    scale_y_continuous(limits = y_limits) +
    coord_cartesian(xlim = x_limits) +
    scale_shape_manual("Model", values = model_shapes, drop = FALSE) +
    scale_fill_manual("Datasets", values = dataset_fills, drop = FALSE) +
    scale_colour_manual("Validation",
                        values = setNames(validation_colour, validation_label)) +
    guides(
      # a neutral mid fill in the shape key, so the shapes are legible
      # there without implying any particular evidence tier
      shape = guide_legend(order = 1,
                           override.aes = list(fill = "gray70", size = 2.6)),
      # a circle in the fill key, so the three fills are compared on one
      # outline rather than on whichever model happens to come first
      fill = guide_legend(order = 2,
                          override.aes = list(shape = 21, size = 2.6,
                                              colour = "black")),
      colour = guide_legend(order = 3,
                            override.aes = list(alpha = 0.85),
                            theme = theme(legend.key.height = unit(15, "pt"),
                                          legend.key.width = unit(14, "pt"))))
}


######################################################
#  Both horizons, three panels each, for one cut.
######################################################

make_figure4 <- function(cut_wanted, figure_file) {
  
  d <- get_rmst_data(cut_wanted)
  
  forests <- lapply(time_points, function(tt) {
    
    # Limits are shared within a pair, so the two arm panels are read on
    # one scale and the two difference panels on another. They are driven
    # by the estimates and the validation reference, not by the full
    # intervals: a single wide interval would otherwise set the scale for
    # every panel and squash everything else into the middle. Intervals
    # that run past the edge are arrowed rather than expanding the axis.
    # get_validation returns NULL at 30 years, where there is nothing to
    # validate against; NULL$value is NULL in R, and c() drops it, so the
    # limits then come from the estimates alone
    arms_limits <- shared_limits(
      d, "rmst", tt, 0,
      extra = c(get_validation(tt, "Crizotinib", 0)$value,
                get_validation(tt, "Alectinib", 0)$value))
    
    diff_limits <- shared_limits(
      d, "irmst", tt, c(0, 0.035), symmetric = TRUE,
      extra = c(get_validation(tt, NULL, 0)$value,
                get_validation(tt, NULL, 0.035)$value))
    
    cat(sprintf("  %2d years: arms %.1f to %.1f, difference %.1f to %.1f\n",
                tt, arms_limits[1], arms_limits[2],
                diff_limits[1], diff_limits[2]))
    
    control <- one_forest(d, "rmst", "Crizotinib", tt, 0,
                          paste0("RMST at ", tt, "-y\nfor Crizotinib"),
                          x_limits = arms_limits,
                          validation = get_validation(tt, "Crizotinib", 0))
    active <- one_forest(d, "rmst", "Alectinib", tt, 0,
                         paste0("RMST at ", tt, "-y\nfor Alectinib"),
                         x_limits = arms_limits,
                         validation = get_validation(tt, "Alectinib", 0))
    diff <- one_forest(d, "irmst", NULL, tt, 0,
                       paste0("\u0394RMST at ", tt, "-y"),
                       x_limits = diff_limits, zero_line = TRUE,
                       validation = get_validation(tt, NULL, 0))
    diff_disc <- one_forest(d, "irmst", NULL, tt, 0.035,
                            paste0("\u0394RMST at ", tt,
                                   "-y,\nwith discounting"),
                            x_limits = diff_limits, zero_line = TRUE,
                            validation = get_validation(tt, NULL, 0.035))
    
    # taken from the 10 year row deliberately: the 30 year panels have no
    # validation line, so their legend would be missing that third key
    forest_legend <- cowplot::get_legend(active)
    
    forest_plots <- plot_grid(control + theme(legend.position = "none"),
                              active + theme(legend.position = "none"),
                              diff + theme(legend.position = "none"),
                              diff_disc + theme(legend.position = "none"),
                              align = "h", rel_widths = c(1, 1, 1, 1),
                              nrow = 1)
    
    list(plots = forest_plots, legend = forest_legend)
  })
  
  # Spacers either side of the legend rather than a margin on the legend
  # itself: legend.box.margin is applied inside the plot the legend came
  # from, and get_legend does not reliably carry it across, which is why
  # the legend kept sitting flush against the panels on one side and the
  # figure edge on the other. Explicit NULL columns in the grid always
  # take the width they are given.
  plot_all <- plot_grid(
    plot_grid(NULL, forests[[1]]$plots, NULL, forests[[2]]$plots,
              rel_heights = c(0.05, 0.5, 0.05, 0.5),
              labels = c("(a)", "", "(b)", ""), label_size = 12,
              label_x = 0, ncol = 1, align = "v"),
    NULL,
    forests[[1]]$legend,
    NULL,
    rel_widths = c(1, 0.03, 0.30, 0.03), nrow = 1)
  
  jpeg(file = figure_file, width = 7.4, height = 6.5, units = "in",
       res = 600, quality = 100)
  print(plot_all)
  dev.off()
  
  cat("wrote", figure_file, "\n")
}


######################################################
#  Both cuts.
######################################################

make_figure4("February 2017", "Figures/Figure_4.jpg")
make_figure4("December 2017", "Figures/Sup_figure_13_forest_dec.jpg")