# Figure 2. Kaplan-Meier curves.
#   (a) overall survival by arm in the ALEX trial, February 2017 cut
#   (b) overall survival for first line crizotinib across data sources
#
# Rewritten to read the reconstructed pseudo-IPD rather than the live
# Flatiron cohort, so it runs without database access. The Flatiron rows
# carry a weight column: the MAIC cohort is 271 rows at a tenth of a
# patient each, giving the published effective sample size of 27.1, so
# survfit has to be given the weights or the risk table will read 271.
#
# Getting the two panels to line up, which took several goes. The label
# used to be the plot title of each panel, positioned with a negative
# hjust, so where it landed depended on that panel's own internal
# layout. Those layouts differ: ggsurvfit sizes each risk-table label
# column to its own widest row label, and panel (a)'s labels are short
# arm names while panel (b)'s are multi-line data source names.
#
# Three things fix it. The label sits in its own row as a separate plot,
# as in Figure 3, so the grid places it and it owes nothing to the panel
# beneath. The row labels are padded to a common width with non-breaking
# spaces, which makes the label column the same width in both. And both
# gtables are then set to the element-wise maximum of their column
# widths as a backstop.
#
# The padding character matters and is the thing that took longest to
# find. Padding with ordinary spaces does nothing at all: a device drops
# trailing whitespace when it measures a string, so a padded label is
# measured at the width of its text. A non-breaking space is measured
# like any other glyph.
#
# cowplot::align_plots does not help here either, for a related reason:
# it works by finding named cells such as "axis-l" in a gtable, and
# ggsurvfit_build returns a composite of the plot stacked on the risk
# table, where those names sit a level down and out of its reach. The
# same nesting is why matching the outer widths alone was not enough.
#
# Inputs, in Data/
#   trial_data.rds                                ALEX, February 2017
#   trial_IPD_OS_PROFILE_1014_Nov_2016.csv        historic trial
#   trial_IPD_OS_flatiron.csv                     RWE, both cohorts

library(dplyr)
library(readr)
library(survival)
library(ggplot2)
library(ggsurvfit)
library(cowplot)

# Every risk-table row label in both panels is padded to this many
# characters, which is what makes the label column the same width in
# both and so lines the two y axes up.
#
# The padding is non-breaking spaces, not ordinary ones. An ordinary
# trailing space is dropped when the device measures a string, so padding
# with spaces changes nothing; a non-breaking space is measured like any
# other glyph. Written as an escape so the file stays plain ASCII.
label_width <- 16
nbsp <- "\u00a0"

pad_lines <- function(x, width = label_width) {
  vapply(x, function(s) {
    lines <- strsplit(s, "\n", fixed = TRUE)[[1]]
    padded <- vapply(lines, function(l)
      if (nchar(l) >= width) l
      else paste0(l, strrep(nbsp, width - nchar(l))), character(1))
    paste(padded, collapse = "\n")
  }, character(1), USE.NAMES = FALSE)
}

# Position of the "(a)" and "(b)" labels within their own rows.
#
# These sit in plots of their own, so the coordinates are simply a
# fraction of that row: x runs 0 at the left edge to 1 at the right, y
# runs 0 at the bottom to 1 at the top. Raise label_x to move both
# labels right; lower label_y to move them down.
#
# The label plots use clip = "off", so text is still drawn where it
# crosses the edge of its row rather than being cut off, and the label
# rows are given more height than the text needs. Both are needed: clip
# alone would let a lowered label run onto the plot beneath, and extra
# height alone would still cut the text once label_y went far enough
# down. If you want the label lower than this allows, raise the two row
# heights in rel_heights near the end of the script rather than pushing
# label_y further.
label_x <- 0.03
label_y <- 0.5
label_size <- 5

# Position of the y axis title, per panel.
#
# The axis title is rotated, and its justification arguments follow the
# text's own frame rather than the page's. So vjust moves the title
# perpendicular to the axis, which is horizontally on the page, and hjust
# moves it along the axis, which is vertically. hjust is therefore left
# at 0.5 in both panels, centring the title against its axis, and is not
# a horizontal control at all.
#
# Two dials move the title horizontally, and they work in opposite
# directions:
#
#   y_title_vjust    larger moves the title toward the axis, to the
#                    right on the page; smaller, including negative,
#                    moves it away to the left. This is the main dial.
#   y_title_margin   margin(r = ) adds a gap on the axis side, so larger
#                    pushes the title left. Use it for fine adjustment
#                    once vjust is roughly right, or set it to 0 and
#                    work with vjust alone.
#
# If the title ends up so far left that it runs off the figure, raise
# the l value in panel_margin below to make room for it.
panel_margin <- margin(t = 5.5, r = 5.5, b = 5.5, l = 17)

y_title_vjust_a <- 5
y_title_vjust_b <- 6

y_title_margin_a <- margin(r = 10)
y_title_margin_b <- margin(r = 0)


# ---- data ---------------------------------------------------------------

trial_data <- readRDS("Data/trial_data.rds")

profile <- read_csv("Data/trial_IPD_OS_PROFILE_1014_Nov_2016.csv",
                    show_col_types = FALSE) %>%
  transmute(time, status, weight = 1, trt = "Crizotinib",
            dataset = "PROFILE-1014")

flatiron <- read_csv("Data/trial_IPD_OS_flatiron.csv", show_col_types = FALSE)

alex_criz <- trial_data %>%
  filter(trt == "Crizotinib") %>%
  transmute(time, status, weight = 1, trt = "Crizotinib", dataset = "ALEX")

rwe_unweighted <- flatiron %>%
  filter(cohort == "flatiron") %>%
  transmute(time, status, weight, trt = "Crizotinib",
            dataset = "Flatiron RWE unweighted")

rwe_maic <- flatiron %>%
  filter(cohort == "flatiron_maic") %>%
  transmute(time, status, weight, trt = "Crizotinib",
            dataset = "Flatiron RWE MAIC weighted")

surv_by_arm <- trial_data %>%
  mutate(trt = factor(trt, levels = c("Alectinib", "Crizotinib"),
                      labels = pad_lines(c("Alectinib", "Crizotinib"))))

surv_by_source <- bind_rows(alex_criz, profile, rwe_unweighted, rwe_maic) %>%
  mutate(dataset = factor(dataset,
                          levels = c("ALEX",
                                     "PROFILE-1014",
                                     "Flatiron RWE unweighted",
                                     "Flatiron RWE MAIC weighted"),
                          labels = pad_lines(c(
                            "ALEX,\nCrizotinib",
                            "PROFILE-1014",
                            "Flatiron RWE\nunweighted",
                            "Flatiron RWE\nMAIC weighted"))))

cat("effective sample sizes entering panel (b)\n")
surv_by_source %>%
  group_by(dataset) %>%
  summarise(rows = n(), effective_n = round(sum(weight), 1),
            deaths = round(sum(weight * status), 1),
            longest = round(max(time), 2), .groups = "drop") %>%
  as.data.frame() %>%
  print(row.names = FALSE)


# ---- the panel label, in a row of its own --------------------------------
#
# As in Figure 3: a plot containing only text, so the grid decides where
# it sits rather than the panel below it.

panel_label <- function(label) {
  ggplot() +
    annotate("text", x = label_x, y = label_y, label = label, hjust = 0,
             vjust = 0.5, size = label_size, fontface = "bold") +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    # clip = "off" so the text is still drawn where it crosses the edge
    # of the panel. Without it, lowering label_y pushes the bottom of the
    # text past y = 0 and that part is simply cut away.
    coord_cartesian(clip = "off") +
    theme_void() +
    theme(plot.margin = unit(c(0.04, 0.1, 0.04, 0.04), "cm"))
}


# ---- panel (a), by arm ---------------------------------------------------

plot_2a <-
  survfit2(Surv(time, status) ~ trt, data = surv_by_arm) %>%
  ggsurvfit(linewidth = 1) +
  add_confidence_interval() +
  add_risktable(
    risktable_height = 0.18,
    risktable_stats = c("{round(n.risk, digits = 1)}"),
    stats_label = c("Number at risk"),
    theme = list(
      theme_risktable_default(),
      theme(
        axis.text.y = element_text(hjust = 0, margin = margin(r = -5, l = -5)),
        plot.caption = element_text(hjust = -5),
        plot.title.position = "plot"
      )
    )
  ) +
  theme(axis.text.x = element_text(hjust = 1)) +
  scale_ggsurvfit() +
  theme_classic() +
  theme(legend.position.inside = c(0.12, 0.6),
        legend.position = "inside",
        legend.text = element_text(size = 8),
        legend.key.spacing.y = unit(0, "pt"),
        legend.box.spacing = unit(-6, "pt"),
        plot.margin = panel_margin,
        axis.title.y = element_text(vjust = y_title_vjust_a,
                                    margin = y_title_margin_a,
                                    hjust = 0.5)) +
  scale_colour_manual(NULL, values = c("#4393c3", "#d6604d")) +
  scale_fill_manual(NULL, values = c("#4393c3", "#d6604d")) +
  scale_y_continuous("Overall Survival", limits = c(0, 1),
                     labels = scales::percent) +
  scale_x_continuous("Time (years)", breaks = seq(0, 3.5, by = 0.25),
                     limits = c(0, 3))

plot_2a_build <- ggsurvfit_build(plot_2a)


# ---- panel (b), by data source -------------------------------------------

plot_2b <-
  survfit2(Surv(time, status) ~ dataset, data = surv_by_source,
           weights = surv_by_source$weight) %>%
  ggsurvfit(linewidth = 1) +
  add_confidence_interval() +
  add_risktable(
    risktable_height = 0.34,
    risktable_stats = c("{round(n.risk, digits = 1)}"),
    stats_label = c("Number at risk"),
    theme = list(
      theme_risktable_default(),
      theme(
        axis.text.y = element_text(hjust = 0, margin = margin(r = -5, l = -5)),
        plot.caption = element_text(hjust = -5),
        plot.title.position = "plot"
      )
    )
  ) +
  scale_ggsurvfit() +
  theme_classic() +
  theme(legend.position.inside = c(0.115, 0.35),
        legend.position = "inside",
        legend.text = element_text(size = 8),
        legend.spacing.y = unit(0, "pt"),
        legend.spacing.x = unit(0, "pt"),
        legend.key.spacing.y = unit(1.2, "pt"),
        legend.box.spacing = unit(-6, "pt"),
        plot.margin = panel_margin,
        axis.title.y = element_text(vjust = y_title_vjust_b,
                                    margin = y_title_margin_b,
                                    hjust = 0.5)) +
  scale_x_continuous("Time (years)", breaks = seq(0, 5.5, by = 0.5)) +
  scale_y_continuous("Overall Survival", limits = c(0, 1),
                     labels = scales::percent) +
  scale_colour_manual(NULL,
                      values = c("#d6604d", "#7CAE00", "#C77CFF", "deepskyblue3")) +
  scale_fill_manual(NULL,
                    values = c("#d6604d", "#7CAE00", "#C77CFF", "deepskyblue3"))

plot_2b_build <- ggsurvfit_build(plot_2b)


# ---- assemble -------------------------------------------------------------
#
# The two panels are lined up by setting both gtables to the element-wise
# maximum of their column widths. Every column then measures the same in
# both, including whichever one holds the risk-table row labels, so the
# plotting areas start and end at the same place.
#
# Two earlier attempts did not work, for reasons worth recording. Padding
# the row labels to a common character count failed because most
# graphics devices do not measure trailing whitespace, so a padded label
# is measured at the width of the text alone. align_plots failed because
# it works by finding named cells such as "axis-l" in a gtable, and
# ggsurvfit_build returns a composite of the plot stacked on the risk
# table, where those names are nested a level down and out of its reach.

as_gtable <- function(p) {
  if (inherits(p, "gtable")) return(p)
  if (inherits(p, "patchwork")) return(patchwork::patchworkGrob(p))
  ggplot2::ggplotGrob(p)
}

match_widths <- function(a, b) {
  ga <- as_gtable(a)
  gb <- as_gtable(b)
  if (length(ga$widths) != length(gb$widths)) {
    warning("the two panels have ", length(ga$widths), " and ",
            length(gb$widths), " width columns, so they cannot be matched ",
            "column by column; they are left as they are and will not line up")
    return(list(ga, gb))
  }
  common <- grid::unit.pmax(ga$widths, gb$widths)
  ga$widths <- common
  gb$widths <- common
  list(ga, gb)
}

aligned <- match_widths(plot_2a_build, plot_2b_build)

# panel (b) gets more of the height: its risk table has four rows to
# panel (a)'s two, and at the old split its bottom row ran off the page
figure_2 <- plot_grid(
  panel_label("(a)"), aligned[[1]],
  panel_label("(b)"), aligned[[2]],
  # the label rows get more height than they need for the text alone, so
  # that lowering label_y moves the label within its own row rather than
  # straight onto the plot beneath it
  ncol = 1, rel_heights = c(0.09, 1, 0.09, 1.35))

print(figure_2)

dir.create("Figures", showWarnings = FALSE)
jpeg(file = "Figures/Figure_2.jpg",
     width = 7.5,
     height = 10.0,
     units = 'in',
     res = 600,
     quality = 100)
print(figure_2)
dev.off()