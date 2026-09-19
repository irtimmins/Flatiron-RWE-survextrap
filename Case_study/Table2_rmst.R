# Table 2: restricted mean survival and difference in restricted mean
# survival, as a formatted Word table.
#
# The February 2017 cut only, which is the company's original base case
# in NICE TA536. The December cut was fitted too and is in
# base_model_all.rds, but it is not shown here.
#
# Three panels, each a heading row with its models indented beneath.
#
#   Validation            the Kaplan-Meier estimate from the trial's
#                         final overall survival analysis, which nothing
#                         in the table below was fitted to
#   Manufacturer models   the two survival models used in the appraisal
#   survextrap models     by the evidence each was given
#
# The survextrap rows are named by model type alone, PH, Non-PH and
# separate arms, since the panel heading already says they are
# survextrap and repeating it in every row wastes the width.
#
# Eight value columns: 10 and 30 years, each giving control, active, the
# difference, and the difference again discounted at 3.5 per cent.
# Discounting is shown only on the difference because that is where it
# can change a decision; the arm level figures stay undiscounted so the
# survival scale is readable.
#
# Thirty years is the lifetime horizon NICE TA536 itself used, which is
# what Fit_exp_and_hybrid.R reproduces the published life years against,
# so this table matches the appraisal directly rather than using a
# horizon of its own.
#
# The validation row has no thirty year entry. The final overall
# survival analysis has 10.6 years of follow-up, so a thirty year
# extrapolation is checked against nothing at all. That cell stays
# visibly empty rather than being filled or the row dropped.
#
# Widths are set for A4 landscape: 1.45 and 1.50 inches for the two
# label columns and 0.92 for each of the eight value columns comes to
# 10.31 against 10.69 usable. Two labels are longer than their column
# and are broken across two lines here rather than left to wrap where
# they happen to fall.

library(dplyr)
library(tidyr)
library(flextable)
library(officer)

TABLE_FONT_SIZE <- 7
INDENT <- 14

round1 <- function(x) sprintf("%.2f", round(x, digits = 2))
fmt <- function(v, lo, hi) paste0(round1(v), " (", round1(lo), ",", round1(hi), ")")

horizons <- c(10, 30)
cut_wanted <- "February 2017"

# consistent with the figures, where the same curve carries this
# description; broken explicitly rather than left to wrap where it falls
validation_source <- "ALEX trial,\nApr 2025 data cut\n(final OS)"

piecewise_label <- "Piecewise: Kaplan-Meier\n+ exponential tail"
flatiron_label <- "ALEX + PROFILE-1014\n+ Flatiron"


# ---- long format from each source ---------------------------------------

keep_disc <- function(d)
  filter(d, disc == 0 | (disc > 0 & trt_label == "contrast"))

trt_label_of <- function(trt)
  case_when(trt == "Crizotinib" ~ "control",
            trt == "Alectinib" ~ "active",
            is.na(trt) ~ "contrast")


survextrap_rows <- function(cut_label_wanted) {
  f <- "Data/base_model_all.rds"
  if (!file.exists(f)) {
    warning(f, " not found; no survextrap rows")
    return(tibble())
  }
  d <- readRDS(f)
  
  if (!all(horizons %in% d$t))
    warning("base_model_all.rds does not contain t = ",
            paste(setdiff(horizons, unique(d$t)), collapse = ", "),
            "; check get_rmst_survextrap() was run with times = c(10, 30)")
  
  # the file carries an internal key (data_cut) and a display label
  # (cut_label); older versions may have only one, or neither
  if ("cut_label" %in% names(d)) {
    # already what this function expects
  } else if ("data_cut" %in% names(d)) {
    warning("base_model_all.rds has data_cut but no cut_label; mapping ",
            "feb_2017/dec_2017 to their display names")
    d <- mutate(d, cut_label = case_when(
      data_cut == "feb_2017" ~ "February 2017",
      data_cut == "dec_2017" ~ "December 2017",
      TRUE ~ NA_character_))
  } else {
    if (cut_label_wanted != "February 2017") {
      warning("base_model_all.rds has no cut_label or data_cut column, so ",
              "it is assumed to be February 2017 only; no survextrap rows ",
              "for ", cut_label_wanted)
      return(tibble())
    }
    d <- mutate(d, cut_label = "February 2017")
  }
  if (!"disc_rate" %in% names(d)) d <- mutate(d, disc_rate = 0)
  
  d %>%
    filter(cut_label == cut_label_wanted,
           variable %in% c("rmst", "irmst"),
           t %in% horizons,
           datasets != "trial_and_all_unweighted",
           df == 6, hsd_rate == 1, (hrsd_rate == 10 | is.na(hrsd_rate))) %>%
    transmute(
      # model type alone: the panel heading already says survextrap
      model_label = case_when(model == "PH" ~ "PH",
                              model == "NON-PH" ~ "Non-PH",
                              model == "Separate_arms" ~ "Separate arms"),
      data_label = case_when(
        datasets == "trial_only" ~ "ALEX only",
        datasets == "trial_and_historic" ~ "ALEX + PROFILE-1014",
        datasets == "trial_and_all_maic" ~ flatiron_label),
      horizon = t, disc = disc_rate,
      trt_label = trt_label_of(trt),
      value = fmt(median, lower, upper)) %>%
    keep_disc()
}


comparator_rows <- function(cut_label) {
  readRDS("Data/exp_hybrid_rmst.rds") %>%
    filter(data_cut == cut_label, time %in% horizons) %>%
    transmute(
      model_label = ifelse(model == "exponential", "Exponential",
                           piecewise_label),
      data_label = "ALEX only",
      horizon = time, disc = disc_rate,
      trt_label = trt_label_of(trt),
      value = fmt(value, lower, upper)) %>%
    keep_disc()
}


validation_rows <- function() {
  f <- "Data/validation_alex_apr_2025.rds"
  if (!file.exists(f)) {
    warning(f, " not found; no validation row")
    return(tibble())
  }
  readRDS(f) %>%
    filter(t == 10) %>%
    transmute(model_label = "Kaplan-Meier",
              data_label = validation_source,
              horizon = t, disc = disc_rate,
              trt_label = trt_label_of(trt),
              value = fmt(rmst, lower, upper)) %>%
    keep_disc()
}


# ---- one block, wide -----------------------------------------------------

col_key <- function(horizon, disc, trt_label)
  paste0("h", horizon, "_", trt_label, ifelse(disc > 0, "_disc", ""))

value_cols <- c("h10_control", "h10_active", "h10_contrast", "h10_contrast_disc",
                "h30_control", "h30_active", "h30_contrast", "h30_contrast_disc")

# Within the survextrap block the data source is the primary sort, so
# the three models sit together for each set of evidence and the effect
# of adding evidence can be read down the column.
data_order <- c("ALEX only", "ALEX + PROFILE-1014", flatiron_label)
survextrap_order <- c("PH", "Non-PH", "Separate arms")
company_order <- c("Exponential", piecewise_label)

to_wide <- function(long) {
  if (!nrow(long)) return(tibble())
  wide <- long %>%
    mutate(col = col_key(horizon, disc, trt_label)) %>%
    select(model_label, data_label, col, value) %>%
    pivot_wider(names_from = col, values_from = value)
  for (cn in value_cols) if (!cn %in% names(wide)) wide[[cn]] <- NA_character_
  wide %>%
    mutate(across(all_of(value_cols), ~ifelse(is.na(.x), "-", .x))) %>%
    select(model_label, data_label, all_of(value_cols))
}

valid <- to_wide(validation_rows())

company <- to_wide(comparator_rows(cut_wanted)) %>%
  arrange(match(model_label, company_order))

survex <- to_wide(survextrap_rows(cut_wanted)) %>%
  arrange(match(data_label, data_order), match(model_label, survextrap_order))

cat("validation:          ", nrow(valid), "row\n")
cat("manufacturer models: ", nrow(company), "rows\n")
cat("survextrap models:   ", nrow(survex), "rows\n")


# ---- stitch into one table ------------------------------------------------

col_names <- c("Model", "Data sources",
               "Control", "Active", "Difference", "Difference\n(Discounted)",
               "Control ", "Active ", "Difference ", "Difference\n(Discounted) ")

named <- function(d) { names(d) <- col_names; as.data.frame(d) }

blank_row <- function(first_cell = "") {
  r <- as.data.frame(as.list(rep("", length(col_names))),
                     stringsAsFactors = FALSE, check.names = FALSE)
  names(r) <- col_names
  r[[1]] <- first_cell
  r
}

group_row <- blank_row("")
group_row[["Control"]]  <- "10 years"
group_row[["Control "]] <- "30 years"

titles <- as.data.frame(as.list(col_names), stringsAsFactors = FALSE,
                        check.names = FALSE)
names(titles) <- col_names

combined <- rbind(group_row, titles,
                  blank_row("Validation"), named(valid),
                  blank_row("Manufacturer models (NICE TA536)"), named(company),
                  blank_row("survextrap models"), named(survex))

n_v <- nrow(valid); n_c <- nrow(company); n_s <- nrow(survex)

row_group <- 1
row_tit   <- 2
row_pan_v <- 3
rows_v    <- row_pan_v + seq_len(n_v)
row_pan_c <- row_pan_v + n_v + 1
rows_c    <- row_pan_c + seq_len(n_c)
row_pan_s <- row_pan_c + n_c + 1
rows_s    <- row_pan_s + seq_len(n_s)

panel_rows <- c(row_pan_v, row_pan_c, row_pan_s)
model_rows <- c(rows_v, rows_c, rows_s)

ft <- flextable(combined)
ft <- delete_part(ft, part = "header")

ft <- merge_at(ft, i = row_group, j = 1:2, part = "body")
ft <- merge_at(ft, i = row_group, j = 3:6, part = "body")
ft <- merge_at(ft, i = row_group, j = 7:10, part = "body")
for (r in panel_rows)
  ft <- merge_at(ft, i = r, j = seq_along(col_names), part = "body")

ft <- bold(ft, i = c(row_group, row_tit, panel_rows), part = "body")
ft <- align(ft, part = "body", align = "center")
ft <- align(ft, j = 1:2, align = "left", part = "body")
ft <- align(ft, i = panel_rows, align = "left", part = "body")

# Black rules only for the frame and the column titles. Everything
# inside is dotted or faint, so the panels read as parts of one table.
rule <- fp_border(color = "black", width = 1)
ft <- border_outer(ft, border = rule, part = "body")
ft <- hline(ft, i = row_group, border = rule, part = "body")
ft <- hline(ft, i = row_tit, border = rule, part = "body")

dotted <- fp_border(color = "grey45", width = 0.9, style = "dotted")
if (n_v) ft <- hline(ft, i = max(rows_v), border = dotted, part = "body")
if (n_c) ft <- hline(ft, i = max(rows_c), border = dotted, part = "body")

faint <- fp_border(color = "grey85", width = 0.5)
for (blk in list(rows_c, rows_s))
  if (length(blk) > 1)
    ft <- hline(ft, i = blk[-length(blk)], border = faint, part = "body")

vrule <- fp_border(color = "grey60", width = 0.75, style = "dotted")
structured <- setdiff(seq_len(nrow(combined)), panel_rows)
ft <- vline(ft, i = structured, j = c(2, 6), border = vrule, part = "body")

ft <- fontsize(ft, size = TABLE_FONT_SIZE, part = "all")
ft <- padding(ft, padding.top = 1, padding.bottom = 1,
              padding.left = 2, padding.right = 2, part = "all")

# The indent goes after the padding call above, not before it. That call
# uses part = "all" and so sets padding.left on every cell; an indent set
# earlier is simply overwritten by it.
ft <- padding(ft, i = model_rows, j = 1, padding.left = INDENT,
              part = "body")

ft <- autofit(ft)
ft <- width(ft, j = "Model", width = 1.40)
ft <- width(ft, j = "Data sources", width = 1.56)

# Control and Active narrowed more than the two difference columns, and
# the discounted difference left at its previous width. The widest
# value seen in testing, a separate arms model at 30 years, needed about
# 0.93in; narrowing control/active below that risks wrapping a cell like
# it onto two lines occasionally. Worth checking the rendered table for
# any such cell once this runs, since it is a real per-cell risk rather
# than a rounding margin.
ctrl_active_cols <- col_names[c(3, 4, 7, 8)]
diff_cols <- col_names[c(5, 9)]
diff_disc_cols <- col_names[c(6, 10)]

ft <- width(ft, j = ctrl_active_cols, width = 0.86)
ft <- width(ft, j = diff_cols, width = 0.89)
ft <- width(ft, j = diff_disc_cols, width = 0.92)
ft <- set_table_properties(ft, layout = "fixed", align = "left")

caption_text <- paste(
  "Table 2. Restricted mean survival (years) at 10 and 30 years, by model and",
  "by the evidence used, for the ALEX February 2017 data cut, the company's",
  "original base case in NICE TA536. Flatiron real world data are weighted by",
  "matching adjusted indirect comparison throughout. Cells give the estimate",
  "with a 95 per cent interval.")

footnote_text <- paste(
  "Control is crizotinib, active is alectinib, difference is alectinib minus",
  "crizotinib. Arm level figures are undiscounted; the difference is shown both",
  "undiscounted and discounted at 3.5 per cent. The exponential and piecewise",
  "models are the two survival models used in NICE TA536; thirty years is the",
  "lifetime horizon used in that appraisal, matched here for direct comparison.",
  "survextrap models are the base case specification, six degrees of freedom,",
  "sigma prior Gamma(2,1) and, for the non-proportional hazards model, tau",
  "prior Gamma(2,10). Separate arms models take external data on the control",
  "arm only and give implausible extrapolations for that arm, as shown in",
  "Supplementary Figure 6. The validation row is the Kaplan-Meier estimate from",
  "the ALEX final overall survival analysis, data cut 28 April 2025,",
  "reconstructed from the published curves; it has no thirty year entry because",
  "follow-up reaches 10.6 years. Intervals for the exponential and piecewise",
  "models are bootstrap percentile intervals resampling patients, and so are",
  "wider than the parameter only intervals a technology appraisal would usually",
  "report.")

dir.create("Tables", showWarnings = FALSE)

doc <- read_docx()
doc <- body_set_default_section(doc, prop_section(
  page_size = page_size(width = 11.69, height = 8.27, orient = "landscape"),
  page_margins = page_mar(top = 0.6, bottom = 0.6, left = 0.5, right = 0.5)))
doc <- body_add_par(doc, caption_text, style = "Normal")
doc <- body_add_flextable(doc, ft)
doc <- body_add_fpar(doc, fpar(ftext(footnote_text,
                                     prop = fp_text(font.size = 8))))
print(doc, target = "Tables/Table2_rmst.docx")

cat("\nwritten to Tables/Table2_rmst.docx\n")
print(as.data.frame(combined), row.names = FALSE)