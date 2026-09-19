# Table 1: patient characteristics across the three data sources.
#
# Six columns under three source headings: the Flatiron real world
# cohort unweighted and MAIC weighted, the two ALEX arms, and the
# PROFILE-1014 crizotinib arm. Characteristics are grouped, with the
# group name on its own row and its categories indented beneath, so the
# table reads down as a list of characteristics rather than across as a
# grid of numbers.
#
# The MAIC weighted column is not a count of people. Its weights were
# derived to reproduce the ALEX crizotinib arm's characteristics, which
# is why its percentages match that column exactly for the matched
# variables, and its total of 27.1 is an effective sample size rather
# than a number of patients. Counts are shown to one decimal place for
# that reason.
#
# Two things in the source numbers that this script reports rather than
# silently accepts, printed when it runs.
#
# PROFILE-1014's ECOG categories sum to 171 against 172 patients, so one
# patient is unaccounted for, presumably missing status.
#
# PROFILE-1014's CNS metastasis "No" count was not supplied and is left
# as a dash. It is derivable as 172 minus 45, but filling it in here
# would be this script inventing a number rather than reporting one.
#
# The layout follows Table 2: all rows are body rows, so the source
# headings and column titles sit where they are wanted, with black rules
# only at the frame and under the titles and everything inside dotted or
# faint.

library(dplyr)
library(tibble)
library(flextable)
library(officer)

TABLE_FONT_SIZE <- 9


# ---- the numbers ---------------------------------------------------------
#
# kind is "count" for the patient total, "group" for a characteristic
# heading with no numbers of its own, and "level" for an indented
# category beneath one.

tbl <- tribble(
  ~kind,   ~label,                     ~fl_unw,     ~fl_maic,      ~alex_criz,  ~alex_alec,  ~profile,
  "count", "Patients, N",              "195",       "27.1",        "151",       "152",       "172",
  
  "group", "Age, years",               "",          "",            "",          "",          "",
  "level", "Mean",                     "60",        "54",          "54",        "56",        "-",
  "level", "Median",                   "61",        "53",          "54",        "58",        "52",
  "level", "Range",                    "24-84",     "24-84",       "18-91",     "25-88",     "22-76",
  
  "group", "Sex",                      "",          "",            "",          "",          "",
  "level", "Male",                     "101 (52%)", "11.5 (42%)",  "64 (42%)",  "68 (45%)",  "68 (40%)",
  "level", "Female",                   "94 (48%)",  "15.6 (58%)",  "87 (58%)",  "84 (55%)",  "104 (60%)",
  
  "group", "Race",                     "",          "",            "",          "",          "",
  "level", "Asian",                    "9 (5%)",    "12.4 (46%)",  "69 (46%)",  "69 (45%)",  "77 (45%)",
  "level", "Non-Asian",                "186 (95%)", "14.7 (54%)",  "82 (54%)",  "83 (55%)",  "95 (55%)",
  
  "group", "Smoking status",           "",          "",            "",          "",          "",
  "level", "Ever smoker",              "92 (47%)",  "9.5 (35%)",   "53 (35%)",  "60 (39%)",  "66 (38%)",
  "level", "Never smoker",             "103 (53%)", "17.6 (65%)",  "98 (65%)",  "92 (61%)",  "106 (62%)",
  
  "group", "ECOG performance status",  "",          "",            "",          "",          "",
  "level", "0 or 1",                   "169 (87%)", "25.3 (93%)",  "141 (93%)", "142 (93%)", "161 (94%)",
  "level", "2",                        "26 (13%)",  "1.8 (7%)",    "10 (7%)",   "10 (7%)",   "10 (6%)",
  
  "group", "Extent of disease",        "",          "",            "",          "",          "",
  "level", "Locally advanced",         "1 (1%)",    "0.1 (0%)",    "6 (4%)",    "4 (3%)",    "4 (2%)",
  "level", "Metastatic",               "194 (99%)", "27.1 (100%)", "145 (96%)", "148 (97%)", "168 (98%)",
  
  "group", "CNS/brain metastasis",     "",          "",            "",          "",          "",
  "level", "Yes",                      "30 (15%)",  "10.4 (38%)",  "58 (38%)",  "64 (42%)",  "45 (26%)",
  "level", "No",                       "165 (85%)", "16.7 (62%)",  "93 (62%)",  "88 (58%)",  "-")


# ---- do the categories add up --------------------------------------------
#
# Every pair of categories under a heading should sum to that column's
# patient total. Anything that does not is printed, rather than left for
# a reader to notice.

totals <- c(fl_unw = 195, fl_maic = 27.1, alex_criz = 151,
            alex_alec = 152, profile = 172)

lead_number <- function(x) suppressWarnings(as.numeric(sub(" .*$", "", x)))

cat("categories that do not sum to the column total\n")
issues <- 0
groups <- which(tbl$kind == "group")
for (g in groups) {
  end <- c(groups, nrow(tbl) + 1)
  end <- min(end[end > g]) - 1
  block <- tbl[(g + 1):end, ]
  if (tbl$label[g] == "Age, years") next          # not categories
  for (col in names(totals)) {
    vals <- lead_number(block[[col]])
    if (all(is.na(vals))) next
    s <- sum(vals, na.rm = TRUE)
    if (abs(s - totals[[col]]) > 0.15) {
      cat(sprintf("  %-24s %-10s sums to %6.1f, total is %6.1f\n",
                  tbl$label[g], col, s, totals[[col]]))
      issues <- issues + 1
    }
  }
}
if (issues == 0) cat("  none\n")
cat("\nblank cells are shown as a dash and are not counted above\n\n")


# ---- build the table ------------------------------------------------------

col_names <- c("Characteristic",
               "Crizotinib,\nunweighted",
               "Crizotinib,\nMAIC-weighted*",
               "Crizotinib",
               "Alectinib",
               "Crizotinib ")     # trailing space: names must be unique

body <- tbl %>% select(-kind) %>% as.data.frame()
names(body) <- col_names

blank_row <- function() {
  r <- as.data.frame(as.list(rep("", length(col_names))),
                     stringsAsFactors = FALSE, check.names = FALSE)
  names(r) <- col_names
  r
}

# the three source headings, each spanning its own columns
source_row <- blank_row()
source_row[[2]] <- "Flatiron RWE"
source_row[[4]] <- "ALEX"
source_row[[6]] <- "PROFILE-1014"

titles <- as.data.frame(as.list(col_names), stringsAsFactors = FALSE,
                        check.names = FALSE)
names(titles) <- col_names

combined <- rbind(source_row, titles, body)

# row positions, after the two header rows
offset <- 2
row_source <- 1
row_titles <- 2
rows_group <- which(tbl$kind == "group") + offset
rows_level <- which(tbl$kind == "level") + offset
rows_count <- which(tbl$kind == "count") + offset

ft <- flextable(combined)
ft <- delete_part(ft, part = "header")

ft <- merge_at(ft, i = row_source, j = 2:3, part = "body")
ft <- merge_at(ft, i = row_source, j = 4:5, part = "body")

ft <- bold(ft, i = c(row_source, row_titles), part = "body")
ft <- bold(ft, i = c(rows_count, rows_group), j = 1, part = "body")

ft <- align(ft, part = "body", align = "center")
ft <- align(ft, j = 1, align = "left", part = "body")

rule <- fp_border(color = "black", width = 1)
ft <- border_outer(ft, border = rule, part = "body")
ft <- hline(ft, i = row_source, border = rule, part = "body")
ft <- hline(ft, i = row_titles, border = rule, part = "body")

# a dotted rule under the patient total and under each block of
# categories, so the characteristics read as separate sections
dotted <- fp_border(color = "grey45", width = 0.9, style = "dotted")
ft <- hline(ft, i = rows_count, border = dotted, part = "body")
for (g in rows_group[-1]) ft <- hline(ft, i = g - 1, border = dotted,
                                      part = "body")

vrule <- fp_border(color = "grey60", width = 0.75, style = "dotted")
ft <- vline(ft, j = c(1, 3, 5), border = vrule, part = "body")

ft <- fontsize(ft, size = TABLE_FONT_SIZE, part = "all")
ft <- padding(ft, padding.top = 1.5, padding.bottom = 1.5,
              padding.left = 3, padding.right = 3, part = "all")

# The indent, applied after the padding call above and not before it.
# That call uses part = "all", so it sets padding.left on every cell in
# the table; an indent set earlier is simply overwritten by it, which is
# what happened here at first and why the categories sat flush with
# their headings.
ft <- padding(ft, i = rows_level, j = 1, padding.left = 16, part = "body")

ft <- autofit(ft)
ft <- width(ft, j = 1, width = 1.85)
ft <- width(ft, j = 2:6, width = 1.05)
ft <- set_table_properties(ft, layout = "fixed", align = "left")

caption_text <- paste(
  "Table 1. Patient characteristics in the Flatiron real world cohort, the",
  "ALEX trial and the PROFILE-1014 trial. Counts are given with the",
  "percentage of that column's patients.")

footnote_text <- paste(
  "* The MAIC-weighted column is a weighted version of the same 195 Flatiron",
  "patients, not a separate group. Its weights were derived to reproduce the",
  "characteristics of the ALEX crizotinib arm, which is why its percentages",
  "match that column for the matched variables, and its total of 27.1 is an",
  "effective sample size rather than a number of people; counts are therefore",
  "shown to one decimal place.",
  "** Mean age was not reported for PROFILE-1014.",
  "PROFILE-1014's ECOG categories sum to 171 of 172 patients, so one patient's",
  "status appears not to have been reported. Its count of patients without CNS",
  "or brain metastasis was not available and is shown as a dash.")

dir.create("Tables", showWarnings = FALSE)

doc <- read_docx()
doc <- body_set_default_section(doc, prop_section(
  page_size = page_size(width = 8.27, height = 11.69, orient = "portrait"),
  page_margins = page_mar(top = 0.8, bottom = 0.8, left = 0.6, right = 0.6)))
doc <- body_add_par(doc, caption_text, style = "Normal")
doc <- body_add_flextable(doc, ft)
doc <- body_add_fpar(doc, fpar(ftext(footnote_text,
                                     prop = fp_text(font.size = 8))))
print(doc, target = "Tables/Table1_characteristics.docx")

cat("written to Tables/Table1_characteristics.docx\n")
print(as.data.frame(combined), row.names = FALSE)