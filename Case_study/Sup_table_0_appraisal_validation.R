# Supplementary Table 1: results reported for the company's overall
# survival model in NICE TA536, alongside our reconstruction.
#
# Three blocks, one per data source and cut:
#   ALEX, Feb 2017     life years, goodness of fit, exponential model
#                      predictions, observed Kaplan-Meier
#   ALEX, Dec 2017     life years, goodness of fit
#   PROFILE-1014       observed Kaplan-Meier
#
# Published values are typed in below next to the document, table and page
# they came from. Our values are never typed in: life years are read from
# the output of Fit_exp_and_hybrid.R and everything else is computed here
# from the reconstructed patient level data.
#
# Every source is written out in full, on the heading row of the panel it
# applies to rather than in a column, which is what lets the table fit on a
# portrait page. The two committee papers packs are listed with their web
# addresses under the table, so a reader can go from any panel to the
# document without a key.
#
# Run after Fit_exp_and_hybrid.R.

library(dplyr)
library(tidyr)
library(tibble)
library(survival)
library(flexsurv)
library(flextable)
library(officer)

TABLE_FONT_SIZE <- 8

cfg <- list(
  cuts = list(
    feb_2017 = list(label = "February 2017", short = "Feb 2017",
                    file = "Data/trial_data.rds", cut_months = 18),
    dec_2017 = list(label = "December 2017", short = "Dec 2017",
                    file = "Data/trial_data_dec_2017.rds", cut_months = 25)),
  life_years = "Data/exp_hybrid_life_years.rds",
  profile    = "Data/trial_IPD_OS_PROFILE_1014_Nov_2016.csv",
  horizon    = 30,
  disc_rate  = 0.035,
  # Treat a Gompertz fit with a negative shape as an exponential, which is
  # what the company's reported values imply they did.
  gompertz_as_company = TRUE,
  out = "Tables/Sup_table_1_TA536.docx")

pack1 <- "Committee papers 1,"
pack2 <- "Committee papers 2,"

dash <- "\u2013"
f2  <- function(x) ifelse(is.na(x), dash, sprintf("%.2f", x))
fc  <- function(x) ifelse(is.na(x), "-", sprintf("%.2f", x))
pct <- function(x, dp = 0)
  ifelse(is.na(x), dash, sprintf(paste0("%.", dp, "f%%"), 100 * x))
surv_risk <- function(s, n)
  ifelse(is.na(s), dash, sprintf("%.0f%% (%d)", 100 * s, as.integer(n)))


# ---- reconstructed data -------------------------------------------------

# Reads either a .rds or a .csv, and supplies the arm for single arm files
# that do not carry one.
read_ipd <- function(file, arm = NULL) {
  if (!file.exists(file)) stop(file, " not found")
  d <- if (grepl("\\.csv$", file, ignore.case = TRUE))
    as_tibble(read.csv(file, stringsAsFactors = FALSE)) else
      as_tibble(readRDS(file))
  if (!all(c("time", "status") %in% names(d)))
    stop(file, " needs time and status columns, has ",
         paste(names(d), collapse = ", "))
  if (!"trt" %in% names(d)) {
    if (is.null(arm)) stop(file, " has no trt column and no arm was supplied")
    d$trt <- arm
  }
  d <- d %>% mutate(trt = as.character(trt)) %>%
    select(all_of(c("time", "status", "trt")))
  if (!is.null(arm)) d <- filter(d, trt == arm)
  if (nrow(d) == 0) stop(file, ": no rows for arm ", arm)
  if (max(d$time) > 15)
    stop(file, ": longest time is ", round(max(d$time), 1),
         ", which looks like months or days rather than years")
  d
}

trial <- lapply(cfg$cuts, function(cc) read_ipd(cc$file))
cut_label <- sapply(cfg$cuts, `[[`, "label")
cut_short <- sapply(cfg$cuts, `[[`, "short")

arm_rate <- function(k, arm) {
  d <- filter(trial[[k]], trt == arm)
  sum(d$status) / sum(d$time)
}


# ---- life years from Fit_exp_and_hybrid.R -------------------------------

ly <- readRDS(cfg$life_years)
need <- c("data_cut", "model", "quantity", "value")
if (!all(need %in% names(ly)))
  stop(cfg$life_years, " needs columns ", paste(need, collapse = ", "))
if ("disc_rate" %in% names(ly))
  ly <- filter(ly, abs(disc_rate - cfg$disc_rate) < 1e-9)

our_ly <- function(cut_lab, model_name) {
  v <- filter(ly, data_cut == cut_lab, model == model_name)
  sapply(c("Alectinib", "Crizotinib", "increment"), function(q) {
    x <- v$value[v$quantity == q]
    if (length(x) == 1) x else NA_real_
  })
}

# One row per model, taking the figures from the ERG documents wherever the
# ERG reported them, which for February 2017 is both models.
published_ly <- tribble(
  ~item, ~cut, ~model, ~alec, ~criz, ~inc,
  "Exponential", "feb_2017", "exponential", 5.17, 4.25, 0.93,
  "Kaplan-Meier to 18 months, exponential tail", "feb_2017", "hybrid",
  5.14, 4.32, 0.83,
  "Kaplan-Meier to 25 months, exponential tail", "dec_2017", "hybrid",
  5.46, 4.42, 1.05)

ly_source <- list(
  feb_2017 = paste("exponential from", pack1, "ERG report, Table 40, p.151;",
                   "piecewise from", pack2,
                   "ERG review of the company's ACD response, Table 6, p.11"),
  dec_2017 = paste(pack2, "Company ACD response, Appendix 3, Table 8, p.11"))

published_ly <- published_ly %>%
  rowwise() %>%
  mutate(ours = list(our_ly(cut_label[[cut]], model))) %>%
  ungroup()

cat("life years, published against ours\n")
for (i in seq_len(nrow(published_ly))) {
  o <- published_ly$ours[[i]]
  cat(sprintf("  %-45s %s  published %.2f %.2f %.2f  ours %s %s %s\n",
              published_ly$item[i], cut_short[[published_ly$cut[i]]],
              published_ly$alec[i], published_ly$criz[i], published_ly$inc[i],
              fc(o[1]), fc(o[2]), fc(o[3])))
}


# ---- goodness of fit ------------------------------------------------------
#
# The company's rank is the order of AIC plus BIC within an arm. That rule
# is checked against the twelve ranks they published for December 2017
# before it is used for anything. Gamma is the generalised gamma: the
# company's BIC minus AIC is 3(log n - 2), so three parameters.

dists <- c(exp = "Exponential", weibull = "Weibull", lnorm = "Log-normal",
           gengamma = "Gamma", llogis = "Log-logistic",
           gompertz = "Gompertz")

published_gof <- tribble(
  ~cut, ~dist, ~aic_a, ~bic_a, ~aic_c, ~bic_c, ~rank_a, ~rank_c,
  "feb_2017", "exp",      246.59, 249.61, 234.24, 237.26, NA, NA,
  "feb_2017", "weibull",  247.98, 254.03, 232.71, 238.74, NA, NA,
  "feb_2017", "lnorm",    247.97, 254.02, 230.88, 236.91, NA, NA,
  "feb_2017", "gengamma", 249.79, 258.86, 232.79, 241.84, NA, NA,
  "feb_2017", "llogis",   247.91, 253.96, 232.10, 238.13, NA, NA,
  "feb_2017", "gompertz", 248.59, 254.63, 234.72, 240.76, NA, NA,
  "dec_2017", "exp",      289.07, 292.10, 271.32, 274.34, 1, 2,
  "dec_2017", "weibull",  289.36, 295.41, 272.18, 278.21, 4, 5,
  "dec_2017", "lnorm",    288.56, 294.61, 268.16, 274.19, 2, 1,
  "dec_2017", "gengamma", 290.50, 299.57, 269.40, 278.45, 6, 4,
  "dec_2017", "llogis",   288.85, 294.89, 270.46, 276.49, 3, 3,
  "dec_2017", "gompertz", 291.07, 297.12, 273.32, 279.35, 5, 6)

pub_gof <- bind_rows(
  published_gof %>% transmute(cut, dist, trt = "Alectinib",
                              aic = aic_a, bic = bic_a, rank_pub = rank_a),
  published_gof %>% transmute(cut, dist, trt = "Crizotinib",
                              aic = aic_c, bic = bic_c, rank_pub = rank_c)) %>%
  group_by(cut, trt) %>%
  mutate(d_aic = aic - min(aic),
         rank_rule = rank(aic + bic, ties.method = "min")) %>%
  ungroup()

check <- filter(pub_gof, !is.na(rank_pub))
if (!all(check$rank_pub == check$rank_rule)) {
  print(as.data.frame(check))
  stop("AIC plus BIC order does not reproduce the company's December ranks")
}
cat("\nAIC plus BIC order reproduces all", nrow(check),
    "ranks the company published for December 2017\n")

fit_arm <- function(d) {
  n <- nrow(d)
  out <- bind_rows(lapply(names(dists), function(dd) {
    f <- tryCatch(flexsurvreg(Surv(time, status) ~ 1, data = d, dist = dd),
                  error = function(e) NULL)
    if (is.null(f))
      return(tibble(dist = dd, loglik = NA_real_, npar = NA_real_,
                    shape = NA_real_))
    tibble(dist = dd, loglik = f$loglik, npar = f$npars,
           shape = if (dd == "gompertz") unname(f$res["shape", "est"])
           else NA_real_)
  }))
  
  g <- out$dist == "gompertz"
  if (any(g) && !is.na(out$shape[g])) {
    cat(sprintf("    gompertz shape %.3f, log likelihood %.2f",
                out$shape[g], out$loglik[g]))
    if (cfg$gompertz_as_company && out$shape[g] < 0) {
      out$loglik[g] <- out$loglik[out$dist == "exp"]
      cat(", negative so set to the exponential's")
    }
    cat("\n")
  }
  
  out %>%
    mutate(aic = -2 * loglik + 2 * npar,
           bic = -2 * loglik + npar * log(n),
           d_aic = aic - min(aic, na.rm = TRUE),
           rank_rule = rank(aic + bic, ties.method = "min", na.last = "keep"))
}

cat("\nfitting parametric models to the reconstructed data\n")
our_gof <- bind_rows(lapply(names(cfg$cuts), function(k) {
  bind_rows(lapply(c("Alectinib", "Crizotinib"), function(a) {
    cat("  ", cut_short[[k]], a, "\n")
    fit_arm(filter(trial[[k]], trt == a)) %>% mutate(cut = k, trt = a)
  }))
}))

gof_cell <- function(tab, k, arm, dd) {
  r <- tab[tab$cut == k & tab$trt == arm & tab$dist == dd, ]
  if (nrow(r) != 1 || is.na(r$d_aic)) return(dash)
  sprintf("%.2f (%d)", r$d_aic, as.integer(r$rank_rule))
}


# ---- landmark survival ----------------------------------------------------

landmarks <- c(12, 24, 30, 36, 48)

# Company submission Document B, Table 18 for alectinib and Table 17 for
# crizotinib, both read from their exponential base case.
published_exp <- tribble(
  ~months, ~alec, ~criz,
  12, 0.85, 0.82,
  24, 0.73, 0.67,
  30, 0.67, 0.60,
  36, 0.62, 0.54,
  48, NA,   0.44)

# The company's own reading of the ALEX curve, alectinib arm only, so this
# tests our digitisation of the figure they were reading.
published_alex_km <- tribble(
  ~months, ~surv, ~n_risk,
  12, 0.84, 119,
  24, 0.73, 38)

published_profile <- tribble(
  ~months, ~surv, ~n_risk,
  12, 0.83, 138,
  24, 0.65, 101,
  30, 0.62, 89,
  36, 0.58, 77,
  48, 0.55, 40)

km_at <- function(d, months, by_arm = TRUE) {
  fit <- if (by_arm) survfit(Surv(time, status) ~ trt, data = d)
  else survfit(Surv(time, status) ~ 1, data = d)
  s <- summary(fit, times = months / 12, extend = TRUE)
  tibble(trt = if (by_arm) sub("^trt=", "", as.character(s$strata))
         else unique(d$trt),
         months = round(12 * s$time), surv = s$surv, n_risk = s$n.risk)
}

alex_km <- km_at(trial$feb_2017, published_alex_km$months)
our_alex <- function(arm, m) {
  r <- filter(alex_km, trt == arm, months == m)
  if (nrow(r) != 1) return(dash)
  surv_risk(r$surv, r$n_risk)
}

profile_n <- NA
our_profile <- tibble(months = landmarks, surv = NA_real_, n_risk = NA_real_)
if (file.exists(cfg$profile)) {
  prof <- read_ipd(cfg$profile, arm = "Crizotinib")
  profile_n <- nrow(prof)
  our_profile <- km_at(prof, landmarks, by_arm = FALSE)
  if (profile_n != 172)
    cat("\nPROFILE-1014 reconstruction has", profile_n, "patients at time",
        "zero, against the 172 randomised to crizotinib that the company's",
        "at-risk numbers imply, so the two are reading different curves\n")
} else {
  cat("\n", cfg$profile, "not found; the PROFILE-1014 block will show",
      "published values only\n")
}


# ---- assemble the rows ----------------------------------------------------
#
# pa, pc, pd are published alectinib, crizotinib and contrast; oa, oc, od
# are ours. Each panel opens with a heading row that spans the table and
# carries both the panel title and where its published values came from.
# Keeping the sources there rather than in a column is what lets the table
# fit the width of a portrait page.

# panel_id is not shown; it tells flextable which rows belong together so a
# panel is not split across a page break.
panel_id <- 0

row_b <- function(text, source)
  tibble(type = "banner", panel_id = panel_id, item = text, source = source,
         pa = "", pc = "", pd = "", oa = "", oc = "", od = "")
row_v <- function(item, pub, ours)
  tibble(type = "value", panel_id = panel_id, item = item, source = "",
         pa = pub[1], pc = pub[2], pd = pub[3],
         oa = ours[1], oc = ours[2], od = ours[3])

panel <- function(text, cut_text, source) {
  panel_id <<- panel_id + 1
  row_b(sprintf("(%s) %s (%s)", letters[panel_id], text, cut_text), source)
}

rows <- list()
add <- function(x) rows[[length(rows) + 1]] <<- x

for (k in names(cfg$cuts)) {
  cut_text <- paste0("ALEX trial, ", cut_short[[k]], " data cut")
  
  add(panel(sprintf("Life years to %d years, discounted at %.1f per cent",
                    cfg$horizon, 100 * cfg$disc_rate),
            cut_text, ly_source[[k]]))
  for (i in which(published_ly$cut == k)) {
    r <- published_ly[i, ]
    add(row_v(r$item, f2(c(r$alec, r$criz, r$inc)), f2(r$ours[[1]])))
  }
  
  add(panel("Overall survival distributions, difference in AIC (rank)",
            cut_text,
            if (k == "feb_2017")
              paste(pack1, "Company submission, Document B, Table 16, p.58")
            else
              paste(pack2, "Company ACD response, Appendix 1, Table 4,",
                    "pp.4-5")))
  for (dd in names(dists))
    add(row_v(dists[[dd]],
              c(gof_cell(pub_gof, k, "Alectinib", dd),
                gof_cell(pub_gof, k, "Crizotinib", dd), ""),
              c(gof_cell(our_gof, k, "Alectinib", dd),
                gof_cell(our_gof, k, "Crizotinib", dd), "")))
  
  if (k == "feb_2017") {
    add(panel("Survival predicted by the exponential model", cut_text,
              paste(pack1, "Company submission, Document B, Table 18, p.65",
                    "for alectinib and Table 17, p.63 for crizotinib;",
                    "20 years from p.64")))
    for (i in seq_len(nrow(published_exp))) {
      m <- published_exp$months[i]
      add(row_v(paste(m, "months"),
                c(pct(published_exp$alec[i]), pct(published_exp$criz[i]), ""),
                c(pct(exp(-arm_rate(k, "Alectinib") * m / 12)),
                  pct(exp(-arm_rate(k, "Crizotinib") * m / 12)), "")))
    }
    add(row_v("20 years", c(pct(0.039, 1), dash, ""),
              c(pct(exp(-arm_rate(k, "Alectinib") * 20), 1),
                pct(exp(-arm_rate(k, "Crizotinib") * 20), 1), "")))
    
    add(panel("Observed Kaplan-Meier survival (number at risk)", cut_text,
              paste(pack1, "Company submission, Document B, Table 18, p.65,",
                    "which gives the alectinib arm only")))
    for (i in seq_len(nrow(published_alex_km))) {
      m <- published_alex_km$months[i]
      add(row_v(paste(m, "months"),
                c(surv_risk(published_alex_km$surv[i],
                            published_alex_km$n_risk[i]), dash, ""),
                c(our_alex("Alectinib", m), our_alex("Crizotinib", m), "")))
    }
  }
}

add(panel("Observed Kaplan-Meier survival (number at risk)",
          "PROFILE-1014 trial, Nov 2016 data cut",
          paste(pack1, "Company submission, Document B, Table 17, p.63;",
                "clarification response B9, p.58")))
for (i in seq_len(nrow(published_profile))) {
  m <- published_profile$months[i]
  o <- our_profile[our_profile$months == m, ]
  add(row_v(paste(m, "months"),
            c("", surv_risk(published_profile$surv[i],
                            published_profile$n_risk[i]), ""),
            c("", if (nrow(o) == 1) surv_risk(o$surv, o$n_risk) else dash, "")))
}

tab <- bind_rows(rows)


# ---- flextable ------------------------------------------------------------

keys <- c("item", "pa", "pc", "pd", "oa", "oc", "od")
ft <- flextable(tab, col_keys = keys)
ft <- set_header_labels(ft, item = "", pa = "Alectinib", pc = "Crizotinib",
                        pd = "Difference", oa = "Alectinib",
                        oc = "Crizotinib", od = "Difference")
ft <- add_header_row(ft, values = c("", "Reported in NICE TA536",
                                    "This study"),
                     colwidths = c(1, 3, 3))
ft <- merge_v(ft, j = "item", part = "header")

row_ban <- which(tab$type == "banner")
for (i in row_ban) ft <- merge_at(ft, i = i, j = 1:7, part = "body")

# The heading row holds the panel title in bold and its source behind it in
# smaller regular type, on the same line, so the source wraps into the full
# width of the table rather than into a column of its own.
for (i in row_ban) {
  ft <- compose(ft, i = i, j = "item", part = "body",
                value = as_paragraph(
                  as_chunk(tab$item[i],
                           props = fp_text(font.size = TABLE_FONT_SIZE,
                                           font.family = "Arial",
                                           bold = TRUE)),
                  as_chunk(paste0("   Source: ", tab$source[i]),
                           props = fp_text(font.size = TABLE_FONT_SIZE - 1.5,
                                           font.family = "Arial"))))
}

ft <- theme_booktabs(ft)
ft <- font(ft, fontname = "Arial", part = "all")
ft <- fontsize(ft, size = TABLE_FONT_SIZE, part = "all")
ft <- bold(ft, part = "header")
ft <- bg(ft, i = row_ban, bg = "grey92", part = "body")
ft <- align(ft, j = 2:7, align = "center", part = "all")
ft <- align(ft, i = row_ban, align = "left", part = "body")
ft <- valign(ft, valign = "top", part = "body")

rule  <- fp_border(color = "black", width = 1)
faint <- fp_border(color = "grey85", width = 0.5)
ft <- hline_top(ft, border = rule, part = "header")
ft <- hline_bottom(ft, border = rule, part = "header")
ft <- hline_bottom(ft, border = rule, part = "body")
ft <- hline(ft, i = setdiff(seq_len(nrow(tab)), c(row_ban - 1, nrow(tab))),
            border = faint, part = "body")
ft <- hline(ft, i = row_ban, border = rule, part = "body")
ft <- hline(ft, i = setdiff(row_ban - 1, 0), border = rule, part = "body")

vrule <- fp_border(color = "grey60", width = 0.75, style = "dotted")
ft <- vline(ft, i = setdiff(seq_len(nrow(tab)), row_ban), j = c(1, 4),
            border = vrule, part = "body")
ft <- vline(ft, j = c(1, 4), border = vrule, part = "header")

ft <- padding(ft, padding.top = 1, padding.bottom = 1,
              padding.left = 2, padding.right = 2, part = "all")
ft <- width(ft, j = "item", width = 1.75)
ft <- width(ft, j = c("pa", "pc", "pd", "oa", "oc", "od"), width = 0.87)
ft <- set_table_properties(ft, layout = "fixed", align = "left")

# Keep each panel whole on a page, and repeat the header on any page the
# table runs onto.
ft <- paginate(ft, group = "panel_id", hdr_ftr = TRUE)


# ---- caption, notes and the document list --------------------------------

caption_text <- paste(
  "Supplementary Table 1. Results reported for the company's overall survival",
  "model in NICE TA536, alongside our reconstruction of that model from",
  "digitised Kaplan-Meier curves. Difference is alectinib minus crizotinib.")

notes <- c(
  paste0(
    "Life years are the discounted area under the modelled survival curve to ",
    cfg$horizon, " years, in weekly cycles with a half-cycle correction, the",
    " convention used in the appraisal. In both piecewise models the",
    " exponential tail is fitted to the whole arm and starts from the",
    " Kaplan-Meier value at the cut point; the documents do not state the",
    " second of these, and it is the convention under which the life years",
    " match."),
  paste(
    "Difference in AIC is the difference from the best fitting distribution",
    "within an arm, and the rank is the order of AIC plus BIC within that arm,",
    "the rule that reproduces all twelve ranks the company published for",
    "December 2017. Absolute AIC values are not comparable between the two",
    "halves of the table, because the company fitted the trial's own patient",
    "data on a different time scale. Gamma is the generalised gamma. Several",
    "of the company's Gompertz fits sit exactly two AIC above the exponential,",
    "meaning the shape parameter was not estimated and the model reduced to an",
    "exponential; where our fitted shape is negative we have done the same."),
  paste0(
    "The company digitised the PROFILE-1014 curve with Engauge Digitizer 9.8,",
    " and their numbers at risk correspond to the 172 patients randomised to",
    " crizotinib",
    if (!is.na(profile_n)) paste0(", against ", profile_n, " in ours") else "",
    ", so the two are reading different curves and the numbers at risk in that",
    " block are not comparable."))

documents <- c(
  paste("Committee papers 1, covering the pre-meeting briefing, the company",
        "submission (Document B), the clarification questions and responses,",
        "and the Evidence Review Group report:",
        "https://www.nice.org.uk/guidance/ta536/documents/committee-papers"),
  paste("Committee papers 2, covering the company's response to the appraisal",
        "consultation document and its appendices, the Evidence Review Group's",
        "review of that response, and the addendum to it:",
        "https://www.nice.org.uk/guidance/ta536/evidence/committee-papers-pdf-4913116525"),
  paste("Final guidance, NICE technology appraisal guidance TA536:",
        "https://www.nice.org.uk/guidance/ta536"))


# ---- write ----------------------------------------------------------------

dir.create(dirname(cfg$out), showWarnings = FALSE, recursive = TRUE)

small <- fp_text(font.size = 7, font.family = "Arial")

doc <- read_docx()
doc <- body_set_default_section(doc, prop_section(
  page_size = page_size(width = 8.27, height = 11.69, orient = "portrait"),
  page_margins = page_mar(top = 0.6, bottom = 0.6, left = 0.5, right = 0.5)))
doc <- body_add_par(doc, caption_text, style = "Normal")
doc <- body_add_flextable(doc, ft)
for (nt in notes) doc <- body_add_fpar(doc, fpar(ftext(nt, prop = small)))
doc <- body_add_fpar(doc, fpar(ftext("Appraisal documents",
                                     prop = update(small, bold = TRUE))))
for (d in documents) doc <- body_add_fpar(doc, fpar(ftext(d, prop = small)))
print(doc, target = cfg$out)
cat("\nwritten to", cfg$out, "\n")