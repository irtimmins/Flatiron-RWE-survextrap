# Background mortality for the additive hazards model, from England and
# Wales life tables.
#
# survextrap can treat the overall hazard as the sum of a known
# background hazard and an excess hazard for the disease, h(t) = hb(t) +
# hc(t), with only the excess estimated. The constraint that these
# patients cannot do better than the matched general population is
# informative exactly where the extrapolation is least supported.
#
# Rates come from the Human Mortality Database, as in Sweeting et al,
# who use HMDHFDplus for the same purpose. England and Wales rather than
# the trial's own international population, because the decision context
# is a NICE appraisal; that mismatch is worth stating in the methods
# rather than leaving implicit.
#
# Sweeting et al match rates to each patient's attained age, sex and
# calendar year. That is not possible here: the cohort is reconstructed
# pseudo-IPD from digitised curves and carries no ages. Rather than
# invent them, the background hazard is built at the cohort level, by
# weighting national rates to the age and sex mix reported in Table 1.
#
# How that averaging is done matters. Taking the mean of the individual
# hazards at each time would be wrong twice over. Mortality is convex in
# age, so the mean hazard across an age distribution is not the hazard
# of the cohort. And a cohort's composition changes as it ages: the
# older members die first, so the survivors are younger than simply
# ageing the baseline distribution suggests, and their mean hazard is
# lower than that of the original mix.
#
# Both are avoided by averaging survival rather than hazard. Each age
# and sex stratum gets its own background survival curve, those are
# averaged with the stratum's weight, and the marginal hazard is the
# derivative of the result:
#
#   Sbar(t) = sum_i w_i exp(-integral of rate_i)
#   hbar(t) = -d log Sbar(t) / dt
#
# which is what the additive model needs, and is survivor weighted by
# construction.
#
# The calendar year is held fixed at randomisation rather than advancing
# with follow-up, so no future improvement in mortality is assumed. That
# is conservative for survival: real rates at a given age will likely be
# lower in twenty years than they were in 2014.
#
# Inputs
#   Human Mortality Database credentials, in the environment as
#   HMD_USER and HMD_PASS, or set below. Registration is free at
#   mortality.org. The download happens once and is cached.
# Outputs
#   Data/life_table_ew.csv   age, sex, rate for the reference year
#   Data/backhaz.rds         time and hazard, for survextrap's backhaz

library(dplyr)
library(tidyr)
library(ggplot2)

cfg <- list(
  country = "GBRTENW",      # England and Wales, total population
  cache = "Data/life_table_ew.csv",
  out = "Data/backhaz.rds",
  
  # The year the cohort was randomised, so the rates are the ones those
  # patients actually faced. ALEX randomised from August 2014 to
  # February 2015, so 2014 covers most of it; worth confirming against
  # the trial publication and switching to 2015 if that fits better.
  # This also sidesteps the pandemic years, which would otherwise need
  # excluding: their excess mortality is not something to project across
  # thirty years.
  ref_year = 2014,
  
  # The cohort as reported in Table 1, pooled across the two ALEX arms:
  # median age 54 and 58 by arm, range 18 to 91, 42 to 45 per cent male.
  # A normal truncated to the observed range is a reasonable summary
  # given only a median and a range. The standard deviation is a
  # judgement, set so the range spans about four of them, and is the one
  # quantity here with no direct support in Table 1.
  age_median = 56,
  age_sd = 11,
  age_range = c(18, 91),
  prop_male = 0.435,
  
  horizon = 30,
  step = 0.5,
  n_cohort = 20000)

dir.create("Data", showWarnings = FALSE)


######################################################
#  The life table, downloaded once then cached.
######################################################

if (file.exists(cfg$cache)) {
  
  life <- read.csv(cfg$cache, stringsAsFactors = FALSE)
  cat("life table read from", cfg$cache, "\n")
  
} else {
  
  if (!requireNamespace("HMDHFDplus", quietly = TRUE))
    stop("HMDHFDplus is needed to download the life tables: ",
         "install.packages(\"HMDHFDplus\")")
  
  user <- Sys.getenv("HMD_USER")
  pass <- Sys.getenv("HMD_PASS")
  
  if (user == "" || pass == "")
    stop("set HMD_USER and HMD_PASS in the environment, or edit them in ",
         "here.\n  Registration is free at mortality.org. For example, in ",
         "~/.Renviron:\n    HMD_USER=you@example.com\n    HMD_PASS=yourpassword")
  
  cat("downloading", cfg$country, "life tables from the HMD\n")
  
  male <- HMDHFDplus::readHMDweb(CNTRY = cfg$country, item = "mltper_1x1",
                                 username = user, password = pass)
  female <- HMDHFDplus::readHMDweb(CNTRY = cfg$country, item = "fltper_1x1",
                                   username = user, password = pass)
  
  years <- sort(unique(male$Year))
  cat("years available:", min(years), "to", max(years), "\n")
  
  if (!cfg$ref_year %in% years)
    stop("the reference year ", cfg$ref_year, " is not in the downloaded ",
         "table, which covers ", min(years), " to ", max(years))
  
  # mx is the central death rate, which is the hazard on an annual scale.
  # qx, the probability of dying within the year, would need converting
  # as -log(1 - qx); the two agree closely except at very old ages.
  life <- bind_rows(
    male %>% filter(Year == cfg$ref_year) %>%
      transmute(age = Age, sex = "male", rate = mx),
    female %>% filter(Year == cfg$ref_year) %>%
      transmute(age = Age, sex = "female", rate = mx))
  
  write.csv(life, cfg$cache, row.names = FALSE)
  cat("saved to", cfg$cache, "\n")
}

life <- life %>% filter(is.finite(rate), rate > 0)

cat("ages", min(life$age), "to", max(life$age), "for",
    paste(sort(unique(life$sex)), collapse = " and "), "\n")

oldest_needed <- max(cfg$age_range) + cfg$horizon
if (max(life$age) < oldest_needed)
  cat("note: the table stops at age", max(life$age), "but the oldest",
      "patient reaches", oldest_needed, "by the end of the horizon;",
      "rates are held at the oldest tabulated age beyond that\n")


######################################################
#  The cohort, as a weighting scheme.
######################################################
#
# No individual here corresponds to a real patient. This is a set of
# weights for the national rates, not a dataset.

set.seed(2024)

ages <- rnorm(cfg$n_cohort, cfg$age_median, cfg$age_sd)
ages <- round(ages[ages >= cfg$age_range[1] & ages <= cfg$age_range[2]])

cohort <- tibble(
  age = ages,
  sex = ifelse(runif(length(ages)) < cfg$prop_male, "male", "female"))

cat("\nweighting cohort:", nrow(cohort), "people, median age",
    median(cohort$age), ",", round(100 * mean(cohort$sex == "male")),
    "per cent male\n")


######################################################
#  The background hazard over time.
######################################################
#
# One background survival curve per age and sex stratum, averaged with
# the stratum's weight, then differentiated. Working in strata rather
# than over individuals is only for speed: ages are whole numbers, so a
# cohort of twenty thousand collapses to a couple of hundred distinct
# combinations carrying counts.

rate_at <- function(age, sex) {
  age <- pmin(pmax(age, min(life$age)), max(life$age))
  life$rate[match(paste(age, sex), paste(life$age, life$sex))]
}

strata <- cohort %>% count(age, sex, name = "w") %>% mutate(w = w / sum(w))

cat("\n", nrow(strata), " age and sex strata\n", sep = "")

# a fine grid for the integration, coarsened afterwards to the intervals
# survextrap is given
fine_step <- 0.05
fine_t <- seq(0, cfg$horizon, by = fine_step)

# rate for each stratum at each time, as a matrix of strata by times
rate_mat <- vapply(fine_t, function(tt)
  rate_at(floor(strata$age + tt), strata$sex), numeric(nrow(strata)))

if (any(is.na(rate_mat)))
  stop("some age and sex combinations were not found; check the table ",
       "covers every single year of age for both sexes")

# cumulative hazard for each stratum, then its survival
cum_haz <- t(apply(rate_mat * fine_step, 1, cumsum)) - rate_mat * fine_step
surv_mat <- exp(-cum_haz)

# the cohort's marginal background survival, and the hazard that implies
surv_bar <- as.vector(strata$w %*% surv_mat)

haz_bar <- -diff(log(surv_bar)) / fine_step
haz_t <- fine_t[-length(fine_t)] + fine_step / 2

# onto the intervals survextrap is given, each taking the mean of the
# fine grid points falling inside it
times <- seq(0, cfg$horizon - cfg$step, by = cfg$step)

backhaz <- tibble(
  time = times,
  hazard = vapply(times, function(tt) {
    inside <- haz_t >= tt & haz_t < tt + cfg$step
    if (!any(inside)) NA_real_ else mean(haz_bar[inside])
  }, numeric(1)))

if (any(is.na(backhaz$hazard)))
  stop("some intervals had no fine grid points in them; reduce step or ",
       "fine_step")

# how much the correct aggregation differs from simply averaging the
# hazards, which is the thing it is avoiding
naive <- vapply(times, function(tt)
  sum(strata$w * rate_at(floor(strata$age + tt), strata$sex)), numeric(1))

cat("\nmarginal hazard against the mean of individual hazards\n")
print(as.data.frame(tibble(
  time = times, marginal = backhaz$hazard, mean_of_hazards = naive) %>%
    filter(time %in% c(0, 5, 10, 15, 20, 25)) %>%
    transmute(time,
              marginal = round(marginal, 4),
              mean_of_hazards = round(mean_of_hazards, 4),
              difference = sprintf("%+.1f%%",
                                   100 * (marginal / mean_of_hazards - 1)))),
  row.names = FALSE)

cat("\n  The two agree at time zero, where no one has died yet, and part\n")
cat("  company later as the higher mortality strata drop out. The\n")
cat("  marginal one is what the additive model needs.\n")

cat("\nbackground hazard over the horizon\n")
print(as.data.frame(backhaz %>%
                      filter(time %in% c(0, 5, 10, 15, 20, 25)) %>%
                      transmute(time,
                                median_age = round(median(cohort$age) + time),
                                hazard = round(hazard, 4),
                                annual_risk = sprintf("%.1f%%", 100 * (1 - exp(-hazard))))),
      row.names = FALSE)

saveRDS(backhaz, cfg$out)
cat("\nwritten to", cfg$out, "\n")


######################################################
#  Is it plausible next to the disease hazard.
######################################################
#
# An additive model cannot have a negative excess hazard, so background
# mortality must sit well below the observed hazard early on. If it does
# not, the cohort's ages or the life table are wrong.

if (file.exists("Data/trial_data_apr_2025.rds")) {
  alex <- readRDS("Data/trial_data_apr_2025.rds")
  crude <- sum(alex$status) / sum(alex$time)
  share <- 100 * backhaz$hazard[1] / crude
  cat("\ncrude overall hazard in the final analysis:", round(crude, 4),
      "per year\n")
  cat("background at time 0:", round(backhaz$hazard[1], 4),
      sprintf(", which is %.1f per cent of it\n", share))
  cat("background at 20 years:",
      round(backhaz$hazard[backhaz$time == 20], 4), "\n")
  if (share > 25)
    cat("\nThat share looks high. Background mortality should be a small\n",
        "part of the total in a cohort with advanced cancer; check the\n",
        "cohort ages and the reference year before using this.\n", sep = "")
}

fig <- ggplot(backhaz, aes(x = time, y = hazard)) +
  geom_step(linewidth = 0.8, colour = "#2C6FB3") +
  scale_y_continuous("Background mortality hazard (per year)") +
  xlab("Time (years)") +
  theme_classic()

print(fig)

dir.create("Figures", showWarnings = FALSE)
jpeg(file = "Figures/background_mortality.jpg",
     width = 5.5, height = 3.5, units = "in", res = 300, quality = 95)
print(fig)
dev.off()