# Fitting the survextrap grid on the cluster.
#
# Eighteen primary scenarios: three treatment effect models by three
# evidence sets by two data cuts. Four sensitivities, all on the fullest
# evidence under non-proportional hazards, February cut, each changing
# one thing so its effect is readable.
#
# Priors come from Data/prior_sets.rds, written by Check_prior_packages.R
# after every set was checked against the landmarks under every spline
# used here. Nothing is retyped.

library(dplyr)
library(tidyr)
library(purrr)
library(survival)
library(survextrap)
library(rslurm)

source("Functions/Survextrap_model.R")


#  Settings

fit_method <- "mcmc"          # "opt" for a quick test run

data_dir <- file.path(getwd(), "Data")
store_models <- file.path(getwd(), "models")
dir.create(store_models, showWarnings = FALSE)

backhaz_file <- file.path(data_dir, "backhaz.rds")

# Knots. With no external data the boundary sits at the end of the trial,
# since nothing informs a knot beyond it; the December cut runs further,
# so its knots start later. The sensitivity sets add a knot at ten and
# one inside the data.
add_knots_feb_trial <- c(2)
add_knots_feb_ext   <- c(2, 5)
add_knots_dec_trial <- c(3)
add_knots_dec_ext   <- c(3, 5)
add_knots_feb_10y   <- c(2, 5, 10)
add_knots_feb_3     <- c(2, 3, 5)

prior_sets <- readRDS(file.path(data_dir, "prior_sets.rds"))
chosen <- prior_sets %>% filter(set == "chosen")
tighter <- prior_sets %>% filter(set == "tighter")

if (!nrow(chosen)) stop("prior_sets.rds has no chosen set")

cuts <- tibble(
  data_cut  = c("feb_2017", "dec_2017"),
  cut_label = c("February 2017", "December 2017"),
  suffix    = c("", "_dec_2017"))

external_path <- function(datasets, suffix) {
  case_when(
    datasets == "trial_only" ~ NA_character_,
    datasets == "trial_and_historic" ~
      paste0(data_dir, "/historic_trial_aggregate", suffix, ".rds"),
    datasets == "trial_and_all_maic" ~
      paste0(data_dir, "/external_data_maic_weighted", suffix, ".rds"),
    datasets == "trial_and_all_unweighted" ~
      paste0(data_dir, "/external_data_unweighted", suffix, ".rds"))
}


#  The eighteen primary scenarios

primary <- expand_grid(
  data_cut = cuts$data_cut,
  model = c("PH", "NON-PH", "Separate_arms"),
  datasets = c("trial_only", "trial_and_historic",
               "trial_and_all_maic")) %>%
  left_join(cuts, by = "data_cut") %>%
  mutate(
    scenario = "primary",
    df = 5,
    add_knots = case_when(
      data_cut == "feb_2017" & datasets == "trial_only" ~
        "add_knots_feb_trial",
      data_cut == "feb_2017" ~ "add_knots_feb_ext",
      data_cut == "dec_2017" & datasets == "trial_only" ~
        "add_knots_dec_trial",
      TRUE ~ "add_knots_dec_ext"),
    hsd_rate = chosen$hsd_rate, eta_med = chosen$eta_med,
    eta_up = chosen$eta_up, tau_rate = chosen$tau_rate)


#  Four sensitivities, one change each
#
# The reference is non-PH on the fullest evidence, February cut. Each row
# below differs from it in one respect, so the comparison is readable
# without a factorial expansion.

ref <- primary %>%
  filter(data_cut == "feb_2017", model == "NON-PH",
         datasets == "trial_and_all_maic")

sens <- bind_rows(
  # how much does the tail assumption drive the answer
  ref %>% mutate(scenario = "knots_10y", add_knots = "add_knots_feb_10y"),
  # does more flexibility inside the data change it
  ref %>% mutate(scenario = "knots_df8", df = 8,
                 add_knots = "add_knots_feb_3"),
  # the tighter of the two prior sets that held under every spline
  ref %>% mutate(scenario = "priors_tighter",
                 hsd_rate = tighter$hsd_rate, eta_med = tighter$eta_med,
                 eta_up = tighter$eta_up, tau_rate = tighter$tau_rate),
  # how much does the MAIC weighting matter
  ref %>% mutate(scenario = "unweighted",
                 datasets = "trial_and_all_unweighted"))

if (!nrow(tighter))
  sens <- sens %>% filter(scenario != "priors_tighter")


#  The grid

base_scenarios <- bind_rows(primary, sens) %>%
  mutate(
    fit_method = fit_method,
    backhaz_file = backhaz_file,
    trial_file = paste0(data_dir, "/trial_data", suffix, ".rds"),
    external_file = external_path(datasets, suffix)) %>%
  mutate(
    store_file = paste0(store_models, "/base_model_", row_number(), ".rds"),
    hazard_survival_file = paste0(store_models, "/base_model_",
                                  row_number(), "_hs.rds"),
    rmst_file = paste0(store_models, "/base_model_",
                       row_number(), "_rmst.rds"),
    diagnostics_file = paste0(store_models, "/base_model_",
                              row_number(), "_diag.rds"))

cat(nrow(base_scenarios), "scenarios:",
    sum(base_scenarios$scenario == "primary"), "primary,",
    sum(base_scenarios$scenario != "primary"), "sensitivity\n\n")
print(as.data.frame(base_scenarios %>%
                      count(scenario, model, datasets, name = "n")), row.names = FALSE)

missing <- unique(c(base_scenarios$trial_file,
                    na.omit(base_scenarios$external_file), backhaz_file))
missing <- missing[!file.exists(missing)]
if (length(missing)) stop("missing: ", paste(missing, collapse = ", "))


#  Fit with slurm

objects_attach <- c("survextrap_mem",
                    "add_knots_feb_trial", "add_knots_feb_ext",
                    "add_knots_dec_trial", "add_knots_dec_ext",
                    "add_knots_feb_10y", "add_knots_feb_3",
                    "takes_disc", "summarise_samples")

package_attach <- c("dplyr", "tidyr", "readr",
                    "survextrap", "rstan", "survival", "posterior")

check_status <- paste0("sacct -S ", as.character(Sys.Date() - 20),
                       " -u ", Sys.info()["user"],
                       " --format=JobID,Jobname,partition,state,elapsed,ncpus -X")

fit_base_slurm <- slurm_apply(
  fit_model,
  base_scenarios %>%
    select(model, datasets, df, add_knots, hsd_rate, eta_med, eta_up,
           tau_rate, fit_method, trial_file, external_file, backhaz_file,
           store_file),
  jobname = "fit_base",
  nodes = 22,
  cpus_per_node = 4,
  submit = TRUE,
  global_objects = objects_attach,
  pkgs = package_attach,
  slurm_options = list(time = '04:00:00',
                       "mem-per-cpu" = '16G'))

system(check_status)


#  Survival, hazard, rmst and diagnostics
#
# Wait for the fitting jobs to finish before submitting these.

run_over <- function(f, in_col, out_col, jobname, time) {
  slurm_apply(
    f,
    base_scenarios %>%
      select(all_of(c(in_col, out_col))) %>%
      rename(model_file = !!in_col, store_file = !!out_col),
    jobname = jobname,
    nodes = 22,
    cpus_per_node = 4,
    submit = TRUE,
    global_objects = objects_attach,
    pkgs = package_attach,
    slurm_options = list(time = time, "mem-per-cpu" = '16G'))
}

sh_base_slurm <- run_over(get_survival_and_hazard_survextrap,
                          "store_file", "hazard_survival_file",
                          "sh_base", '00:30:00')

rmst_base_slurm <- run_over(get_rmst_survextrap,
                            "store_file", "rmst_file",
                            "rmst_base", '00:30:00')

diag_base_slurm <- run_over(get_diagnostics_survextrap,
                            "store_file", "diagnostics_file",
                            "diag_base", '00:30:00')

system(check_status)


#  Combine
#
# Wait for the three jobs above to finish.

read_all <- function(col) {
  bind_rows(lapply(seq_len(nrow(base_scenarios)), function(i) {
    f <- base_scenarios[[col]][i]
    if (!file.exists(f)) return(NULL)
    readRDS(f) %>% mutate(store_file = base_scenarios$store_file[i])
  }))
}

base_results <- bind_rows(read_all("hazard_survival_file"),
                          read_all("rmst_file")) %>%
  left_join(base_scenarios, by = "store_file")

saveRDS(base_results, paste0(data_dir, "/base_model_all.rds"))

diagnostics <- read_all("diagnostics_file") %>%
  left_join(base_scenarios, by = "store_file")

saveRDS(diagnostics, paste0(data_dir, "/base_model_diagnostics.rds"))


#  Check the run
#
# A slurm array can lose individual tasks while the rest complete, so
# this is worth checking before anything is read into a table or figure.

cat("\nscenarios expected:", nrow(base_scenarios),
    " with results:", length(unique(base_results$store_file)), "\n")

lost <- setdiff(base_scenarios$store_file, unique(base_results$store_file))
if (length(lost)) {
  cat("missing results for", length(lost), "scenarios:\n")
  print(as.data.frame(base_scenarios %>%
                        filter(store_file %in% lost) %>%
                        select(scenario, data_cut, model, datasets)), row.names = FALSE)
}

# rhat above 1.01 means the chains have not mixed, and any divergent
# transitions mean part of the posterior was not explored. Either way
# that model's intervals should not be reported.
flagged <- diagnostics %>%
  filter(max_rhat > 1.01 | divergences > 0 | is.na(max_rhat))

cat("\nmodels with poor convergence:", nrow(flagged),
    "of", nrow(diagnostics), "\n")

if (nrow(flagged)) {
  print(as.data.frame(flagged %>%
                        select(scenario, data_cut, model, datasets, arm,
                               max_rhat, min_ess, divergences) %>%
                        arrange(desc(max_rhat))), row.names = FALSE)
}