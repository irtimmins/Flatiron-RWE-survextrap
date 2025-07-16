

library(populationmodels)
library(rwdcohort)


rwdcohort::rwdcohort_be_noisy()

cohort1 <- create_cohort_rwd_flatiron_nsclc(
  force = F,
  flatiron_cohort = "nsclc",
  directory = "Flatiron_data",
  output_format = "cohort"
) %>%
  add_patient_time_at_risk(
    patient_id_column = "patientid",
    index_date_column = "lot_1_startdate",
    entry_date_column = "cohort_entry_date",
    event_date_column = "dateofdeath_imputed_date",
    data_cut_off_date = as.Date("2017-02-09"),
    last_known_date_column = "last_confirmed_activity_date",
    censoring_strategy = "administrative_cutoff",
    prefix = "os1") %>%
  rwdcohort::derive_custom(
    is_stage_iiib_or_iv = (groupstage %in% c("Stage IIIB", "Stage IV", "Stage IVA", "Stage IVB", "Stage IIIB")) | (ismetastatic == TRUE),
    # Ecog 0/1/2
    ecog_0_1 = ecog %in% c(0, 1),
    ecog_0_1_2 = ecog %in% c(0, 1, 2),
    # Receiving the appropriate treatments
    is_alex_regimen = lot_1_linename == "Crizotinib",
    # Receive the treatments after the date that treatment became SoC
    tx_start_after_2011 = lot_1_startdate >= as.Date("2011-08-26") & lot_1_startdate < as.Date("2017-02-09"),
    age_is_known = !is.na(age),
    sex_is_known = sex %in% c("Male", "Female"),
    race_is_known = !is.na(race), 
    smoking_status_is_known = smoking_status %in% c("Ever Smoker", "Never Smoker"),
    ecog_0_1_2_is_known = !is.na(ecog_0_1_2),
    os_is_known = !is.na(os_time_from_index),
    os1_is_known = !is.na(os1_time_from_index),
    race_asian = race == "Asian",
    key_covariates_known = age_is_known & sex_is_known &  race_is_known &
      smoking_status_is_known  & ecog_0_1_2_is_known
  )

# cohort1 %>%  pull_cohort() %>% nrow()
# Stage 4.
# test1 <- cohort1 %>% pull_cohort()
# summary(test1$ismetastatic)

cohort2 <- cohort1 %>% rwdcohort::include_custom("is_stage_iiib_or_iv", "Diagnosed with stage IIIB or IV")
test2 <- cohort2 %>% pull_cohort()
cohort2 %>% pull_cohort() %>% nrow()
summary(as.factor(test2$groupstage))
summary(as.factor(test2$is_stage_iv))
summary(as.factor(test2$ismetastatic))
# sum(test$age >= 18)
#cohort3 <- cohort2 %>% rwdcohort::include_custom(col_incl = "is_advanced", description_incl = "Advanced NSCLC")
cohort3 <- cohort2 %>% 
  derive_age_at_event(
    event_col = "lot_1_startdate",
    def_birth_date_col = "birthyear_imputed_date"
  ) %>%
  rwdcohort::derive_custom("initiatied 1L treatment"= !is.na(lot_1_linesetting)) %>%
  rwdcohort::derive_custom("age_18_at_lot_1_startdate"= age_at_lot_1_startdate >= 18) %>%
  rwdcohort::include_custom("initiatied 1L treatment", "1L initiation") %>%
  rwdcohort::include_custom("age_18_at_lot_1_startdate", "age 18 at 1L initiation")

#test <- cohort3 %>% pull_cohort()
#summary(test$age_at_lot_1_startdate)
#summary(test$age_18_at_lot_1_startdate)

cohort4  <- cohort3 %>%  rwdcohort::include_custom("tx_start_after_2011", "1L start after 2011-08-26") %>%
   rwdcohort::include_custom("is_alex_regimen", "ALEX Regimen")

cohort5 <- cohort4 %>% rwdcohort::include_custom("is_alk_positive", "ALK Positive") 
cohort6  <- cohort5 %>%  rwdcohort::include_custom("ecog_0_1_2", "ECOG 0, 1 or 2")
cohort7  <- cohort6 %>%  rwdcohort::include_custom("key_covariates_known", "Key covariates known")
cohort8  <- cohort7 %>%  rwdcohort::include_custom("os1_is_known", "Overall survival known") 


create_cohort_rwd_flatiron_nsclc <- function (force, flatiron_cohort = c("nsclc", "nsclc_cg"), directory, 
          output_format = "cohort", override_schema = NULL, ...) 
{
  additional_args <- list(...)
  flatiron_cohort <- match.arg(flatiron_cohort)
  init_cohort <- rwdcohort::start_rwdcohort(data_provider = "flatiron", 
                                            cohort = flatiron_cohort, override_schema = override_schema)
  init_cohort$version <- get_latest_version(init_cohort)
  version_date <- as.Date(init_cohort$version, "%Y%m%d")
  name_mapping <- list(nsclc = list(name = "advnsclcedm", source_name = "EDM", 
                                    dataset = "flatiron_edm"), nsclc_cg = list(name = "cgdb", 
                                                                               source_name = "CGDB", dataset = "flatiron_cgdb"))
  init_cohort$name <- as.character(glue::glue("rwd_flatiron_nsclc_{name_mapping[[flatiron_cohort]]$name}_base_1l"))
  init_cohort$data_sources <- tibble::tibble(name = as.character(glue::glue("Flatiron {name_mapping[[flatiron_cohort]]$source_name}")), 
                                             type = "RWD", link = populationmodels:::get_dataset_description_link(name_mapping[[flatiron_cohort]][["dataset"]]), 
                                             aliases = as.character(glue::glue("Flatiron {name_mapping[[flatiron_cohort]]$source_name} ({format(version_date, '%Y %b %d')})")), 
                                             publication_name = NA_character_, publication_year = as.integer(lubridate::year(version_date)))
  init_cohort$output_format <- output_format
  create_cohort_fn <- function(init_cohort) {
    regimen_list <- populationmodels:::get_regimen()
    cohort <- rwdcohort::derive_custom(rwdcohort::define_index_date(rwdcohort::add_patient_lineoftherapy(init_cohort, 
                                                                                                         1), "lot_1_startdate"), obs_window_start = as.Date(-Inf), 
                                       obs_window_end = lot_1_startdate, obs_window_end_biomarkers = lot_1_startdate + 
                                         28, fu_window_start = lot_1_startdate + 1, fu_window_end = as.Date(Inf), 
    )
    cohort <- populationmodels:::add_patient_lung_characteristics(rwdcohort = cohort, 
                                               index_date = "lot_1_startdate", obs_window_start = "obs_window_start", 
                                               obs_window_end = "obs_window_end", obs_window_end_biomarkers = "obs_window_end_biomarkers", 
                                               fu_window_start = "fu_window_start", fu_window_end = "fu_window_end", 
                                               cohort_entry = "cohort_entry_date")
    cohort <- rwdcohort::derive_custom(rwdcohort::add_patient_lineoftherapy(rwdcohort::add_patient_lineoftherapy(cohort, 
                                                                                                                 2), 3), is_advanced = !is.na(advanceddiagnosisdate), 
                                       is_treated = !is.na(lot_1_startdate), is_treated_2l = !is.na(lot_2_startdate))
    cohort <- rwdcohort::add_patient_time_at_risk(rwdcohort::add_patient_lastobserveddate(rwdcohort::add_patient_treatment_settings(cohort)), 
                                                  patient_id_column = "patientid", index_date_column = "lot_1_startdate", 
                                                  entry_date_column = "cohort_entry_date", event_date_column = "dateofdeath_imputed_date", 
                                                  last_known_date_column = "last_confirmed_activity_date", 
                                                  censoring_strategy = "last_confirmed_date", prefix = "os")

    if (cohort$cohort_name == "nsclc_cg") {
      if ("discard_tw" %in% names(additional_args)) {
        discard_tw = additional_args[["discard_tw"]]
        message(glue::glue("A discard time window after index date was provided: {paste(discard_tw)} days"))
      }
      else {
        discard_tw = 14
        message(glue::glue("The default discard time window after index date will be used: {paste(discard_tw)} days"))
      }
      if ("last_cn_to_death_cuttoff_time" %in% names(additional_args)) {
        last_cn_to_death_cuttoff_time = additional_args[["last_cn_to_death_cuttoff_time"]]
        message(glue::glue("A cut-off time between 'last clinic date' and 'death' was provided: {paste(last_cn_to_death_cuttoff_time)} days"))
      }
      else {
        last_cn_to_death_cuttoff_time = 90
        message(glue::glue("A default cut-off time between 'last clinic date' and 'death' will be used: {paste(last_cn_to_death_cuttoff_time)} days"))
      }
      cohort <- extract_pfs_endpoint(cohort, index_lot_column = "lot_1_startdate", 
                                     next_lot_column = "lot_2_startdate", discard_tw = discard_tw, 
                                     last_cn_to_death_cuttoff_time = last_cn_to_death_cuttoff_time)
    }
    cohort$cohort$cohort_columns_n <- ncol(cohort$cohort$cohort_tbl)

    if (init_cohort$output_format == "list") {
      cohort <- list(cohort)
    }
    return(cohort)
  }
  cohort <- get_or_create(init_cohort = init_cohort, force = force, 
                          directory = directory, create_function = create_cohort_fn, 
                          output_format = output_format)
  return(cohort)
}
