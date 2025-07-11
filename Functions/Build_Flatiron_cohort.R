
###################################################
# Function for creating 1L Flatiron RWE dataset
###################################################

create_flatiron_data <- function(treatment = NULL,
                                 data_cut_off_date = as.Date(Inf),
                                 censoring_strategy = "last_confirmed_date",
                                 SoC_date = NULL){
  
  rwdcohort::rwdcohort_be_noisy()
  
  base_cohort <- create_cohort_rwd_flatiron_nsclc_base_1l(
    force = F,
    flatiron_cohort = "nsclc",
    directory = "Flatiron_data",
    output_format = "cohort"
  )  %>%
    add_patient_time_at_risk(
      patient_id_column = "patientid",
      index_date_column = "lot_1_startdate",
      entry_date_column = "cohort_entry_date",
      event_date_column = "dateofdeath_imputed_date",
      data_cut_off_date = data_cut_off_date,
      last_known_date_column = "last_confirmed_activity_date",
      censoring_strategy = censoring_strategy,
      prefix = "os1")

  filtered_cohort <- base_cohort |>
    derive_age_at_event(
      event_col = "lot_1_startdate",
      def_birth_date_col = "birthyear_imputed_date"
    ) %>%
    rwdcohort::derive_custom(
      is_stage_iiib_or_iv = (groupstage %in% c("Stage IIIB", "Stage IV", "Stage IVA", "Stage IVB", "Stage IIIB")) | (ismetastatic == TRUE),
      # Ecog 0/1/2
      ecog_0_1 = ecog %in% c(0, 1),
      ecog_0_1_2 = ecog %in% c(0, 1, 2),
      # Receiving the appropriate treatments
      is_alex_regimen = lot_1_linename == treatment,
      # Receive the treatments after the date that treatment became SoC
      tx_start_after_2011 = lot_1_startdate >= SoC_date,
      age_is_known = !is.na(age),
      sex_is_known = sex %in% c("Male", "Female"),
      race_is_known = !is.na(race), 
      smoking_status_is_known = smoking_status %in% c("Ever Smoker", "Never Smoker"),
      ecog_0_1_2_is_known = !is.na(ecog_0_1_2),
      os_is_known = !is.na(os_time_from_index),
      os1_is_known = !is.na(os1_time_from_index),
      race_asian = race == "Asian",
      key_covariates_known = age_is_known & sex_is_known & race_is_known &
        smoking_status_is_known  & ecog_0_1_2_is_known
    ) |>
    rwdcohort::include_custom("is_stage_iiib_or_iv", "Diagnosed with stage IIIB or IV") |>
    rwdcohort::include_custom("tx_start_after_2011", "1L start after 2011-08-26") |>
    rwdcohort::include_custom("is_alk_positive", "ALK Positive") |>
    rwdcohort::include_custom("is_alex_regimen", "ALEX Regimen") |>
    rwdcohort::include_custom("ecog_0_1_2", "ECOG 0, 1 or 2") |>
    rwdcohort::include_custom("key_covariates_known", "Key covariates known") |>
    rwdcohort::include_custom("os1_is_known", "Overall survival known") 
  
  cohort_data <- filtered_cohort %>%
     pull_cohort() %>%
     mutate(time = os1_time_from_index,
            status = os1_event_status) %>%
    mutate(time = time/12)
  
  return(cohort_data)
  
}











