######################################################
#  ...
######################################################
# test
rwd_data

rwdcohort::rwdcohort_be_noisy()

# create_flatiron_data(treatment = "Crizotinib",
#                      data_cut_off_date = as.Date("2017-02-09"),
#                      censoring_strategy = "administrative_cutoff",
#                      SoC_date = as.Date("2011-08-26"))

base_cohort <- create_cohort_rwd_flatiron_nsclc_base_1l(
  force = F,
  flatiron_cohort = "nsclc",
  directory = "Flatiron_data",
  output_format = "cohort") 

test_base <- base_cohort %>% pull_cohort()

base_cohort_plus_alk <- base_cohort %>%
  connect_refresh() %>%
  derive_custom(
    `60_days_before_1l` = lot_1_startdate - 60,
    `60_days_after_1l` = lot_1_startdate + 60
  ) %>%
  add_patient_biomarker(biomarker_name = "ALK",
                        def_obs_window_start = "60_days_before_1l",
                        def_obs_window_end = "60_days_after_1l",
                        prefix = "new60",
                        result_reduce_strategy = c("results_exist")) %>%
  derive_age_at_event(
    event_col = "lot_1_startdate",
    def_birth_date_col = "birthyear_imputed_date"
  ) %>%
  derive_custom(
    is_stage_iiib_or_iv = (groupstage %in% c("Stage IIIB", "Stage IV", "Stage IVA", "Stage IVB", "Stage IIIB")) | (ismetastatic == TRUE),
    # Ecog 0/1/2
    ecog_0_1 = ecog %in% c(0, 1),
    ecog_0_1_2 = ecog %in% c(0, 1, 2),
    # Receiving the appropriate treatments
    is_alex_regimen = lot_1_linename == "Crizotinib",
    # Receive the treatments after the date that treatment became SoC
    tx_start_after_2011 = lot_1_startdate >= as.Date("2011-08-26"),
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
  ) %>%
  rwdcohort::include_custom("is_stage_iiib_or_iv", "Diagnosed with stage IIIB or IV")  |>
  rwdcohort::include_custom("new60_alk_status_entire_obs_window_any_positive", "ALK Positive") |>
  rwdcohort::include_custom("is_alex_regimen", "ALEX Regimen") |>
  rwdcohort::include_custom("tx_start_after_2011", "1L start after 2011-08-26") |>
  rwdcohort::include_custom("os1_is_known", "Overall survival known") |>
  rwdcohort::include_custom("ecog_0_1_2", "ECOG 0, 1 or 2") |>
  rwdcohort::include_custom("key_covariates_known", "Key covariates known") 
  #   add_patient_time_at_risk(
#     patient_id_column = "patientid",
#     index_date_column = "lot_1_startdate",
#     entry_date_column = "cohort_entry_date",
#     event_date_column = "dateofdeath_imputed_date",
#     data_cut_off_date = data_cut_off_date,
#     last_known_date_column = "last_confirmed_activity_date",
#     censoring_strategy = censoring_strategy,
#     prefix = "os1")

test <- base_cohort_plus_alk %>% pull_cohort()
names(test)
summary(test$new60_alk_status_entire_obs_window_any_positive)
summary(test$alk_status_entire_obs_window_any_positive)

all_cohort$diagnosisdate

filtered_cohort <- base_cohort %>%
  connect_refresh() %>%
  derive_custom(
  #  `anytime_before_1l` = -Inf,
  #  `60_days_after_1l` = lot_1_startdate + 60
    `anytime_before_1l` = diagnosisdate-90,
    `60_days_after_1l` = diagnosisdate+365) %>%
  #   `anytime_before_1l` = -Inf,
  #  `60_days_after_1l` = Inf) %>%
  add_patient_biomarker(biomarker_name = "ALK",
                        def_obs_window_start = "anytime_before_1l",
                        def_obs_window_end = "60_days_after_1l",
                        result_reduce_strategy = "results_exist",
                        prefix = "new60") %>%
  add_patient_time_at_risk(
    patient_id_column = "patientid",
    index_date_column = "lot_1_startdate",
    entry_date_column = "cohort_entry_date",
    event_date_column = "dateofdeath_imputed_date",
    data_cut_off_date = as.Date("2017-02-09"),
    last_known_date_column = "last_confirmed_activity_date",
    censoring_strategy = "administrative_cutoff",
    prefix = "os1") %>%
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
    is_alex_regimen = lot_1_linename == "Crizotinib",
    # Receive the treatments after the date that treatment became SoC
    tx_start_after_2011 = lot_1_startdate >= as.Date("2011-08-26"),
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
  rwdcohort::include_custom("is_stage_iiib_or_iv", "Diagnosed with stage IIIB or IV") #  |>
  rwdcohort::include_custom("new60_alk_status_entire_obs_window_any_positive", "ALK Positive") |>
  rwdcohort::include_custom("is_alex_regimen", "ALEX Regimen") |>
  rwdcohort::include_custom("tx_start_after_2011", "1L start after 2011-08-26") |>
  rwdcohort::include_custom("os1_is_known", "Overall survival known") |>
  rwdcohort::include_custom("ecog_0_1_2", "ECOG 0, 1 or 2") |>
  rwdcohort::include_custom("key_covariates_known", "Key covariates known") 

names(all_cohort)
all_cohort <- filtered_cohort %>% pull_cohort()
summary(all_cohort$obs_window_start)
summary(all_cohort$obs_window_end_biomarkers)
summary(all_cohort$alk_status_entire_obs_window_any_positive)
summary(all_cohort$is_alk_positive)
summary(all_cohort$diagnosisdate)
summary(as.factor(all_cohort$is_alk_confirmed_negative))

summary(as.factor(all_cohort$new60_alk_status_entire_obs_window_any_positive))

summary(as.factor(all_cohort$new60_alk_status_entire_obs_window_results_exist))



# all_cohort$win

library(ggplot2)
# all_cohort %>%
#   mutate(year_diagnosis = format(diagnosisdate, "%Y")) %>%
#   mutate(year_1l = format(lot_1_startdate, "%Y")) %>%
#   ggplot(aes(x =year_diagnosis , y = year_1l ))+
#   theme_classic()+
#   geom_point()

test1 <- all_cohort %>%
  mutate(year_diagnosis = format(diagnosisdate, "%Y")) %>%
  mutate(year_diagnosis = as.integer(year_diagnosis)) %>%
  mutate(alk_status_known = new60_alk_status_entire_obs_window_results_exist) %>%
  group_by(year_diagnosis) %>%
  summarise( count = n(),
    proportion_test = 100*sum(alk_status_known)/n()) %>%
  mutate(proportion_test = round(proportion_test, 1)) %>%
  filter(year_diagnosis >= 2011, year_diagnosis <= 2017)#%>%
  View()
  
  test1
  
  
  test2 <- all_cohort %>%
    mutate(year_diagnosis = format(diagnosisdate, "%Y")) %>%
    mutate(year_diagnosis = as.integer(year_diagnosis)) %>%
    mutate(alk_status_known = new60_alk_status_entire_obs_window_results_exist) %>%
    group_by(year_diagnosis) %>%
    summarise( test = sum(alk_status_known),
               count = n(),
               proportion_test = 100*sum(alk_status_known)/n()) %>%
    mutate(proportion_test = round(proportion_test, 1)) %>%
    filter(year_diagnosis >= 2011, year_diagnosis <= 2017)#%>%
  
  test2
  View()

  
  test2
  