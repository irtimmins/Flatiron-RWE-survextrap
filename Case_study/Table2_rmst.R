

base_model_all <- readRDS("/projects/aa/klvq491/Flatiron_ansclc/models/base_model_all.rds")

#base_rmst <-  
  
round1 <- function(x)(
  sprintf("%.2f", round(x, digits = 2)))

survextrap_results <- base_model_all %>%
  filter(variable %in% c("rmst" ,"irmst")) %>%
  filter(t %in% c(5,20)) %>%
  mutate_at(c("median", "lower", "upper"),round1) %>%
  mutate(value = paste0(median, " (", lower, ", ", upper,")")) %>%
  select(variable, trt, t, value, model, datasets, df, hsd_rate, hrsd_rate, store_file) %>%
  mutate(trt_label = case_when(trt == "Crizotinib" ~ "control",
                               trt == "Alectinib" ~ "active",
                               is.na(trt) ~ "contrast"),
           .before = "model") %>%
  mutate(t_label = case_when(t == 5 ~ "short",
                               t == 20 ~ "long"),
         .before = "model") %>%
  group_by(store_file) %>%
  select(-c(t,trt)) %>%
  pivot_wider(names_from = c("variable","t_label","trt_label"), values_from = c("value")) %>%
  # mutate(model_type = "Survextrap", .before = "model") %>%
  mutate(model_label = case_when(model == "PH" ~ "PH",
                                 model == "NON-PH" ~ "Non-PH",
                                 model == "Separate_arms" ~ "Separate arms"),
         .before = "model")  %>%
  mutate(data_label = case_when(datasets == "trial_only" ~ "ALEX trial only",
                                datasets == "trial_and_historic" ~ "ALEX+PROFILE-1014",
                                datasets == "trial_and_all" ~ "ALEX+PROFILE-1014+Flatiron RWE"),
         .after = model_label) %>%
  mutate(model_label = paste0("Survextrap, ", model_label)) %>%
  mutate(sigma_label = paste0("Gamma(2,", hsd_rate, ")"), .before = df) %>%
  mutate(tau_label = case_when(model == "NON-PH" ~ paste0("Gamma(2,", hrsd_rate, ")"),
                               model != "NON-PH" ~ NA),
         .after = sigma_label) #%>%
  # ungroup() %>%
  # select(-c(store_file,model, hsd_rate, hrsd_rate)) %>%


#View(base_rmst)


primary_results <- 
  survextrap_results %>%
  filter(df == 6, hsd_rate == 1, hrsd_rate == 10 | is.na(hrsd_rate)) %>%
  arrange(model) %>%
  arrange(data_label) %>%
  ungroup(store_file) %>%
  select(-c(store_file, hsd_rate, hrsd_rate, model, datasets))

primary_results %>% View()




hybrid_rmst <- readRDS("/projects/aa/klvq491/Flatiron_ansclc/models/hybrid_model_rmst.rds")
exp_rmst <- readRDS("/projects/aa/klvq491/Flatiron_ansclc/models/exp_model_rmst.rds")


hybrid_rmst <- hybrid_rmst %>%
  mutate(model_label = "Piecewise, exponential tail")

exp_rmst <- exp_rmst %>%
  mutate(model_label = "Exponential")

nice_ta536_results <- exp_rmst %>%
  bind_rows(hybrid_rmst) %>%
  filter(variable %in% c("rmst" ,"irmst")) %>%
  rename(t = time) %>%
  filter(t %in% c(5,20)) %>%
  mutate_at(c("value", "lower", "upper"),round1) %>%
  mutate(value = paste0(value, " (", lower, ", ", upper,")")) %>%
  select(-c(se, lower, upper)) %>%
  mutate(trt_label = case_when(trt == "Crizotinib" ~ "control",
                               trt == "Alectinib" ~ "active",
                               is.na(trt) ~ "contrast")) %>%
  mutate(t_label = case_when(t == 5 ~ "short",
                             t == 20 ~ "long")) %>%
  select(-c(t,trt)) %>%
  mutate(data_label = "ALEX trial only", .after = "model_label") %>%
  pivot_wider(names_from = c("variable","t_label","trt_label"), values_from = c("value")) 
  



###
trial_data <- readRDS("Data/trial_IPD_OS_ALEX_Nov_2019.rds")

control_valid <- survfit(Surv(time, status) ~ 1, data = trial_data %>% filter(trt == "Crizotinib"))
km_summary_control <- summary(control_valid, rmean=5)
RMST_control <- km_summary_control$table["rmean"] %>% unname()
RMST_se_control <- km_summary_control$table["se(rmean)"] %>% unname()


active_valid <- survfit(Surv(time, status) ~ 1, data = trial_data %>% filter(trt == "Alectinib"))
km_summary_active <- summary(active_valid, rmean=5)
RMST_active <- km_summary_active$table["rmean"] %>% unname()
RMST_se_active <- km_summary_active$table["se(rmean)"] %>% unname()

RMST_difference <- RMST_active-RMST_control
RMST_se_difference <- sqrt(RMST_se_control^2+RMST_se_active^2)

validation_results <- 
  tibble(value = c(RMST_control, RMST_active, RMST_difference),
         se = c(RMST_se_control, RMST_se_active, RMST_se_difference),
         trt_label = c("control", "active", "contrast"),
         variable = c("rmst", "rmst", "irmst")) %>%
  mutate(model_label = "Validation, Kaplan-Meier",
         data_label = "ALEX trial, final OS at 5 years") %>%
  mutate(upper = value+1.96*se,
         lower = value-1.96*se) %>%
  mutate_at(c("value", "lower", "upper"),round1) %>%
  mutate(value = paste0(value, " (", lower, ", ", upper,")")) %>%
  select(-c(se, lower, upper)) %>%
  select(model_label, data_label, trt_label, value, variable) %>%
  mutate(t = 5) %>%
  mutate(t_label = case_when(t == 5 ~ "short",
                           t == 20 ~ "long")) %>%
  select(-t) %>%
  pivot_wider(names_from = c("variable","t_label","trt_label"), values_from = c("value")) 


table2 <- primary_results %>%
  bind_rows(nice_ta536_results) %>%
  bind_rows(validation_results) %>%
  select(model_label, data_label, rmst_short_control, rmst_short_active, irmst_short_contrast,
         rmst_long_control, rmst_long_active, irmst_long_contrast)
View(table2)
library(readr)
write_csv(table2, "Tables/Table2_rmst.csv")  

