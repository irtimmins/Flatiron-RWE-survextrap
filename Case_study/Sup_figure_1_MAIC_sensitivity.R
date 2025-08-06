
######################################################
#  Load packages.
######################################################

library(populationmodels)
library(rwdcohort)
library(ggplot2)
library(survival)
library(survminer)
library(maic)
library(azci)
library(bshazard)
library(dplyr)
library(tidyr)
library(ggplot2)
library(survival)
library(survminer)
library(survextrap)
library(ggh4x)
library(condsurv)
library(cowplot)
library(flexsurv)

######################################################
#  Derive Flatiron cohort for Crizotinib
######################################################

# Source dependencies.
source("Functions/Build_Flatiron_cohort.R")
source("Functions/Derive_conditional_survival.R")
source("Functions/Hybrid_model.R")

######################################################
#  Derive Flatiron cohort for Crizotinib
######################################################

# 262 - 39 ALK confirmed negative, 206 ALK confirmed positive. 17 unknown.
# Filter down to 206 ALK confirmed positive.

rwd_data <- create_flatiron_data(treatment = "Crizotinib",
                                 data_cut_off_date = as.Date("2017-02-09"),
                                 censoring_strategy = "administrative_cutoff",
                                 SoC_date = as.Date("2011-08-26"))


summary(rwd_data$age_at_lot_1_startdate)


######################################################
#  Apply MAIC weighting to Flatiron cohort
######################################################

ALEX_table1 <-
  tibble::tribble(
    ~characteristic, ~mean, ~sd, ~n,
    "age", 53.8, 13.5, 151,
    "sex_male", 0.424, NA, 64,
    "sex_female", 0.576, NA, 87,
    "race_asian", 0.457, NA,  69,
    "race_non_asian", 0.543, NA, 82,
    "ecog_0_1", 0.934, NA, 141,
    "ecog_2", 0.066, NA, 10,
    "smoking_status_ever", 0.351, NA, 53,
    "smoking_status_never", 0.649, NA, 98,
    "brain_mets", 0.384, NA, 58
  )

# 

ALEX_variables1 =  c(
  "age"
)

ALEX_variables2 =  c(
  "age",
  "sex_male"
)

ALEX_variables3 =  c(
  "age",
  "sex_male",
  "ecog_0_1"
)

ALEX_variables4 =  c(
  "age",
  "sex_male",
  "ecog_0_1",
  "smoking_status_ever"
)

ALEX_variables5 =  c(
  "age",
  "sex_male",
  "ecog_0_1",
  "smoking_status_ever",
  "brain_mets"
)

ALEX_variables6 =  c(
  "age",
  "sex_male",
  "ecog_0_1",
  "smoking_status_ever",   
  "brain_mets",
  "race_asian"
)

# Cohort with MAIC weights.

i <- 1
trial_data <- readRDS("Data/trial_IPD_OS_ALEX_Feb_2017_v4.rds")
historic_trial_data <- readRDS("Data/trial_IPD_OS_PROFILE_1014_Nov_2016.rds")

for(i in 1:6){

  rwd_data_weighted <- create_weighted_cohort(
    cohort_data = rwd_data,
    reference_table1 = ALEX_table1,
    match_variables = get(paste0("ALEX_variables", i)),
  )[["cohort_data"]]  %>%
    rename(weight = patient_weight)
  
  ESS <- sum(rwd_data_weighted$weight)^2/sum(rwd_data_weighted$weight^2)
  assign(paste0("ESS", i), ESS)
  
  
  #head(rwd_data_weighted$weight)
  # Cohort without MAIC weights.
  
  rwd_data_unweighted <- rwd_data_weighted %>%
    mutate(weight = 1,
           effective_sample_size = n())
  
  ######################################################
  #  Plot RWD KM curves alongside trial data.
  ######################################################
  
  
  
  km_trial_data  <- survfit(Surv(time, status) ~ trt, data=trial_data)
  km_trial_plot <- ggsurvplot(km_trial_data, data=trial_data)["data.survplot"] [[1]] %>%
    select(time, surv, trt) %>%
    mutate(dataset = "ALEX") %>%
    as_tibble() %>%
    add_row(time = 0, surv = 1, trt = "Alectinib", dataset = "ALEX") %>%
    add_row(time = 0, surv = 1, trt = "Crizotinib", dataset = "ALEX") %>%
    mutate(trt = as.factor(trt))
  
  km_historic_trial_data  <- survfit(Surv(time, status) ~ 1, data=historic_trial_data)
  km_historic_trial_plot <- ggsurvplot(km_historic_trial_data, data=historic_trial_data)["data.survplot"][[1]] %>%
    mutate(trt = "Crizotinib",
           dataset = "PROFILE-1014") %>%
    select(time, surv, trt, dataset) %>%
    as_tibble()
  
  km_rwd_data_weighted <- survfit(Surv(time, status) ~ 1, data=rwd_data_weighted, 
                                  weights = rwd_data_weighted$weight)
  km_rwd_plot <- ggsurvplot(km_rwd_data_weighted, data=rwd_data_weighted)["data.survplot"][[1]] %>%
    mutate(trt = "Crizotinib",
           dataset = "Flatiron RWE MAIC Weighted") %>%
    select(time, surv, trt, dataset)
  
  km_rwd_data_unweighted <- survfit(Surv(time, status) ~ 1, data=rwd_data_unweighted)
  km_rwd_plot_unweighted <- ggsurvplot(km_rwd_data_unweighted, data=rwd_data_unweighted)["data.survplot"][[1]] %>%
    mutate(trt = "Crizotinib",
           dataset = "Flatiron RWE Unweighted") %>%
    select(time, surv, trt, dataset)
  
  if(i == 1){
    
    comb_maic <- bind_rows(km_trial_plot,
                           km_historic_trial_plot,
                           km_rwd_plot,
                           km_rwd_plot_unweighted) %>%
      mutate(maic = i,
             ESS = ESS)
    
    comb_maic_weights <- rwd_data_weighted %>%
      mutate(maic = i)
    
  } else{
    
    temp_maic <- bind_rows(km_trial_plot,
                           km_historic_trial_plot,
                           km_rwd_plot,
                           km_rwd_plot_unweighted)  %>%
      mutate(maic = i,
             ESS = ESS)
    
    comb_maic <- bind_rows(comb_maic,  
                           temp_maic)
    
    temp_maic_weights <-  rwd_data_weighted %>%
      mutate(maic = i)
    
    comb_maic_weights <- bind_rows(comb_maic_weights,  
                                   temp_maic_weights)
    
    
  }
  
  
   rwd_data_weighted %>%
    ggplot()+
    geom_histogram(aes(x = weight), colour = "gray50", binwidth = 0.01)+
    #  geom_histogram(aes(x = weight, colour = race_asian, fill = race_asian), binwidth = 0.01)+
    theme_classic()+
    scale_y_continuous("Frequency")+
    scale_x_continuous("Weight")
  
  
  

}


#summary(as.factor(comb_maic$maic))
maic_ESS <- tibble(x = 4, y = 0.8, maic = 1:6, label = 0,
       dataset = "ALEX")
for(i in 1:nrow(maic_ESS)){
  maic_ESS$label[i] <- paste0("ESS=", round(get(paste0("ESS", i)), digits = 1))
}
maic_ESS <- maic_ESS %>% 
  mutate(maic = factor(maic, 
         levels = 1:6,
         labels = paste0("MAIC adjusted for\n",
                         c("Age",
                           "Age, Sex",
                           "Age, Sex, ECOG",
                           "Age, Sex, ECOG, \nSmoking status",
                           "Age, Sex, ECOG, \nSmoking status, \nBrain mets.",
                           "Age, Sex, ECOG, \nSmoking status,\nBrain mets., Asian race"))))
maic_ESS

plot_maic_a <- 
comb_maic %>%
  mutate(dataset = factor(dataset, 
                          levels = c("ALEX", 
                                     "PROFILE-1014", 
                                     "Flatiron RWE MAIC Weighted", 
                                     "Flatiron RWE Unweighted"),
                          labels = c("ALEX", 
                                     "PROFILE-1014", 
                                     "Flatiron RWE \nMAIC weighted", 
                                     "Flatiron RWE \nunweighted"))) %>%
  filter(dataset != "PROFILE-1014") %>%
  filter(trt == "Crizotinib") %>%
  mutate(maic = factor(maic, 
                       levels = 1:6,
                       labels = paste0("MAIC adjusted for\n",
                                       c("Age",
                                  "Age, Sex",
                                  "Age, Sex, ECOG",
                                  "Age, Sex, ECOG, \nSmoking status",
                                  "Age, Sex, ECOG, \nSmoking status, \nBrain mets.",
                                  "Age, Sex, ECOG, \nSmoking status,\nBrain mets., Asian race")))) %>%
  ggplot(aes(x = time, y = surv, colour = dataset))+
  theme_classic()+
  theme(#legend.position.inside = c(0.66,0.93),
        #legend.position = "inside",
        legend.text=element_text(size=9),
        legend.spacing.y = unit(6, "pt"),
        legend.spacing.x = unit(0, "pt"),
        legend.key.spacing.y = unit(6, "pt"),
        legend.box.spacing = unit(-10, "pt"),
        legend.title = element_text(size=12,
                                    margin = margin(l = 0, r = 0, b = 0, t = 0)),
        strip.text.x = element_text(size = 8)) +
  geom_text(data = maic_ESS,
            aes(x = x, y = y, label = label),
            x = 3.5, y = 0.9, size = 3, fontface = "bold", colour = "black")+
  geom_step(linewidth = 0.8)+
  scale_colour_manual("Dataset", 
                      values = c("#d6604d", "#7CAE00", "#00BFC4", "#C77CFF" ))+
  scale_y_continuous("Overall Survival", 
                     limits = c(0,1), labels = scales::percent)+
  scale_x_continuous("Time (years)",  breaks = 0:5, limits = c(0,5.5))+
  facet_wrap(~maic, ncol = 3)
#ggtitle("Crizotinib, Overall Survival, Feb 2017 Datacut")
#plot_maic_a




plot_maic_b <- comb_maic_weights %>%
  mutate(maic = factor(maic, 
                       levels = 1:6,
                       labels = paste0("MAIC adjusted for\n",c("Age",
                                                               "Age, Sex",
                                                               "Age, Sex, ECOG",
                                                               "Age, Sex, ECOG, \nSmoking status",
                                                               "Age, Sex, ECOG, \nSmoking status, \nBrain mets.",
                                                               "Age, Sex, ECOG, \nSmoking status,\nBrain mets., Asian race")))) %>%
  ggplot()+
  geom_histogram(aes(x = weight), colour = "gray50", binwidth = 0.01)+
  #  geom_histogram(aes(x = weight, colour = race_asian, fill = race_asian), binwidth = 0.01)+
  theme_classic()+
  theme(#legend.position.inside = c(0.66,0.93),
    #legend.position = "inside",
    legend.text=element_text(size=9),
    legend.spacing.y = unit(6, "pt"),
    legend.spacing.x = unit(0, "pt"),
    legend.key.spacing.y = unit(6, "pt"),
    legend.box.spacing = unit(-10, "pt"),
    legend.title = element_text(size=12,
                                margin = margin(l = 0, r = 0, b = 0, t = 0)),
    strip.text.x = element_text(size = 8)) +
  scale_y_continuous("Frequency")+
  scale_x_continuous("Weight")+
 facet_wrap(~maic, ncol = 3)




figure_maic <- plot_grid(plot_maic_a+
                            theme(  plot.title = element_text(hjust = -0.1,
                                                              size=14, face="bold"))+
                            labs(title = "(a)"),
                         plot_maic_b+
                            theme(  plot.title = element_text(hjust = -0.1,
                                                              size=14, face="bold"))+
                            labs(title = "(b)"),
                          align = "v",
                          rel_heights  =c(0.5,0.5),
                          ncol = 1)
figure_maic


tiff(file = "Figures/Sup_figure_1_maic.tiff",   
     width = 7.6, 
     height = 9.5,
     units = 'in',  
     res = 300, 
     compression = "lzw")
print(figure_maic)
dev.off()


