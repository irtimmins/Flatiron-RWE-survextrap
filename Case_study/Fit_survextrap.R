
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
library(purrr)
library(rslurm)

######################################################
#  Load helper functions.
######################################################

source("Functions/Survextrap_model.R")

######################################################
#  Get data
######################################################

trial_data <- readRDS("Data/trial_data.rds")
historic_trial_aggregate <- readRDS("Data/historic_trial_aggregate_no_overlap.rds")
external_data_aggregate <- readRDS("Data/external_data_no_overlap.rds")

store_models <- "/projects/aa/klvq491/Flatiron_ansclc/models/"

######################################################
#  Specify scenarios.
######################################################
##?survextrap()

add_knots1 <- c(2, 3, 5)
#add_knots2 <- c(2, 2.5, 3.5, 5)
#add_knots3 <- c(2, 2.5, 3.5, 4, 4.5, 5)

base_scenarios <- expand_grid(
  model = c("PH", "NON-PH", "Separate_arms"),
  datasets = c("trial_only", "trial_and_historic", "trial_and_all"),
  df = c(3,6,10),
  hsd_rate = c(1,5,10),
  hrsd_rate = c(1,5,10),
  add_knots = paste0("add_knots", 1),
  fit_method = "mcmc") %>%
  mutate(hrsd_rate = if_else(model == "NON-PH", hrsd_rate, NA)) %>%
  distinct() %>%
  mutate(
    store_file = paste0(store_models, "base_model_", row_number(), ".rds"),
    hazard_survival_file = paste0(store_models, "base_model_", row_number(), "_hs.rds"),
    rmst_file = paste0(store_models, "base_model_", row_number(), "_rmst.rds"))

#View(base_scenarios)
#summary(as.factor(base_scenarios$hrsd_rate))
#summary(as.factor(base_scenarios$model))

# knots_sensitivity <- expand_grid(
#   model = "PH",
#   datasets = c("trial_only", "trial_and_historic", "trial_and_all"),
#   df = c(3,6,10)) %>%
#   mutate(
#     hsd_rate = 5,
#     hrsd_rate = if_else(model == "NON-PH", 5, NA)) %>%
#   mutate(
#     store_file = paste0(store_wd, "knots_model_", row_number(), ".rds"),
#     hazard_survival_file = paste0(store_wd, "knots_model_", row_number(), "_hs.rds"),
#     rmst_file = paste0(store_wd, "knots_model_", row_number(), "_rmst.rds"))


######################################################
# Fit models, using pmap or slurm.
######################################################

# using pmap.

pmap(base_scenarios %>% 
       select(-hazard_survival_file, -rmst_file),
     fit_model)

# pmap(knots_sensitivity %>% 
#        select(-hazard_survival_file, -rmst_file),
#      fit_model)

######################################################

# or using slurm:

user <- Sys.info()["user"]
check_status <- paste0("sacct -S ", as.character(Sys.Date()-20),
                       " -u ", user ,
                       " --format=JobID,Jobname,partition,state,elapsed,ncpus -X")

objects_attach <- c("trial_data",
                    "historic_trial_aggregate",
                    "external_data_aggregate",
                    "survextrap_mem",
                    paste0("add_knots", 1:1))

package_attach <- c("dplyr", "tidyr", "readr",
                   "survextrap", "rstan", "survival",
                    "posterior")

fit_base_slurm <- slurm_apply(
  fit_model, 
  base_scenarios %>% 
    select(-hazard_survival_file, -rmst_file), 
  jobname = "fit_base",
  nodes = 20, 
  cpus_per_node = 4, 
  submit = T,
  global_objects = objects_attach,
  pkgs = package_attach,
  slurm_options = list(time='00:30:00',
                       partition='core',
                       "mem-per-cpu"= '16G'))

# fit_sensistivity_slurm <- slurm_apply(
#   fit_model, 
#   knots_sensitivity %>% 
#     select(-hazard_survival_file, -rmst_file), 
#   jobname = "fit_knots_sensitivity",
#   nodes = 4, 
#   cpus_per_node = 4, 
#   submit = T,
#   global_objects = objects_attach,
#   pkgs = package_attach,
#   slurm_options = list(time='01:00:00',
#                        partition='core',
#                        "mem-per-cpu"= '16G'))
# check results.

system(check_status)

results_0 <- readRDS("_rslurm_fit_base/results_0.RDS")
# results_0 <- readRDS("_rslurm_fit_knots_sensitivity/results_0.RDS")
# results_0

######################################################
#  Get survival, hazard and rmst.
######################################################

pmap(base_scenarios %>% 
       select(store_file, hazard_survival_file) %>%
       rename(model_file = store_file, store_file =  hazard_survival_file),
     get_survival_and_hazard_survextrap)

pmap(base_scenarios %>% 
       select(store_file, rmst_file) %>%
       rename(model_file = store_file, store_file =  rmst_file),
     get_rmst_survextrap)

sh_base_slurm <- slurm_apply(
  get_survival_and_hazard_survextrap, 
  base_scenarios %>%  
    select(store_file, hazard_survival_file) %>%
    rename(model_file = store_file, store_file =  hazard_survival_file), 
  jobname = "sh_base",
  nodes = 10, 
  cpus_per_node = 4, 
  submit = T,
  global_objects = objects_attach,
  pkgs = package_attach,
  slurm_options = list(time='00:30:00',
                       partition='core',
                       "mem-per-cpu"= '16G'))


rmst_base_slurm <- slurm_apply(
  get_rmst_survextrap, 
  base_scenarios %>% 
    select(store_file, rmst_file) %>%
    rename(model_file = store_file, store_file =  rmst_file), 
  jobname = "rmst_base",
  nodes = 10, 
  cpus_per_node = 4, 
  submit = T,
  global_objects = objects_attach,
  pkgs = package_attach,
  slurm_options = list(time='00:30:00',
                       partition='core',
                       "mem-per-cpu"= '16G'))

# pmap(knots_sensitivity %>% 
#        select(store_file, hazard_survival_file) %>%
#        rename(model_file = store_file, store_file =  hazard_survival_file),
#      get_survival_and_hazard_survextrap)
# 
# pmap(knots_sensitivity %>% 
#        select(store_file, rmst_file) %>%
#        rename(model_file = store_file, store_file =  rmst_file),
#      get_rmst_survextrap)
#
# test <- readRDS("/scratch/klvq491/case_study_nice_ta536/base_model_1_rmst.rds")
# test2 <- readRDS("/scratch/klvq491/case_study_nice_ta536/base_model_7_rmst.rds")
# 
# View(test)
# View(test2)


######################################################
#  Plot results.
######################################################

# read in results from base scenarios.

base_scenarios

for(i in 1:nrow(base_scenarios)){

  temp_hs <- readRDS(base_scenarios$hazard_survival_file[i]) %>%
    mutate(store_file = base_scenarios$store_file[i])
  
   temp_rmst <- readRDS(base_scenarios$rmst_file[i])  %>%
     mutate(store_file = base_scenarios$store_file[i])
  
  if(i == 1){
    
    base_results <- temp_hs %>% 
      bind_rows(temp_rmst) 
    
  } else {
    
    base_results <-  base_results %>%
      bind_rows(temp_hs)   %>%
      bind_rows(temp_rmst)
    
    
  }
}

base_results <- base_results %>%
  left_join(base_scenarios, by = "store_file")


saveRDS(base_results, 
        paste0(store_models, "base_model_all.rds"))


############################################################
# Survival/hazard plot figure.
############################################################

# Specify models to plot.
df_value <- 6
hsd_rate_value <- 1
hrsd_rate_value <- 10

km_data_all <- readRDS("Data/km_data_all.rds")

for(i in 1:3){
  
 # i <- 1
  
  model_type <- c("PH", "NON-PH", "Separate_arms")[i]
  
  figure_file <- c("Figures/Sup_figure_5_ph.tiff",
                   "Figures/Figure_3_nonph.tiff",
                   "Figures/Sup_figure_6_sep_arms.tiff")[i]
  
  #model_type <- "PH"

  figure3_a <- base_results  %>%
    filter(variable == "survival", add_knots == "add_knots1") %>%
    filter_model(model = model_type, df = 6, hrsd_rate = hrsd_rate_value) %>%
    filter(hsd_rate == hsd_rate_value) %>%
    filter(t > 0) %>%
    mutate(datasets = factor(datasets, 
                             levels = c("trial_only", 
                                        "trial_and_historic",
                                        "trial_and_all"
                             ),
                             labels = c("ALEX trial data only", 
                                        "ALEX + PROFILE-1014",
                                        "ALEX + PROFILE-1014 + Flatiron RWE"
                             ))) %>%  
    ggplot()+
    theme_classic()+
    theme(legend.position = "bottom",
          legend.title = element_text(margin = margin(l = unit(4.0, 'cm'), r = unit(12.0, 'cm'))),
          legend.key.spacing.x =  unit(0.3, 'cm'),
          legend.box.spacing = unit(0, "inch"),
          plot.margin = unit(c(0.1,0,0,0), "cm"),
          legend.spacing = unit(c(0,0,0,0), "cm"),
          legend.spacing.x = unit(0, "mm"),
          legend.spacing.y = unit(0, "mm"))+
    geom_ribbon(aes(x = t, y = median, ymin = lower, ymax = upper, colour = trt, fill = trt),
                alpha = 0.25, colour =  NA)+
    geom_line(aes(x = t, y = median, colour = trt))+
    geom_line(data = km_data_all %>% 
                filter(dataset == "ALEX trial") %>% 
                select(-c(dataset,model)),
              aes(x = time, y = surv, colour = trt) )+
    scale_colour_discrete("Treatment")+
    scale_fill_discrete("Treatment")+
    scale_y_continuous(limits = c(0,1),labels = scales::percent)+
    ylab("Overall Survival")+
    xlab("Time (years)")+
    facet_wrap(~datasets, ncol = 1)
  
  figure3_a
  
  figure3_b <- base_results  %>%
    filter(variable == "hazard", add_knots == "add_knots1") %>%
    filter_model(model = model_type, df = 6, hrsd_rate = hrsd_rate_value) %>%
    filter(hsd_rate == hsd_rate_value) %>%
    filter(t > 0) %>%
    mutate(datasets = factor(datasets, 
                             levels = c("trial_only", 
                                        "trial_and_historic",
                                        "trial_and_all"
                             ),
                             labels = c("ALEX trial data only", 
                                        "ALEX + PROFILE-1014",
                                        "ALEX + PROFILE-1014 + Flatiron RWE"
                             ))) %>%  
    filter(t > 0) %>%
    ggplot()+
    theme_classic()+
    theme(legend.box = "horizontal",
          legend.box.spacing = unit(0, "inch"),
          legend.title = element_text(margin = margin(b = 0.1)),
          plot.margin = unit(c(0.1,0,0,0), "cm"),
          legend.spacing = unit(c(0,0,0,0), "cm"),
          legend.spacing.x = unit(0, "mm"),
          legend.spacing.y = unit(0, "mm"))+
    geom_ribbon(aes(x = t, y = median, ymin = lower, ymax = upper, colour = trt, fill = trt),
                alpha = 0.25, colour =  NA)+
    geom_line(aes(x = t, y = median, colour = trt))+
    scale_colour_discrete("Treatment")+
    scale_fill_discrete("Treatment")+
    scale_y_continuous("Hazard", limits = c(0, 0.45))+
    xlab("Time (years)")+
    facet_wrap(~datasets, ncol = 1)
  
  figure3_b
  
  plot_legend <- get_legend_35(figure3_a)
  
  plot_figure_3_no_legend <- plot_grid(
    figure3_a +
      theme(legend.position="none",
            plot.title = element_text(hjust = -0.1,
                                      vjust = -0.02,
                                      size=14, face="bold"))+
      labs(title = "(a)"),
    NULL,
    figure3_b+
      theme(legend.position="none",
            plot.title = element_text(hjust = -0.1,
                                      vjust = -0.02,
                                      size=14, face="bold"))+
      labs(title = "(b)"), 
    rel_widths = c(0.51,0.01, 0.5),
    ncol = 3)
  
  plot_figure_3 <- plot_grid(plot_figure_3_no_legend,
                             NULL,
                             plot_legend, 
                             rel_heights = c(0.9,-0.01, 0.1),
                             ncol = 1)
  
  plot_figure_3
  
  tiff(file = figure_file,   
       width = 7.2, 
       height = 6.0,
       units = 'in',  
       res = 300, 
       compression = "lzw")
  print(plot_figure_3)
  dev.off()
  
  
}


############################################################
# Forest plot figure.
############################################################


for(i in c(1,2)){

  time_point <- c(5,20)[i]
  

  vlines_rmst <- list(seq(from = 2, to = 4, by = 0.5),
                      seq(from = 2.5, to = 12.5, by = 2.5))
  vlines_irmst <- list(seq(from = -0.5, to = 1, by = 0.5),
                       seq(from = -5, to = 2.5, by = 2.5))
  alpha_lines <- 0.4
  
  low_increment = c(0.1, 0.4)
  
  low_x_axis <- base_results %>%
    filter(variable == "rmst") %>%
    filter(t == time_point) %>%
    pull(lower) %>%
    min() - low_increment[i]
  
  high_x_axis <- base_results %>%
    filter(variable == "rmst") %>%
    filter(t == time_point) %>%
    pull(upper) %>%
    max() + 0.1
  
  control_forest <- base_results %>%
    filter(variable == "rmst") %>%
    filter(trt == "Crizotinib") %>%
    filter(t > 0) %>%
    mutate(datasets = factor(datasets, 
                             levels = c("trial_only", 
                                        "trial_and_historic",
                                        "trial_and_all"
                             ),
                             labels = c("ALEX trial data only", 
                                        "ALEX + PROFILE-1014",
                                        "ALEX + PROFILE-1014 + Flatiron RWE"))) %>%  
    mutate(model = factor(model, 
                          levels = c("PH", 
                                     "NON-PH",
                                     "Separate_arms"
                          ),
                          labels = c("PH", 
                                     "Non-PH",
                                     "Separate arms"))) %>%  
    filter(df == df_value) %>%
    filter(hrsd_rate == hrsd_rate_value | is.na(hrsd_rate)) %>%
    filter(hsd_rate == hsd_rate_value) %>%
    filter(t == time_point) %>%
    arrange(datasets) %>%
    #group_by(trt) %>%
    mutate(group_id = -row_number()) %>%
    # mutate(index = paste0("group", row_number())) %>%
    ggplot()+
    theme_classic()+
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.line.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title.y=element_blank(),
         axis.text.x = element_text( size = 8), #, angle = 45, vjust = 0.5, hjust=0.5),
   #      axis.title.x = element_text( face="bold",size = 10),
         axis.title.x = element_text( size = 10),
          legend.background = element_blank(),
          legend.box.background = element_rect(colour = "black"),
          legend.key.spacing.y = unit(3, "pt"),
          legend.text = element_text(size=8)) + 
    geom_vline(xintercept = vlines_rmst[[i]],
               colour = "gray70",
               alpha = alpha_lines)+
    geom_point(aes( x = median, y=group_id, colour = datasets, shape =  model),
               alpha = 1,
               stroke = 1,
               size = 2)+
    geom_linerange(aes(xmin = lower, xmax = upper, y = group_id, colour = datasets))+
    scale_x_continuous( paste0("RMST at ", time_point, "-years\n", "for Crizotinib"),
                        limits = c(low_x_axis, high_x_axis))+
    scale_shape_discrete("Model")+
    scale_colour_discrete("Datasets")+
    guides(                              
      shape = guide_legend(override.aes=list(colour = "gray60",
                                             fill = "gray60")))
  
  active_forest <- base_results %>%
    filter(variable == "rmst") %>%
    filter(trt == "Alectinib") %>%
    filter(t > 0) %>%
    mutate(datasets = factor(datasets, 
                             levels = c("trial_only", 
                                        "trial_and_historic",
                                        "trial_and_all"
                             ),
                             labels = c("ALEX trial data only", 
                                        "ALEX + PROFILE-1014",
                                        "ALEX + PROFILE-1014 + Flatiron RWE"))) %>%  
    mutate(model = factor(model, 
                          levels = c("PH", 
                                     "NON-PH",
                                     "Separate_arms"
                          ),
                          labels = c("PH", 
                                     "Non-PH",
                                     "Separate arms"))) %>%  
    filter(df == df_value) %>%
    filter(hrsd_rate == hrsd_rate_value | is.na(hrsd_rate)) %>%
    filter(hsd_rate == hsd_rate_value) %>%
    filter(t == time_point) %>%
    arrange(datasets) %>%
    #group_by(trt) %>%
    mutate(group_id = -row_number()) %>%
    # mutate(index = paste0("group", row_number())) %>%
    ggplot()+
    theme_classic()+
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.line.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title.y=element_blank(),
          axis.text.x = element_text( size = 8), #, angle = 45, vjust = 0.5, hjust=0.5),
      #    axis.title.x = element_text( face="bold",size = 10),
         axis.title.x = element_text(size = 10),
          legend.background = element_blank(),
          legend.box.background = element_rect(colour = "black"),
          legend.key.spacing.y = unit(3, "pt"),
          legend.text = element_text(size=8)) + 
    geom_vline(xintercept = vlines_rmst[[i]],
               colour = "gray70",
               alpha = alpha_lines)+
    geom_point(aes( x = median, y=group_id, colour = datasets, shape =  model),
               alpha = 1,
               stroke = 1,
               size = 2)+
    geom_linerange(aes(xmin=  lower, xmax = upper, y = group_id, colour = datasets))+
    scale_x_continuous( paste0("RMST at ", time_point, "-years\n", "for Alectinib"),
                        limits = c(low_x_axis, high_x_axis))+
    scale_shape_discrete("Model")+
    scale_colour_discrete("Datasets")+
    guides(                              
      shape = guide_legend(override.aes=list(colour = "gray60",
                                             fill = "gray60")))
  
  difference_forest <- base_results %>%
    filter(variable == "irmst") %>%
    filter(t > 0) %>%
    mutate(datasets = factor(datasets, 
                             levels = c("trial_only", 
                                        "trial_and_historic",
                                        "trial_and_all"
                             ),
                             labels = c("ALEX trial data only", 
                                        "ALEX + PROFILE-1014",
                                        "ALEX + PROFILE-1014 + Flatiron RWE"))) %>%  
    mutate(model = factor(model, 
                          levels = c("PH", 
                                     "NON-PH",
                                     "Separate_arms"
                          ),
                          labels = c("PH", 
                                     "Non-PH",
                                     "Separate arms"))) %>%  
    filter(df == df_value) %>%
    filter(hrsd_rate == hrsd_rate_value | is.na(hrsd_rate)) %>%
    filter(hsd_rate == hsd_rate_value) %>%
    filter(t == time_point) %>%
    arrange(datasets) %>%
    #group_by(trt) %>%
    mutate(group_id = -row_number()) %>%
    # mutate(index = paste0("group", row_number())) %>%
    ggplot()+
    theme_classic()+
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.line.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title.y=element_blank(),
          axis.text.x = element_text( size = 8), #, angle = 45, vjust = 0.5, hjust=0.5),
     #     axis.title.x = element_text( face="bold",size = 10),
          axis.title.x = element_text(size = 10),
          legend.background = element_blank(),
          legend.box.background = element_rect(colour = "black"),
          legend.key.spacing.y = unit(3, "pt"),
          legend.text = element_text(size=8)) + 
    geom_vline(xintercept = vlines_irmst[[i]],
               colour = "gray70",
               alpha = alpha_lines)+
    geom_vline(xintercept = 0,
               colour = "gray20",
               alpha = alpha_lines)+
    geom_point(aes( x = median, y=group_id, colour = datasets, shape =  model),
               alpha = 1,
               stroke = 1,
               size = 2)+
    geom_linerange(aes(xmin=  lower, xmax = upper, y = group_id, colour = datasets))+
    scale_x_continuous( paste0("Difference in \n RMST at ", time_point, "-years\n") )+
    scale_shape_discrete("Model")+
    scale_colour_discrete("Datasets")+
    guides(                              
      shape = guide_legend(override.aes=list(colour = "gray60",
                                             fill = "gray60")))
  
  
  forest_legend <- cowplot::get_legend(
    # create some space to the left of the legend
    active_forest + theme(legend.box.margin = margin(0, 12, 0, 0))
  )

  forest_plots <- plot_grid(control_forest+
                              theme(legend.position="none",
                                    plot.title = element_text(size=10, face="bold")),#+
                          #    labs(title = plot_label),
                            active_forest+
                              theme(legend.position="none",
                                    plot.title = element_text(size=10, face="bold")),#+
                       #       labs(title = " "),
                            difference_forest+
                              theme(legend.position="none",
                                    plot.title = element_text(size=10, face="bold")),#+
                            #  labs(title = " "),
                            align = "h",
                            rel_widths=c(1,1,1),
                            nrow = 1)
  print(i)
  
  assign(paste0("forest",i), forest_plots)
}

# forest
# forest1
# forest2
# 
# plot_grid(forest_plots,
#           forest_legend,
#           rel_widths=c(1, 0.6),
#           axis = "l"
# )

plot_all <- ggdraw(
  plot_grid(plot_grid(NULL,forest1, NULL, forest2, rel_heights = c(0.06, 0.5, 0.06, 0.5),  
                      labels = c("(a)", "", "(b)", ""),
                      label_size = 12,
                      label_x = 0, ncol = 1,align = "v"),
            plot_grid(NULL, plot_grid(forest_legend , NULL, ncol = 1, rel_heights = c(1, 100)), ncol=1),
            rel_widths=c(1, 0.65)))

plot_all

tiff(file = "Figures/Figure_4.tiff",   
     width = 7.5, 
     height = 5.0,
     units = 'in',  
     res = 300, 
     compression = "lzw")
print(plot_all)
dev.off()


############################################################
# Sensitivity analysis for Non-PH Model, on Sigma prior.
############################################################


sensitivity_figure1_a <- base_results  %>%
  filter(variable == "survival") %>%
  filter(model == "NON-PH") %>%
  filter(df == 6) %>%
  filter(hrsd_rate == 10) %>% 
  filter(datasets == "trial_and_all") %>%
  filter(t > 0) %>%
  mutate(prior = paste0("'\u03c3~'*`G`*`amma`*`(`*", 
                        2, "*`,`*" , 
                        hsd_rate, "*`)`")) %>%
  mutate(prior = factor(hsd_rate, labels = prior, levels = hsd_rate)) %>%
  ggplot()+
  theme_classic()+
  theme(legend.position = "bottom",
        legend.title = element_text(margin = margin(l = unit(4.0, 'cm'), r = unit(12.0, 'cm'))),
        legend.key.spacing.x =  unit(0.3, 'cm'),
        legend.box.spacing = unit(0, "inch"),
        plot.margin = unit(c(0.1,0,0,0), "cm"),
        legend.spacing = unit(c(0,0,0,0), "cm"),
        legend.spacing.x = unit(0, "mm"),
        legend.spacing.y = unit(0, "mm"))+
  geom_ribbon(aes(x = t, y = median, ymin = lower, ymax = upper, colour = trt, fill = trt),
              alpha = 0.25, colour =  NA)+
  geom_line(aes(x = t, y = median, colour = trt))+
  geom_line(data = km_data_all %>% 
              filter(dataset == "ALEX trial") %>% 
              select(-c(dataset,model)),
            aes(x = time, y = surv, colour = trt) )+
  scale_colour_discrete("Treatment")+
  scale_fill_discrete("Treatment")+
  scale_y_continuous(limits = c(0,1),labels = scales::percent)+
  ylab("Overall Survival")+
  xlab("Time (years)")+
  facet_wrap(~prior, ncol = 1, labeller = label_parsed)

sensitivity_figure1_a

sensitivity_figure1_b <- base_results  %>%
  filter(variable == "hazard") %>%
  filter(model == "NON-PH") %>%
  filter(df == 6) %>%
  filter(hrsd_rate == 10) %>%
  filter(datasets == "trial_and_all") %>%
  filter(t > 0) %>%
  mutate(prior = paste0("'\u03c3~'*`G`*`amma`*`(`*", 
                        2, "*`,`*" , 
                        hsd_rate, "*`)`")) %>%
  mutate(prior = factor(hsd_rate, labels = prior, levels = hsd_rate)) %>%
  ggplot()+
  theme_classic()+
  theme(legend.box = "horizontal",
        legend.box.spacing = unit(0, "inch"),
        legend.title = element_text(margin = margin(b = 0.1)),
        plot.margin = unit(c(0.1,0,0,0), "cm"),
        legend.spacing = unit(c(0,0,0,0), "cm"),
        legend.spacing.x = unit(0, "mm"),
        legend.spacing.y = unit(0, "mm"))+
  geom_ribbon(aes(x = t, y = median, ymin = lower, ymax = upper, colour = trt, fill = trt),
              alpha = 0.25, colour =  NA)+
  geom_line(aes(x = t, y = median, colour = trt))+
  scale_colour_discrete("Treatment")+
  scale_fill_discrete("Treatment")+
  scale_y_continuous("Hazard", limits = c(0, 0.45))+
  xlab("Time (years)")+
  facet_wrap(~prior, ncol = 1, labeller = label_parsed)

sensitivity_figure1_b

plot_legend <- get_legend_35(sensitivity_figure1_a)

sensitivity_figure1_no_legend <- plot_grid(
  sensitivity_figure1_a +
    theme(legend.position="none",
          plot.title = element_text(hjust = -0.1,
                                    vjust = -0.02,
                                    size=14, face="bold"))+
    labs(title = "(a)"),
  NULL,
  sensitivity_figure1_b+
    theme(legend.position="none",
          plot.title = element_text(hjust = -0.1,
                                    vjust = -0.02,
                                    size=14, face="bold"))+
    labs(title = "(b)"), 
  rel_widths = c(0.51,0.01, 0.5),
  ncol = 3)

sensitivity_figure1 <- plot_grid(sensitivity_figure1_no_legend,
                                 NULL,
                                 plot_legend, 
                                 rel_heights = c(0.9,-0.01, 0.1),
                                 ncol = 1)

sensitivity_figure1

tiff(file = "Figures/Sup_figure_7_sigma.tiff",   
     width = 7.2, 
     height = 6.0,
     units = 'in',  
     res = 300, 
     compression = "lzw")
print(sensitivity_figure1)
dev.off()


############################################################
# Sensitivity analysis for Non-PH Model, on Tau prior.
############################################################


sensitivity_figure2_a <- base_results  %>%
  filter(variable == "survival") %>%
  filter(model == "NON-PH") %>%
  filter(df == 6) %>%
  filter(hsd_rate == 1) %>%
  filter(datasets == "trial_and_all") %>%
  filter(t > 0) %>%
  mutate(prior = paste0("'\u03C4~'*`G`*`amma`*`(`*", 
                        2, "*`,`*" , 
                        hrsd_rate, "*`)`")) %>%
  mutate(prior = factor(hrsd_rate, labels = prior, levels = hrsd_rate)) %>%
  ggplot()+
  theme_classic()+
  theme(legend.position = "bottom",
        legend.title = element_text(margin = margin(l = unit(4.0, 'cm'), r = unit(12.0, 'cm'))),
        legend.key.spacing.x =  unit(0.3, 'cm'),
        legend.box.spacing = unit(0, "inch"),
        plot.margin = unit(c(0.1,0,0,0), "cm"),
        legend.spacing = unit(c(0,0,0,0), "cm"),
        legend.spacing.x = unit(0, "mm"),
        legend.spacing.y = unit(0, "mm"))+
  geom_ribbon(aes(x = t, y = median, ymin = lower, ymax = upper, colour = trt, fill = trt),
              alpha = 0.25, colour =  NA)+
  geom_line(aes(x = t, y = median, colour = trt))+
  geom_line(data = km_data_all %>% 
              filter(dataset == "ALEX trial") %>% 
              select(-c(dataset,model)),
            aes(x = time, y = surv, colour = trt) )+
  scale_colour_discrete("Treatment")+
  scale_fill_discrete("Treatment")+
  scale_y_continuous(limits = c(0,1),labels = scales::percent)+
  ylab("Overall Survival")+
  xlab("Time (years)")+
  facet_wrap(~prior, ncol = 1, labeller = label_parsed)

sensitivity_figure2_a

sensitivity_figure2_b <- base_results  %>%
  filter(variable == "hazard") %>%
  filter(model == "NON-PH") %>%
  filter(df == 6) %>%
  filter(hsd_rate == 1) %>%
  filter(datasets == "trial_and_all") %>%
  filter(t > 0) %>%
  mutate(prior = paste0("'\u03C4~'*`G`*`amma`*`(`*", 
                        2, "*`,`*" , 
                        hrsd_rate, "*`)`")) %>%
  mutate(prior = factor(hrsd_rate, labels = prior, levels = hrsd_rate)) %>%
  ggplot()+
  theme_classic()+
  theme(legend.box = "horizontal",
        legend.box.spacing = unit(0, "inch"),
        legend.title = element_text(margin = margin(b = 0.1)),
        plot.margin = unit(c(0.1,0,0,0), "cm"),
        legend.spacing = unit(c(0,0,0,0), "cm"),
        legend.spacing.x = unit(0, "mm"),
        legend.spacing.y = unit(0, "mm"))+
  geom_ribbon(aes(x = t, y = median, ymin = lower, ymax = upper, colour = trt, fill = trt),
              alpha = 0.25, colour =  NA)+
  geom_line(aes(x = t, y = median, colour = trt))+
  scale_colour_discrete("Treatment")+
  scale_fill_discrete("Treatment")+
  scale_y_continuous("Hazard", limits = c(0, 0.45))+
  xlab("Time (years)")+
  facet_wrap(~prior, ncol = 1, labeller = label_parsed)

sensitivity_figure2_b

plot_legend <- get_legend_35(sensitivity_figure2_a)

sensitivity_figure2_no_legend <- plot_grid(
  sensitivity_figure2_a +
    theme(legend.position="none",
          plot.title = element_text(hjust = -0.1,
                                    vjust = -0.02,
                                    size=14, face="bold"))+
    labs(title = "(a)"),
  NULL,
  sensitivity_figure2_b+
    theme(legend.position="none",
          plot.title = element_text(hjust = -0.1,
                                    vjust = -0.02,
                                    size=14, face="bold"))+
    labs(title = "(b)"), 
  rel_widths = c(0.51,0.01, 0.5),
  ncol = 3)

sensitivity_figure2 <- plot_grid(sensitivity_figure2_no_legend,
                           NULL,
                           plot_legend, 
                           rel_heights = c(0.9,-0.01, 0.1),
                           ncol = 1)

sensitivity_figure2

tiff(file = "Figures/Sup_figure_8_tau.tiff",   
     width = 7.2, 
     height = 6.0,
     units = 'in',  
     res = 300, 
     compression = "lzw")
print(sensitivity_figure2)
dev.off()



############################################################
# Sensitivity analysis on knot positions.
############################################################


sensitivity_figure3_a <-   base_results  %>%
  filter(variable == "survival") %>%
  filter(model == "NON-PH") %>%
  filter(datasets == "trial_and_all") %>%
  filter(hrsd_rate == 10) %>%
  filter(hsd_rate == 1) %>%
  filter(t > 0) %>%
  mutate(df_model = paste0("df = ", df)) %>%
  mutate(df_model = factor(df, levels = df, labels = df_model)) %>%
  ggplot()+
  theme_classic()+
  theme(legend.position = "bottom",
        legend.title = element_text(margin = margin(l = unit(4.0, 'cm'), r = unit(12.0, 'cm'))),
        legend.key.spacing.x =  unit(0.3, 'cm'),
        legend.box.spacing = unit(0, "inch"),
        plot.margin = unit(c(0.1,0,0,0), "cm"),
        legend.spacing = unit(c(0,0,0,0), "cm"),
        legend.spacing.x = unit(0, "mm"),
        legend.spacing.y = unit(0, "mm"))+
  geom_ribbon(aes(x = t, y = median, ymin = lower, ymax = upper, colour = trt, fill = trt),
              alpha = 0.25, colour =  NA)+
  geom_line(aes(x = t, y = median, colour = trt))+
  geom_line(data = km_data_all %>% 
              filter(dataset == "ALEX trial") %>% 
              select(-c(dataset,model)),
            aes(x = time, y = surv, colour = trt) )+
  scale_colour_discrete("Treatment")+
  scale_fill_discrete("Treatment")+
  scale_y_continuous(limits = c(0,1),labels = scales::percent)+
  ylab("Overall Survival")+
  xlab("Time (years)")+
  facet_wrap(~df_model, ncol = 1)

sensitivity_figure3_a

sensitivity_figure3_b <- base_results  %>%
  filter(variable == "hazard") %>%
  filter(model == "NON-PH") %>%
  filter(datasets == "trial_and_all") %>%
  filter(hrsd_rate == 10) %>%
  filter(hsd_rate == 1) %>%
  filter(t > 0) %>%
  mutate(df_model = paste0("df = ", df)) %>%
  mutate(df_model = factor(df, levels = df, labels = df_model)) %>%
  ggplot()+
  theme_classic()+
  theme(legend.box = "horizontal",
        legend.box.spacing = unit(0, "inch"),
        legend.title = element_text(margin = margin(b = 0.1)),
        plot.margin = unit(c(0.1,0,0,0), "cm"),
        legend.spacing = unit(c(0,0,0,0), "cm"),
        legend.spacing.x = unit(0, "mm"),
        legend.spacing.y = unit(0, "mm"))+
  geom_ribbon(aes(x = t, y = median, ymin = lower, ymax = upper, colour = trt, fill = trt),
              alpha = 0.25, colour =  NA)+
  geom_line(aes(x = t, y = median, colour = trt))+
  scale_colour_discrete("Treatment")+
  scale_fill_discrete("Treatment")+
  scale_y_continuous("Hazard", c(0, 0.45))+
  xlab("Time (years)")+
  facet_wrap(~df_model, ncol = 1)

sensitivity_figure3_b

plot_legend <- get_legend_35(sensitivity_figure3_a)

sensitivity_figure3_no_legend <- plot_grid(
  sensitivity_figure3_a +
    theme(legend.position="none",
          plot.title = element_text(hjust = -0.1,
                                    vjust = -0.02,
                                    size=14, face="bold"))+
    labs(title = "(a)"),
  NULL,
  sensitivity_figure3_b+
    theme(legend.position="none",
          plot.title = element_text(hjust = -0.1,
                                    vjust = -0.02,
                                    size=14, face="bold"))+
    labs(title = "(b)"), 
  rel_widths = c(0.51,0.01, 0.5),
  ncol = 3)

sensitivity_figure3 <- plot_grid(sensitivity_figure3_no_legend,
                                 NULL,
                                 plot_legend, 
                                 rel_heights = c(0.9,-0.01, 0.1),
                                 ncol = 1)

sensitivity_figure3

tiff(file = "Figures/Sup_figure_9_df.tiff",   
     width = 7.2, 
     height = 6.0,
     units = 'in',  
     res = 300, 
     compression = "lzw")
print(sensitivity_figure3)
dev.off()




