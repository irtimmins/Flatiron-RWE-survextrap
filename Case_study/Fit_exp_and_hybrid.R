
library(dplyr)
library(ggplot2)
library(survival)
library(flexsurv)

####################################################
# Hybrid model.
####################################################

# load helper functions.
source("Functions/Hybrid_model.R")
source("Functions/Survextrap_model.R")


trial_data <- readRDS("Data/trial_IPD_OS_ALEX_Feb_2017_v4.rds")
km_data_all <- readRDS("Data/km_data_all.rds")

control <- trial_data %>%
  filter(trt == "Crizotinib") %>%
  select(time, status)

active <- trial_data %>%
  filter(trt == "Alectinib") %>%
  select(time, status)

set.seed(1212)

hybrid_control <- fit_hybrid_model(data = control, cut_point = 1.5)
hybrid_control_survival <- hybrid_model_survival(hybrid_control, 
                                                 t = seq(from = 0, to = 20, length.out = 1e3))

hybrid_active <- fit_hybrid_model(data = active, cut_point = 1.5)
hybrid_active_survival <- hybrid_model_survival(hybrid_active , 
                                                t = seq(from = 0, to = 20, length.out = 1e3))
sum(active$status[active$time > 1.5])
sum(control$status[control$time > 1.5])

hybrid_active_survival <- hybrid_active_survival %>%
  mutate(trt = "Alectinib") %>%
  select(trt, time, value, se, lower, upper)

hybrid_control_survival <- hybrid_control_survival %>%
  mutate(trt = "Crizotinib")  %>%
  select(trt, time, value, se, lower, upper)
#
hybrid_survival_plot <- 
  hybrid_active_survival %>%
  bind_rows(hybrid_control_survival) %>%
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
  geom_line(data = km_data_all %>% 
              filter(dataset == "ALEX trial") %>% 
              select(-c(dataset,model)),
            aes(x = time, y = surv, colour = trt),
            alpha = 0.5)+
  geom_ribbon(aes(x = time, ymin = lower, ymax = upper, fill = trt),
              alpha= 0.25, colour = NA)+
  geom_vline(xintercept = 1.5, alpha = 0.4)+
  geom_step(aes(x = time, y = value, colour = trt))+
  scale_colour_discrete("Treatment")+
  scale_fill_discrete("Treatment")+
  scale_y_continuous(limits = c(0,1.01),labels = scales::percent)+
  ylab("Overall Survival")+
  xlab("Time (years)")

hybrid_survival_plot
#max(hybrid_active_survival$upper)
#max(hybrid_control_survival$upper)




###################################################
# RMST results.
####################################################

exp_control <-  flexsurvreg(Surv(time, status) ~ 1, data = control,  dist = "exp")
exp_active <- flexsurvreg(Surv(time, status) ~ 1, data = active,  dist = "exp")
head(hybrid_control_survival)

exp_control_survival <- summary(exp_control, 
                                t = seq(from = 0, to = 20, length.out = 1e3),
                                se = T)[[1]] %>% 
  as_tibble() %>%
  mutate(variable = "Survival", trt = "Crizotinib") %>%
  rename(value = est, lower = lcl, upper = ucl) %>%
  select(variable, trt, time, value, se, lower, upper)

exp_active_survival <- summary(exp_active, 
                                t = seq(from = 0, to = 20, length.out = 1e3),
                                se = T)[[1]] %>% 
  as_tibble() %>%
  mutate(variable = "Survival", trt = "Alectinib") %>%
  rename(value = est, lower = lcl, upper = ucl) %>%
  select(variable, trt, time, value, se, lower, upper)

exp_survival_plot <- exp_active_survival %>%
  bind_rows(exp_control_survival) %>%
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
  geom_line(data = km_data_all %>% 
              filter(dataset == "ALEX trial") %>% 
              select(-c(dataset,model)),
            aes(x = time, y = surv, colour = trt),
            alpha = 0.5)+
  geom_ribbon(aes(x = time, ymin = lower, ymax = upper, fill = trt),
              alpha= 0.25, colour = NA)+
  geom_vline(xintercept = 1.5, alpha = 0.4)+
  geom_step(aes(x = time, y = value, colour = trt))+
  scale_colour_discrete("Treatment")+
  scale_fill_discrete("Treatment")+
  scale_y_continuous(limits = c(0,1.01),labels = scales::percent)+
  ylab("Overall Survival")+
  xlab("Time (years)")

exp_survival_plot


plot_legend <- get_legend_35(exp_survival_plot)

plot_figure_no_legend <- plot_grid(
  exp_survival_plot +
    theme(legend.position="none",
          plot.title = element_text(hjust = -0.1,
                                    vjust = -0.02,
                                    size=14, face="bold"))+
    labs(title = "(a)"),
  NULL,
   hybrid_survival_plot+
    theme(legend.position="none",
          plot.title = element_text(hjust = -0.1,
                                    vjust = -0.02,
                                    size=14, face="bold"))+
    labs(title = "(b)"), 
  rel_widths = c(0.51,0.01, 0.5),
  ncol = 3)

plot_figure <- plot_grid(plot_figure_no_legend,
                           NULL,
                           plot_legend, 
                           rel_heights = c(0.9,00.005, 0.1),
                           ncol = 1)

plot_figure 

tiff(file = "Figures/Sup_figure_4_exp_hybrid.tiff",   
     width = 7.2, 
     height = 3.7,
     units = 'in',  
     res = 300, 
     compression = "lzw")
print(plot_figure)
dev.off()


#######################################

set.seed(1212)

#test
hybrid_control_rmst <- hybrid_model_rmst(hybrid_control, t = c(5, 10, 15, 20))  %>%
  mutate(trt = "Crizotinib") %>%
  mutate(variable = "rmst")

hybrid_active_rmst <- hybrid_model_rmst(hybrid_active, t = c(5, 10, 15, 20)) %>%
  mutate(trt = "Alectinib") %>%
  mutate(variable = "rmst")

rmst_hybrid <- hybrid_active_rmst %>%
  bind_rows(hybrid_control_rmst) %>%
  relocate(c("variable", "trt"), .before =  "time")


irmst_hybrid <- hybrid_control_rmst %>%
  bind_rows(hybrid_active_rmst) %>%
  select(-c(lower, upper)) %>%
  pivot_wider(values_from = c("value","se"), names_from = "trt" ) %>%
  mutate(value = value_Alectinib - value_Crizotinib,
         se = sqrt(se_Alectinib^2+se_Crizotinib^2)) %>%
  select(time, value, se) %>%
  mutate(lower = value - 1.96*se,
         upper = value + 1.96*se) %>%
  mutate(trt = NA) %>%
  mutate(variable = "irmst") %>%
  relocate(c(variable, trt), .before = time)
  
hybrid_all <- rmst_hybrid %>%
  bind_rows(irmst_hybrid)
  
hybrid_all

saveRDS(hybrid_all, "/projects/aa/klvq491/Flatiron_ansclc/models/hybrid_model_rmst.rds")

# test <- readRDS("/projects/aa/klvq491/Flatiron_ansclc/models/base_model_1_rmst.rds")
# test

test2 <- readRDS("/projects/aa/klvq491/Flatiron_ansclc/models/base_model_1_hs.rds")
test2


exp_control_rmst <- summary(exp_control, 
                            t = c(5, 10, 15, 20),
                            type = "rmst",
                            se = T)[[1]] %>% 
  as_tibble() %>%
  mutate(variable = "rmst", trt = "Crizotinib") %>%
  rename(value = est, lower = lcl, upper = ucl) %>%
  select(variable, trt, time, value, se, lower, upper)

exp_active_rmst <- summary(exp_active, 
                               t = c(5, 10, 15, 20),
                               type = "rmst",
                               se = T)[[1]] %>% 
  as_tibble() %>%
  mutate(variable = "rmst", trt = "Alectinib") %>%
  rename(value = est, lower = lcl, upper = ucl) %>%
  select(variable, trt, time, value, se, lower, upper)

exp_irmst <- exp_control_rmst %>%
    bind_rows(exp_active_rmst) %>%
    select(-c(lower, upper)) %>%
    pivot_wider(values_from = c("value","se"), names_from = "trt" ) %>%
  mutate(value = value_Alectinib - value_Crizotinib,
         se = sqrt(se_Alectinib^2+se_Crizotinib^2)) %>%
  select(time, value, se) %>%
  mutate(lower = value - 1.96*se,
         upper = value + 1.96*se) %>%
  mutate(trt = NA) %>%
  mutate(variable = "irmst") %>%
  relocate(c(variable, trt), .before = time)

exp_all <- exp_active_rmst %>%
  bind_rows(exp_control_rmst) %>%
  bind_rows(exp_irmst)

saveRDS(exp_all, "/projects/aa/klvq491/Flatiron_ansclc/models/exp_model_rmst.rds")

# exp_all
# hybrid_all




