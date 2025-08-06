
library(ggsurvfit)
library(patchwork)

new_surv_data0 <- trial_data  %>%
  mutate(trt = factor(trt, 
                          levels = c("Alectinib",
                                     "Crizotinib"),
                          labels = c("Alectinib         ", 
                                     "Crizotinib         "))) 

new_surv_data1 <- trial_data %>%
  filter(trt == "Crizotinib") %>%
  mutate(weight = 1) %>%
  mutate(trt = "Crizotinib",
         dataset = "ALEX")   

new_surv_data2 <- historic_trial_data %>%
  mutate(weight = 1) %>%
  mutate(trt = "Crizotinib",
         dataset = "PROFILE-1014")  

new_surv_data3 <- rwd_data_weighted %>%
  mutate(trt = "Crizotinib",
         dataset = "Flatiron RWE MAIC weighted") %>%
  select(time,status,weight,trt,dataset)  

new_surv_data4 <- rwd_data_unweighted %>%
  mutate(trt = "Crizotinib",
         dataset = "Flatiron RWE unweighted")  %>%
  select(time,status,weight,trt,dataset)  

new_surv_data_all <- bind_rows(new_surv_data1,
                               new_surv_data2,
                               new_surv_data3,
                               new_surv_data4) %>%
  mutate(dataset = factor(dataset, 
                          levels = c("ALEX", 
                                     "PROFILE-1014", 
                                     "Flatiron RWE unweighted",
                                     "Flatiron RWE MAIC weighted"),
                          labels = c("ALEX,\nCrizotinib", 
                                     "PROFILE-1014", 
                                     "Flatiron RWE\nunweighted",
                                     "Flatiron RWE\nMAIC weighted"))) #%>%

plot_2a <- 
  survfit2(Surv(time, status) ~ trt, data = new_surv_data0) %>%
  ggsurvfit(linewidth = 1) +
  add_confidence_interval() +
  add_risktable(
    risktable_height = 0.18,
    risktable_stats = c("{round(n.risk, digits = 1)}"),
    stats_label = c("Number at risk"),
    theme = list(
      theme_risktable_default(), 
      theme(
        axis.text.y = element_text(hjust = 0, margin = margin(r = -5, l = -5)),
        plot.caption = element_text(hjust = -5),
        plot.title.position = "plot"
      )
    )
  ) + 
  theme(axis.text.x = element_text(hjust = 1))+
  scale_ggsurvfit()+
  theme_classic()+
  scale_x_continuous(breaks = seq(from = 0, to = 2.5, by = 0.25))+
  theme_classic()+
   theme(legend.position.inside = c(0.12,0.6),
          legend.position = "inside",
          legend.text=element_text(size=8),
         legend.key.spacing.y = unit(0, "pt"),
          legend.box.spacing = unit(-6, "pt"),
         axis.title.y = element_text(vjust=6, hjust = 0.5)) +
  scale_colour_manual(NULL, values = c("#4393c3","#d6604d"))+
  scale_fill_manual(NULL, values = c("#4393c3","#d6604d"))+
  scale_y_continuous("Overall Survival", 
                     limits = c(0,1), labels = scales::percent)+
  scale_x_continuous("Time (years)", breaks = seq(from = 0, to = 3.5, by = 0.25), limits = c(0,3))+
  theme(  plot.title = element_text(hjust = -0.2,
                                    size=14, face="bold"))+
  labs(title = "(a)")

plot_2a_build <-  ggsurvfit_build(plot_2a)

plot_2a_build

plot_2b <- survfit2(Surv(time, status) ~ dataset, data=new_surv_data_all, weights = new_surv_data_all$weight) %>%
  ggsurvfit(linewidth = 1) +
  add_confidence_interval() +
  # add_risktable(risktable_stats = c("{round(n.risk, digits = 1)}"),
  #                                   stats_label = c("Number at risk"),
  #               risktable_height = 0.3
  #             ) +
  add_risktable(
    risktable_height = 0.34,
    risktable_stats = c("{round(n.risk, digits = 1)}"),
    stats_label = c("Number at risk"),
     theme = list(
       theme_risktable_default(), 
       theme(
         axis.text.y = element_text(hjust = 0, margin = margin(r = -5, l = -5)),
         plot.caption = element_text(hjust = -5),
         plot.title.position = "plot"
       )
     )
  ) + 
  scale_ggsurvfit()+
  theme_classic()+
  theme(legend.position.inside = c(0.115,0.35),
    legend.position = "inside",
    legend.text=element_text(size=8),
    legend.spacing.y = unit(0, "pt"),
    legend.spacing.x = unit(0, "pt"),
    legend.key.spacing.y = unit(1.2, "pt"),
    legend.box.spacing = unit(-6, "pt"),
    axis.title.y = element_text(vjust=6, hjust = 0.5))+
  scale_x_continuous("Time (years)", breaks = seq(from = 0, to = 5.5, by = 0.5))+
  scale_y_continuous("Overall Survival", 
                     limits = c(0,1), labels = scales::percent)+
  scale_colour_manual(NULL, 
                      values = c("#d6604d", "#7CAE00","#C77CFF", "deepskyblue3" ))+
  scale_fill_manual(NULL, 
                    values = c("#d6604d", "#7CAE00", "#C77CFF","deepskyblue3" ))+
  theme(  plot.title = element_text(hjust = -0.2,
                                    size=14, face="bold"))+
  labs(title = "(b)")

plot_2b_build <-  ggsurvfit_build(plot_2b)

figure_2 <- wrap_plots(plot_2a_build, plot_2b_build, ncol = 1, heights = c(1,1.18))

figure_2

tiff(file = "Figures/Figure_2.tiff",   
     width = 7.5, 
     height = 9.0,
     units = 'in',  
     res = 300, 
     compression = "lzw")
print(figure_2)
dev.off()

