
######################################################
# Log-cumulative hazard
######################################################


plot_sup_3a <-km_all %>%
  mutate(dataset = factor(dataset, 
                          levels = c("ALEX", "PROFILE-1014", 
                                     "Flatiron RWE MAIC Weighted", "Flatiron RWE Unweighted"))) %>%
  filter(dataset == "ALEX") %>%
  ggplot(aes(x = time, y = log(-log(surv)), colour = trt))+
  theme_classic()+
  theme(legend.position.inside = c(0.6,0.3),
        legend.position = "inside",
        legend.text=element_text(size=9),
        legend.key.spacing.y = unit(-8, "pt"),
        legend.box.spacing = unit(-10, "pt"),
        legend.title = element_text(size=12,
                                    margin = margin(l = 0, r = 0, b = 0, t = 0))) +
  geom_step(linewidth = 0.8)+
  scale_colour_manual("Treatment", values = c("#4393c3","#d6604d"))+
  scale_y_continuous("Log Cumulative Hazard")+
  scale_x_continuous("Time (years)",breaks = 0:5, limits = c(0,3))


plot_sup_3b <- km_all %>%
  mutate(dataset = factor(dataset, 
                          levels = c("ALEX", 
                                     "PROFILE-1014", 
                                     "Flatiron RWE MAIC Weighted", 
                                     "Flatiron RWE Unweighted"),
                          labels = c("ALEX", 
                                     "PROFILE-1014", 
                                     "Flatiron RWE MAIC weighted", 
                                     "Flatiron RWE unweighted"))) %>%
  filter(trt == "Crizotinib") %>%
  ggplot(aes(x = time, y = log(-log(surv)), colour = dataset))+
  theme_classic()+
  theme(legend.position.inside = c(0.6,0.3),
        legend.position = "inside",
        legend.text=element_text(size=9),
        legend.spacing.y = unit(-6, "pt"),
        legend.spacing.x = unit(0, "pt"),
        legend.key.spacing.y = unit(-8, "pt"),
        legend.box.spacing = unit(-10, "pt"),
        legend.title = element_text(size=12,
                                    margin = margin(l = 0, r = 0, b = 0, t = 0))) +
  geom_step(linewidth = 0.8)+
  scale_colour_manual("Dataset", 
                      values = c("#d6604d", "#7CAE00", "#00BFC4", "#C77CFF" ))+
  scale_y_continuous("Log Cumulative Hazard")+
  scale_x_continuous("Time (years)",  breaks = 0:5, limits = c(0,5.5))

plot_sup_3b

figure_sup_3 <- plot_grid(plot_sup_3a+
                            theme(  plot.title = element_text(hjust = -0.1,
                                                              size=14, face="bold"))+
                            labs(title = "(a)"),
                          plot_sup_3b+
                            theme(  plot.title = element_text(hjust = -0.1,
                                                              size=14, face="bold"))+
                            labs(title = "(b)"),
                          align = "h",
                          rel_widths  =c(0.5,0.5),
                          ncol = 2)
figure_sup_3


tiff(file = "Figures/Sup_figure_3.tiff",   
     width = 7.6, 
     height = 3.7,
     units = 'in',  
     res = 300, 
     compression = "lzw")
print(figure_sup_3)
dev.off()



######################################################
# Analyse empirical hazards for trial and RWD
######################################################

trial_smooth_bs_control <- bshazard(Surv(time, status) ~ 1, 
                                    data = trial_data  %>% filter(trt == "Crizotinib"))

trial_smooth_control <- trial_smooth_bs_control[c("time", "hazard")] %>%
  as_tibble() %>%
  mutate(data = "ALEX", trt = "Crizotinib")

trial_smooth_bs_active <- bshazard(Surv(time, status) ~ 1, 
                                   data = trial_data  %>% filter(trt == "Alectinib"))

trial_smooth_active <- trial_smooth_bs_active[c("time", "hazard")] %>%
  as_tibble() %>%
  mutate(data = "ALEX", trt = "Alectinib")

historic_trial_bs <- bshazard(Surv(time, status) ~ 1, 
                              data = historic_trial_data)
historic_smooth <- historic_trial_bs[c("time", "hazard")] %>%
  as_tibble() %>%
  mutate(data = "PROFILE-1014", trt = "Crizotinib")

rwd_bs <- bshazard(Surv(time, status) ~ 1, 
                   data = rwd_data)

rwd_smooth <- rwd_bs[c("time", "hazard")] %>%
  as_tibble() %>%
  mutate(data = "Flatiron RWE", trt = "Crizotinib")

sup_figure3 <- bind_rows(trial_smooth_control,
                         trial_smooth_active,
                         historic_smooth,
                         rwd_smooth) %>%
  ggplot(aes(x = time, y = hazard, colour = data, linetype = trt))+
  theme_classic()+
  scale_colour_discrete("Dataset")+
  scale_linetype_discrete("Treatment")+
  geom_line()+
  geom_vline(xintercept = max(trial_data$time)) +
  ylim(c(0,0.5))+
  ylab("Hazard")+
  xlab("Time (years)")#+
#ggtitle("Empirical hazards")

tiff(file = "Figures/Sup_figure_3.tiff",   
     width = 5.8, 
     height = 3.4,
     units = 'in',  
     res = 300, 
     compression = "lzw")
print(sup_figure3)
dev.off()

