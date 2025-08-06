######################################################
# MAIC weights
######################################################

sup_figure1 <- rwd_data_weighted %>%
  ggplot()+
  geom_histogram(aes(x = weight), colour = "gray50", binwidth = 0.01)+
  #  geom_histogram(aes(x = weight, colour = race_asian, fill = race_asian), binwidth = 0.01)+
  theme_classic()+
  scale_y_continuous("Frequency")+
  scale_x_continuous("Weight")


tiff(file = "Figures/Sup_Figure_1_test.tiff",   
     width = 5.8, 
     height = 2.4,
     units = 'in',  
     res = 300, 
     compression = "lzw")
print(sup_figure1)
dev.off()
