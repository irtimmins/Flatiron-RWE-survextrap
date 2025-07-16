#?standsurv()

test_data <- readRDS("Data/trial_IPD_OS_ALEX_Feb_2017.rds")
test_data <- test_data %>%
  filter(trt == "Alectinib") %>%
  select(time, status)


test_cetux <- cetux %>%
  mutate(time = years, status = d) %>%
  filter(treat == "Control") %>%
  select(time, status)


# Cut off point at 18 months (in years).
cut_point <- 1.5

hybrid_model_test <- fit_hybrid_model(data = test_data, 
                                      cut_point = 1.5)

hybrid_model_test

test_survival <- hybrid_model_survival(hybrid_model_test, t = seq(from = 0, to = 10, length.out = 100))

test_survival %>%
  ggplot(aes(x = time, y = value))+
  theme_classic()+
  geom_ribbon(aes(x = time, ymin = lower, ymax = upper),
              alpha= 0.2)+
  geom_line()+
  geom_vline(xintercept = 1.5)+
  scale_y_continuous(limits = c(0,1))

test_survival_log_scale <- hybrid_model_survival_log_scale(hybrid_model_test, t = seq(from = 0, to = 10, length.out = 100))
test_survival_log_scale %>%
  ggplot(aes(x = time, y = value))+
  theme_classic()+
  geom_ribbon(aes(x = time, ymin = lower, ymax = upper),
              alpha= 0.2)+
  geom_line()+
  geom_vline(xintercept = 1.5)+
  scale_y_continuous(limits = c(0,1))



test_rmst <- hybrid_model_rmst(hybrid_model_test, t = 20)
test_rmst
test_rmst_boot <- hybrid_model_rmst_bootstrap(data = test_data, cut_point = 1.5, t= 20, N = 1e4)
summary(test_rmst_boot)
sd(test_rmst_boot)
test_rmst_boot_semi <- hybrid_model_rmst_semi_bootstrap(data = test_data, cut_point = 1.5, t = 20, N = 1e4)
summary(test_rmst_boot_semi)
head(test_rmst_boot_semi)
sd(test_rmst_boot_semi)


hist(test_rmst_boot_semi, breaks = 500)
hist(log(test_rmst_boot_semi), breaks = 500)
# View(test_mst_boot_semi)


test_rmst %>%
  filter(time < 1.8, time >= 1.3) %>%
  ggplot(aes(x = time, y = value))+
  theme_classic()+
  geom_ribbon(aes(x = time, ymin = lower, ymax = upper),
              alpha= 0.2)+
  geom_line()+
  geom_vline(xintercept = 1.5)+
  scale_x_continuous(limits = c(1.3,1.8)) 
#+
#  scale_y_continuous(limits = c(0,1))



hybrid_model_test_cetux <- fit_hybrid_model(data = test_cetux, 
                                            cut_point = 1.5)

test_cetux_rmst <- hybrid_model_rmst(hybrid_model_test_cetux, t = 1.55)
test_cetux_rmst

test_cetux_rmst_boot <- hybrid_model_rmst_bootstrap(data = test_cetux, cut_point = 1.5, t= 1.55, N = 1e4)
summary(test_cetux_rmst_boot)
sd(test_cetux_rmst_boot)
hist(test_cetux_rmst_boot, breaks = 500, xlim = c(3,7))

test_cetux_rmst_boot_semi <- hybrid_model_rmst_semi_bootstrap(data = test_cetux, cut_point = 1.5, t = 1.55, N = 1e4)
summary(test_cetux_rmst_boot_semi)
head(test_cetux_rmst_boot_semi)
sd(test_cetux_rmst_boot_semi)
## View(test_cetux_rmst_boot_semi)
hist(test_cetux_rmst_boot_semi, breaks = 500, xlim = c(3,7))

# test_survival_cetux <- 
# ............
test_survival_cetux <- hybrid_model_survival(hybrid_model_test_cetux, t = seq(from = 0, to = 10, length.out = 1e3))


test_survival_cetux %>%
  ggplot(aes(x = time, y = value))+
  theme_classic()+
  geom_ribbon(aes(x = time, ymin = lower, ymax = upper),
              alpha= 0.2)+
  geom_step()+
  geom_vline(xintercept = 1.5)+
  scale_y_continuous(limits = c(0,1))



test_cetux_survival_boot <- hybrid_model_survival_bootstrap(data = test_cetux, cut_point = 1.5, t= 3, N = 1e2)
sd(test_cetux_survival_boot)

hybrid_model_survival(hybrid_model = hybrid_model_test_cetux, t = 3)

