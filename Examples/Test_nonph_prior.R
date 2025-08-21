

library(tibble)
library(dplyr)
library(survextrap)

trial_data <- readRDS("Data/trial_IPD_OS_ALEX_Feb_2017_v4.rds")

prior_hr(p_normal(0, 2.5))

#colons


nd_mod <- survextrap(Surv(time, status) ~ trt, data=trial_data, fit_method="opt")
sp <- nd_mod$mspline

#mspline <- mspline_spec(Surv(years, d) ~ 1, data=cetux, df=6, add_knots=20)

mspline <- mspline_spec(Surv(time, status) ~ trt, df = 10, data=trial_data, add_knots =  c(2,3,5))

prior_hscale <- p_meansurv(median=5, upper=10, mspline=mspline)
#prior_haz_const(mspline, prior_hscale = prior_hscale)

#set.seed(1)
prior_hsd <- p_gamma(2, 1)
prior_haz_sd(mspline = mspline,
             prior_hsd = prior_hsd,
             prior_hscale = prior_hscale)

prior_loghr <- p_hr(median=1, upper=10)
prior_hr(prior_loghr)

set.seed(1)
prior_hrsd <- p_gamma(2, 5)

prior_hr_sd(mspline = mspline,                              
            prior_hsd = prior_hsd,
            prior_loghr = prior_loghr,
            prior_hscale = prior_hscale,
            prior_hrsd = prior_hrsd,
            formula = ~treat,
            nonprop = ~treat,
            newdata = data.frame(treat=1), 
            newdata0 = data.frame(treat=0),
            nsim = 10000)

#?prior_hr_sd




trial_data
############

mspline <-  mspline_spec(Surv(time,status) ~ trt, 
                             data = trial_data, 
                             df=10)

mspline_rwe <-  mspline_spec(Surv(time,status) ~ trt, 
                         data = trial_data, 
                         df=10, 
                         add_knots = c(2, 3, 5))

prior_hscale <- p_meansurv(median=5, upper=20, mspline=mspline)
prior_hscale_rwe <- p_meansurv(median=5, upper=20, mspline=mspline_rwe)
prior_hsd <- p_gamma(2, 1)
prior_loghr <- p_hr(median=1, upper=10)
prior_hrsd <- p_gamma(2, 5)



#trial_data
prior_hr_sd(mspline = mspline,                              
            prior_hsd = prior_hsd,
            prior_loghr = prior_loghr,
            prior_hscale = prior_hscale,
            prior_hrsd = prior_hrsd,
            formula = ~trt,
            nonprop = ~trt,
            newdata = tibble(trt="Alectinib") %>% mutate(trt=factor(trt, levels= c("Crizotinib", "Alectinib"))), 
            newdata0 = tibble(trt="Crizotinib")  %>% mutate(trt= factor(trt, levels= c("Crizotinib", "Alectinib"))),
            nsim = 1000)


test_mod <- survextrap(Surv(time,status) ~ trt, 
                               data = trial_data , 
                               mspline = mspline,
                               prior_hscale = prior_hscale,
                               prior_loghr = prior_loghr,
                               prior_hsd = prior_hsd,
                               fit_method = "mcmc")


plot(test_mod, tmax = 20)
rmst(test_mod, t = 20)
rmst(test_mod, t = 5)

test_mod_nonph <- survextrap(Surv(time,status) ~ trt, 
                           data = trial_data, 
                           mspline = mspline,
                           nonprop = T,
                           prior_hscale = prior_hscale,
                           prior_loghr = prior_loghr,
                           prior_hsd = prior_hsd,
                           prior_hrsd = prior_hrsd, 
                           fit_method = "mcmc")

plot(test_mod_nonph, tmax = 20)
rmst(test_mod_nonph, t = 20)
rmst(test_mod_nonph, t = 5)

test_mod_nonph_rwe <- survextrap(Surv(time,status) ~ trt, 
                             data = trial_data, 
                             external = external_data_unweighted %>%
                               filter(start >= 2.5),
                             mspline = mspline_rwe,
                             nonprop = T,
                             prior_hscale = prior_hscale_rwe,
                             prior_loghr = prior_loghr,
                             prior_hsd = prior_hsd,
                             prior_hrsd = prior_hrsd,
                             fit_method = "mcmc")

plot(test_mod_nonph_rwe, tmax = 20)
rmst(test_mod_nonph_rwe, t = 20)
rmst(test_mod_nonph_rwe, t = 5)
