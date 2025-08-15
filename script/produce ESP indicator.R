# Develop ESP indicator for juvenile snow crab energetic condition:
#annual mean % DWt of hepatopancreas

# Erin Fedewa

#NOTE: In 2025, ESP indicator contribution was transitioned to a model-based output
  #in order to correct for julian day effects 

# load ----
library(tidyverse)
library(lubridate)
library(rstan)
library(brms)
library(bayesplot)
library(marginaleffects)
library(emmeans)
library(MARSS)
library(corrplot)
library(factoextra)
library(patchwork)
library(modelr)
library(broom.mixed)
library(pROC)
library(ggthemes)
library(tidybayes)
library(RColorBrewer)
library(knitr)
library(loo)
library(sjPlot)
library(hrbrthemes)
library(bayestestR)
source("./script/stan_utils.R")

sc_condition <- read.csv("./output/condition_master.csv")

#######################################
#Calculate observed mean condition (ESP indicator prior to 2025)
sc_condition %>%
  filter(lme != "NA", #one crab collected outside the sampling design
         lme != "NBS",
         !vial_id %in% c("2019-65","2019-67","2019-68","2019-71","2019-66"),
         maturity != 1) %>%
  group_by(year) %>%
  summarise(avg_condition = mean(Perc_DWT, na.rm=T)) -> cond

#Plots 
missing <- data.frame(year = 2020)

sc_condition %>%
  filter(lme != "NA", #one crab collected outside the sampling design
         !vial_id %in% c("2019-65","2019-67","2019-68","2019-71","2019-66"),
         maturity != 1) %>%
  group_by(year, lme) %>%
  summarise(sd = sd(Perc_DWT, na.rm = TRUE),
            cond = mean(Perc_DWT, na.rm=T)) %>%
  bind_rows(missing) %>%
  arrange(year) %>%
  filter(lme != "NA") %>%
  mutate(lme = recode(lme, EBS = "Eastern Bering Sea",
                      NBS = "Northern Bering Sea")) -> plot

#plot annual mean as bar plot 
plot %>% 
  ggplot(aes(as.factor(year), cond, fill=year)) +
  geom_bar(stat="identity") + 
  geom_errorbar(aes(ymin = cond-sd, ymax = cond+sd), width = 0.2, color="gray45") +
  facet_wrap(~lme) +
  theme_bw() +
  theme_ipsum(axis_title_just = "cc", axis_title_size = 11, axis_text_size =10) +
  labs(y = "Snow Crab Condition (% DWT)", x = "") +
  theme(axis.text=element_text(size=12)) +
  theme(legend.position = "none") +
  geom_hline(yintercept=22.6, color="red") + #adding prelim threshold from starvation lab experiment
  geom_rect(aes(xmin=0, xmax=Inf, ymin=22.6, ymax=22.6 + 2.9), fill="red", alpha = 0.05) + 
  geom_rect(aes(xmin=0, xmax=Inf, ymin=22.6-2.9, ymax=22.6), fill="red", alpha = 0.05) +
  geom_vline(data = subset(plot, lme == "Eastern Bering Sea"), aes(xintercept = 1.5), linetype="dashed") +
  geom_text(data = subset(plot, lme == "Eastern Bering Sea"), aes(x = .7, y=44, label = "Mid-collapse"),
            size = 2.5, color = "#D55E00") +
  geom_text(data = subset(plot, lme == "Eastern Bering Sea"), aes(x = 3, y=44, label = "Post-collapse"),
            size = 2.5, color = "#084594") 
ggsave("./figures/observed_condition.png", height = 6, width = 7, units = "in", dpi = 300)


############################################################################
#Improvements for 2025- deriving annual means from a model controlling for seasonality 
  #and crab size (see Fedewa et al 2025 for further detail)

#data wrangling- EBS dataset  
sc_condition %>%
  mutate(julian=yday(parse_date_time(start_date, "ymd", "US/Alaska"))) %>%
  filter(lme == "EBS", 
         !vial_id %in% c("2019-65","2019-67","2019-68","2019-71","2019-66"), 
         maturity != 1,
         Perc_DWT >= 0) %>%
  mutate(year = as.factor(year),
         sex = as.factor(sex),
         region = as.factor(sample_region),
         station = as.factor(station_id),
         julian = as.numeric(julian),
         perc_dwt = as.numeric(Perc_DWT)) -> ebs.dat 

#EBS ANNUAL MEANS
ebs_annual_final_formula <-  bf(perc_dwt | trunc(lb = 0) ~ s(cw, k = 3) + s(julian, k = 3) +
                                  year + (1 | region))  

ebs_annual_final <- brm(ebs_annual_final_formula,
                        data = ebs.dat,
                        family = gaussian,
                        cores = 4, chains = 4, iter = 10000, warmup = 1000,
                        save_pars = save_pars(all = TRUE), seed = 3,
                        control = list(adapt_delta = 0.999, max_treedepth = 14))

#Save model output 
saveRDS(ebs_annual_final, file = "./output/ebs_annual_final.rds")
ebs_annual_final <- readRDS("./output/ebs_annual_final.rds")

#MCMC convergence diagnostics 
check_hmc_diagnostics(ebs_annual_final$fit)
neff_lowest(ebs_annual_final$fit)
rhat_highest(ebs_annual_final$fit)
summary(ebs_annual_final) #dramatically lower condition in 2019
bayes_R2(ebs_annual_final) #r2 = .52 

#Diagnostic Plots
plot(ebs_annual_final, ask = FALSE)
plot(conditional_effects(ebs_annual_final), ask = FALSE)
mcmc_plot(ebs_annual_final, prob = 0.95)
mcmc_neff(neff_ratio(ebs_annual_final)) #Effective sample size: All ratios > 0.1

#Posterior Predictive Check Plots:
pp_check(ebs_annual_final)
pp_check(ebs_annual_final, type = "ecdf_overlay")
pp_check(ebs_annual_final, type = "stat", stat = "mean")
pp_check(ebs_annual_final, type = "stat", stat = "min")
pp_check(ebs_annual_final, type = "stat", stat = "max")

#-----------------------------------------------------------------------------------
#Extract conditional effect of year 

conditional_effects(ebs_annual_final, effect = "year")

ce1s_1 <- conditional_effects(ebs_annual_final, effect = "year", re_formula = NA,
                              probs = c(0.025, 0.975)) 
ce1s_1$year %>%
  dplyr::select(year, estimate__, lower__, upper__) %>%
  rename(annual_mean = estimate__, 
         lower_CI = lower__,
         upper_CI = upper__) -> year_ebs

#plot model-derived annual mean as bar plot 
year_ebs %>% 
  ggplot(aes(as.factor(year), annual_mean)) +
  geom_bar(stat="identity", fill="#a6bddb") + 
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI), width = 0.2, color="gray45") +
  theme_bw() +
  theme_ipsum(axis_title_just = "cc", axis_title_size = 11, axis_text_size =10) +
  labs(y = "Snow Crab Condition (% DWT)", x = "") +
  theme(axis.text=element_text(size=12)) +
  theme(legend.position = "none") +
  geom_hline(yintercept=22.6, color="red") + #adding prelim threshold from starvation lab experiment
  geom_rect(aes(xmin=0, xmax=Inf, ymin=22.6, ymax=22.6 + 2.9), fill="red", alpha = 0.05) + 
  geom_rect(aes(xmin=0, xmax=Inf, ymin=22.6-2.9, ymax=22.6), fill="red", alpha = 0.05) 
ggsave("./figures/model_condition.png", height = 6, width = 7, units = "in", dpi = 300)

#Write output
missing <- data.frame(year = 2020)

year_ebs %>%
  mutate(year = as.numeric(as.character(year))) %>%
  bind_rows(missing) %>%
  arrange(year) %>%
  write.csv(file="./output/opilio_condition.csv")

