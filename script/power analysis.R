#Perform a power analysis to inform standard bottom trawl hepatopancreas
  #sampling efforts, i.e. ensuring enough samples collected to detect year to 
  #year differences in energetic condition of snow crab population

# load ----
library(tidyverse)
library(pwr)
library(effsize)
library(gt)
library(sjPlot)
library(marginaleffects)
library(patchwork)

sc_condition <- read.csv("./output/condition_master.csv")

######################################
#quick exploration of sample sizes by year/region of 2019-2025 collections 

# Sample sizes by year 
sc_condition %>%
  group_by(year,lme) %>%
  count() %>%
  filter(lme != "NA") 

#EBS sample sizes by region
sc_condition %>%
  filter(lme == "EBS") %>%
  group_by(year,sample_region) %>%
  count() %>%
  pivot_wider(names_from = year, values_from = n) 

#EBS # of stations sampled per year
sc_condition %>%
  filter(lme == "EBS") %>%
  group_by(year) %>%
  summarise(no_stations = n_distinct(station_id)) 

#EBS # of stations sampled per region each year
sc_condition %>%
  filter(lme == "EBS") %>%
  group_by(year, sample_region) %>%
  summarise(no_stations = n_distinct(station_id)) %>%
  ggplot() +
  geom_bar(aes(x=as.factor(sample_region), y=no_stations), stat='identity') +
  facet_wrap(~year) +
  theme_bw() +
  labs(x= "Region", y = "Number of stations sampled") 

# EBS # of stations where max of 5 males or 5 females are sampled
sc_condition %>%
  filter(lme == "EBS") %>%
  group_by(year, sample_region, sex, station_id) %>%
  count() %>%
  filter(case_when(sex ==1 ~ n >= 5,
                   sex == 2 ~ n >= 5)) %>%
  group_by(year, sample_region) %>%
  count() %>%
  ggplot() +
  geom_bar(aes(x=as.factor(sample_region), y=n), stat='identity') +
  facet_wrap(~year) +
  theme_bw() +
  labs(x= "Region", y = "Number of stations where max # of males or females was sampled") 

#Now let's plot this as a % of the number of stations sampled within each region
sc_condition %>%
  filter(lme == "EBS") %>%
  group_by(year, sample_region, sex, station_id) %>%
  count() %>%
  filter(case_when(sex ==1 ~ n >= 5,
                   sex == 2 ~ n >= 5)) %>%
  group_by(year, sample_region) %>%
  count() %>%
  rename(n_max = n) -> n_max
  
sc_condition %>%
  filter(lme == "EBS") %>%
  group_by(year, sample_region) %>%
  summarise(n_stations = n_distinct(station_id)) %>%
  right_join(n_max) %>%
  mutate(perc = (n_max/n_stations)*100) %>%
  ggplot() +
  geom_bar(aes(x=as.factor(sample_region), y=perc), stat='identity') +
  facet_wrap(~year) +
  theme_bw() +
  labs(x= "Region", y = "% of stations sampled where max # \n of males or females was collected") 

#######################################
#Power analysis step 1: Estimate the size of the effect between year and hepato WWT:DWT in 
  #current dataset (Cohen's D)

# A custom function to calculate Cohen's d for two independent samples 
effect_size <- function(u1, u2, sd1, sd2) {
  abs(u1 - u2) / sqrt((sd1^2 + sd2^2) / 2)
}

#calculate summary stats from data
sc_condition %>%
  filter(lme == "EBS", 
         !vial_id %in% c("2019-65","2019-67","2019-68","2019-71","2019-66"),
         maturity != 1) %>%
  group_by(year) %>%
  summarise(avg_condition = mean(Perc_DWT, na.rm=T),
            sd_condition = sd(Perc_DWT, na.rm=T),
            sum_samples = n()) 

#to calculate, let's use 2019 and 2025 means/sd's
effect_size(17.3, 30.7, 8.8, 5.20)
#by using 2019, we estimate a huge effect size, which isn't surprising

#Let's estimate effect size from two more similar years
effect_size(30.8, 30.7, 4.98, 5.20)
#0.02

#and two more years
effect_size(34.1, 32.6, 5.05, 8.24 )

#Using this effect size and know sample sizes, let's estimate our power for two years
pwr.t.test(n = 98,
  d = effect_size(17.3, 30.7, 8.8, 5.20),
  sig.level = 0.05,
  type = "two.sample")
#100% probability of detecting a ~50% reduction in energetic condition

#Now let's half our samples in these two years to see how power changes
sc_condition %>%
  filter(lme == "EBS", 
         !vial_id %in% c("2019-65","2019-67","2019-68","2019-71","2019-66"),
         maturity != 1) %>%
  group_by(year) %>%
  slice_sample(n=75) %>%
  summarise(avg_condition = mean(Perc_DWT, na.rm=T),
            sd_condition = sd(Perc_DWT, na.rm=T),
            sum_samples = n()) 

pwr.t.test(n = 75,
           d = effect_size(18.1, 32.0, 10.2, 6.24),
           sig.level = 0.05,
           type = "two.sample")
#100% probability of detecting a true effect

#And what about years that are more similar
pwr.t.test(n = 250,
           d = effect_size(32.0, 30.9, 4.46, 4.63),
           sig.level = 0.05,
           type = "two.sample")
#25% probability of detecting a true effect

#Now let's compare result of lm looking at annual differences in energetic conditon
  #when sample sizes are cut in half each year 
sc_condition %>%
  filter(lme == "EBS", 
         !vial_id %in% c("2019-65","2019-67","2019-68","2019-71","2019-66"),
         maturity != 1) -> dat

m1 <- lm(Perc_DWT ~ as.factor(year), data=dat) 
summary(m1)
plot_predictions(m1, by="year") -> m1_plot


#Now if we cut samples in half each year can we still detect a significant effect of year?
sc_condition %>%
  filter(lme == "EBS", 
         !vial_id %in% c("2019-65","2019-67","2019-68","2019-71","2019-66"),
         maturity != 1) %>%
  group_by(year) %>%
  slice_sample(n=75) -> split_dat

m2 <- lm(Perc_DWT ~ as.factor(year), data=split_dat) 
summary(m2)
plot_predictions(m2, by="year") -> m2_plot

m1_plot + m2_plot

#Seems that because we were able to detect such a large effect in 2019 despite only
  #75 samples the dataset seems pretty robust to reducing sample size 

#Step 2: Calculate the required sample size within a year given the size of the 
  #effect you have, the level of alpha (0.5) and the amount of power (0.9)
#Here we're solving for sample size, n, per group (i.e. within year), where k
#is our 6 year timeseries 

#large effect size
large <- pwr.anova.test(k = 6, n = NULL, f =.4, sig.level = 0.05, power = .9)
plot.power.htest(large)

#medium effect size
medium <- pwr.anova.test(k = 6, n = NULL, f =.25, sig.level = 0.05, power = .9)
plot.power.htest(medium)

#small effect size
small <- pwr.anova.test(k = 6, n = NULL, f =.1, sig.level = 0.05, power = .9)
plot.power.htest(small)
#so we'd not be able to detect small effect size


