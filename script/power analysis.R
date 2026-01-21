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
library(ggthemes)
library(hrbrthemes)
library(crabpack)
library(lubridate)
library(magrittr)
library(akgfmaps)

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

#and now an effect size representative of half the difference between 2019 and 2021
effect_size(32.6, 25, 8, 8)

#Using this effect size and know sample sizes, let's estimate our power for two years
pwr.t.test(n = 120,
  d = effect_size(32.6, 25, 8, 8),
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

#############################################################
#Additional requests from data working group:

#Evaluate evidence for sex-specific differences in energetic condition
  #i.e. can we just collect females for monitoring condition?

#plot annual mean, both sexes
sc_condition %>%
  filter(lme == "EBS", 
         !vial_id %in% c("2019-65","2019-67","2019-68","2019-71","2019-66"),
         maturity != 1) %>%
  group_by(year, lme) %>%
  summarise(sd = sd(Perc_DWT, na.rm = TRUE),
            cond = mean(Perc_DWT, na.rm=T)) %>%
  arrange(year) %>%
  filter(lme != "NA") %>%
  ggplot(aes(as.factor(year), cond, fill=year)) +
  geom_bar(stat="identity") + 
  geom_errorbar(aes(ymin = cond-sd, ymax = cond+sd), width = 0.2, color="gray45") +
  theme_bw() +
  #theme_ipsum(axis_title_just = "cc", axis_title_size = 11, axis_text_size =10) +
  labs(y = "Snow Crab Condition (% DWT)", x = "") +
  theme(axis.text=element_text(size=12)) +
  theme(legend.position = "none") +
  geom_hline(yintercept=22.6, color="red") + #adding prelim threshold from starvation lab experiment
  ggtitle("EBS Both Sexes")
 
#males only
sc_condition %>%
  filter(lme == "EBS", 
         !vial_id %in% c("2019-65","2019-67","2019-68","2019-71","2019-66"),
         maturity != 1,
         sex == 1) %>%
  group_by(year, lme) %>%
  summarise(sd = sd(Perc_DWT, na.rm = TRUE),
            cond = mean(Perc_DWT, na.rm=T)) %>%
  arrange(year) %>%
  filter(lme != "NA") %>%
  ggplot(aes(as.factor(year), cond, fill=year)) +
  geom_bar(stat="identity") + 
  geom_errorbar(aes(ymin = cond-sd, ymax = cond+sd), width = 0.2, color="gray45") +
  theme_bw() +
  #theme_ipsum(axis_title_just = "cc", axis_title_size = 11, axis_text_size =10) +
  labs(y = "Snow Crab Condition (% DWT)", x = "") +
  theme(axis.text=element_text(size=12)) +
  theme(legend.position = "none") +
  geom_hline(yintercept=22.6, color="red") + #adding prelim threshold from starvation lab experiment
  ggtitle("EBS Males only")

#females only
sc_condition %>%
  filter(lme == "EBS", 
         !vial_id %in% c("2019-65","2019-67","2019-68","2019-71","2019-66"),
         maturity != 1,
         sex == 2) %>%
  group_by(year, lme) %>%
  summarise(sd = sd(Perc_DWT, na.rm = TRUE),
            cond = mean(Perc_DWT, na.rm=T)) %>%
  arrange(year) %>%
  filter(lme != "NA") %>%
  ggplot(aes(as.factor(year), cond, fill=year)) +
  geom_bar(stat="identity") + 
  geom_errorbar(aes(ymin = cond-sd, ymax = cond+sd), width = 0.2, color="gray45") +
  theme_bw() +
  #theme_ipsum(axis_title_just = "cc", axis_title_size = 11, axis_text_size =10) +
  labs(y = "Snow Crab Condition (% DWT)", x = "") +
  theme(axis.text=element_text(size=12)) +
  theme(legend.position = "none") +
  geom_hline(yintercept=22.6, color="red") + #adding prelim threshold from starvation lab experiment
  ggtitle("EBS Females only")

#Evaluate evidence for size signal in energetic condition
  
#bin data, plot by %DWT
sc_condition %>%
  filter(lme == "EBS", 
         !vial_id %in% c("2019-65","2019-67","2019-68","2019-71","2019-66"),
         maturity != 1) %>%
  mutate(size_bin = cut(cw, breaks=c(35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110))) %>%
  group_by(size_bin, sex, year) %>%
  summarise(cond = mean(Perc_DWT, na.rm=T),
            sample_size = n()) %>%
  filter(size_bin != "NA") -> plot_dat

#male plot
plot_dat %>% 
  filter(sex ==1) %>%
  ggplot(aes(as.factor(size_bin), cond, fill=sample_size)) +
  geom_col(stat="identity") +
  geom_hline(yintercept=22.6, color="red") +
facet_wrap(~year)

#female plot
plot_dat %>% 
  filter(sex ==2) %>%
  ggplot(aes(as.factor(size_bin), cond, fill=sample_size)) +
  geom_col(stat="identity") +
  geom_hline(yintercept=22.6, color="red") +
  facet_wrap(~year)

###########################
#hypothetical sampling effort of random uniform sampling design

#if we used size criteria and proposed goal of 180 samples per year in the EBS,
  #what would spatial coverage look like if we sampled 1 crab per station without
  #strata and per region sampling goals?

#pull data from crabpack
## Pull haul and snow crab specimen data ----
snow <- get_specimen_data(species = "SNOW",
                          region = "EBS",
                          channel = "KOD")

haul <- snow$haul

#Filter data for condition sampling criteria- note this is very imperfect b/c 
  #we can't assign male maturity to specimen-level data. We'll use a size cutoff 
  #and shell condition 2 only to attempt to eliminate some of the mature males
snow$specimen %>%
  left_join(haul) %>%
  filter(YEAR %in% 2017:2025,
         SEX == 1 & SIZE >= 70 & SIZE <= 95 & SHELL_CONDITION < 3 | 
           SEX == 2 & SIZE >= 45 & CLUTCH_SIZE == 0 ) %>%
  group_by(YEAR, STATION_ID, HAUL, SEX, START_DATE, MID_LATITUDE, MID_LONGITUDE) %>%
  summarise(condition_samples = n()) -> dat

#Max number of total samples/stations sampled per sex if only 1 sample collected
  #per station
dat %>% 
  filter(SEX == 2) %>%
  group_by(YEAR) %>%
  summarize(n = n())

dat %>% 
  filter(SEX == 1) %>%
  group_by(YEAR) %>%
  summarize(n = n())

#total sample size by year
dat %>%
  group_by(YEAR) %>%
  summarize(n = n())

#splice data for 90 females per year. limiting to years when we have 
  #condition data
dat %>%
  filter(SEX == 2,
         YEAR %in% 2019:2025) %>%
  group_by(YEAR) %>%
  arrange(START_DATE) %>%
  group_by(YEAR) %>%
  slice_head(n = 60) -> female_dat

#splice data for 90 males per year
dat %>%
  filter(SEX == 1,
         YEAR %in% 2019:2025) %>%
  group_by(YEAR) %>%
  arrange(START_DATE) %>%
  group_by(YEAR) %>%
  slice_head(n = 60) -> male_dat

#combine for final dataset to plot
female_dat %>%
  bind_rows(male_dat) %>%
  group_by(YEAR, STATION_ID, MID_LATITUDE, MID_LONGITUDE) %>%
  summarise(n = sum(condition_samples))-> map_data

#plot maps of sampling effort
 
## SET COORDINATE REFERENCE SYSTEMS (CRS) --------------------------------------
in.crs <- "+proj=longlat +datum=NAD83" #CRS is in lat/lon
map.crs <- "EPSG:3338" # final crs for mapping/plotting: Alaska Albers

## LOAD SHELLFISH ASSESSMENT PROGRAM GEODATABASE -------------------------------
survey_gdb <- "./data/SAP_layers" 
survey_strata <- terra::vect(survey_gdb, layer = "EBS.NBS_surveyarea")
#EBS/NBS Boundary line
boundary <- st_read(layer = "EBS_NBS_divide", survey_gdb) 

## LOAD ALASKA REGION LAYERS (FROM AKGFMAPS R package) -----------------------------------
ebs_layers <- akgfmaps::get_base_layers(select.region = "ebs", set.crs = "EPSG:3338")
ebs_survey_areas <- ebs_layers$survey.area
ebs_survey_areas$survey_name <- c("Eastern Bering Sea", "Northern Bering Sea")

#Transform survey crab data into spatial data frame
map_data %>% 
  # Convert lat/long to an sf object
  st_as_sf(coords = c("MID_LONGITUDE", "MID_LATITUDE"), crs = st_crs(4326)) %>%
  #st_as_sf needs crs of the original coordinates- need to transform to Alaska Albers
  st_transform(crs = st_crs(3338)) -> crab_dat

#Transform condition crab data into spatial data frame
sc_condition %>%
  filter(lme == "EBS", 
         !vial_id %in% c("2019-65","2019-67","2019-68","2019-71","2019-66"),
         maturity != 1) %>%
  rename(YEAR=year, MID_LATITUDE=mid_latitude, MID_LONGITUDE=mid_longitude) %>%
  group_by(YEAR, MID_LATITUDE, MID_LONGITUDE) %>%
  summarise(n_crab=n()) %>%
  # Convert lat/long to an sf object
  st_as_sf(coords = c("MID_LONGITUDE", "MID_LATITUDE"), crs = st_crs(4326)) %>%
  #st_as_sf needs crs of the original coordinates- need to transform to Alaska Albers
  st_transform(crs = st_crs(3338)) -> condition_dat

#map
ggplot() +
  geom_sf(data = ebs_layers$survey.grid, fill=NA, color=alpha("grey80"))+
  geom_sf(data = ebs_survey_areas, fill = NA) +
  geom_sf(data = ebs_layers$akland, fill = "grey80", color = "black") +
  #add hypothetical crab sampling
  geom_sf(data=crab_dat, color = "grey30", alpha = .6) +
  #add existing random stratified sampling
  geom_sf(data=condition_dat, aes(size = n_crab), color = "blue", alpha = .6) +
  geom_sf(data= boundary, linewidth = 1, color = "grey40") +
  scale_x_continuous(limits = ebs_layers$plot.boundary$x,
                     breaks = ebs_layers$lon.breaks) +
  scale_y_continuous(limits = ebs_layers$plot.boundary$y,
                     breaks = ebs_layers$lat.breaks) +
  scale_size_continuous(range = c(1,4)) +
  theme_bw() +
  facet_wrap(~YEAR) +
  labs(x="", y="", size = expression(paste("Snow crab \n samples"))) +
  theme(legend.position="bottom") +
  guides(size = guide_legend(theme = theme(
    legend.title = element_text(size = 9)))) +
  theme(plot.margin = margin(0,-5,0,-5)) +
  theme(axis.text=element_text(size=8))

#Now let's try uniform random sampling without a sample cap, and instead with reduced
  #sampling frequency of only sampling odd-numbered hauls

#female dataset subset for odd hauls only 
dat %>%
  filter(SEX == 2) %>%
  group_by(YEAR) %>%
  filter(HAUL %% 2 == 1) -> female_dat_odd

#number of female samples collected
female_dat_odd %>%
  group_by(YEAR) %>%
  summarize(n = n())

#number of odd stations with >1 female that meets criteria 
female_dat_odd %>%
  group_by(YEAR) %>%
  filter(condition_samples > 1) %>%
  summarize(n = n())

#male dataset subset for odd hauls only 
dat %>%
  filter(SEX == 1) %>%
  group_by(YEAR) %>%
  filter(HAUL %% 2 == 1) -> male_dat_odd

#number of male samples collected
male_dat_odd %>%
  group_by(YEAR) %>%
  summarize(n = n())

#number of odd stations with >1 male that meets criteria 
male_dat_odd %>%
  group_by(YEAR) %>%
  filter(condition_samples > 1) %>%
  summarize(n = n())

#total sample size by year
dat %>%
  group_by(YEAR) %>%
  filter(HAUL %% 2 == 1) %>%
  group_by(YEAR) %>%
  summarize(n = n())

#Using effect size from 7% reduction in energetic condition and these sample sizes, 
  #let's estimate our power 
pwr.t.test(n = 75,
           d = effect_size(32.6, 25, 8, 8),
           sig.level = 0.05,
           type = "two.sample")
#99% probability of detecting a ~7% reduction in energetic condition even with 
  #lowest sample size

#combine for final dataset to plot map
female_dat_odd %>%
  bind_rows(male_dat_odd) %>%
  group_by(YEAR, STATION_ID, MID_LATITUDE, MID_LONGITUDE) %>%
  summarise(n = sum(condition_samples))-> map_data_odd

#and plot spatial extent of sampling effort

#Transform survey crab data into spatial data frame
map_data_odd %>% 
  # Convert lat/long to an sf object
  st_as_sf(coords = c("MID_LONGITUDE", "MID_LATITUDE"), crs = st_crs(4326)) %>%
  #st_as_sf needs crs of the original coordinates- need to transform to Alaska Albers
  st_transform(crs = st_crs(3338)) -> crab_dat

#map
ggplot() +
  geom_sf(data = ebs_layers$survey.grid, fill=NA, color=alpha("grey80"))+
  geom_sf(data = ebs_survey_areas, fill = NA) +
  geom_sf(data = ebs_layers$akland, fill = "grey80", color = "black") +
  #add hypothetical crab sampling layer
  geom_sf(data=crab_dat, color = "grey30", alpha = .6) +
  geom_sf(data= boundary, linewidth = 1, color = "grey40") +
  scale_x_continuous(limits = ebs_layers$plot.boundary$x,
                     breaks = ebs_layers$lon.breaks) +
  scale_y_continuous(limits = ebs_layers$plot.boundary$y,
                     breaks = ebs_layers$lat.breaks) +
  scale_size_continuous(range = c(1,4)) +
  theme_bw() +
  facet_wrap(~YEAR) +
  labs(x="", y="", size = expression(paste("Snow crab \n samples"))) +
  theme(legend.position="bottom") +
  guides(size = guide_legend(theme = theme(
    legend.title = element_text(size = 9)))) +
  theme(plot.margin = margin(0,-5,0,-5)) +
  theme(axis.text=element_text(size=8))





