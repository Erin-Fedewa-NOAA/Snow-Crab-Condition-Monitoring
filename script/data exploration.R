#Data exploration plots & data summaries for energetic condition dataset

# Author: EJF

# load ----
library(tidyverse)
library(sf)
library(ggmap)
library(gganimate)
library(viridis)
library(ggridges)
library(RColorBrewer)
library(broom)

condition_master <- read.csv("./output/condition_master.csv")

####################################
#SAMPLE SIZES AND SPATIAL EXTENT

#plot sampling locations by year 
#Basemaps
usa <- raster::getData("GADM", country = c("USA"), level = 1, path = "./data")
can <- raster::getData("GADM", country = c("CAN"), level = 1, path = "./data")

#Sample size by year map
condition_master %>% 
  group_by(year, mid_latitude, mid_longitude) %>%
  summarise(n_crab=n()) %>%
  ggplot() + 
  geom_polygon(data = usa, aes(x = long, y = lat, group = group))+
  geom_point(aes(x = mid_longitude, y = mid_latitude, size=n_crab), color= "light blue")+
  coord_quickmap(xlim = c(-179, -158), ylim = c(53, 66)) +
  theme_bw() +
  facet_wrap(~year)
  ggsave("./figures/data exploration/n_year.png", dpi=300)

# Sample sizes by year plot
condition_master %>%
  group_by(year,lme) %>%
  count() %>%
  filter(lme != "NA") %>% #one crab collected outside the sampling design
  ggplot() +
  geom_bar(aes(x=as.factor(lme), y= n), stat='identity') +
  facet_wrap(~year) +
  theme_bw() +
  labs(x= "", y = "Sample size")

#Sample sizes by maturity
condition_master %>%
  filter(maturity != "NA") %>%
  group_by(lme,year,maturity) %>%
  count() %>%
  ggplot() +
  geom_bar(aes(x=as.factor(maturity), y= n), stat='identity') +
  facet_wrap(~year) +
  theme_bw() +
  labs(x= "Maturity Status", y = "Sample size") #45 mature males in dataset

#Sample sizes by sex
condition_master %>%
  filter(sex != "NA") %>%
  group_by(year,sex) %>%
  count() %>%
  ggplot() +
  geom_bar(aes(x=as.factor(sex), y= n), stat='identity') +
  facet_wrap(~year) +
  theme_bw() +
  labs(x= "Sex", y = "Sample size")

#Size composition sampled by region/yr
condition_master %>%
  mutate(Sex = recode_factor(sex, '1' = "M", '2' = "F")) %>%
  filter(lme != "NA") %>% #one crab collected outside the sampling design
  group_by(year) %>%
  ggplot() +
  geom_density(aes(x=cw, fill=Sex), position = "stack", binwidth = 2) +
  scale_fill_manual(values=c("#00BFC4", "#F8766D")) +
  facet_grid(lme~year) +
  theme_bw() +
  labs(x= "Snow crab carapace width (mm)", y = "Count")

#Pull very large females that appear to be Tanner crab and re-plot
condition_master %>%
  mutate(Sex = recode_factor(sex, '1' = "M", '2' = "F")) %>%
  filter(lme != "NA", #one crab collected outside the sampling design
         !vial_id %in% c("2019-65","2019-67","2019-68","2019-71","2019-66")) %>% #large females likely tanners
  group_by(year) %>%
  ggplot() +
  geom_density(aes(x=cw, fill=Sex), position = "stack", binwidth = 2) +
  scale_fill_manual(values=c("#00BFC4", "#F8766D")) +
  facet_grid(lme~year) +
  theme_bw() +
  labs(x= "Snow crab carapace width (mm)", y = "Count")
ggsave("./figures/data exploration/size_comp.png", dpi=300)
#will need to factor size into models! 

#Size range sampled across years by region
condition_master %>%
  filter(lme != "NA", #one crab collected outside the sampling design
         !vial_id %in% c("2019-65","2019-67","2019-68","2019-71","2019-66")) %>%
  group_by(lme, year, sex) %>%
  summarize(avg_cw = mean(cw, na.rm=T), 
            max_cw = max(cw, na.rm=T), 
            min_cw = min(cw, na.rm=T))


#data wrangling - use this dataset for all further analyses!!
condition_master %>%
  filter(lme != "NA", #one crab collected outside the sampling design
         !vial_id %in% c("2019-65","2019-67","2019-68","2019-71","2019-66","2019-207", "2019-212"),#likely tanners
         maturity != 1) -> new.dat
  
#%DWT by lme and year
new.dat %>%
  ggplot(aes(factor(year), Perc_DWT)) +
  geom_boxplot() +
  facet_wrap(~lme) +
  theme_bw() +
  labs(x= "", y = "% DWT in hepatopancreas")
ggsave("./figures/data exploration/DWT_year.png", dpi=300)




