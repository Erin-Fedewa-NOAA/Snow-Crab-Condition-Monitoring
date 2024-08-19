#GOALS: 
# 1) Append haul data to snow crab biometrics datasheet
#2) Add a column for maturity using clutch codes/chela heights
    #Distribution based cutlines for males published in 2019 tech memo
#3) Calculate snow crab CPUE at each station
#4) Join haul data with hepato DWT:WWT data 

#Follow ups: Change chela NBS vrs EBS script section to recognize 01 vrs 02 string

# Author: EJF

# load ----
library(tidyverse)

#############################
#Append maturity to biometrics data
bio_dat <- read.csv("./data/2019_2024 biometrics data.csv")

#Determine male maturity via distribution-based cutline method/clutch codes
bio_dat %>%
  mutate(CW = as.numeric(CW),
      maturity = case_when((Sex == 2 & CH_CC > 0) ~ 1, #mature female (EBS & NBS)
                              (Sex == 2 & CH_CC == 0) ~ 0, #immature female (EBS & NBS)
                              #Define conditions for EBS immature/mature male using cutlines 
                              (Sex == 1 & Cruise %in% c(201901,202101,202201,202301, 202401) & 
                                 log(CH_CC) < -2.20640 + 1.13523 * log(CW)) | #immature male EBS
                                      (Sex == 1 & CW < 50 & Cruise %in% c(201901,202101,202201,202301)) |
                                      (Sex == 1 & CH_CC == 0) ~ 0, #in 2024, datasheets read "imm" for males because
                                          #maturity was confirmed with maturity app. Coded as 0 
                              (Sex == 1 & Cruise %in% c(201901,202101,202201,202301,202401) &
                                 log(CH_CC) >= -2.20640 + 1.13523 * log(CW)) ~ 1, #mature male EBS
                              #NBS male cutlines
                              (Sex == 1 & Cruise %in% c(201902,202102,202202,202302) & 
                                 log(CH_CC) < -1.916947 + 1.070620 * log(CW))| (Sex == 1 & CW < 40 &
                                  Cruise %in% c(201902,202102,202202,202302)) ~ 0, #immature male NBS
                           (Sex == 1 & Cruise %in% c(201902,202102,202202,202302) &
                              log(CH_CC) >= -1.916947 + 1.070620 * log(CW)) ~ 1)) -> cond_mat #mature male NBS
                                                            

#############
#Append EBS & NBS haul data 
ebs_haul <- read.csv("./data/crabhaul_opilio.csv")
nbs_haul <- read.csv("./data/haul_opilio_nbs.csv")
  
#Combine ebs and nbs haul files 
ebs_haul %>%
  bind_rows(nbs_haul) %>% 
  rename_with(tolower) %>%
  filter(cruise %in% c(201901, 201902, 202101, 202102, 202201, 202202, 202301, 202302, 202401),
         haul_type==3) %>%
  select(vessel, cruise, haul, mid_latitude, mid_longitude, gis_station, 
         bottom_depth, gear_temperature, start_date) %>%
  distinct() -> snow_haul

#Join haul and biometric datasets 
cond_mat %>% 
  rename_with(tolower) %>%
  left_join(snow_haul, by = c("cruise", "vessel", "haul")) -> mat_haul

#Add in sampling regions associated with each station
  #read in lookup table
regions <- read.csv("./data/regions_lookup.csv")

#join
mat_haul %>%
  left_join(regions, by="gis_station") %>%
#calculate additional WWT:DWT/FA metrics
  mutate(DWT_WWT = hepato_dwt/hepato_wwt,
       Perc_DWT = DWT_WWT*100,
       WWT_DWT = hepato_wwt/hepato_dwt) %>%
  mutate(year = as.numeric(str_extract(cruise, "\\d{4}"))) %>%
  write_csv(file="./output/condition_master.csv")






