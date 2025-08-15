#GOALS: 
# 1) Append haul-level data to snow crab biometrics datasheet
#2) Add a column for maturity using clutch codes/chela heights
#3) Join haul data with hepato DWT:WWT data 


# Author: EJF

# load ----
library(tidyverse)
library(crabpack)
library(stringr)

#############################
#Append maturity to biometrics data
bio_dat <- read.csv("./data/2019_2025 biometrics data.csv")

#Determine male maturity via distribution-based cutline method/clutch codes
bio_dat %>%
  mutate(CW = as.numeric(CW),
      maturity = case_when((Sex == 2 & CH_CC > 0) ~ 1, #mature female (EBS & NBS)
                              (Sex == 2 & CH_CC == 0) ~ 0, #immature female (EBS & NBS)
                              #Define conditions for EBS immature/mature male using cutlines 
                              (Sex == 1 & str_detect(Cruise, "01$") & 
                                 log(CH_CC) < -2.360406 + 1.176558 * log(CW)) | #immature male EBS
                                      (Sex == 1 & CW < 50 & str_detect(Cruise, "01$")) |
                                      (Sex == 1 & CH_CC == 0) ~ 0, #in 2024, datasheets read "imm" for males because
                                          #maturity was confirmed with maturity app. Hard coded as 0 
                              (Sex == 1 & Cruise %in% c(201901,202101,202201,202301,202401) &
                                 log(CH_CC) >= -2.360406 + 1.176558 * log(CW)) ~ 1, #mature male EBS
                              #NBS male cutlines
                              (Sex == 1 & str_detect(Cruise, "02$") & 
                                 log(CH_CC) < -1.916947 + 1.070620 * log(CW)) | (Sex == 1 & CW < 40 &
                                      str_detect(Cruise, "02$")) ~ 0, #immature male NBS
                           (Sex == 1 & str_detect(Cruise, "02$") &
                              log(CH_CC) >= -1.916947 + 1.070620 * log(CW)) ~ 1)) %>% #mature male NBS
                  rename_with(tolower) -> cond_mat 
                                                            
#############
## Pull haul and snow crab specimen data from crabpack ----
snow <- get_specimen_data(species = "SNOW",
                          region = "EBS",
                          channel = "KOD")

snow_nbs <- get_specimen_data(species = "SNOW",
                              region = "NBS",
                              channel = "KOD")

ebs_haul <- snow$haul
nbs_haul <- snow_nbs$haul

#Join haul and biometric datasets 
ebs_haul %>%
  bind_rows(nbs_haul) %>% 
  rename_with(tolower) %>%
  right_join(cond_mat, by = c("cruise", "vessel", "haul")) -> mat_haul
#Note the NAs for 202502 because NBS survey not done yet at run time! 

#Add in sampling regions associated with each station
  #read in lookup table
regions <- read.csv("./data/regions_lookup.csv") 

#join and write csv
mat_haul %>%
  left_join(regions, by="station_id") %>%
#calculate additional WWT:DWT/FA metrics
  mutate(DWT_WWT = hepato_dwt/hepato_wwt,
       Perc_DWT = DWT_WWT*100,
       WWT_DWT = hepato_wwt/hepato_dwt) %>%
  write_csv(file="./output/condition_master.csv")






