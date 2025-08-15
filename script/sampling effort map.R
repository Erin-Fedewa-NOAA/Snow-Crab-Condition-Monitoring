#Script to create map of sampling effort overlayed with ice extent
#1) ERA 5 ice cover processing to add ice extent boundary to map
#2) Map of sampling locations

#2026 to do: add a map with sampling locations overlaid on summer bottom temps

# load ----
library(tidyverse)
library(tidync)
library(lubridate)
library(magrittr)
library(akgfmaps)

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

#LOAD CRAB CONDITION DATA ----------------------------
condition_master <- read.csv("./output/condition_master.csv")

###########################################################
#To download and process ice cover data from ERA 5 monthly averaged data
# 1) Navigate here (will need to login): https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels-monthly-means?tab=overview
# 2) Click on "Download data" tab
# 3) Click on the product type you’d like. We have been using “Monthly averaged reanalysis”
# 4) For ice cover, click on the “Other” drop down, and then check the “Sea-ice cover” box. 
# 5) Select years you’d like the data to cover. Initial pull has been divided into two time stanzas, as
#a larger request can result in an error with processing 
# 6) Select months you’d like the data to cover (Jan-Apr in this case) 
# 7) Select geographical area for the data. For the data below, the region has been set to:
#North 65°, West -178°, South 56°, and East -165°. 
# 8) Select NetCDF(experimental) as the data format and unarchived as download format, and submit form to query/download data. 

################################################################
### Process Ice Data
tidync("./Data/ERA5_ice_2000_2025.nc") %>% 
  hyper_tibble() %>% 
  separate(valid_time, into = c("Year", "Month", "Time"), sep = "-") %>%
  select(-Time) %>%
  mutate(Year = as.numeric(Year),
         Month = as.numeric(Month)) %>%
  #Filter for sea ice extent threshold 
  #We'll use 15% as our threshold for sea ice extent from estimates of sea ice 
  #concentration (e.g. concentration >= 0.15 is ice-covered)- per NSIDC
  filter(siconc >= 0.15) %>%
  filter(Year > 2018, Month == 03, Year != 2020) %>%
  rename(year = Year) %>%
  #transform sea ice data into spatial data frame
  sf::st_as_sf(coords = c("longitude", "latitude"), crs = in.crs) %>%
  sf::st_transform(sf::st_crs(map.crs)) %>%
  vect(.) %>%
  mask(., survey_strata) %>% #this bounds ice extent data to survey region
  sf::st_as_sf() -> ice_extent  

#Transform crab data into spatial data frame
condition_master %>% 
  filter(cruise != 202502) %>%   #eliminating for now since we don't have lat/lons
  group_by(year, mid_latitude, mid_longitude) %>%
  summarise(n_crab=n()) %>%
  # Convert lat/long to an sf object
  st_as_sf(coords = c("mid_longitude", "mid_latitude"), crs = st_crs(4326)) %>%
  #st_as_sf needs crs of the original coordinates- need to transform to Alaska Albers
  st_transform(crs = st_crs(3338)) -> crab_dat

#Add ice extent layer and crab data to map
ggplot() +
  geom_sf(data = ebs_layers$survey.grid, fill=NA, color=alpha("grey80"))+
  geom_sf(data = ebs_survey_areas, fill = NA) +
  geom_sf(data = ebs_layers$akland, fill = "grey80", color = "black") +
  #add ice extent
  geom_sf(data=ice_extent , aes(), color = "#9ECAE1", alpha = 0.25 ) +
  #add crab sampling
  geom_sf(data=crab_dat, aes(size = n_crab), color = "grey30", alpha = .6) +
  geom_sf(data= boundary, linewidth = 1, color = "grey40") +
  scale_x_continuous(limits = ebs_layers$plot.boundary$x,
                     breaks = ebs_layers$lon.breaks) +
  scale_y_continuous(limits = ebs_layers$plot.boundary$y,
                     breaks = ebs_layers$lat.breaks) +
  scale_size_continuous(range = c(1,4)) +
  theme_bw() +
  facet_wrap(~year) +
  labs(x="", y="", size = expression(paste("Snow crab \n samples"))) +
  theme(legend.position="bottom") +
  guides(size = guide_legend(theme = theme(
    legend.title = element_text(size = 9)))) +
  theme(plot.margin = margin(0,-5,0,-5)) +
  theme(axis.text=element_text(size=8)) -> map
ggsave(map, file = "./figures/sample_site_map.png", dpi=300, width = 7.5, height = 7, units = "in")
