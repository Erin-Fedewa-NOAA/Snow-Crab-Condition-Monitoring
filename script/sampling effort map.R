#Script to create map of sampling effort overlayed with ice extent
#1) ERA 5 ice cover processing to add ice extent boundary to map
#2) Map of sampling locations

#Authors: E. Fedewa, E. Ryznar (ERA 5 ice data/mapping)

#To download and process ice cover data from ERA 5 monthly averaged data
# 1) Navigate here (will need to login): https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels-monthly-means?tab=overview
# 2) Click on "Download data" tab
# 3) Click on the product type you’d like. We have been using “Monthly averaged reanalysis”
# 4) For ice cover, click on the “Other” drop down, and then check the “Sea-ice cover” box. 
# 5) Select years you’d like the data to cover. Sometimes including the most recent year with earlier years 
#    in your selection results in an error, which can be resolved by downloading recent year data separately.
# 6) Select months you’d like the data to cover (Jan-Apr in this case) 
# 7) Select geographical area for the data. Subsetting to a specific region speeds up the data querying/download.
#    For the data below, the region has been North 66°, West -175°, South 52°, and East -155°. 
# 8) Select NetCDF(experimental) as the data format, and submit form to query/download data. 

# **NOTE** 
# For requests which return a mixture of ERA5 and ERA5T data  
# (such as for data from the 1st of the month), instantaneous variables (e.g temperature) come from ERA5T (which has 'experiment version'  of 5) 
# while accumulated variables (fluxes, precipitation) come from both datasets with the following structure:
# 00-06 UTC on 1 day of the month from ERA5 (expver 1)
# 07-23 UTC on 1 day of the month (and the following dates up to 5 day from present) from ERA5T (expver 5)

# When these data are converted to netCDF a new dimension is created called expver containing 1 and 5. 
#Moreover, a single time coordinate is used which covers the entire requested period.

### LOAD PACKAGES -------------------------------------------------------------------------------------------------------

library(tidyverse)
library(tidync)
library(sf)
library(terra)
library(akgfmaps)
library(ggridges)
library(patchwork)
library(hrbrthemes)
library(ggtext)
library(ggpubr)

#install.packages("remotes")
#remotes::install_github("afsc-gap-products/akgfmaps")

### PROCESS ICE DATA ----------------------------------------------------------
data <- tidync("./data/ERA5_ice_1975_2024.nc") %>% 
  hyper_tibble() %>% 
  separate(date, into = c("Year", "Month", "Time"), sep = c(4,6)) %>%
  select(-Time) %>%
  mutate(Year = as.numeric(Year),
         Month = as.numeric(Month))

## SET COORDINATE REFERENCE SYSTEMS (CRS) --------------------------------------
in.crs <- "+proj=longlat +datum=NAD83" #CRS is in lat/lon
map.crs <- "EPSG:3338" # final crs for mapping/plotting: Alaska Albers

## LOAD SHELLFISH ASSESSMENT PROGRAM GEODATABASE -------------------------------
survey_gdb <- "./data/SAP_layers" 
survey_strata <- terra::vect(survey_gdb, layer = "EBS.NBS_surveyarea")
#EBS/NBS Boundary line
boundary <- st_read(layer = "EBS_NBS_divide", survey_gdb) 

## LOAD CONDITION AND SURVEY DATA ---------------------------------------------
condition_master <- read.csv("./output/condition_master.csv")
ebs_haul <- read.csv("./data/crabhaul_opilio.csv")
nbs_haul <- read.csv("./data/crabhaul_opilio_nbs.csv")
ebs_strata <- read.csv("./data/crabstrata_opilio.csv")
nbs_strata <- read.csv("./data/crabstrata_opilio_nbs.csv")

## LOAD ALASKA REGION LAYERS (FROM AKGFMAPS R package) -----------------------------------
ebs_layers <- akgfmaps::get_base_layers(select.region = "ebs", set.crs = "EPSG:3338")
ebs_survey_areas <- ebs_layers$survey.area
ebs_survey_areas$survey_name <- c("Eastern Bering Sea", "Northern Bering Sea")

#Survey areas plot
ggplot() +
  geom_sf(data = ebs_layers$akland) +
  geom_sf(data = ebs_survey_areas, mapping = aes(fill = survey_name)) +
  scale_x_continuous(limits = ebs_layers$plot.boundary$x,
                     breaks = ebs_layers$lon.breaks) +
  scale_y_continuous(limits = ebs_layers$plot.boundary$y,
                     breaks = ebs_layers$lat.breaks) +
  scale_fill_viridis_d(name = "Survey") +
  theme_bw()

#EBS/NBS shelf Survey grid plot
ggplot() +
  geom_sf(data = ebs_layers$akland) +
  geom_sf(data = ebs_layers$survey.grid, fill = NA) +
  geom_sf(data= boundary, linewidth = 2) +
  scale_x_continuous(limits = ebs_layers$plot.boundary$x,
                     breaks = ebs_layers$lon.breaks) +
  scale_y_continuous(limits = ebs_layers$plot.boundary$y,
                     breaks = ebs_layers$lat.breaks) +
  theme_bw()

###  MAP  -----------------------------------------------------------

#Filter for sea ice extent threshold 
#We'll use 15% as our threshold for sea ice extent from estimates of sea ice 
#concentration (e.g. concentration >= 0.15 is ice-covered)- per NSIDC
data %>%
  filter(Year > 2018, Month == 03, Year != 2020) %>%
  filter(siconc >= 0.15) %>%
  #there are some ice-covered areas in the AI that really drive the southern 
  #bound of the hull so we'll eliminate these few points manually 
  filter(latitude > 55.5) %>%
  rename(year = Year) %>%
  #transform sea ice data into spatial data frame
  sf::st_as_sf(coords = c("longitude", "latitude"), crs = in.crs) %>%
  sf::st_transform(sf::st_crs(map.crs)) %>%
  vect(.) %>%
  mask(., survey_strata) %>% #this bounds ice extent data to survey region
  sf::st_as_sf() -> ice_extent  

#Transform crab data into spatial data frame
condition_master %>% 
  group_by(year, mid_latitude, mid_longitude) %>%
  summarise(n_crab=n()) %>%
  # Convert lat/long to an sf object
  st_as_sf(coords = c("mid_longitude", "mid_latitude"), crs = st_crs(4326)) %>%
  #st_as_sf needs crs of the original coordinates- need to transform to Alaska Albers
  st_transform(crs = st_crs(3338)) -> crab_dat

#Use the spatial data frame to generate a convex hull around the data extent
ice_hull  <-  st_simplify(st_buffer(st_convex_hull(st_union(st_geometry(ice_extent))), 
                                    dist = 15000), dTolerance = 5000)

#Add hull, ice extent layer and crab data to map
ggplot() +
  geom_sf(data = ebs_layers$survey.grid, fill=NA, color=alpha("grey80"))+
  geom_sf(data = ebs_survey_areas, fill = NA) +
  geom_sf(data = ebs_layers$akland, fill = "grey80", color = "black") +
  #add hull for sea ice extent- though we won't use this approach here....
  #geom_sf(data = ice_hull,
  #fill = NA,
  #color = alpha("red", 0.85),
  #linewidth = 1) +
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
