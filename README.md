# Snow-Crab-Condition-Monitoring
This dataset includes immature snow crab sampled for hepatopancreas annually on the NOAA eastern and northern Bering Sea bottom trawl surveys. Water content of the hepatopancreas has been validated as a rapid field metric for energetic condition in juvenile snow crab, therfore collecting samples annually enables current year monitoring of energetic condition. This metric is reported annually in the Snow Crab Ecosystem and Socioeconomic Profile, and serves as a critical indicator for stock health, recruitment success, and survival to the terminal molt. 

# METADATA:

Bibliographic
-------------

| Published     | 02/02/2023   |
| ------------- | ------------- |
| Keywords      | snow crab |
|               |   Bering Sea |
|               |  fatty acids |
|               | energetic condition |
|               | marine heatwave |

 

Coverage
--------

### Temporal

| Begin    | 2019-06-01 |
| ------------- | ------|
| End   | 2024-08-31 |

 
 

### Spatial

| LME     |                     |
| ------------- | ------------- |
|                | Eastern Bering Sea |
|               |   Northern Bering Sea |

![Rplot](https://github.com/Erin-Fedewa-NOAA/Snow-Crab-Condition/blob/main/figures/data%20exploration/n_year.png?raw=true)



Attributes
----------
One master datasets has been produced for further modeling via the "append haul" script. "condition_master.csv" includes all snow crab biometric, energetic condition and haul level data, with data attributes listed below. 

| Name    |    Description   |   Unit    |
| ------- | ---------------- | ---------- |
| `year`     |        Year of specimen collection | numeric
|  `vial_id`    |    Unique ID for hepatopancreas sample by year. AKK=Alaska Knight samples, V=Vesteraalen samples   |   numeric/text
| `cruise` | Cruise ID for Bering Sea bottom trawl surveys. YYYY-01 indicates EBS surveys, YYYY-02 indicates NBS surveys. See RACEBASE or AKFIN for additional NOAA cruise metadata  |   ID Code
| `gis_station`   | Alpha-numeric designation for the station established in the design of AFSC standardized surveys | numeric/text   
| `area_swept`   |   Unit of effort for AFSC bottom trawl surveys: computed by distance towed*mean net width   | numeric, in ha
|  `cpue`   |   Station-level snow crab density, calculated as CPUE"   |   numeric, crab/nmi^2
|  `vessel`  |     ID number of the vessel used to collect data for that haul associated with vessel name    |   numeric
|  `haul`      |  Uniquely identifies a sampling event (haul) within an AFSC cruise. It is a sequential number, in chronological order of occurrence |  numeric
|  `sex` | Sex of specimen sampled. 1=Male, 2-Female   |  numeric
|  `cw`  |  Carapace width of specimen sampled | numeric, in mm
|  `ch_cc`   |   Chela height (males) or clutch code (females) used to determine maturity. 0 clutch code for females signifies 000, or immature   | numeric, in mm
 | `crab_wgt`  |   Whole crab weight of specimen sampled (prior to dissection). Only individuals with no missing limbs weighed. Weights were eliminated in 2024 | numeric, in g
|  `hepato_wwt`    | Total wet weight of hepatopancreas (subtracted from weight of weigh boat)  | numeric, in g
|  `hepato_dwt`    |Total dry weight of hepatopancreas after drying at 70C for 48hrs (subtracted from weight of weigh boat)  | numeric, in g
|  `maturity`  |  Maturity of specimen sampled. 0=Immature, 1=Mature"  | numeric                                                                                        
| `mid_latitude`       |   Latitude of specimen collection. Designates latitude at start of haul for AFSC standardized surveys    | numeric
|  `mid_longitude`    | Longitude of specimen collection. Designates longitude at start of haul for AFSC standardized surveys | decimal degree
 | `bottom_depth`    |    Bottom depth at station for AFSC standardized surveys  | numeric, in m
|  `gear_temperature`   |    Bottom temperature at sampling station | degree C
 | `start_date`     |        Date of sampling | date, month/day/year
 | `sample_region`     |        Pre-defined spatial sampling strata for EBS and NBS bottom trawl survey hepatopancreas collections | numeric
|  `lme` |   Large marine ecosystem: EBS=eastern Bering Sea, NBS=northern Bering Sea   |  text
|  `DWT_WWT` |   hepatopancreas dry weight/hepatopancreas wet weight   |  numeric, in g
|  `Perc_DWT` |   hepatopancreas dry weight/ hepatopancreas wet weight x 100; percentage hepatopancreas dry weight    |  percentage
|  `WWT_DWT` |   hepatopancreas wet weight/dry weight    |  numeric, in g

:::

 
   # NOAA License
This repository is a scientific product and is not official communication of the National Oceanic and Atmospheric Administration, or the United States Department of Commerce. All NOAA GitHub project code is provided on an ‘as is’ basis and the user assumes responsibility for its use. Any claims against the Department of Commerce or Department of Commerce bureaus stemming from the use of this GitHub project will be governed by all applicable Federal law. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by the Department of Commerce. The Department of Commerce seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by DOC or the United States Government.

Software code created by U.S. Government employees is not subject to copyright in the United States (17 U.S.C. §105). The United States/Department of Commerce reserve all rights to seek and obtain copyright protection in countries other than the United States for Software authored in its entirety by the Department of Commerce. To this end, the Department of Commerce hereby grants to Recipient a royalty-free, nonexclusive license to use, copy, and create derivative works of the Software outside of the United States.
