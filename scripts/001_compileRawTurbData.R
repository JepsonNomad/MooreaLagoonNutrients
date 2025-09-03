## This script collates Turbinaria %N data from UGA and UCSB labs
## into one long-format dataset for analysis in the lagoon nutrients project

## Author: Christian John
## Date: 13 May 2025


#### Load packages ----
library(tidyverse) # Data wrangling

#### Define parameters ----
## Working with rainy season data only
samplingSeason = "WET"

if(samplingSeason == "ALL"){
  samplingDates = c("2016-January", "2016-May", "2016-July", 
                    "2020-August", "2021-May", "2022-April",
                    "2023-April")
}else if(samplingSeason == "WET"){
  samplingDates = c("2016-May", "2021-May", "2022-April", "2023-April")
}else{
  break("Can't understand sampling instructions! Did you choose a sampling season?")
}

#### Import data ----
## University of Georgia lab results
turb_UGA = read_csv("~/Documents/Burkepile_Lab/Data/ATI/EDI_turbinaria/data/clean/turbinaria_CHN_UGA.csv")

## UC Santa Barbara lab results
turb_UCSB = read_csv("~/Documents/Burkepile_Lab/Data/ATI/EDI_turbinaria/data/clean/turbinaria_CHN_UCSB.csv") %>%
  ## The 2016 May data were for validation only at UCSB and the full 2016 May
  ## sampling was processed at UGA. Remove these 20 validation points.
  filter(FieldSeason != "2016-May")

#### Data wrangling ----
## Bind the two labs' datasets by row. Select relevant columns
## and retain only data from field seasons of interest.
## Condense data 
turb_compiled = plyr::rbind.fill(turb_UGA, turb_UCSB) %>%
  select(Site, Date, Percent_N, FieldSeason) %>%
  ## Select just relevant seasonal data
  filter(FieldSeason %in% samplingDates) %>%
  arrange(FieldSeason, Site) %>%
  ## Site 080 in 2021 had two Turbinaria measurements on the same day
  ## Take the mean of these two so that there's one measurement per
  ## site per season:
  group_by(Site, FieldSeason) %>%
  summarize(across(everything(), mean)) %>%
  ungroup()

#### Write output ----
## Save the result to a tidy .csv file
write_csv(turb_compiled,
          "data/generated/turb_N_compiled.csv")


