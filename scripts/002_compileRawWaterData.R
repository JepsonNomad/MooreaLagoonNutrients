## This script compiles water column chemistry data (nutrients + fDOM)
## into one long-format dataset for analysis in the lagoon nutrients project

## Author: Christian John
## Date: 13 May 2025


#### Load packages ----
library(tidyverse) 


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
## Water column nutrients
water_nutrients = read_csv("~/Documents/Burkepile_Lab/Data/ATI/EDI_water_nutrients/data/clean/lagoon_nutrients_long.csv") %>%
  filter(FieldSeason %in% samplingDates) %>%
  ## Sites 004 and 061 in 2022 had extreme N+N values and are likely errors
  ## They were over 7 standard deviations above the mean.
  ## Remove them:
  mutate(Nitrite_plus_Nitrate = case_when(Site == "004" & Year == 2022 ~ NA,
                                          Site == "061" & Year == 2022 ~ NA,
                                          TRUE ~ Nitrite_plus_Nitrate)) %>%
  ## Site 043 in 2022 had an extreme Phosphate value that is likely an error
  ## It was approx 12 standard deviations above the mean.
  ## Remove it:
  mutate(Phosphate = case_when(Phosphate > 0.6 ~ NA,
                               TRUE ~ Phosphate))

## Water column fDOM
water_fdom = read_csv("~/Documents/Burkepile_Lab/Data/ATI/EDI_fDOM/data/clean/watercolumn_fDOM.csv") %>%
  filter(FieldSeason %in% samplingDates)


#### Data wrangling ----
## Combine the water column nutrient and fdom data
water_chemistry = full_join(water_nutrients,
                            water_fdom,
                            by = c("Site", "FieldSeason",
                                   "Date",
                                   "Year","Month","Day"))

## Check for duplicates
sum(duplicated(paste0(water_chemistry$Site, water_chemistry$Year)))


#### Write output ----
## Save the result to a tidy .csv file
write_csv(water_chemistry,"data/generated/water_N_compiled.csv")
