library(tidyverse)
library(phyloseq)

#### Turbinaria %N ----
read_csv("data/generated/turb_N_compiled.csv") %>%
  group_by(FieldSeason) %>%
  summarize(nSites = length(unique(Site)))


#### Water nutrients ----
read_csv("~/Documents/Burkepile_Lab/Data/ATI/EDI_water_nutrients/data/clean/lagoon_nutrients_long.csv") %>%
  group_by(Year, FieldSeason) %>%
  summarize(nSites = length(unique(Site)))


#### fDOM ----
read_csv("~/Documents/Burkepile_Lab/Data/ATI/EDI_fDOM/data/clean/watercolumn_fDOM.csv") %>%
  group_by(Year, FieldSeason) %>%
  summarize(nSites = length(unique(Site)))


#### Microbes ----
r2021 <- 
  readr::read_rds("~/Documents/Burkepile_Lab/Data/ATI/EDI_microbes/data/clean/ati-2021-physeq-rare.RDS")
r2022 <- 
  readr::read_rds("~/Documents/Burkepile_Lab/Data/ATI/EDI_microbes/data/clean/ati-2022-physeq-rare.RDS")
r2023 <- 
  readr::read_rds("~/Documents/Burkepile_Lab/Data/ATI/EDI_microbes/data/clean/ati-2023-physeq-rare.RDS")

c2021 <- microViz::samdat_tbl(r2021) %>% pull(Site) %>% unique() %>% length()
c2022 <- microViz::samdat_tbl(r2022) %>% pull(Site) %>% unique() %>% length()
c2023 <- microViz::samdat_tbl(r2023) %>% pull(Site) %>% unique() %>% length()

c(c2021,c2022,c2023)
