## This script generates Figure 4
## Impacts of recent precipitation on lagoon nutrients

## Author: Christian John
## Date: 13 May 2025


#### Load Libraries ----
library(sf)
library(tidyverse)
library(biscale)
library(cowplot)
library(lmerTest)
library(RcppRoll)
library(visreg)
source("scripts/000_reefCols.R")


#### Import data ----
## Turbinaria data
## (one measure per site per year)
dataTurbInSitu <- read_csv("data/generated/turb_N_imputed.csv")

## Detrended Turbinaria (accounting for distance from shore)
## (one measure per site)
dataTurbDetrended <- read_csv("data/generated/detrendedNutrients.csv") %>%
  rename(Nsite = Percent_N) %>%
  select(Site, Nsite, detrended)

## Water data
dataWater<-read_csv("data/generated/water_N_compiled.csv")

## Site metadata
dataSitemeta <- read_csv("~/Documents/Burkepile_Lab/Data/ATI/EDI_site_locations/data/clean/sites_metadata_DRAFT.csv")

## MeteoFrance rain data
dataRainAll<-read_csv("../../MeteoFrance/data/CJ_2025-01-06/MooreaPrecipitation.csv")
rainStations<-read_csv("../../MeteoFrance/data/CJ_2025-01-06/stationMetadata.csv") %>%
  select(-geometry)

## Moorea outline
moo <- st_read("~/Documents/Burkepile_Lab/Moorea_mapping_quicklayers/data/Moorea_LOC_ILE.shp") %>%
  st_transform(4326)

dataRainLocs = dataRainAll %>%
  select(NUM_POSTE,NOM_USUEL,LAT,LON) %>%
  group_by(NUM_POSTE,NOM_USUEL) %>%
  summarize(LAT = mean(LAT),
            LON = mean(LON)) %>%
  st_as_sf(coords = c("LON","LAT")) %>%
  st_set_crs(4326)

ggplot(dataRainLocs) +
  geom_sf(data  = moo) +
  geom_sf(aes(col = NOM_USUEL))

## > Supplementary materials s1 ----
## Count sites in each year for supplement
read_csv("data/generated/turb_N_compiled.csv") %>%
  group_by(FieldSeason) %>%
  summarize(nSites = n())
# 1 2016-May       171
# 2 2021-May       193
# 3 2022-April     190
# 4 2023-April     198
dataWater %>%
  group_by(FieldSeason) %>%
  summarize(nSites = n())
# FieldSeason nSites
# <chr>        <int>
#   1 2021-May       195
# 2 2022-April     191
# 3 2023-April     200


#### Data wrangling ----
## Wrangle precipitation data
dataRainSelect = dataRainAll %>%
  filter(NUM_POSTE %in% rainStations$NUM_POSTE) %>%
  left_join(rainStations) %>%
  mutate(Date = ymd(AAAAMMJJ)) %>%
  arrange(Date, NUM_POSTE) 

sampleDates = unique(dataWater$Date)
unique(dataRainSelect$NUM_POSTE)
unique(dataRainSelect$NOM_USUEL)
rainSummary = dataRainSelect %>%
  mutate(Date = ymd(AAAAMMJJ)) %>%
  arrange(Date, NUM_POSTE) %>%
  filter(Date > as.Date("2021-01-01"),
         Date < as.Date("2024-01-01")) %>%
  group_by(Date) %>%
  summarize(nStat = n(),
            Rain = mean(RR))
rainSummary

rainSum = rainSummary %>%
  mutate(recentRain2 = roll_sum(x = Rain, 2, align = "right", fill = NA),
         recentRain3 = roll_sum(x = Rain, 3, align = "right", fill = NA),
         recentRain5 = roll_sum(x = Rain, 5, align = "right", fill = NA),
         recentRain10 = roll_sum(x = Rain, 10, align = "right", fill = NA),
         recentRain30 = roll_sum(x = Rain, 30, align = "right", fill = NA),
         recentRain60 = roll_sum(x = Rain, 60, align = "right", fill = NA))

rainSum %>%
  ggplot() +
  geom_vline(xintercept = sampleDates, col = "grey") +
  # geom_line(aes(x = Date, y = Rain)) +
  geom_line(aes(x = Date, y = recentRain3), col = "dodgerblue1", lwd = 1.5) +
  geom_line(aes(x = Date, y = recentRain10), col = "dodgerblue2", lwd = 1.25) +
  geom_line(aes(x = Date, y = recentRain30), col = "dodgerblue3", lwd = 1) +
  geom_line(aes(x = Date, y = recentRain60), col = "dodgerblue4", lwd = 0.75) +
  ylab("Recent rainfall\n(3, 10, 30, and 60-day cumulative, mm)") +
  CJsBasics::BasicTheme

write_csv(x = rainSum,
          file = "data/generated/rainSummaries.csv")

## > Combine datasets ----
## Fuse water and turbinaria data
data_raw <- left_join(dataWater, dataTurbInSitu) %>%
  left_join(dataTurbDetrended) %>%
  left_join(dataSitemeta)

## Filter to just the sites and times with water column data
data_filtered <- data_raw %>%
  filter(complete.cases(Phosphate, Silicate, Nitrite_plus_Nitrate)) %>%
  select(-Ammonia) %>% ## Remove ammonia which had QAQC issues
  left_join(rainSum, by = "Date") %>%
  mutate(Habitat = factor(Habitat,
                          levels = names(reefColors)))

#### Direct comparisons ----
data_filtered %>%
  group_by(Site, Distance_from_shore, Cleared_area) %>%
  summarize(meanN = mean(Nsite)) %>%
  ggplot(aes(x = Distance_from_shore,
             y = meanN)) +
  geom_point(col = "black", size = 2.5) +
  geom_point(aes(col = Cleared_area), size = 2) +
  stat_smooth(method = "lm") +
  scale_color_viridis_c() +
  ylab("Mean Turbinaria %N") +
  CJsBasics::BasicTheme

data_filtered %>%
  ggplot(aes(x = Distance_from_shore,
             y = Percent_N)) +
  geom_point(col = "black", size = 2.5) +
  geom_point(aes(col = Cleared_area), size = 2) +
  # stat_smooth(method = "lm") + 
  stat_smooth(method = "lm",
              aes(group = Year),
              col = "black") + # Pattern is consistent across years
  scale_color_viridis_c() +
  ylab("Instantaneous Turbinaria %N") +
  CJsBasics::BasicTheme

summary(lm(Nitrite_plus_Nitrate ~ Distance_from_shore + recentRain3,
           data = data_filtered))
## Recent rain elevates water column N across the lagoon.

#### Exploring over gradients of land clearing ----
lmerH2O_3 = lmer(Nitrite_plus_Nitrate ~ Cleared_area + recentRain3 + (1|Site),
                 data = data_filtered)
## Rain over the last couple days elevates water column N across the lagoon.
lmerH2O_30 = lmer(Nitrite_plus_Nitrate ~ Cleared_area + recentRain30 + (1|Site),
                  data = data_filtered)
## Rain over the last month reduces water column N
lmerTurb_3 = lmer(Percent_N ~ Cleared_area + recentRain3 + (1|Site),
                  data = data_filtered)
## Rain over the last couple days reduces Turb N, which is higher near watersheds with a lot of cleared area
lmerTurb_30 = lmer(Percent_N ~ Cleared_area + recentRain30 + (1|Site),
                   data = data_filtered)
## Rain over the last month increases Turb N

bH2O_3 = summary(lmerH2O_3)$coefficients[3,1] %>% signif(2)
pH2O_3 = summary(lmerH2O_3)$coefficients[3,5] %>% signif(2)
bH2O_30 = summary(lmerH2O_30)$coefficients[3,1] %>% signif(2)
pH2O_30 = summary(lmerH2O_30)$coefficients[3,5] %>% signif(2)
bTurb_3 = summary(lmerTurb_3)$coefficients[3,1] %>% signif(2)
pTurb_3 = summary(lmerTurb_3)$coefficients[3,5] %>% signif(2)
bTurb_30 = summary(lmerTurb_30)$coefficients[3,1] %>% signif(2)
pTurb_30 = summary(lmerTurb_30)$coefficients[3,5] %>% signif(2)

#### Investigating effects of timescale ----
## For many time-scales (i.e. 1 day ago through and including 60 days ago),
## Calculate recent cumulative rain
timeWindow = c(1:60)
rainSumList = lapply(timeWindow,
                     FUN = function(x){
                       rainSummary %>%
                         mutate(recentRain = roll_sum(x = Rain, 
                                                      x, 
                                                      align = "right", 
                                                      fill = NA),
                                timeScale = x)
                     })
rainSumWindow = do.call("rbind",
                        rainSumList)

## For each time-scale, construct a linear mixed effects model with 
## nutrients as a response to cumulative precipitation over that time 
## scale and cleared land. Extract the slope coefficient and standard 
## error of the slope for each linear model.
rainComp = lapply(
  timeWindow,
  FUN = function(x){
    ## Create a data.frame with lagoon nutrients, 
    ## recent precip, and cleared land
    lmDat = data_filtered %>%
      select(Site, Date, Cleared_area,
             Nitrite_plus_Nitrate, Percent_N) %>%
      left_join(rainSumWindow %>%
                  filter(timeScale == x), 
                by = "Date")
    
    ## Build model for Turbinaria N
    ## Fit linear mixed effects model
    timeScaleLMTurb = lmer(Percent_N ~ log(recentRain + 1) +
                             Cleared_area +
                             (1|Site),
                           data = lmDat)
    ## Can also do this with a simple lm and not account for site;
    ## These results are almost identical, with marginally
    ## wider confidence around mid-to-long time-scales
    # timeScaleLMTurb = lm(Percent_N ~ log(recentRain + 1) +
    #                       Cleared_area,
    #                     data = lmDat)
    ## Extract the coefficients and standard errors
    turbrainFx = summary(timeScaleLMTurb)$coefficients[2,1]
    turbrainFxSE = summary(timeScaleLMTurb)$coefficients[2,2]
    turbclearedFx = summary(timeScaleLMTurb)$coefficients[3,1]
    turbclearedFxSE = summary(timeScaleLMTurb)$coefficients[3,2]

    ## Repeat for water column nutrients    
    ## Fit linear mixed effects model
    timeScaleLMH20 = lmer(Nitrite_plus_Nitrate ~ log(recentRain + 1) +
                            Cleared_area +
                            (1|Site),
                          data = lmDat)
    ## Can also do this with a lm and not account for site
    # timeScaleLMH20 = lm(Nitrite_plus_Nitrate ~ log(recentRain + 1) +
    #                       Cleared_area,
    #                     data = lmDat)
    ## Extract the coefficients and standard errors
    h20rainFx = summary(timeScaleLMH20)$coefficients[2,1]
    h20rainFxSE = summary(timeScaleLMH20)$coefficients[2,2]
    h20clearedFx = summary(timeScaleLMH20)$coefficients[3,1]
    h20clearedFxSE = summary(timeScaleLMH20)$coefficients[3,2]
    
    ## Compile effect of cumulative recent precip on nutrients 
    outDF = data.frame("timeScale" = x,
                       "turbrainFx" = turbrainFx,
                       "turbrainFxSE" = turbrainFxSE,
                       "turbclearedFx" = turbclearedFx,
                       "turbclearedFxSE" = turbclearedFx,
                       "h20rainFx" = h20rainFx,
                       "h20rainFxSE" = h20rainFxSE,
                       "h20clearedFx" = h20clearedFx,
                       "h20clearedFxSE" = h20clearedFx)
    return(outDF)
  })

rainEffects = do.call("rbind",
                      rainComp)

#### Figure 4b ----
## Plot effect of rain on lagoon nutrients
## over many time-scales.
rainTurb = rainEffects %>%
  ggplot(aes(x = timeScale)) +
  geom_vline(aes(xintercept = 3),
             lty = 1,
             lwd = 2,
             col = "grey90") +
  geom_vline(aes(xintercept = 30),
             lty = 1,
             lwd = 2,
             col = "grey90") +
  geom_hline(aes(yintercept = 0),
             lty = 2,
             col = "grey40") +
  geom_ribbon(aes(ymin = turbrainFx - turbrainFxSE,
                  ymax = turbrainFx + turbrainFxSE),
              col = "green4",
              fill = "green4",
              alpha = 0.2,
              lwd = 0.1) +
  geom_line(aes(y = turbrainFx),
            col = "green4",
            lwd = 2) +
  coord_flip() +
  ylab(expression(atop("Effect of cumulative", 
                       "precip. on "*italic(Turbinaria)*" %N"))) +
  xlab("") +
  CJsBasics::BasicTheme
rainTurb

rainH20 = rainEffects %>%
  ggplot(aes(x = timeScale)) +
  geom_vline(aes(xintercept = 3),
             lty = 1,
             lwd = 2,
             col = "grey90") +
  geom_vline(aes(xintercept = 30),
             lty = 1,
             lwd = 2,
             col = "grey90") +
  geom_hline(aes(yintercept = 0),
             lty = 2,
             col = "grey40") +
  geom_ribbon(aes(ymin = h20rainFx - h20rainFxSE,
                  ymax = h20rainFx + h20rainFxSE),
              col = "dodgerblue2",
              fill = "dodgerblue2",
              alpha = 0.2,
              lwd = 0.1) +
  geom_line(aes(y = h20rainFx),
            col = "dodgerblue2",
            lwd = 2) +
  coord_flip() +
  ylab(expression(atop("Effect of cumulative",
                       "precip. on water col. N"))) +
  xlab("Timescale (days)") +
  scale_y_continuous(breaks = c(0,-0.1,-0.2)) +
  CJsBasics::BasicTheme
rainH20

#### ... Compile Fig 4a ----
precipEffectsPlot = cowplot::plot_grid(rainH20,
                                       rainTurb,
                                       nrow = 1)
precipEffectsPlot

#### Figure 4b ----
## Contrast plots examining effect of cleared area on lagoon nutrients
## at different scales of precipitation inputs
## Generate contrast plots for 3-day and 30-day time-scales
## for both Turbinaria N and water column N
regressionPlot1 = visreg(fit = lmerTurb_3, 
                         xvar = "Cleared_area",
                         # scale = "response",
                         type = "contrast",
                         partial=T,
                         line = list(col = "green4"),
                         point = list(col = "grey80"),
                         gg = T) +
  geom_hline(aes(yintercept = 0),
             lty = 2,
             lwd = 0.5,
             col = "grey60") +
  ylab(expr(paste(italic("Turbinaria"), " %N"))) +
  xlab(expression(paste("Cleared area (",km^2,")"))) +
  ggtitle("3 days") +
  CJsBasics::BasicTheme
regressionPlot2 = visreg(fit = lmerTurb_30, 
                         xvar = "Cleared_area",
                         # scale = "response",
                         type = "contrast",
                         partial=T,
                         line = list(col = "green4"),
                         point = list(col = "grey80"),
                         gg = T) +
  geom_hline(aes(yintercept = 0),
             lty = 2,
             lwd = 0.5,
             col = "grey60") +
  ylab(expr(paste(italic("Turbinaria"), " %N"))) +
  xlab(expression(paste("Cleared area (",km^2,")"))) +
  ggtitle("30 days") +
  CJsBasics::BasicTheme
regressionPlot3 = visreg(fit = lmerH2O_3, 
                         xvar = "Cleared_area",
                         # scale = "response",
                         type = "contrast",
                         partial=T,
                         line = list(col = "dodgerblue2"),
                         point = list(col = "grey80"),
                         gg = T) +
  geom_hline(aes(yintercept = 0),
             lty = 2,
             lwd = 0.5,
             col = "grey60") +
  ylab("Water column N") +
  xlab(expression(paste("Cleared area (",km^2,")"))) +
  ggtitle("3 days") +
  CJsBasics::BasicTheme
regressionPlot4 = visreg(fit = lmerH2O_30, 
                         xvar = "Cleared_area",
                         # scale = "response",
                         type = "contrast",
                         partial=T,
                         line = list(col = "dodgerblue2"),
                         point = list(col = "grey80"),
                         gg = T) +
  geom_hline(aes(yintercept = 0),
             lty = 2,
             lwd = 0.5,
             col = "grey60") +
  ylab("Water column N") +
  xlab(expression(paste("Cleared area (",km^2,")"))) +
  ggtitle("30 days") +
  CJsBasics::BasicTheme

## ... Compile Fig 4 b ----
visRegPlots = cowplot::plot_grid(regressionPlot3,
                                 regressionPlot4,
                                 regressionPlot1,
                                 regressionPlot2,
                                 nrow = 2,
                                 ncol = 2,
                                 labels = "",
                                 align = "hv")
visRegPlots

#### Compile Fig 4 a and b ----
Figure4 = cowplot::plot_grid(precipEffectsPlot,
                             visRegPlots,
                             nrow = 1,
                             labels = c("a","b"))

#### Export ----
ggsave("plots/figure4.jpg",
       Figure4,
       width = 11, height = 6, units = "in", dpi = 300)
ggsave("plots/figure4.png",
       Figure4,
       width = 11, height = 6, units = "in", dpi = 300)
ggsave("plots/figure4.tiff",
       Figure4,
       width = 11, height = 6, units = "in", dpi = 300)
