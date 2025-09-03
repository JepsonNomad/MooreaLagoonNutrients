## This script generates Supplementary materials S2
## Comparison of Turbinaria %N analyzed at different labs

## Author: Christian John
## Date: 13 May 2025


#### Load packages ----
library(tidyverse)

#### Import data ----
ucsb = read_csv("~/Documents/Burkepile_Lab/Data/ATI/EDI_turbinaria/data/clean/turbinaria_CHN_UCSB.csv") %>%
  filter(FieldSeason == "2016-May") %>%
  rename(Nucsb = Percent_N) %>%
  select(Site, Nucsb)
geor = read_csv("~/Documents/Burkepile_Lab/Data/ATI/EDI_turbinaria/data/clean/turbinaria_CHN_UGA.csv") %>%
  filter(FieldSeason == "2016-May") %>%
  rename(Nuga = Percent_N) %>%
  select(Site, Nuga)

#### Data wrangling ----
## Fuse UGA and UCSB data
xcorrdat = geor %>%
  full_join(ucsb) %>%
  filter(complete.cases(.))

#### Data analysis ----
## Test correlation among UCSB and UGA data
xcorrstats = cor.test(xcorrdat$Nucsb, xcorrdat$Nuga)
xcorrstats$statistic
xcorrstats$p.value
xcorrstats$estimate

#### Visualize results ----
## Plot correlation and annotate statistics
comparisonPlot = xcorrdat %>%
  ggplot(aes(x = Nuga, y = Nucsb)) +
  geom_abline() +
  geom_point(col = "grey20") +
  stat_smooth(method = "lm", se = F, col = "black") +
  annotate("text", x = 0.8, y = 0.43, 
                label = paste0("r = ", round(xcorrstats$estimate, 2)),
                size = 4.5, label.size = 0, fill = "transparent") +
  annotate("text", x = 0.8, y = 0.4, 
                label = paste0("p = ", round(xcorrstats$p.value, 6)), 
                size = 4.5, label.size = 0, fill = "transparent") +
  xlim(min(c(xcorrdat$Nuga, xcorrdat$Nucsb))-0.05,
       max(c(xcorrdat$Nuga, xcorrdat$Nucsb))+0.05) +
  ylim(min(c(xcorrdat$Nuga, xcorrdat$Nucsb))-0.05,
       max(c(xcorrdat$Nuga, xcorrdat$Nucsb))+0.05)  +
  xlab(expression(paste(italic(Turbinaria), " %N (UGA)"))) +
  ylab(expression(paste(italic(Turbinaria), " %N (UCSB)"))) +
  CJsBasics::BasicTheme +
  theme(plot.margin = margin(0.15,0.15,0.15,0.15, unit = "in"))

#### Export ----
ggsave("plots/supplement/s2_xcorr_TurbN.jpg",
       comparisonPlot, width = 4, height = 4, units = "in", dpi = 300)
