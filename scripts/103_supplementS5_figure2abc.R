## This script generates panels a-c for Supplementary materials S3
## Repeating Fig 2 results with non-interpolated data

## Author: Christian John
## Date: 14 May 2025

#### Import packages ----
library(tidyverse) # Data wrangling
library(cowplot) # Plot wrangling
library(lme4)
library(lmerTest)
library(emmeans)
library(multcomp)
library(multcompView)
library(GGally)
source("scripts/000_reefCols.R")

select = dplyr::select
filter = dplyr::filter

set.seed(11091989) # Berlin wall comes down


#### Import data ----
## ... Site locations ----
ati = read_csv("~/Documents/Burkepile_Lab/Data/ATI/EDI_site_locations/data/clean/sites_metadata_DRAFT.csv") %>%
  select(Site, Latitude, Longitude, Habitat, Island_shore) %>% 
  mutate(Habitat = factor(Habitat, levels = c("Bay","Fringing reef",
                                              "Mid lagoon", "Back reef",
                                              "Reef crest", "Reef pass")),
         Island_shore = case_when(Island_shore == "N" ~ "North shore",
                                  Island_shore == "SW" ~ "Southwest shore",
                                  Island_shore == "SE" ~ "Southeast shore"))

## ... Raw water column nutrient data ----
water_N = read_csv("data/generated/water_N_compiled.csv")
## Check the distribution of water N values
water_N %>% left_join(ati) %>%
  select(Site, Year, Nitrite_plus_Nitrate) %>%
  rename("N" = Nitrite_plus_Nitrate) %>%
  mutate(N_scal = scale(N)[,1]) %>%
  ggplot(aes(x = N_scal)) +
  geom_histogram(bins=100)

## ... Raw turbinaria nutrient data
sitesWithMissingN = read_csv("data/generated/turbSitesWithMissingN.csv")$Site
turb_N = read_csv("data/generated/turb_N_compiled.csv") %>% 
  mutate(Year = year(lubridate::ym(FieldSeason))) %>%
  select(Site, Year, FieldSeason, Percent_N) %>%
  ## Note: remove sites that had interpolation ----
  filter(!(Site %in% sitesWithMissingN))

## ... Compile turb and water data ----
compiled_N = water_N %>%
  full_join(turb_N, by = c("Site", "FieldSeason", "Year")) %>%
  left_join(ati, by = "Site")


## Analysis of habitat effects on lagoon nutrients ----
#### Figure 2a ----
## ... lmer of Turb N ----
## Turbinaria %N significantly varied by habitat type (conditional R2 = 0.38), with the highest %N at fringing reef sites and in bays and their passes, and lower at mid-lagoon and back reef sites (Fig 2a). 
turb_diffs_lmer = lmer(Percent_N ~ Habitat + (1|Year), data = compiled_N)
summary(turb_diffs_lmer)
MuMIn::r.squaredGLMM(turb_diffs_lmer)
em_turb_lmer = emmeans(turb_diffs_lmer, "Habitat")
contr_turb_lmer = contrast(em_turb_lmer, "pairwise")
plot(turb_diffs_lmer)
turb_lmer_cld = cld(em_turb_lmer, Letters = letters, reversed=T) %>%
  dplyr::select(Habitat, emmean, .group) %>%
  mutate(label = gsub(" ","",.group))

## ... Plot ----
Fig2a = compiled_N %>%
  group_by(Site, Habitat) %>%
  summarize(N = mean(Percent_N, na.rm = T)) %>% 
  ggplot() +
  geom_boxplot(aes(x = Habitat, 
                   y = N, 
                   col = Habitat),
               outliers = F,
               show.legend = F) +
  geom_jitter(aes(x = Habitat,
                  y = N,
                  col = Habitat),
              shape = 1,
              cex = 1.5,
              stroke = 0.5,
              width = 0.2,
              show.legend = F) +
  geom_label(data = turb_lmer_cld,
             aes(x = Habitat, label = label),
             y = 0.3725,
             label.size = 0, 
             size = 4.5,
             alpha = 0) +
  scale_x_discrete(labels = c("Bay","Frng. reef", "Mid lagoon", "Back reef", "Reef pass")) +
  scale_y_continuous(limits = c(0.365,NA)) +
  scale_color_manual(values = reefColors) +
  ylab(expression(paste(italic("Turbinaria"), " %N"))) +
  xlab("") +
  CJsBasics::BasicTheme +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title.y = element_text())
Fig2a


#### Figure 2b ----
## ... lmer of H2O N ----
## Water column chemistry also varied by habitat (conditional R2 = 0.45, Supplementary materials S5). Water column nitrites + nitrates showed a similar pattern of variation, with highest water column N in bays, followed by fringing reef, reef pass, and finally back reef and reef crest sites (Fig 2b). 
H2O_diffs_lmer = lmer(Nitrite_plus_Nitrate ~ Habitat + (1|Year), data = compiled_N)
summary(H2O_diffs_lmer)
MuMIn::r.squaredGLMM(H2O_diffs_lmer)
em_H2O_lmer = emmeans(H2O_diffs_lmer, "Habitat")
contr_H2O_lmer = contrast(em_H2O_lmer, "pairwise")
plot(H2O_diffs_lmer)
H2O_lmer_cld = cld(em_H2O_lmer, Letters = letters, reversed=T) %>%
  dplyr::select(Habitat, emmean, .group) %>%
  mutate(label = gsub(" ","",.group))
## ... Plot ----
Fig2b = compiled_N %>%
  group_by(Site, Habitat) %>%
  summarize(N = mean(Nitrite_plus_Nitrate, na.rm = T)) %>%
  ggplot() +
  geom_boxplot(aes(x = Habitat, 
                   y = N, 
                   col = Habitat),
               outliers = F,
               show.legend = F) +
  geom_jitter(aes(x = Habitat,
                  y = N,
                  col = Habitat),
              shape = 1,
              cex = 1.5,
              stroke = 0.5,
              width = 0.2,
              show.legend = F) +
  geom_label(data = H2O_lmer_cld,
             aes(x = Habitat, label = label),
             y = 0.095,
             label.size = 0, 
             size = 4.5,
             alpha = 0) +
  scale_x_discrete(labels = c("Bay","Frng. reef", "Mid lagoon", "Back reef", "Reef pass")) +
  scale_y_continuous(limits = c(0.075,NA)) +
  scale_color_manual(values = reefColors) +
  ylab("Water col. N (µM)") +
  xlab("") +
  CJsBasics::BasicTheme +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title.y = element_text())
Fig2b


## Figure 2ab ----
## Combine figure 2 a and b panels
Fig2ab = cowplot::plot_grid(Fig2a,
                            Fig2b,
                            ncol = 2,
                            align = "hv")


## Tests for interannual consistency ----
turb_diffs_lm = lm(Percent_N ~ Habitat*as.factor(Year),
                   data = compiled_N)
summary(turb_diffs_lm)
summary(aov(turb_diffs_lm))
## Note, interaction between habitat and year is not significant.
H2O_diffs_lm = lm(Nitrite_plus_Nitrate ~ Habitat*as.factor(Year),
                  data = compiled_N)
summary(H2O_diffs_lm)
summary(aov(H2O_diffs_lm))
## Note, interaction between habitat and year is not significant.

h2o_turb_relationships = lm(Percent_N ~ Nitrite_plus_Nitrate*as.factor(Year), 
                            data = compiled_N)
summary(h2o_turb_relationships)
summary(aov(h2o_turb_relationships))

## Compare Turbinaria %N and Water column [N] ----
annualTurbVsWater = compiled_N %>%
  ggplot() +
  geom_point(aes(y = Percent_N,
                 x = Nitrite_plus_Nitrate,
                 col = Habitat),
             pch = 1) +
  stat_smooth(aes(y = Percent_N,
                  x = Nitrite_plus_Nitrate),
              col = "black",
              method = "lm",
              # se = F,
              lwd = 1.25)  +
  facet_wrap(~Year) +
  ylab(expression(paste(italic("Turbinaria"), " %N"))) +
  xlab("Water col. Nitrite + Nitrate (µM)") +
  CJsBasics::BasicTheme +
  scale_color_manual(values = reefColors)
annualTurbVsWater

## Ask: Does Turb N predict water col N?
n_comparison1 = lm(Nitrite_plus_Nitrate ~ Percent_N + as.factor(Year), 
                   data = compiled_N)
summary(n_comparison1)

## Ask if Turbinaria N is important with partial pooling by year
n_comparison2 = lmer(Nitrite_plus_Nitrate ~ Percent_N + (1|Year), 
                     data = compiled_N %>% mutate(Year= as.factor(Year)))
summary(n_comparison2)
MuMIn::r.squaredGLMM(n_comparison2)

## Check for a significant interaction by year
n_comparison3 = lm(Nitrite_plus_Nitrate ~ Percent_N*Year, 
                   data = compiled_N %>% mutate(Year= as.factor(Year)))
summary(n_comparison3)
summary(aov(n_comparison3))
## Note, interaction is not significant

## Ask if water col N predicts Turb N
n_comparison4 = lmer(Percent_N ~ Nitrite_plus_Nitrate + (1|Year), 
                     data = compiled_N%>%mutate(Year= as.factor(Year)))
summary(n_comparison4)
MuMIn::r.squaredGLMM(n_comparison4)

# In general, Turbinaria %N was positively associated with water column N (Fig 2c; p < 0.001; R2 = 0.09)
longTermQuestion = compiled_N %>%
  ## Use only years in which both water column and Turbinaria data were collected
  filter(Year %in% c(2021:2023)) %>%
  ## Idenfify mean N by site
  group_by(Site, Habitat) %>%
  summarize(meanTurbN = mean(Percent_N, 
                             na.rm. = F),
            meanWatrN = mean(Nitrite_plus_Nitrate, 
                             na.rm. = F)) %>%
  ungroup() %>%
  filter(complete.cases(.))
longTermLM = lm(meanTurbN ~ meanWatrN, data = longTermQuestion)
summary(longTermLM)
summary(aov(longTermLM))

#### Figure 2c ----
Fig2c = longTermQuestion %>%
  ggplot() +
  geom_point(aes(x = meanWatrN,
                 y = meanTurbN,
                 col = Habitat),
             shape = 1,
             cex = 1.5,
             stroke = 0.75) +
  stat_smooth(aes(x = meanWatrN,
                  y = meanTurbN),
              col = "black",
              method = "lm",
              lwd = 1.25)  +
  xlab("Water col. N (µM)") +
  ylab(expression(paste(italic("Turbinaria"), " %N"))) +
  scale_color_manual(values = reefColors) +
  CJsBasics::BasicTheme +
  theme(legend.position = "none",
        legend.spacing.y = unit(0,"in"),
        legend.key.height = unit(0.05,"in"),
        panel.background = element_blank(),
        plot.background = element_blank())
Fig2c


#### Figure 2abc combined ----
turb_water_compare = plot_grid(Fig2ab,
                               Fig2c,
                               ncol = 1,
                               # align = "v",
                               # axis = "trbl",
                               rel_heights = c(0.8,1),
                               labels = NULL)
turb_water_compare
ggsave("plots/supplement/s3_figure2_noImputedVals/turb_water_compare.jpg",
       turb_water_compare,
       width = 3.55, height = 5.25, units = "in", dpi = 600)


