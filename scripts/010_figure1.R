## This script generates all panels for Figure 1.

## Author: Christian John
## Date: 13 May 2025

#### Load packages ----
library(tidyverse)
library(sf)
library(terra)
library(lme4)
library(lmerTest)
library(emmeans)
library(multcomp)
library(multcompView)
source("scripts/000_reefCols.R")

#### Set seed ----
set.seed(64) # Number of genera in Sciuridae according to ASM

#### Import data ----
## ... Bathymetry data (north shore) ----
bathy_df = read_table("/Users/christianjohn/Documents/Burkepile_Lab/Data/data_MCRLTER/knb-lter-mcr.1036.6/MCR_LTER_Bathymetry_Point_Data_North_Shore_version1.0_20111123.csv", col_names = F) %>%
  set_names("x","y","z")
bathy_rast = rast(bathy_df, type="xyz",
                  crs = "epsg:32706")
plot(bathy_rast)

## ... Bathymetry data (island-wide) ----
lagoon = rast("~/Documents/Burkepile_Lab/Moorea_mapping_quicklayers/data/lagoon_3m.tif") %>%
  classify(rcl=matrix(data=c(-99,0,NA,
                             0,5,1),
                      ncol=3,byrow=T))
lagoon_sf = as.polygons(lagoon) %>%
  st_as_sf()

## ... Moorea outline ----
moo = read_sf("~/Documents/Burkepile_Lab/Moorea_mapping_quicklayers/data/Moorea_LOC_ILE.shp") %>%
  st_transform(st_crs(bathy_rast))

## ... LTER sites ----
lterPoly = read_sf("/Users/christianjohn/Documents/Burkepile_Lab/Data/data_MCRLTER/MCR_sites/MCR_LTER_Sites_polygons.kml") %>%
  filter(Name != "Cook's Bay") %>%
  st_transform(st_crs(bathy_rast))

## ... Reef crest line ----
crest = read_sf("~/Documents/Burkepile_Lab/Moorea_mapping_quicklayers/data/crest.shp") %>%
  st_transform(st_crs(bathy_rast))

## ... Sampling site locations ----
ati = read_csv("~/Documents/Burkepile_Lab/Data/ATI/EDI_site_locations/data/clean/sites_metadata_DRAFT.csv") %>%
  st_as_sf(coords = c("Longitude", "Latitude"),
           crs = 4326) %>%
  st_transform(st_crs(bathy_rast))
lter_Sites = ati %>%
  filter(grepl("LTER_", x = Site))

## ... DEM ----
## DEM and hillshade
dem = rast("~/Documents/Burkepile_Lab/Data/Moorea_DSM/Moorea_DSM_full.tif") %>%
  aggregate(4) %>%
  project(crs(bathy_rast))
hill = rast("~/Documents/Burkepile_Lab/Moorea_mapping_quicklayers/data/hillshade.tif") %>%
  project(crs(bathy_rast))

## ... Land use ----
lulc = rast("~/Documents/Burkepile_Lab/Data/data_MCRLTER/knb-lter-mcr.6005.2/consensusLandUse.tif") %>%
  aggregate(4, fun="modal") %>%
  project(crs(bathy_rast),method = "near")
lulc_attr = readxl::read_xlsx(path = "~/Documents/Burkepile_Lab/Data/data_MCRLTER/knb-lter-mcr.6005.2/metadata_landuse.xlsx",sheet = "Attribute Codes") %>%
  dplyr::select(Code, Definition)
## Prune raster layers to ROI and aggregate
lulc_mask = mask(lulc, vect(moo))
## Convert rasters to data.frames that are ggplot-friendly
lulc_df = as.data.frame(lulc_mask, xy=T) %>%
  tibble() %>%
  left_join(lulc_attr %>%
              rename(classnm = Code)) %>%
  mutate(Definition = factor(Definition,
                             levels = 
                               sort(lulc_attr$Definition)))

## ... Turbinaria N data ----
## Import Turbinaria data and do a quick cleaning operation;
## rename habitats to match Lagoon Nutrients naming conventions
## and aggregate dataset to match Lagoon Nutrients island-wide sampling protocol
turb = read_csv("/Users/christianjohn/Documents/Burkepile_Lab/Data/data_MCRLTER/knb-lter-mcr.20.21/MCR_LTER_Macroalgal_CHN_2005_to_2022_20230713.csv") %>%
  filter(Genus == "Turbinaria") %>%
  mutate(Habitat = case_when(Habitat == "Back Reef" ~ "Back reef",
                             Habitat == "Reef Crest" ~ "Reef crest",
                             Habitat == "Fringe" ~ "Fringing reef")) %>%
  mutate(Habitat = factor(Habitat, levels = c("Fringing reef",
                                              "Back reef",
                                              "Reef crest"))) %>%
  filter(is.na(Comment)) %>% ## Get rid of 4 insane observations that are clearly flukes
  filter(!is.na(N)) %>%
  group_by(Year, Site, Habitat, Genus) %>%
  ## This summary is to match the sampling resolution of the ATI data; 
  ## i.e. one estimate per site per year that comes from a composite of all 
  ## Turb samples at that site in that year.
  summarize(N = mean(N)) %>%
  mutate(Year = as.factor(Year))
turb


#### Data wrangling ----
## ... Create a lagoon depth gradient ----
myline = st_as_sf(data.frame(x = c(-149.837863,
                                   -149.840104),
                             y = c(-17.502665,
                                   -17.472586)),
                  coords = c("x","y"),
                  crs = 4326) %>%
  st_transform(st_crs(dem)) %>%
  mutate(line = "line") %>%
  group_by(line) %>%
  summarise() %>% # union points into lines using our created lineid
  st_cast("LINESTRING")
myline_vect = vect(myline)

## ... Extract bathymetry and topography from transect ----
mybathy = extract(bathy_rast, myline_vect, xy=T) %>%
  tibble() %>%
  filter(!is.na(z)) %>%
  mutate(product = "Bathymetry",
         z = z*-1)
myelev = extract(dem, myline_vect, xy=T) %>%
  tibble() %>%
  filter(y <= (min(mybathy$y))) %>%
  rename(z = Moorea_DSM_full) %>%
  mutate(product = "Topography")
mybathy
myelev
## Join the bathymetry and topography data
surface_transect = rbind(mybathy, myelev)

## Figure 1a ----
## Create a data.frame to add custom tick marks to x axis
shoreDistDF = data.frame(x = c(-250,0,250,500,750,1000),
                         y = rep(0, 6),
                         lab = c("-250","","250","500","750","1000"),
                         laby = rep(-1,6))
## Subset the surface transect to relevant region
surface_transect_selected = surface_transect %>%
  mutate(distFromShore = y - min(mybathy$y)) %>%
  filter(distFromShore > -300,
         distFromShore < 1150) 
## Create a set of colors to pair reef types with position on surface_transect
## This is a little hacky but it works
depthGradCols = c(reefColors[rep("Fringing reef",4)],
                  "grey80","grey80", "grey80","grey80","grey80",
                  "grey80","grey80","grey80","grey80","grey80",
                  rep(reefColors["Mid lagoon"], 7),
                  colorRampPalette(c(reefColors["Mid lagoon"], reefColors["Back reef"]))(6),
                  rep(reefColors["Back reef"], 6),
                  colorRampPalette(c(reefColors["Back reef"], reefColors["Reef crest"]))(4),
                  rep(reefColors["Reef crest"],1),
                  "grey80","grey80","grey80","grey80","grey80",
                  "grey80","grey80","grey80","grey80")

depthGrad = surface_transect_selected %>%
  ggplot() +
  geom_hline(yintercept = 0,lty=2,linewidth=0.5,
             col = "grey40") +
  geom_ribbon(aes(x = distFromShore, ymin = z, ymax = 0,
                  fill = product),
              alpha = 0.25,
              linewidth = 0) +
  geom_line(aes(x = distFromShore, y = z),
            data = surface_transect_selected %>% filter(distFromShore < 0),
            col = "tan4",
            linewidth = 1) +
  ggnewscale::new_scale_color() +
  geom_line(aes(x = distFromShore, y = z,
                col = distFromShore),
            data = surface_transect_selected %>% filter(distFromShore >= 0),
            linewidth = 1) +
  geom_point(data = shoreDistDF,
             aes(x = x, y = y),
             pch = "|",
             cex = 2,
             col = "grey40") +
  geom_label(data = shoreDistDF,
             aes(x = x, y = laby, label = lab),
             cex = 3.5,
             label.size = NA,
             alpha = 0,
             col = "grey40") +
  xlab("") +
  ylab("Surface height (m)") +
  scale_y_continuous(breaks = seq(from = -10, to = 30, by = 5)) +
  scale_fill_manual(values = c("skyblue3","tan4")) +
  scale_color_gradientn(colors = depthGradCols) +
  CJsBasics::BasicTheme +
  theme(legend.position = "none",
        plot.margin = margin(0.05,0.25,0.05,0.25,"in"),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(color = "grey40",
                                    size = 12,
                                    vjust = -1),
        axis.text.y = element_text(color = "grey40",
                                   size = 10),
        axis.ticks.y = element_line(color = "grey40"))
depthGrad
ggsave("plots/depthGradient.png",
       depthGrad,
       width = 6, height = 6, units = "in")


#### Figure 1b ----
lterMap = ggplot() +
  geom_sf(data = moo) +
  geom_raster(data = lulc_df,
              aes(x = x, y = y, fill = Definition),
              show.legend = FALSE) +
  scale_fill_manual(values = c("tan",
                               "tan",
                               "tan",
                               "tan",
                               "darkgreen",
                               "tan",
                               "tan",
                               "tan",
                               "tan",
                               "lightblue")) +
  ggnewscale::new_scale_fill() +
  geom_raster(data = as.data.frame(hill, xy = T),
              aes(x = x, y = y, fill = hillshade),
              alpha = 0.35,
              show.legend = F) +
  scale_fill_gradient(low = "black", high = "white") +
  geom_sf(data = lagoon_sf %>% st_transform(st_crs(moo)),
          fill = "grey90",
          col = "grey90",
          lwd = 0) +
  geom_sf(data = crest,
          fill = "transparent",
          col = "grey40",
          lwd = 0.25) +
  geom_sf(data = moo,
          fill = "transparent",
          col = "grey40",
          lwd = 0.25) +
  geom_sf(data = lterPoly,
          fill = "transparent",
          col = "black",
          lwd = 0.5) +
  geom_sf(data = lter_Sites,
          col = "black",
          size = 5,
          show.legend = F) +
  geom_sf(data = lter_Sites,
          col = "white",
          size = 4,
          show.legend = F) +
  geom_sf(data = lter_Sites,
          aes(col = Habitat),
          size = 4,
          alpha = 0.75,
          show.legend = F) +
  scale_color_manual(values = reefColors) +
  ggnewscale::new_scale_color() +
  geom_line(data = surface_transect_selected %>%
              mutate(coordID = 1:nrow(.)) %>%
              filter(coordID %in% c(1,max(coordID))),
            aes(x = x, y = y),
            col = "black",
            linewidth = 3,
            linetype = 1) +
  geom_line(data = surface_transect_selected %>%
              filter(distFromShore < 10) %>%
              filter(distFromShore %in% c(min(distFromShore),
                                          max(distFromShore))),
            aes(x = x, y = y),
            col = "tan4",
            linewidth = 2,
            linetype = 1) +
  geom_line(data = surface_transect_selected %>%
              filter(distFromShore >= 0) %>%
              filter(distFromShore %in% c(min(distFromShore),
                                          max(distFromShore))) %>%
              right_join(tibble(x = seq(min((.)$x), max((.)$x)))) %>%
              mutate(y = if_else(is.na(y),
                                 approx(x, y, x)$y,
                                 y)) %>%
              mutate(distFromShore = if_else(is.na(distFromShore),
                                             approx(x, distFromShore, x)$y,
                                             distFromShore)) %>%
              arrange(distFromShore),
            aes(x = x, y = y,
                col = distFromShore),
            show.legend = F,
            linewidth = 2,
            linetype = 1) +
  geom_text(aes(x = 192002,y = 8065700),
            label = "5km",
            size = 7) +
  geom_segment(aes(x = 189502, xend = 194502,
                   y = 8066170, yend = 8066170),
               col = "black",
               lwd = 1.5) +
  scale_color_gradientn(colors = c(reefColors[rep("Fringing reef",2)],
                                   "grey80","grey80", "grey80",
                                   "grey80","grey80","grey80",
                                   reefColors[c(rep("Mid lagoon",5),
                                                rep("Back reef",4),
                                                rep("Reef crest",1))],
                                   "grey80","grey80","grey80","grey80")) +
  theme_void()
lterMap
ggsave("plots/lter_MooreaMap.png",
       lterMap,
       width = 8, height = 6, units = "in", dpi = 600)


#### Turbinaria N by habitat ----
turb %>%
  group_by(Year, Site, Habitat) %>%
  summarize(N = mean(N, na.rm = T)) %>%
  ggplot() +
  geom_boxplot(aes(x = Habitat, y = N)) +
  geom_point(aes(x = Habitat, y = N)) +
  facet_wrap(~Site) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 7)) +
  ylab("Turbinaria %N") +
  CJsBasics::BasicTheme

## ... Check differences among years ----
## Check two ways: 
## ...... First, control for year using a mixed effects model ----
turb_lmer_data = turb %>%
  dplyr::select(Year, Site, Habitat, Genus, N) %>%
  filter(!is.na(N))
turb_diffs_lmer = lmer(N ~ Habitat + (1|Year) + (1|Site), data = turb_lmer_data)
summary(turb_diffs_lmer)
em_lmer = emmeans(turb_diffs_lmer, "Habitat")
contr_lmer = contrast(em_lmer, "pairwise")
## Across the long-term MCR LTER dataset, Turbinaria %N was highest at fringing sites (p < 0.001) but did not significantly differ among back reef and reef crest sites (p = 0.3).
contr_lmer
plot(turb_diffs_lmer)
## Create plottable comparisons using the emmeans results
turb_lmer_cld = cld(em_lmer, Letters = letters, reversed=T) %>%
  dplyr::select(Habitat, emmean, .group) %>%
  mutate(label = gsub(" ","",.group))

## ...... Then, test for consistency among years using a fixed interaction ----
turb_diffs_lm = lm(N ~ Habitat*Year, data = turb)
summary(turb_diffs_lm)
summary(aov(turb_diffs_lm))
## Interaction not significant, indicating that pattern among habitats
## is consistent across years
em_lm = emmeans(turb_diffs_lm, "Habitat")
contrast(em_lm, method = "pairwise")
turb_lm_cld = cld(em_lm, Letters = letters, reversed=T) %>%
  dplyr::select(Habitat, emmean, .group) %>%
  mutate(label = gsub(" ","",.group))
turb_lm_complabels = turb_lm_cld %>%
  left_join(turb %>%
              group_by(Habitat) %>%
              summarize(Nse = plotrix::std.error(N, na.rm = T),
                        N75 = quantile(N, 0.75, na.rm = T),
                        Nmean = mean(N, na.rm = T))) %>% 
  mutate(labx = rep(2013,3), laby = seq(from = 1.4,1.2,length.out = 3),
         text = paste0(Habitat," ",round(Nmean,2),"±",round(Nse,2)," (",label,")"))
turb_lm_complabels

#### Figure 1c ----
## Plot overall mean Turb N for each site, boxes for each habitat
## Include stat labels for lmer comparing habitats with year as random intercept
turb_lmer_boxplots = turb_lmer_data %>%
  group_by(Site, Habitat) %>%
  summarize(N = mean(N, na.rm = T)) %>%
  ggplot() +
  geom_boxplot(aes(x = Habitat, 
                   y = N, 
                   col = Habitat),
               outliers = F,
               show.legend = F) +
  geom_jitter(aes(x = Habitat,
                  y = N,
                  col = Habitat),
              size = 4,
              show.legend = F,
              width = 0.1,
              alpha = 0.75) +
  geom_label(data = turb_lmer_cld,
             aes(x = Habitat, label = label),
             y = 0.56,
             label.size = 0, 
             size = 6,
             alpha = 0) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 7)) +
  scale_y_continuous(limits = c(0.55,NA)) +
  scale_color_manual(values = reefColors) +
  ylab(expression(paste(italic("Turbinaria"), " %N"))) +
  xlab("") +
  CJsBasics::BasicTheme +
  theme(panel.background = element_blank(),
        plot.background = element_blank(),
        strip.text = element_text(size = 18),
        axis.title = element_text(size = 18),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 14))
turb_lmer_boxplots
ggsave("plots/lter_turbN_boxplots_lmer.png",
       turb_lmer_boxplots,
       width = 3, height = 3.2, units = "in")

## Plot time series of Turbinaria %N
## Comparison labels here come from emmeans of lmer 
# turb_timeseries = turb %>%
#   filter(!is.na(N)) %>%
#   group_by(Site, Habitat,Year) %>%
#   summarize(N = mean(N)) %>%
#   ggplot() +
#   geom_line(aes(x = as.numeric(as.character(Year)), y = N,
#                 group = paste0(Site,Habitat),
#                 col = Habitat),
#             lwd = 1) +
#   geom_text(data = turb_diffs_complabels,
#             mapping = aes(x = labx-6,
#                           y = laby,
#                           label = as.character(text),
#                           col = Habitat),
#             fontface = "bold",
#             hjust = 0,
#             size = 4,
#             show.legend = F) +
#   xlab("Year") +
#   ylab(expression(paste(italic("Turbinaria"), " %N"))) +
#   CJsBasics::BasicTheme +
#   scale_color_manual(values = reefColors) +
#   theme(legend.position = "none",
#         legend.background = element_blank(),
#         legend.key.height = unit(0.15,"in"),
#         panel.background = element_blank(),
#         plot.background = element_blank())
# turb_timeseries

#### Figure 1d ----
## Plot impacts of land clearing on Turbinaria N
landuse_stats_df = turb %>% 
  group_by(Site, Habitat) %>% 
  summarize(N = mean(N, na.rm = T)) %>%
  left_join(lter_Sites %>%
              st_drop_geometry() %>%
              separate(Site, into = c("LTER","Site","Hab")) %>%
              mutate(Site = paste0("LTER ", Site)) %>%
              dplyr::select(Site, Habitat, Distance_from_shore,
                            Cleared_area, Watershed_area))

## Generate summary of association between Turb N and land clearing
summary(lm(N ~ Watershed_area,
           data = landuse_stats_df))
summary(lm(N ~ Cleared_area,
           data = landuse_stats_df))
clearedP = summary(lm(N ~ Cleared_area,
                      data = landuse_stats_df))$coefficients[2,4]
clearedR2 = summary(lm(N ~ Cleared_area,
                       data = landuse_stats_df))$r.squared
summary(lm(N ~ Distance_from_shore,
           data = landuse_stats_df))
mymodel = lm(N ~ Cleared_area + Distance_from_shore,
             data = landuse_stats_df)
## Nutrient enrichment was positively associated with land clearing, and 63% of the variability of overall mean Turbinaria %N was explained by distance from shore and cleared area in a simple linear model (p < 0.001). 
summary(mymodel)
summary(aov(mymodel))
modelR2 = summary(mymodel)$r.squared %>%
  round(2)

## Within the level of individual habitat types, land clearing and distance from shore explained over 90% of the variability of Turbinaria %N (p = 0.01) for fringing reef sites, but these factors were not significantly associated with Turbinaria %N for back reef or reef crest sites (p = 0.7 and 0.8, respectively).
lm(N ~ Cleared_area + Distance_from_shore,
   data = landuse_stats_df %>% filter(Habitat == "Fringing reef")) %>% summary()
lm(N ~ Cleared_area + Distance_from_shore,
   data = landuse_stats_df %>% filter(Habitat == "Back reef")) %>% summary()
lm(N ~ Cleared_area + Distance_from_shore,
   data = landuse_stats_df %>% filter(Habitat == "Reef crest")) %>% summary()

## Plot association
lter_clearedArea = landuse_stats_df %>%
  ggplot(aes(x = Cleared_area, 
             y = N)) +
  stat_smooth(method = "lm",
              col = "black") +
  geom_point(aes(col = Habitat),
             size = 4,
             alpha = 0.75) +
  scale_color_manual("",
                     values = reefColors) +
  xlab(expr(Cleared~area~(km^2))) +
  ylab(expr(paste(italic("Turbinaria"), " %N"))) +
  CJsBasics::BasicTheme +
  theme(legend.position = c(0.7,0.175),
        legend.text = element_text(size = 14),
        legend.key.height = unit(0.1,
                                 "in"),
        legend.key.width = unit(0.05,
                                "in"),
        legend.background = element_rect(fill = "transparent",
                                         color = "transparent"),
        panel.background = element_rect(fill = "transparent",
                                        color = "transparent"),
        panel.border = element_rect(fill = "transparent",
                                    color = "transparent"),
        plot.background = element_rect(fill = "transparent",
                                       color = "transparent"),
        strip.text = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.title.x = element_text(size = 18,
                                    vjust = 1.5),
        axis.text.x = element_text(size = 16,
                                   vjust = 0.5),
        axis.text.y = element_text(size = 14)) +
  guides(col = guide_legend(label.position = "left",
                            label.hjust = 1))
lter_clearedArea
ggsave("plots/lter_turbN_clearedArea.png",
       lter_clearedArea,
       width = 3, height = 3.2, units = "in", dpi = 300)


## Supplementary materials S7 ----
## A function to compare habitat types within the site level
compfx = function(mySite){
  sitediffs = turb %>%
    filter(Site == mySite) %>%
    group_by(Year, Habitat) %>%
    summarize(N = mean(N, na.rm = T))  %>%
    filter(!is.nan(N))
  site_diffs_lm = lmer(N ~ Habitat + (1|Year), data = sitediffs)

  em_site_lm = emmeans(site_diffs_lm, "Habitat")
  em_site_lm_cld = cld(em_site_lm, Letters = letters, reversed=T) %>%
    dplyr::select(Habitat, emmean, .group) %>%
    mutate(label = gsub(" ","",.group))
  complabels = em_site_lm_cld %>%
    left_join(turb %>%
                group_by(Habitat) %>%
                summarize(Nse = plotrix::std.error(N, na.rm = T),
                          N75 = quantile(N, 0.75, na.rm = T),
                          Nmean = mean(N, na.rm = T))) %>% 
    mutate(labx = rep(2013,3), laby = seq(from = 1.4,1.2,length.out = 3),
           Site = mySite,
           text = paste0(Habitat," ",round(Nmean,2),"±",round(Nse,2)," (",label,")"))
  rownames(complabels) <- NULL
  return(complabels)
}

turb_sitelabels = do.call("rbind",
                     lapply(unique(turb$Site),
                            compfx))
turb_sitelabels

turb_sitewise_boxplots = turb %>%
  group_by(Year, Site, Habitat) %>%
  summarize(N = mean(N, na.rm = T)) %>%
  ggplot() +
  geom_boxplot(aes(x = Habitat, 
                   y = N, 
                   col = Habitat),
               show.legend = F) +
  geom_jitter(aes(x = Habitat,
                  y = N,
                  col = Habitat),
              show.legend = F,
              width = 0.1) +
  geom_label(data = turb_sitelabels,
             aes(x = Habitat, y = N75, label = label),
             label.size = 0,
             size = 5,
             nudge_x = 0.3,
             nudge_y = 0.2,
             alpha = 0) +
  facet_wrap(~Site,strip.position = "right") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 7)) +
  scale_color_manual(values = reefColors) +
  ylab(expression(paste(italic("Turbinaria"), " %N"))) +
  xlab("") +
  CJsBasics::BasicTheme +
  theme(panel.background = element_blank(),
        plot.background = element_blank(),
        strip.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))
turb_sitewise_boxplots
ggsave("plots/supplement/s7_lter_turbN_sitewise_boxplots.png",
       turb_sitewise_boxplots,
       width = 8, height = 3, units = "in")

