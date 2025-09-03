## This script generates panels d and e for Supplementary materials S3
## Repeating Fig 2 results with non-interpolated data

## Author: Christian John
## Date: 14 May 2025


#### Load packages ----
library(tidyverse) # Data wrangling
library(cowplot) # Plot wrangling
library(terra) # Raster data
library(sf) # Simple features
library(gstat) # Kriging operations
library(emmeans) # Post hoc comparisons
library(multcomp) # Plot labels for stats
source("scripts/000_reefCols.R")


set.seed(1959) # Homage to Santa Rosalia

#### Import data ----
## ... Nutrient data ----
## Turb N
sitesWithMissingN = read_csv("data/generated/turbSitesWithMissingN.csv")$Site
turb_N = read_csv("data/generated/turb_N_compiled.csv") %>%
  ## Note: remove sites that had interpolation ----
  filter(!(Site %in% sitesWithMissingN))
imputedSites = read_csv("data/generated/turb_N_imputed.csv") %>%
  pull(Site) %>%
  unique()

## Water N
water_N = read_csv("data/generated/water_N_compiled.csv")

## ... Raster data ----
## Hillshade
## r is a hillshade layer for plotting the final bivariate color map 
r = rast("data/spatialLayers/hillshade.tif") %>%
  project("epsg:32706")

## Land use
lulc = rast("~/Documents/Burkepile_Lab/Manuscripts_coauthorship/202x_Neumann_landuse/consensusLandUse.tif") %>%
  aggregate(4, fun="modal") %>%
  project(crs(r),method = "near")
lulc_attr = readxl::read_xlsx(path = "/Users/christianjohn/Documents/Burkepile_Lab/Manuscripts_coauthorship/202x_Neumann_landuse/EDI_metadata/metadata_landuse.xlsx",sheet = "Attribute Codes") %>%
  dplyr::select(Code, Definition)

## Lagoon bathymetry
lagoon = rast("~/Documents/Burkepile_Lab/Moorea_mapping_quicklayers/data/lagoon_3m.tif") %>%
  classify(rcl=matrix(data=c(-99,0,NA,
                             0,5,1),
                      ncol=3,byrow=T))
lagoon_sf = as.polygons(lagoon) %>%
  st_as_sf()

## ... Vector data ----
## Reef crest
reefcrest = read_sf("~/Documents/Burkepile_Lab/Moorea_mapping_quicklayers/data/crest.shp") %>%
  st_transform(crs(r))

## Moorea perimeter
Moorea_perimeter = st_read("data/spatialLayers/Moorea_LOC_ILE.shp") %>%
  st_transform(crs(r))

## Site metadata
## Sampling locations: This is used for identifying sampling location coordinates
## and distance from shore
samplingLocations = read_csv(
  "~/Documents/Burkepile_Lab/Data/ATI/EDI_site_locations/data/clean/sites_metadata_DRAFT.csv") %>%
  st_as_sf(coords = c("Longitude","Latitude")) %>%
  st_set_crs(4326) %>%
  st_transform(crs(r)) 

sitesLocs = unique(samplingLocations$Site)
sitesData = unique(turb_N$Site)
## Check for sites that have data but no metadata (would be bad)
sitesData[!(sitesData%in%sitesLocs)] ## There are no data sites for which we don't have metadata
sitesLocs[!(sitesLocs%in%sitesData)] ## There are some sites that we don't have data for (not an issue just means we have metadata for other spots in the lagoon)


#### Data wrangling ----
## Begin by prepping spatial data
## ... Raster prep ----
## Clip to region of interest
r = mask(r,vect(Moorea_perimeter))
## Aggregate to coarser resolution if desired
rAgg = aggregate(r,factor = 4)
## Convert to data.frame and scale pixel values as necessary
mapRaster = as.data.frame(r, xy=T)

## Prune land use layers to ROI and aggregate
lulc_mask = mask(lulc, vect(Moorea_perimeter))
## Convert rasters to data.frames that are ggplot-friendly
lulc_df = as.data.frame(lulc_mask, xy=T) %>%
  tibble() %>%
  left_join(lulc_attr %>%
              rename(classnm = Code)) %>%
  mutate(Definition = factor(Definition,
                             levels = 
                               sort(lulc_attr$Definition))) %>%
  filter(Definition != "Water")

## ... Prep Turbinaria data for analysis ----
## Create map point vector data
mapData = turb_N %>%
  group_by(Site) %>%
  summarize(Percent_N = mean(Percent_N)) %>%
  dplyr::select(Site, Percent_N) %>%
  left_join(samplingLocations %>%
              st_drop_geometry()) 

#### Calculate residual Turbinaria N ----
## ... Visualize dependence of N on distance from shore ----
water_lm_data = water_N %>%
  group_by(Site) %>%
  summarize(N = mean(Nitrite_plus_Nitrate)) %>%
  left_join(samplingLocations %>%
              st_drop_geometry())

water_detrend_lm = lm(water_lm_data$N ~ water_lm_data$Distance_from_shore)
plot(water_lm_data$N ~ water_lm_data$Distance_from_shore)
abline(water_detrend_lm,col="red2")
summary(water_detrend_lm)

# Mean Turbinaria %N was highest near shorelines and diminished linearly with increasing distance from shore (decreasing by 0.036% per km from shore; Fig 2d, R2 = 0.092, p < 0.001).
LDST_detrend_lm = lm(mapData$Percent_N ~ mapData$Distance_from_shore)
plot(mapData$Percent_N ~ mapData$Distance_from_shore)
abline(LDST_detrend_lm,col="red2")
summary(LDST_detrend_lm)
# Call:
#   lm(formula = mapData$Percent_N ~ mapData$Distance_from_shore)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.174545 -0.039311  0.001511  0.034227  0.229908 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                  6.023e-01  6.748e-03  89.263  < 2e-16 ***
#   mapData$Distance_from_shore -5.258e-05  1.281e-05  -4.105 5.87e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.05925 on 202 degrees of freedom
# Multiple R-squared:  0.077,	Adjusted R-squared:  0.07243 
# F-statistic: 16.85 on 1 and 202 DF,  p-value: 5.866e-05

## ... Calculate residual variability ----
mapData$detrended = residuals(LDST_detrend_lm)
## Save output (used in downstream analyses)
## This save should only be done for the main figure; don't want to overwrite.
# write_csv(x = mapData, file = "data/generated/detrendedNutrients.csv")

## Figure s3d ----
## ... Plot distance trend with residuals ----
Fig2d = mapData %>%
  ggplot() +
  geom_point(aes(x = Distance_from_shore,
                 y = Percent_N,
                 col = detrended)) +
  stat_smooth(aes(x = Distance_from_shore,
                  y = Percent_N),
              method = "lm",
              col = "black",
              lwd = 1) +
  scale_color_gradientn("Residual\nN (%)",
                        colors = rev(c(RColorBrewer::brewer.pal(11,"BrBG")[1],
                                       RColorBrewer::brewer.pal(11,"BrBG")[1],
                                       RColorBrewer::brewer.pal(11,"BrBG")[1],
                                       RColorBrewer::brewer.pal(11,"BrBG")[1],
                                       RColorBrewer::brewer.pal(11,"BrBG")))) +
  xlab("Distance from shore (m)") +
  ylab(expression(paste(italic("Turbinaria"), " %N"))) +
  scale_x_continuous(breaks = c(0,500,1000)) +
  CJsBasics::BasicTheme +
  theme(legend.position = c(0.85,0.8))
Fig2d
# ggsave("plots/supplement/s3_figure2_noImputedVals/distanceTrend_TurbN.jpg",
#        Fig2d,
#        width = 4, height = 4, units = "in", dpi = 300)


#### Kriging ----
## ... Prep kriging data ----
turb_forKrig = samplingLocations %>%
  left_join(mapData) %>%
  filter(!is.na(Percent_N))

## Start by creating a window over which to apply kriging predictions
xrange <- c(st_bbox(turb_forKrig)$xmin - 1000, st_bbox(turb_forKrig)$xmax + 1000)
yrange <- c(st_bbox(turb_forKrig)$ymin - 1000, st_bbox(turb_forKrig)$ymax + 1000)
(xrange["xmax"] - xrange["xmin"])/1000 ## 19.7m intervals on x axis
(yrange["ymax"] - yrange["ymin"])/1000 ## 16.7m intervals on y axis.
## Round up to 20 m
gridInterval = 20
# Create a grid for prediction that extends the range
grid <- expand.grid(x = seq(xrange[1], xrange[2], by = gridInterval), 
                    y = seq(yrange[1], yrange[2], by = gridInterval))
grid_sf <- st_as_sf(grid, coords = c("x", "y"), crs = st_crs(turb_forKrig))

## ... Fit variogram ----
# Fit variogram to the detrended N data
vgm_model <- variogram(detrended ~ 1, data = turb_forKrig)
vgm_model
str(vgm_model)
fit_vgm <- fit.variogram(vgm_model, model = vgm("Sph"))
fit_vgm

#### Run the kriging operation ----
## Apply the kriging procedure to the detrended N data
kriging_result <- krige(detrended ~ 1, turb_forKrig, grid_sf, model = fit_vgm)
## Convert the kriging result to a tibble keeping only values within our roi
krige_df <- kriging_result %>%
  st_filter(lagoon_sf %>% st_transform(st_crs(kriging_result))) %>%
  mutate(x = st_coordinates(.)[,1],
         y = st_coordinates(.)[,2],
         z = var1.pred) %>%
  st_drop_geometry() %>%
  dplyr::select(x,y,z)
## Preview
ggplot(krige_df) +
  geom_raster(aes(x = x, y = y, fill = z))


#### Figure S3e ----
## ... Nitrogen residual map ----
Fig2e = ggplot() +
  geom_raster(data = krige_df,
              aes(x = x,
                  y = y,
                  fill = z), 
              show.legend = FALSE) +
  scale_fill_gradientn("Residual\nN (%)",
                       colors = rev(c(RColorBrewer::brewer.pal(11,"BrBG")[1],
                                      RColorBrewer::brewer.pal(11,"BrBG")[1],
                                      RColorBrewer::brewer.pal(11,"BrBG")[1],
                                      RColorBrewer::brewer.pal(11,"BrBG")))) +
  ggnewscale::new_scale("fill") +
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
  ggnewscale::new_scale("fill") +
  geom_raster(data = mapRaster,
              aes(x = x,
                  y = y,
                  fill = hillshade),
              alpha = 0.35,
              show.legend = F) +
  scale_fill_gradient(low = "black", high = "white") +
  geom_sf(data = Moorea_perimeter,
          fill = "grey90",
          alpha = 0.5,
          lwd = 0.5) +
  geom_sf(data = lagoon_sf, 
          fill = "transparent",
          col = "grey60",
          lwd = 0.25) +
  geom_sf(data = samplingLocations %>% filter(Site %in% imputedSites),
          size = 3.9,
          col = "grey60") +
  geom_sf(data = turb_forKrig,
          col = "black",
          show.legend = F,
          size = 4) +
  scale_color_manual(values = reefColors) +
  ggnewscale::new_scale_color() +
  geom_sf(data = turb_forKrig,
          aes(col = detrended),
          size = 2) +
  geom_text(aes(x = 192002,y = 8065700),
            label = "5km",
            size = 5) +
  geom_segment(aes(x = 189502, xend = 194502,
                   y = 8066170, yend = 8066170),
               col = "black",
               lwd = 1.5) +
  scale_color_gradientn("Residual\nN (%)",
                        colors = rev(c(RColorBrewer::brewer.pal(11,"BrBG")[1],
                                       RColorBrewer::brewer.pal(11,"BrBG")[1],
                                       RColorBrewer::brewer.pal(11,"BrBG")[1],
                                       RColorBrewer::brewer.pal(11,"BrBG")[1],
                                       RColorBrewer::brewer.pal(11,"BrBG")))) +
  theme_void()
residLegend = cowplot::get_legend(Fig2e + 
                                    theme(legend.key.height = unit(0.15,"in"),
                                          legend.title = element_text(size = 14),
                                          legend.text = element_text(size = 14)))

## Compile Figures 2 d and e
Fig2de <- cowplot::ggdraw() +
  cowplot::draw_plot(Fig2d + 
                       theme(legend.position = "none",
                             axis.title = element_text(size = 16),
                             axis.text = element_text(size = 14)),
                     0.005, 0.005, 0.325, 0.4) +
  cowplot::draw_plot(residLegend,
                     0.75, 0.01, 0.25, 0.3) +
  cowplot::draw_plot(Fig2e + 
                       theme(legend.position = "none"),
                     0.05, 0, 0.95, 1)
## Save
ggsave(paste0("plots/supplement/s3_figure2_noImputedVals/map_residualVariation.jpg"),
       Fig2de,
       height = 6, width = 7, units = "in", dpi = 600)

