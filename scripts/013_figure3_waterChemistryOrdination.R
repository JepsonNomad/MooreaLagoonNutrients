## This script generates Figure 3 and relevant supplementary materials
## Variability in water chemistry

## Author: Christian John
## Date: 13 May 2025


#### Load Libraries ----
library(tidyverse)
library(cowplot)
library(sf)
library(RcppRoll)
library(vegan)
library(ggordiplots)
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

## Observed water column nutrients
dataWater<-read_csv("data/generated/water_N_compiled.csv")

## Site metadata 
dataSites <- read_csv("~/Documents/Burkepile_Lab/Data/ATI/EDI_site_locations/data/clean/sites_metadata_DRAFT.csv")

#### Data wrangling ----
## > Combine datasets ----
## Fuse water and turbinaria data
data_raw <- dataWater %>%
  left_join(dataTurbInSitu, by = c("Site", "FieldSeason")) %>%
  left_join(dataTurbDetrended, by = "Site") %>%
  left_join(dataSites, by = "Site")
data_raw 

## Make sure there aren't any inadvertent duplications after all this joining
data_raw %>% 
  mutate(dupes = duplicated(paste0(Site, FieldSeason))) %>%
  pull(dupes) %>%
  sum()

## Filter to just the sites and times with water column data
data_filtered <-data_raw %>%
  filter(complete.cases(Phosphate, Silicate, Nitrite_plus_Nitrate)) %>%
  filter(!is.na(Percent_N)) %>% # remove the data that have no turbinaria data
  select(-Ammonia) ## Remove ammonia which had QAQC issues

## > Scale data ----
## Scale the data by year 
## Using a log(x+1) followed by z score transformation
data_logged <- data_filtered %>%
  group_by(Year) %>%
  mutate_at(vars(c("Phosphate","Silicate","Nitrite_plus_Nitrate",
                   "Ultra.Violet.Humic.like","Tyrosine.like","Visible.Humic.like",
                   "Marine.Humic.like","Tryptophan.like","Lignin.like","M.C",
                   "HIX","BIX","FI")), 
            function(x){log(x+1)}) %>%
  ungroup()

data_scaled = data_logged %>%
  group_by(Year) %>%
  mutate_at(vars(c("Phosphate","Silicate","Nitrite_plus_Nitrate",
                   "Ultra.Violet.Humic.like","Tyrosine.like","Visible.Humic.like",
                   "Marine.Humic.like","Tryptophan.like","Lignin.like","M.C",
                   "HIX","BIX","FI")), 
            list(scale = scale)) %>%
  ungroup()

## Nsite, i.e. overall mean Turbinaria N, 
## is strongly predictive of instantaneous N, 
## i.e. 1-time measure of Turb N at a site
data_scaled %>% 
  ggplot(aes(x = Nsite, 
             y = Percent_N)) + 
  geom_point() + 
  stat_smooth(method = "lm")


#### CJ Main results ----
## two ways to look at this: permanova and envfit. Permanova looks at groups; envfit looks at continuous environmental gradients.

## > Envfit ----
## Explanatory variables

## For this we need two splits of one dataset:
pc_data <- data_scaled %>%
  select(Site,
         Date,
         Phosphate_scale:FI_scale,
         Percent_N,
         Nsite,
         detrended,
         Habitat,
         Watershed_area,
         Cleared_area, 
         Distance_from_shore) %>%
  filter(complete.cases(.))

## First split, the nutrient+fdom community data
nutrientfdom_data_community <- pc_data %>%
  select(Phosphate_scale:FI_scale)

## Second split, the environmental data.
data_envfit = pc_data %>%
  mutate(logTurbinariaN = log(Nsite+1),
         logDetrendN = log(detrended+1),
         logAreakm2 = log(Watershed_area+1),
         logDistance = log(Distance_from_shore+1),
         logClearedArea = log(Cleared_area+1))

## Perform principal components analysis
pc <- prcomp(nutrientfdom_data_community)
print(pc) #axis loadings
biplot(pc) #Note, these loadings are essentially the same as the biplot2 findings from princomp() above, but now calculated using the prcomp() function instead
summary(pc) #sd, variance explained of axes

## Detrended Turbinaria N is not as good a predictor of water chemistry. 
m3 <- envfit(pc ~ 
               logDistance + logClearedArea + logTurbinariaN + detrended,
             data = data_envfit, 
             perm = 999)
m3
m3b <- bioenv(nutrientfdom_data_community ~ 
                logDistance + logClearedArea + logTurbinariaN + detrended,
              data = data_envfit)
m3b
summary(m3b)

## Generalized additive models identify significant associations between continuous environmental parameters and the ordination space (b; bold arrows; r2=0.15, 0.03, 0.15, and 0.09 for Distance from shore, Cleared Area, Turbinaria %N, and detrended Turbinaria %N after accounting for distance from shore respectively, and p<0.01 in all cases).
## Compile loadings from PCA and from explanatory analysis
myspecies_expl = names(m3$vectors$r[m3$vectors$pvals<0.05])
expl.m3 <- as.data.frame(scores(m3, c("vectors"))) %>%
  mutate(species = rownames(.)) %>%
  filter(species %in% myspecies_expl) %>%
  mutate(species = case_when(species == "areakm2" ~ "Watershed area",
                             species == "log(areakm2 + 1)" ~ "Watershed area",
                             species == "Distance" ~ "Dist. from shore",
                             species == "logDistance" ~ "Dist. from shore",
                             species == "log(Distance + 1)" ~ "Distance from shore",
                             species == "Cleared" ~ "Cleared (%)",
                             species == "log(Cleared + 1)" ~ "Cleared (%)",
                             species == "clearedArea" ~ "Cleared area",
                             species == "logClearedArea" ~ "Cleared area",
                             species == "log(clearedArea + 1)" ~ "Cleared area",
                             species == "Turbinaria_N" ~ "Turbinaria N mean",
                             species == "logTurbinariaN" ~ "Turbinaria N mean",
                             species == "log(Turbinaria_N + 1)" ~ "Turbinaria N mean",
                             species == "detrended" ~ "Residual Turbinaria N",
                             species == "Turbinaria_V" ~ "Turbinaria N variance",
                             TRUE ~ species),
         speciesP = m3$vectors$pvals[m3$vectors$pvals<0.05],
         printSpecies = paste0(species,ifelse(speciesP<=0.001,
                                              " ***",
                                              ifelse(speciesP <= 0.01,
                                                     " **",
                                                     ifelse(speciesP <= 0.05,
                                                            " *",
                                                            "")))),
         PClab1 = 15*PC1 + 0.4*sign(PC1),
         PClab2 = 15*PC2 + 0.4*sign(PC2),
         PCarrow1 = 15*PC1,
         PCarrow2 = 15*PC2)


## > Permanova ----
data_envfit
permanova1 = vegan::adonis2(formula = nutrientfdom_data_community ~ Habitat,
                            data = data_envfit,
                            method = "euclidean")
permanova1

h20dist = vegan::vegdist(data_scaled %>%
                           select(c("Phosphate","Silicate","Nitrite_plus_Nitrate",
                                    "Ultra.Violet.Humic.like","Tyrosine.like",
                                    "Visible.Humic.like","Marine.Humic.like",
                                    "Tryptophan.like","Lignin.like","M.C",
                                    "HIX","BIX","FI")), na.rm = T, 
                         method = "euclidean")
str(h20dist)

## PERMANOVA revealed significant differences in the water chemistry community by discrete habitat type (Fig. 3a, p = 0.001).
## Also Figure 3 caption:
## PERMANOVA reveals significant differences in water column chemistry among habitat types (a; R2=0.08, F=12.45, p=0.001; ellipses based on 1 standard deviation, and centroids shown as bold dots)
permanova_distbased = vegan::adonis2(formula = h20dist ~ Habitat,
                                     data = data_scaled, method = "euclidean")
permanova_distbased
pred_distbased = data_scaled %>%
  mutate(logTurbinariaN = log(Nsite+1),
         logAreakm2 = log(Watershed_area+1),
         logDistance = log(Distance_from_shore+1),
         logClearedArea = log(Cleared_area+1)) %>%
  select(logTurbinariaN,
         detrended,
         logDistance,
         logClearedArea)
## Further, overall mean Turbinaria %N (p = 0.001), distance from shore (p = 0.001), residual Turbinaria %N after accounting for distance from shore (p = 0.001), and the area of cleared land in the upstream watershed (p = 0.001) were all significantly associated with the first two water chemistry principal components (Fig. 3b).
m3distbased <- envfit(h20dist ~ 
                        logDistance + logClearedArea + logTurbinariaN + detrended,
                      data = pred_distbased, 
                      perm = 999)
m3distbased

m3bdistbased <- bioenv(h20dist ~ 
                         logDistance + logClearedArea + logTurbinariaN + detrended,
                       data = pred_distbased)
m3bdistbased
summary(m3bdistbased)
library(AICcPermanova)
AllModelsdist <- make_models(vars = names(pred_distbased))
AllModelsdist
NonColineardist <- filter_vif(all_forms = AllModelsdist, 
                              env_data = pred_distbased,
                              ncores = 4,
                              verbose = T)
NonColineardist
Fitteddist <- fit_models(
  all_forms = NonColineardist,
  com_data = h20dist,
  env_data = pred_distbased,
  ncores = 4,
  method = "euclidean",
  verbose = T
)
Fitteddist
Selecteddist <- select_models(Fitteddist)
topModVars = Selecteddist[1,1] %>%
  gsub(pattern = "Distance ~ ", replacement = "") %>%
  gsub(pattern = " ", replacement = "") %>%
  strsplit(split = "\\+") %>%
  unlist()
topModVars
Summary <- akaike_adjusted_rsq(Selecteddist)
Summary
topModel = adonis2(
  h20dist ~ 
    detrended +
    logDistance,
  data = pred_distbased)
topModel

# Residual Turbinaria %N after accounting for distance from shore was also associated with variation along PC2 (p = 0.001)
resid_pc_assocdat = data_envfit %>%
  cbind(pc$x)
lm(PC2 ~ detrended, data = resid_pc_assocdat) %>% summary()

#### Figure 3a ----
pc_centroids = data_envfit %>%
  cbind(pc$x) %>%
  group_by(Habitat) %>%
  summarize(PC1 = mean(PC1),
            PC2 = mean(PC2))
myordplt = gg_ordiplot(
  ord = pc,
  groups = factor(
    data_envfit$Habitat,
    levels = reefLevels[reefLevels%in%(unique(data_envfit$Habitat))]),
  ellipse = T,
  kind = "sd",
  pt.size = 2,
  # show.groups = F,
  hull = F,
  spiders = F)
pcInfo = as.data.frame(pc$rotation) %>%
  mutate("Group" = rownames(.)) %>%
  mutate(Group = gsub(".like_scale","-like",Group),
         Group = gsub(".H"," H",Group),
         Group = gsub(".V"," V",Group),
         Group = gsub("_scale","",Group),
         Group = gsub("M.C","M:C",Group),
         Group = gsub("[_]"," ",Group))
Figure3a = myordplt$plot +
  scale_color_manual(values = setNames(
    paste0(reefColors,"70"),
    names(reefColors))) +
  ggnewscale::new_scale_color() +
  geom_path(data = myordplt$df_ellipse,
            aes(x=x, y=y, group = Group),
            col = "black",
            lwd = 1.95) +
  geom_path(data = myordplt$df_ellipse,
            aes(x=x, y=y, color=Group),
            lwd = 1.5) +
  geom_point(data = pc_centroids,
             aes(x = PC1, y = PC2, col = Habitat),
             col = "black",
             size = 6) +
  geom_point(data = pc_centroids,
             aes(x = PC1, y = PC2, col = Habitat),
             size = 5) +
  scale_color_manual("Habitat", values = reefColors) +
  coord_cartesian(ylim=c(-6,8),
                  xlim=c(-10,4))  +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))
Figure3a

envfitplot = gg_ordiplot(
  ord = pc,
  groups = factor(
    data_envfit$Habitat,
    levels = reefLevels[reefLevels%in%(unique(data_envfit$Habitat))]),
  ellipse = F,
  pt.size = 2,
  hull = F,
  spiders = F,
  plot = F)

#### Figure 3b ----
Figure3b = pcInfo %>% 
  ggplot() +
  geom_segment(data = pcInfo,
               aes(x = 0, xend = PC1*10, y = 0, yend = PC2*10),
               arrow=arrow(length=unit(0.2,"cm")), 
               alpha=0.75, color="black") +
  geom_segment(data=expl.m3,
               aes(x=0,xend=PCarrow1,y=0,yend=PCarrow2),
               arrow = arrow(length = unit(0.3, "cm")),
               lwd = 1.5,
               col = "black") +
  geom_segment(data=expl.m3,
               aes(x=0,xend=PCarrow1,y=0,yend=PCarrow2),
               arrow = arrow(length = unit(0.3, "cm")),
               lwd = 1,
               col = "black") +
  geom_label(data = pcInfo,
             aes(x = (PC1*10 + sign(PC1)*0.5), y = (PC2*10 + sign(PC2)*0.5),
                 label = Group),
             alpha = 0.15,
             col = "transparent") +
  geom_label(data = pcInfo,
             aes(x = (PC1*10 + sign(PC1)*0.5), y = (PC2*10 + sign(PC2)*0.5),
                 label = Group),
             alpha = 0.5,
             col = "transparent") +
  geom_label(data = pcInfo,
             aes(x = (PC1*10 + sign(PC1)*0.5), y = (PC2*10 + sign(PC2)*0.5),
                 label = Group),
             alpha = 0.25,
             col = "black") +
  geom_label(data=expl.m3,
             aes(x=PClab1,y=PClab2,
                 label=printSpecies),
             fontface = "bold",
             size=5, 
             col = "black",
             fill = alpha(c("white"),0.75)) +
  coord_cartesian(ylim=c(-6,8),
                  xlim=c(-10,4))  +
  scale_x_continuous(breaks = c(0),labels="0") +
  
  scale_y_continuous(breaks = c(0),labels="0") +
  xlab("PC1") +
  ylab("PC2") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16))
Figure3b

fig3brepeldata = pcInfo %>%
  select(Group, PC1, PC2)  %>% 
  mutate(type = "chem",
         PC1 = PC1*10,
         PC2 = PC2*10) %>%
  rbind(expl.m3 %>% mutate(Group = printSpecies) %>% 
          select(Group, PCarrow1, PCarrow2) %>% 
          rename("PC1" = PCarrow1,
                 "PC2" = PCarrow2) %>%
          mutate(type = "env")) %>%
  mutate(typealpha = ifelse(type == "chem", 0.5, 1),
         typesize = ifelse(type == "chem", 4, 5),
         typeface = ifelse(type == "chem", 1, 2),
         typewidth = ifelse(type == "chem", 0.25, 1.25))
Figure3brepel = pcInfo %>%
  ggplot() +
  geom_segment(data = fig3brepeldata,
               aes(x = 0, xend = PC1, 
                   y = 0, yend = PC2,
                   col = Group,
                   linewidth = I(typewidth)),
               arrow=arrow(length=unit(0.2,"cm")), color="black") +
  ggrepel::geom_label_repel(data=fig3brepeldata,
                            aes(x=PC1,y=PC2,
                                label=Group,
                                col = Group, 
                                fontface = I(typeface),
                                size = I(typesize)),
                            alpha = 0.75,
                            min.segment.length = 5,
                            max.overlaps = 100,
                            col = "black") +
  coord_cartesian(ylim=c(-6,8),
                  xlim=c(-10,4))  +
  scale_x_continuous(breaks = c(0),labels="0") +
  
  scale_y_continuous(breaks = c(0),labels="0") +
  scale_color_manual(values = c("grey60", "black")) +
  xlab("PC1") +
  ylab("PC2") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16))
Figure3brepel

Figure3_Legend <- cowplot::get_legend(
  envfitplot$plot +
    scale_color_manual("Habitat",
                       values = reefColors) + 
    theme(legend.key=element_blank(),
          legend.key.height = unit(0.15,"in"),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12)))
Figure3a_withLegend = cowplot::ggdraw() +
  cowplot::draw_plot(Figure3a, 0, 0, 1, 1)+
  cowplot::draw_plot(Figure3_Legend,
                     0.725, 0.75, 0.25, 0.25) 

Figure3 = cowplot::plot_grid(Figure3a_withLegend,
                             Figure3brepel,
                             nrow = 1,
                             labels = "auto")
ggsave("plots/figure3.jpg", Figure3, width = 12, height = 6, units = "in", dpi = 300)
ggsave("plots/figure3.png", Figure3, width = 12, height = 6, units = "in", dpi = 300)


#### Water chemistry data by year ----
## .... Plots of figure 3 building blocks ----
## Just the chemistry community loadings first
fig3_arrows = gginnards::delete_layers(Figure3b, idx = c(2,3,4,5,7))
fig3_arrows
ggsave("plots/waterchemistry/figure3_waterchemistry.jpg",
       fig3_arrows,
       width = 6, height = 6, units = "in",dpi = 300)

fig3_environ = gginnards::delete_layers(Figure3b, idx = c(4,5))
fig3_environ
ggsave("plots/waterchemistry/figure3_environmental.jpg",
       fig3_environ,
       width = 6, height = 6, units = "in",dpi = 300)

fig3_points = gginnards::delete_layers(Figure3a, idx = c(2,3,4:5,6,9:12))
fig3_points
ggsave("plots/waterchemistry/figure3_points.jpg",cowplot::ggdraw() +
         cowplot::draw_plot(fig3_points, 0, 0, 1, 1)+
         cowplot::draw_plot(Figure3_Legend, 0.1,0.1,0.25,0.25) ,
       width = 6, height = 6, units = "in",dpi = 300)

fig3_ellipse = gginnards::delete_layers(Figure3a, idx = c(2,9:12))
fig3_ellipse
ggsave("plots/waterchemistry/figure3_ellipse.jpg",cowplot::ggdraw() +
         cowplot::draw_plot(fig3_ellipse, 0, 0, 1, 1)+
         cowplot::draw_plot(Figure3_Legend, 0.1,0.1,0.25,0.25) ,
       width = 6, height = 6, units = "in",dpi = 300)

#### Supplementary materials S10 ----
## Supplemental version of Fig 3; broken down by year
generateYearOrdPlot = function(x = data_scaled, myYear = 2021,
                               writeplot = F){
  yrDatScaled = x %>%
    select(Site,
           Year,
           Date,
           Phosphate_scale:FI_scale,
           Percent_N,
           Nsite,
           detrended,
           Habitat,
           Watershed_area,
           Cleared_area, 
           Distance_from_shore) %>%
    mutate(logAreakm2 = log(Watershed_area+1),
           logDistance = log(Distance_from_shore+1),
           logClearedArea = log(Cleared_area+1),
           logTurbinariaN = log(Nsite+1)) %>% 
    filter(Year == myYear) %>%
    filter(complete.cases(.))
  
  yrPC = prcomp(yrDatScaled %>%
                  select(Phosphate_scale:FI_scale))
  
  ## Check: Does this pattern apply every year?
  yrEnv = envfit(yrPC ~ 
                   logDistance + logClearedArea + logTurbinariaN + detrended,
                 data = yrDatScaled, 
                 perm = 999)
  ## Compile loadings from PCA and from explanatory analysis
  yrmyspecies_expl = names(yrEnv$vectors$r[yrEnv$vectors$pvals<0.05])
  expl.yrEnv <- as.data.frame(
    scores(yrEnv, c("vectors"))) %>%
    mutate(species = rownames(.)) %>%
    filter(species %in% yrmyspecies_expl) %>%
    mutate(species = case_when(species == "areakm2" ~ "Watershed area",
                               species == "log(areakm2 + 1)" ~ "Watershed area",
                               species == "Distance" ~ "Dist. from shore",
                               species == "logDistance" ~ "Dist. from shore",
                               species == "log(Distance + 1)" ~ "Distance from shore",
                               species == "Cleared" ~ "Cleared (%)",
                               species == "log(Cleared + 1)" ~ "Cleared (%)",
                               species == "clearedArea" ~ "Cleared area",
                               species == "logClearedArea" ~ "Cleared area",
                               species == "log(clearedArea + 1)" ~ "Cleared area",
                               species == "Turbinaria_N" ~ "Turbinaria N mean",
                               species == "logTurbinariaN" ~ "Turbinaria N mean",
                               species == "log(Turbinaria_N + 1)" ~ "Turbinaria N mean",
                               species == "detrended" ~ "Residual Turbinaria N",
                               species == "Turbinaria_V" ~ "Turbinaria N variance",
                               TRUE ~ species),
           speciesP = yrEnv$vectors$pvals[yrEnv$vectors$pvals<0.05],
           printSpecies = paste0(species,ifelse(speciesP<=0.001,
                                                " ***",
                                                ifelse(speciesP <= 0.01,
                                                       " **",
                                                       ifelse(speciesP <= 0.05,
                                                              " *",
                                                              "")))),
           PClab1 = 15*PC1 + 0.4*sign(PC1),
           PClab2 = 15*PC2 + 0.4*sign(PC2),
           PCarrow1 = 15*PC1,
           PCarrow2 = 15*PC2)
  
  suppYearORD = gg_ordiplot(ord = yrPC,
                            groups = yrDatScaled$Habitat,
                            ellipse = T,
                            pt.size = 0.5,
                            # show.groups = F,
                            hull = F,
                            spiders = F,
                            plot = F)
  yrPC_centroids = yrDatScaled %>%
    cbind(yrPC$x) %>%
    group_by(Habitat) %>%
    summarize(PC1 = mean(PC1),
              PC2 = mean(PC2))
  
  yrPCInfo = as.data.frame(yrPC$rotation) %>%
    mutate("Group" = rownames(.)) %>%
    mutate(Group = gsub(".like_scale","",Group),
           Group = gsub("_scale","",Group))
  
  
  suppYearPERM = suppYearORD$plot +
    scale_color_manual(values = setNames(
      paste0(reefColors,"70"),
      names(reefColors))) +
    ggnewscale::new_scale_color() +
    geom_segment(data = yrPCInfo,
                 aes(x = 0, xend = PC1*10, y = 0, yend = PC2*10),
                 arrow=arrow(length=unit(0.2,"cm")), 
                 alpha=0.75, color="black") +
    geom_path(data = suppYearORD$df_ellipse,
              aes(x=x, y=y, group = Group),
              col = "black",
              lwd = 1.95) +
    geom_path(data = suppYearORD$df_ellipse,
              aes(x=x, y=y, color=Group),
              lwd = 1.5) +
    geom_label(data = yrPCInfo,
               aes(x = (PC1*10 + sign(PC1)*0.5), y = (PC2*10 + sign(PC2)*0.5),
                   label = Group),
               alpha = 0.15,
               col = "transparent") +
    geom_label(data = yrPCInfo,
               aes(x = (PC1*10 + sign(PC1)*0.5), y = (PC2*10 + sign(PC2)*0.5),
                   label = Group),
               alpha = 0.5,
               col = "transparent") +
    geom_label(data = yrPCInfo,
               aes(x = (PC1*10 + sign(PC1)*0.5), y = (PC2*10 + sign(PC2)*0.5),
                   label = Group),
               alpha = 0.25,
               col = "black") +
    geom_point(data = yrPC_centroids,
               aes(x = PC1, y = PC2, col = Habitat),
               col = "black",
               size = 6) +
    geom_point(data = yrPC_centroids,
               aes(x = PC1, y = PC2, col = Habitat),
               size = 5) +
    scale_color_manual(values = reefColors) +
    # coord_cartesian(xlim = c(-10,5),
    #                 ylim = c(-10,7)) +
    ggtitle(myYear)  +
    theme_bw() +
    theme(legend.position = "none",
          plot.title = element_text(size = 14))
  suppYearPERM
  
  suppYrEnvRepelData = yrPCInfo %>%
    select(Group, PC1, PC2)  %>% 
    mutate(type = "chem",
           PC1 = PC1*10,
           PC2 = PC2*10) %>%
    rbind(expl.yrEnv %>% mutate(Group = printSpecies) %>% 
            select(Group, PCarrow1, PCarrow2) %>% 
            rename("PC1" = PCarrow1,
                   "PC2" = PCarrow2) %>%
            mutate(type = "env")) %>%
    mutate(typealpha = ifelse(type == "chem", 0.5, 1),
           typesize = ifelse(type == "chem", 4, 5),
           typeface = ifelse(type == "chem", 1, 2),
           typewidth = ifelse(type == "chem", 0.25, 1.25))
  
  yrEnvfitplot = gg_ordiplot(ord = yrPC,
                             groups = yrDatScaled$Habitat,
                             ellipse = F,
                             pt.size = 0.5,
                             hull = F,
                             spiders = F,
                             plot = F)
  suppYearENV = yrEnvfitplot$plot +
    geom_segment(data = yrPCInfo,
                 aes(x = 0, xend = PC1*10, y = 0, yend = PC2*10),
                 arrow=arrow(length=unit(0.2,"cm")), 
                 alpha=0.75, color="black") +
    geom_segment(data=expl.yrEnv,
                 aes(x=0,xend=PCarrow1,y=0,yend=PCarrow2),
                 arrow = arrow(length = unit(0.3, "cm")),
                 lwd = 1.5,
                 col = "black") +
    geom_segment(data=expl.yrEnv,
                 aes(x=0,xend=PCarrow1,y=0,yend=PCarrow2),
                 arrow = arrow(length = unit(0.3, "cm")),
                 lwd = 1,
                 col = "black") +
    ggrepel::geom_label_repel(data=suppYrEnvRepelData,
                              aes(x=PC1,y=PC2,
                                  label=Group,
                                  col = Group, 
                                  fontface = I(typeface),
                                  size = I(typesize)),
                              alpha = 0.75,
                              min.segment.length = 5,
                              max.overlaps = 100,
                              col = "black") +
    scale_color_manual(values = reefColors) +
    theme_bw() +
    theme(legend.position = "none")
  
  
  suppYear_a_withLegend <- cowplot::ggdraw() +
    cowplot::draw_plot(suppYearPERM, 0, 0, 1, 1) +
    cowplot::draw_plot(Figure3_Legend,
                       0.1, 0.1, 0.25, 0.25) 
  suppYear_full = cowplot::plot_grid(suppYearPERM,
                                     suppYearENV,
                                     nrow = 1, align = "h", axis = "tb")
  if(writeplot){
    ggsave(paste0("plots/supplement/waterchemisty_ordination_",myYear,".jpg"),suppYear_full,
           width = 12, height = 6, units = "in",dpi = 300)
  }
  return(list(suppYearPERM, suppYearENV))
}

## Within individual years, water column nutrients were consistently positively associated with Turbinaria %N and negatively associated with distance from shore every year, and also positively associated with cleared area of the upstream watershed in 2021 (Supplementary materials S10; p < 0.01 in all cases). 
yrOrd2021 = generateYearOrdPlot(x = data_scaled, myYear = 2021)
yrOrd2022 = generateYearOrdPlot(x = data_scaled, myYear = 2022)
yrOrd2023 = generateYearOrdPlot(x = data_scaled, myYear = 2023)

allYearsOrdPlots = cowplot::plot_grid(yrOrd2021[[1]] + coord_cartesian(), 
                                      yrOrd2021[[2]] + coord_cartesian(),
                                      yrOrd2022[[1]] + coord_cartesian(), 
                                      yrOrd2022[[2]] + coord_cartesian(),
                                      yrOrd2023[[1]] + coord_cartesian(), 
                                      yrOrd2023[[2]] + coord_cartesian(),
                                      align = "hv",
                                      axis = "trbl",
                                      nrow = 3, byrow = T)
ggsave(filename = "plots/supplement/s10_waterchemistry_allYearsOrdPlots.jpg",
       plot = allYearsOrdPlots,
       width = 12, height = 14)
