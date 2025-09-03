## This script generates Supplementary materials s9
## Water column chemistry plots and tables

## Author: Christian John
## Date: 14 May 2025


#### Load packages ----
library(tidyverse) # Data wrangling
library(cowplot) # Plot wrangling
library(sf) # Simple features
library(multcomp) # generate plot letters with emmeans
library(lme4) # fitting linear mixed effects models
library(lmerTest) # lmer significance testing
library(emmeans) # post-hoc comparisons
source("scripts/000_reefCols.R")

set.seed(1992) #R's inception

#### Import data ----
## ... Water data ----
water_N = read_csv("data/generated/water_N_compiled.csv") 

## ... Vector data ----
## Sampling locations: This is used for identifying sampling location coordinates
samplingLocations = read_csv(
  "~/Documents/Burkepile_Lab/Data/ATI/EDI_site_locations/data/clean/sites_metadata_DRAFT.csv") %>%
  mutate(Habitat = factor(Habitat, levels = names(reefColors))) %>%
  st_as_sf(coords = c("Longitude","Latitude")) %>%
  st_set_crs(4326)

sitesLocs = unique(samplingLocations$Site)


#### Prep data ----
## Fuse habitat info with water chemistry data
water_sf = right_join(samplingLocations,
                      water_N)

#### Plot results ----
## A function to generate boxplots with emmeans significance letters
## for a given water chemistry variable
plotWaterChem = function(x = water_sf,
                         H2Ovar = "Phosphate"){
  
  ## Make variables look nicer
  prettyvar = gsub("M.C","M:C",H2Ovar)
  prettyvar = gsub("_"," ",prettyvar)
  prettyvar = gsub("[.]"," ",prettyvar)
  prettyvar = gsub(" like","-like",prettyvar)
  
  ## Create data.frame for individual variable
  vardf = x %>% 
    dplyr::select(Site, all_of(H2Ovar), Year, Habitat, Island_shore) %>%
    st_drop_geometry() %>%
    rename(val = 2) %>%
    dplyr::filter(complete.cases(.))
  
  ## Mixed effects model
  varlmer = lmer(val ~ Habitat + (1|Year), data = vardf)
  summary(varlmer)
  varlm = lm(val ~ Habitat*Year, data = vardf)
  summary(aov(varlm))
  H2Ostats = list(varlmer, varlm)
  
  ## Post-hoc comparisons
  ## Using lmer's so we can compare just among habitat after
  ## accounting for interannual differences
  em_var_lmer = emmeans(varlmer, "Habitat")
  var_lmer_cld = cld(em_var_lmer, Letters = letters, reversed=T) %>%
    dplyr::select(Habitat, emmean, .group) %>%
    mutate(label = gsub(" ","",.group))
  
  H2Oboxplots = vardf %>%
    group_by(Site, Habitat) %>%
    summarize(val = mean(val)) %>%
    ggplot(aes(x = Habitat,
               y = val)) +
    geom_boxplot(aes(col = Habitat),
                 outlier.color = "transparent",
                 show.legend = F) +
    geom_jitter(aes(col = Habitat),
                shape = 1,
                cex = 2,
                stroke = 1,
                width = 0.2,
                show.legend = F) +
    geom_label(data = var_lmer_cld,
               aes(x = Habitat, 
                   y = 0, 
                   label = label),
               label.size = 0, size = 4.5,
               alpha = 0) +
    ylab(prettyvar) +
    xlab("") +
    scale_x_discrete(labels = scales::label_wrap(7)) +
    scale_color_manual(values = reefColors) +
    CJsBasics::BasicTheme + 
    theme(axis.title = element_text(size = 18),
          axis.text = element_text(size = 16),
          strip.text = element_text(size = 18))
  return(list(H2Oboxplots,
              H2Ostats))
}


## ... Apply plotting function ----
## Nutrients
H2O_pho = plotWaterChem(x = water_sf, H2Ovar = "Phosphate")
H2O_sil = plotWaterChem(x = water_sf, H2Ovar = "Silicate")
H2O_nit = plotWaterChem(x = water_sf, H2Ovar = "Nitrite_plus_Nitrate")
## fDOM components
H2O_uvh = plotWaterChem(x = water_sf, H2Ovar = "Ultra.Violet.Humic.like")
H2O_vhl = plotWaterChem(x = water_sf, H2Ovar = "Visible.Humic.like")
H2O_mhl = plotWaterChem(x = water_sf, H2Ovar = "Marine.Humic.like")
H2O_tyr = plotWaterChem(x = water_sf, H2Ovar = "Tyrosine.like")
H2O_try = plotWaterChem(x = water_sf, H2Ovar = "Tryptophan.like")
H2O_lig = plotWaterChem(x = water_sf, H2Ovar = "Lignin.like")
## fDOM indices
H2O_hix = plotWaterChem(x = water_sf, H2Ovar = "HIX")
H2O_bix = plotWaterChem(x = water_sf, H2Ovar = "BIX")
H2O_m.c = plotWaterChem(x = water_sf, H2Ovar = "M.C")
H2O_fix = plotWaterChem(x = water_sf, H2Ovar = "FI")


#### Wrangle plots ----
cowplot::plot_grid(H2O_pho[[1]],
                   H2O_sil[[1]],
                   H2O_nit[[1]],
                   nrow = 1,
                   ncol = 3) %>%
  ggsave(filename = paste0("plots/supplement/s9_nutrients.jpg"),
         width = 12, height = 4, units = "in", dpi = 300)

cowplot::plot_grid(H2O_uvh[[1]],
                   H2O_vhl[[1]],
                   H2O_mhl[[1]],
                   H2O_tyr[[1]],
                   H2O_try[[1]],
                   H2O_lig[[1]],
                   nrow = 2,
                   ncol = 3) %>%
  ggsave(filename = paste0("plots/supplement/s9_fDOM_components.jpg"),
         width = 14, height = 8, units = "in", dpi = 300)

cowplot::plot_grid(H2O_hix[[1]],
                   H2O_bix[[1]],
                   H2O_m.c[[1]],
                   H2O_fix[[1]],
                   nrow = 1,
                   ncol = 4) %>%
  ggsave(filename = paste0("plots/supplement/s9_fDOM_indices.jpg"),
         width = 16, height = 4, units = "in", dpi = 300)


#### Generate stat outputs ----
myAOVsummary = function(myLM = lm_uvhl,
                        longName = "Ultra-violet humic-like"){
  myAOV = aov(myLM)
  myAOVsum = summary(myAOV)
  data.frame("term" = rownames(myAOVsum[[1]]),
             data.frame(myAOVsum[[1]])) %>%
    remove_rownames() %>%
    mutate(term = gsub(" ", "", term),
           Sum.Sq = round(Sum.Sq, 5),
           Mean.Sq = round(Mean.Sq, 5),
           `F` = round(F.value, 3),
           p = round(`Pr..F.`, 3),
           ptags = CJsBasics::pStars(p),
           Parameter = longName) %>%
    mutate(Parameter = ifelse(term == "Habitat", Parameter, "")) %>%
    replace(is.na(.), "") %>%
    dplyr::select(Parameter, term, Df, Sum.Sq, Mean.Sq, `F`, p, ptags) %>%
    rbind(rep("",ncol(.)))
}


aovs_components = rbind(myAOVsummary(H2O_uvh[[2]][[2]], "Ultra-violet humic-like"),
                        myAOVsummary(H2O_vhl[[2]][[2]], "Visible humic-like"),
                        myAOVsummary(H2O_mhl[[2]][[2]], "Marine humic-like"),
                        myAOVsummary(H2O_tyr[[2]][[2]], "Tyrosine-like"),
                        myAOVsummary(H2O_try[[2]][[2]], "Tryptophan-like"),
                        myAOVsummary(H2O_lig[[2]][[2]], "Lignin-like"))

aovs_indices = rbind(myAOVsummary(H2O_hix[[2]][[2]], "HIX"),
                     myAOVsummary(H2O_bix[[2]][[2]], "BIX"),
                     myAOVsummary(H2O_m.c[[2]][[2]], "M:C"),
                     myAOVsummary(H2O_fix[[2]][[2]], "FI"))

aovs_nutrients = rbind(myAOVsummary(H2O_nit[[2]][[2]], "Nitrate"),
                       myAOVsummary(H2O_sil[[2]][[2]], "Silicate"),
                       myAOVsummary(H2O_pho[[2]][[2]], "Phosphate"))

sjPlot::tab_df(aovs_components, 
               col.header = c(names(aovs_components)[1:7],""),
               file = "statOutputs/supplement/s9/aov_fdomComponents.doc")
sjPlot::tab_df(aovs_indices, 
               col.header = c(names(aovs_indices)[1:7],""),
               file = "statOutputs/supplement/s9/aov_fdomIndices.doc")
sjPlot::tab_df(aovs_nutrients, 
               col.header = c(names(aovs_nutrients)[1:7],""),
               file = "statOutputs/supplement/s9/aov_nutrients.doc")

