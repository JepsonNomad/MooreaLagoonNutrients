#### microbial community analyses
#### this script loads & merges phyloseq objects for 2021-2023 lagoon microbial communities
#### and regresses microbial community ordinations against environmental factors

#### Load packages ----
library(tidyverse)
library(phyloseq)
library(microViz)
library(vegan)
library(AICcPermanova)
library(multcompView)
library(cowplot)
source("scripts/000_reefCols.R")

#### Set seed ----
seed = T
if(seed){
  set.seed(140) # Calories in a Yuengling
}


#### Define functions -----
## Convert phyloseq objects to vegan OTU's based on McCurdie and Holmes 2013
## https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0061217#s6
## This function came from their supplementary materials S2
veganotu <- function(physeq) {
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

## A function to create p value asterisk strings from envfit objects
get_ptags = function(x){
  pvals_yr = x$vectors$pvals %>%
    set_names(NULL)
  ptags_yr = CJsBasics::pStars(pvals_yr)
  return(ptags_yr)
}

## A ggplot-friendly Tukey test visualizer
ggTukey <- function(myTukey){
  df_tuk = as.data.frame(myTukey[1])
  colnames(df_tuk)[2:3]<- c("min",
                            "max")
  df_plot = data.frame(id = row.names(df_tuk),
                       min = df_tuk$min,
                       max = df_tuk$max)
  plot_tukey = df_plot %>%
    ggplot() +
    geom_errorbar(aes(x = id,
                      ymin=min,
                      ymax=max),
                  width = 0.2)+
    geom_hline(yintercept=0,
               color="grey60")+
    labs(x = NULL)+
    coord_flip()+
    # scale_y_continuous(breaks = c(-0.3,0,0.3)) +
    CJsBasics::BasicTheme+
    theme(plot.title = element_text(hjust = 1))
  return(plot_tukey)
}

#### Import data ----
## Turbinaria data
## (one measure per site per year)
## dataTurbInSitu <- read_csv("data/generated/turb_N_imputed.csv")

## Detrended Turbinaria (accounting for distance from shore)
## (one measure per site)
dataTurbDetrended <- read_csv("data/generated/detrendedNutrients.csv") %>%
  rename(Turbinaria_N = Percent_N) %>%
  select(Site, Turbinaria_N, detrended)

## Water data
dataWater<-read_csv("data/generated/water_N_compiled.csv") %>%
  group_by(Site) %>%
  summarize(Nitrite_plus_Nitrate = mean(Nitrite_plus_Nitrate),
            Silicate = mean(Silicate),
            Phosphate = mean(Phosphate),
            HIX = mean(HIX),
            BIX = mean(BIX),
            FI = mean(FI),
            M.C = mean(M.C),
            Ultra.Violet.Humic.like = mean(Ultra.Violet.Humic.like),
            Visible.Humic.like = mean(Visible.Humic.like),
            Marine.Humic.like = mean(Marine.Humic.like),
            Lignin.like = mean(Lignin.like),
            Tryptophan.like = mean(Tryptophan.like),
            Tyrosine.like = mean(Tyrosine.like))

## Site metadata
dataSitemeta <- read_csv("~/Documents/Burkepile_Lab/Data/ATI/EDI_site_locations/data/clean/sites_metadata_DRAFT.csv")

sitedata = dataTurbDetrended %>%
  full_join(dataWater) %>%
  full_join(dataSitemeta) %>%
  filter(complete.cases(.))

## > Import microbial data ----
## 2021
r2021 <- 
  readr::read_rds("~/Documents/Burkepile_Lab/Data/ATI/EDI_microbes/data/clean/ati-2021-physeq-rare.RDS")
sample_names(r2021) <- paste0("2021-", samdat_tbl(r2021)$Site)
p2021 = r2021 %>%
  phyloseq::phyloseq(
    otu_table(.),
    tax_table(.),
    sample_data(.)) %>%
  ps_join(sitedata, by = "Site") %>%
  ps_mutate(Year = 2021)

## 2022
r2022 <- 
  readr::read_rds("~/Documents/Burkepile_Lab/Data/ATI/EDI_microbes/data/clean/ati-2022-physeq-rare.RDS")
sample_names(r2022) <- paste0("2022-", samdat_tbl(r2022)$Site)
p2022 = r2022 %>%
  phyloseq::phyloseq(
    otu_table(.),
    tax_table(.),
    sample_data(.)) %>%
  ps_join(sitedata, by = "Site") %>%
  ps_mutate(Year = 2022)

## 2023
r2023 <- 
  readr::read_rds("~/Documents/Burkepile_Lab/Data/ATI/EDI_microbes/data/clean/ati-2023-physeq-rare.RDS")
sample_names(r2023) <- paste0("2023-", samdat_tbl(r2023)$Site)
p2023 = r2023 %>%
  phyloseq::phyloseq(
    otu_table(.),
    tax_table(.),
    sample_data(.)) %>%
  ps_join(sitedata, by = "Site") %>%
  ps_mutate(Year = 2023)

## Fuse the three years of microbial data
microbes_raw = p2021 %>%
  merge_phyloseq(p2022) %>%
  merge_phyloseq(p2023) %>%
  ps_arrange(Site, Year)

## Merge by site so that we have one "community" per site
microbes_merged = microbes_raw %>%
  merge_samples("Site") %>%
  ps_select(-Site, -Habitat, -Island_shore) %>%
  ps_mutate(Site = sample_names(.)) %>%
  ps_join(sitedata %>% select(Site, Habitat, Island_shore), 
          by = "Site") %>%
  ps_filter(complete.cases(samdat_tbl(.)))
## Extract site ID's so we can work with them down the line
merged_siteIDs =  sample_names(microbes_merged)

## Check the merge. Note that none of the metadata in this phyloseq object will 
## be retained.
microbes_merged %>% samdat_tbl() # %>% View()

## Validate to check for consistent sample names
microbes_valid <- phyloseq_validate(microbes_merged, remove_undetected = TRUE)

## Join together the site metadata and environmental data
microbes_prep = microbes_valid %>%
  # ps_mutate(site = merged_siteIDs) %>%
  # ps_select(site) %>% 
  # ps_join(dataTurbi) %>%
  # ps_join(dataWater) %>%
  # ps_join(dataLanduse) %>%
  # ps_join(dataShoredist) %>%
  # ps_join(dataHabitat) %>%
  tax_fix(unknowns = c("endosymbionts", "uncultured",
                       "Unknown_Family", "Unknown_Family Family")) %>%
  ps_filter(!is.na(Turbinaria_N)) #%>%
## Can look at just one habitat, e.g. Fringing reef.
# ps_filter(Habitat == "Fringing reef")

## Order is retained between the following two overrides so should be ok.
rownames(microbes_prep@sam_data) <- sample_names(microbes_prep)

## When feeding into a bray-curtis dissimilarity matrix, no transformation is
## needed but it's best to include the identity transformation just so we 
## don't forget this step in other situations where a (e.g. log) transformation
## is warranted
taxLevel = "Genus"
microbes_transform = microbes_prep %>%
  tax_transform(trans = "identity", rank = taxLevel) 
## Calculate distance
microbes_distance = microbes_transform %>%
  dist_calc("bray") 

## .... Ordinate ----
## Ordinate using the approaches selected in the function call
## Ordinate using PCoA with the bray-curtis dissimilarity matrix
microbes_PCoA = microbes_distance %>%
  ord_calc(method = "PCoA") 


## .... Explanatory variables ----
## Use bioenv to test all environmental variables simultaneously
allPredictors = c("Habitat",
                  "logShoreDist",
                  "logClearedArea",
                  "Turbinaria_N",
                  "detrended",
                  "Nitrite_plus_Nitrate",
                  "Silicate",
                  "Phosphate",
                  "Lignin.like",
                  "Tyrosine.like",
                  "Tryptophan.like",
                  "Visible.Humic.like",
                  "Marine.Humic.like",
                  "Ultra.Violet.Humic.like",
                  "HIX",
                  "BIX",
                  "M.C",
                  "FI")
myPredictors = allPredictors[c(4,6:8,15,16)]

env = samdat_tbl(microbes_transform)[,myPredictors]
# microbes_bioenv = bioenv(veganotu(microbes_transform), env)
# microbes_envfit = envfit(veganotu(microbes_transform)[complete.cases(env),],
#                          env[complete.cases(env),])
# predictors_bioenv = myPredictors[microbes_bioenv$models[[microbes_bioenv$whichbest]]$best]
# sink(file = paste0("statOutputs/allYears_bioenv_",taxLevel,".txt"))
# print(microbes_bioenv)
# print(summary(microbes_bioenv))
# microbes_PCoA %>% 
#   dist_permanova(variables = predictors_bioenv) %>%
#   print()
# print(microbes_envfit)
# sink(file= NULL)


## Use permanova to test for environmental relationships with community diversity
## Test distance matrix as response to environmental variables. 
envi_Complete = env %>%
  filter(complete.cases(.))
bray_Complete = as.matrix(microbes_PCoA@dist)[(1:nrow(env))[complete.cases(env)],
                                              (1:nrow(env))[complete.cases(env)]]

## Microbial communities differed by habitat (Fig. 5a, PERMANOVA, p = 0.001).
adonis2(
  bray_Complete ~ 
    Habitat,
  data = samdat_tbl(microbes_transform)[,"Habitat"])


#### Model selection ----
## Create all possible models from candidate predictors
AllModels <- make_models(vars = myPredictors)
AllModels
## Filter to include only non-colinear models
NonColinear <- filter_vif(all_forms = AllModels, 
                          env_data = envi_Complete,
                          ncores = 4,
                          verbose = T)
NonColinear
## Fit each of the filtered models
Fitted <- fit_models(
  all_forms = NonColinear,
  com_data = bray_Complete,
  env_data = envi_Complete,
  ncores = 4,
  method = "bray",
  verbose = T
)
## Select the best model(s)
Selected <- select_models(Fitted)
topModVars = Selected[1,1] %>%
  gsub(pattern = "Distance ~ ", replacement = "") %>%
  gsub(pattern = " ", replacement = "") %>%
  strsplit(split = "\\+") %>%
  unlist()
topModVars
Summary <- akaike_adjusted_rsq(Selected)
Summary

## Based on the model selection framework, patterns in microbial communities were best explained by Turbinaria %N and water column nitrite + nitrate, silicate, phosphate, and HIX (Table 1; p < 0.01 in all cases). 
topModel = adonis2(
  bray_Complete ~ 
    Turbinaria_N +
    Nitrite_plus_Nitrate + 
    Silicate +
    Phosphate +
    HIX,
  data = envi_Complete)
topModel

#### Table 1 ----
# Generate output txt file
sink(file = paste0("statOutputs/allYears_permanova_PCoA_",taxLevel,".txt"))
## Summarize variability along environmental axes and by habitat
print(topModel)
sink(file= NULL)
topModelTable = as.data.frame(topModel) %>%
  mutate(Variable = rownames(.),
         SumOfSqs = round(SumOfSqs, 2),
         R2 = round(R2, 2),
         `F` = round(`F`, 2)) %>%
  select(Variable, Df, SumOfSqs, R2, `F`, `Pr(>F)`)
write_csv(topModelTable,
          paste0("statOutputs/allYears_permanova_PCoA_",taxLevel,".csv"))

## Run envfit on the rotated coordinates of the microbial community
## In order to visualize environ. vectors in this ordination space
microbes_PCoA_ordination = ord_get(microbes_PCoA)
ordsum = summary(microbes_PCoA_ordination, scaling = 2, axes = 2)
microbes_PCoA_envfit = envfit(
  microbes_PCoA@ord, 
  samdat_tbl(microbes_transform)[,topModVars], 
  na.rm = T,
  perm = 999)
microbes_PCoA_envfit_ptags = get_ptags(microbes_PCoA_envfit)

microbes_PCoA_envfit_df <- as.data.frame(
  scores(microbes_PCoA_envfit, c("vectors"))) %>%
  mutate(species = rownames(.),
         sp_ptag = microbes_PCoA_envfit_ptags) %>%
  filter(sp_ptag != "",
         sp_ptag != ".") %>%
  mutate(species = case_when(species == "logShoreDist" ~ "Dist. from shore",
                             species == "logClearedArea" ~ "Cleared area",
                             species == "detrended" ~ "Detrended Turbinaria %N",
                             species == "Turbinaria_N" ~ "Turbinaria %N",
                             species == "Visible.Humic.like" ~ "Visible humic",
                             species == "Marine.Humic.like" ~ "Marine humic",
                             species == "Ultra.Violet.Humic.like" ~ "Ultra-violet humic",
                             species == "Tyrosine.like" ~ "Tyrosine-like",
                             species == "Tryptophan.like" ~ "Tryptophan-like",
                             species == "Lignin.like" ~ "Lignin-like",
                             TRUE ~ species),
         species_print = paste0(species, sp_ptag),
         xoffset = 0,
         yoffset = 0,
         MDSlab1 = 2*MDS1 + 0.1*sign(MDS1) + xoffset,
         MDSlab2 = 2*MDS2 + 0.1*sign(MDS2) + yoffset,
         MDSarrow1 = 2*MDS1 + xoffset,
         MDSarrow2 = 2*MDS2 + yoffset)

HabitatCentroids = ordsum$sites %>% 
  as.data.frame() %>%
  mutate(Site = rownames(.)) %>%
  # mutate(sample_name = rownames(.)) %>%
  # separate(sample_name, into = c("Year","site"), sep = "-") %>%
  left_join(sitedata) %>%
  select(Site, Habitat, MDS1, MDS2) %>%
  group_by(Habitat) %>%
  summarize(MDS1 = mean(MDS1),
            MDS2 = mean(MDS2))

#### Dispersion test ----
## While mid lagoon, back reef, and reef pass sites had similar microbial communities (indicated by highly overlapping regions in PCoA space), fringing reefs and bays featured divergent microbial communities.
## Ask: Do different habitats have different "spreads" in the ordination?
## Calculate PERMDISP (i.e. betadisper in vegan)
dispersionallYears = betadisper(
  microbes_PCoA@dist,
  samdat_tbl(microbes_PCoA)$Habitat,
  type = "median",
  bias.adjust = T,
  add = FALSE) 

## An analysis of homogeneity of group dispersions showed that fringing reefs had more variable communities than mid-lagoon and back reef sites, and bays and reef passes had more variable communities than backreef sites (p < 0.05 in all cases). 
## Use a Tukey test to see how different groups differ
tukeydisp = dispersionallYears %>% TukeyHSD()
tukeydisp
## Generate multiple comparison letter labels for plotting
mic_cld <- multcompLetters(tukeydisp$group[,4])
mic_complabels = mic_cld$Letters %>%
  as.data.frame() %>%
  rename("Letters" = 1) %>%
  mutate(Habitat = rownames(.),
         Habitat = factor(Habitat, levels = reefLevels)) %>%
  arrange(Habitat) %>%
  mutate(HabitatLabs = paste0(Habitat, " (", Letters, ")"))
## Create a data.frame with the tukey object and add significance stars
tukeyDat = tukeydisp$group %>%
  as.data.frame() %>%
  mutate(sig = CJsBasics::pStars(`p adj`),
         group = rownames(.))
## Plot the group differences
tukeyDispPlot = ggTukey(tukeydisp) +
  geom_text(data = tukeyDat,
            aes(label = sig,
                x = group,
                y = -0.275),
            hjust = 0,
            nudge_x = -0.1) +
  scale_y_continuous(breaks = c(-0.2,0,0.2),
                     lim = c(-0.275,NA)) +
  ylab("Group differences")
tukeyDispPlot
ggsave("plots/microbes/allYears_habitat_permdisp_tukey.jpg",
       tukeyDispPlot,
       height = 4, width = 4, units = "in")


## Get info from PCoA object for plot
biplot_scores = vegan::scores(ord_get(microbes_PCoA), display = "sites")
biplot_eigVals <- vegan::eigenvals((ord_get(microbes_PCoA)))
biplot_explainedVar <- ((biplot_eigVals[c(1,2)] / sum(biplot_eigVals))*100) %>%
  round(digits = 1)

#### Figure 5a ----
## Plot PCoA data
pcoa_loadings_biplot_points = biplot_scores %>%
  as_tibble() %>%
  mutate(Site = samdat_tbl(microbes_PCoA)$Site,
         Habitat = samdat_tbl(microbes_PCoA)$Habitat,
         Habitat = factor(Habitat, levels = reefLevels)) %>%
  ggplot(aes(x = MDS1,
             y = MDS2)) +
  geom_point(aes(col = Habitat),
             shape = 16,
             size = 2,
             show.legend = F) +
  scale_color_manual("",
                     values = setNames(
                       paste0(reefColors,"70"),
                       names(reefColors))) +
  ggnewscale::new_scale_color() +
  stat_chull(aes(col = Habitat,
                 fill = Habitat),
             alpha = 0.1,
             lwd = 0.05,
             show.legend = F) +
  geom_point(data = HabitatCentroids,
             aes(x = MDS1, y = MDS2),
             col = "black",
             size = 5.5,
             show.legend = F) +
  geom_point(data = HabitatCentroids,
             aes(x = MDS1, y = MDS2,
                 col = Habitat),
             size = 4) +
  scale_color_manual("",
                     values = reefColors,
                     labels = mic_complabels$HabitatLabs) +
  scale_fill_manual("",
                    values = reefColors,na.value = "transparent") +
  scale_y_continuous(name = paste0("MDS2 (", biplot_explainedVar[2],"%)"),
                     lim = c(min(c(ordsum[["sites"]][,"MDS2"],
                                   microbes_PCoA_envfit_df$MDSlab2),
                                 na.rm = T) - 0.05,
                             max(c(ordsum[["sites"]][,"MDS2"],
                                   microbes_PCoA_envfit_df$MDSlab2),
                                 na.rm = T) + 0.05)) +
  scale_x_continuous(name = paste0("MDS1 (", biplot_explainedVar[1],"%)"),
                     lim = c(min(c(ordsum[["sites"]][,"MDS1"],
                                   microbes_PCoA_envfit_df$MDSlab1),
                                 na.rm = T) - 0.05,
                             max(c(ordsum[["sites"]][,"MDS1"],
                                   microbes_PCoA_envfit_df$MDSlab1),
                                 na.rm = T) + 0.05)) +
  theme_bw() +
  theme(legend.position = c(0.79,0.85),
        legend.key.height = unit(0.001,"in"),
        legend.key.width = unit(0.001,"in"),
        legend.text = element_text(size = 10,
                                   hjust = 1),
        legend.text.position = "left",
        legend.key.spacing.y = unit(0.001,"in"),
        legend.spacing.x = unit(-0.1,"in"), 
        legend.title = element_blank(),
        plot.title = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.border = element_rect(color="black", linewidth = 1))
pcoa_loadings_biplot_points

## Figure 5a. Create an equivalent biplot but with the axes flipped so that MDS1
## is the vertical axis, to facilitate intercomparison with MDS1 vs distance 
## from shore plot in Fig 5 b.
pcoa_loadings_biplot_points_flipped = biplot_scores %>%
  as_tibble() %>%
  mutate(Site = samdat_tbl(microbes_PCoA)$Site,
         Habitat = samdat_tbl(microbes_PCoA)$Habitat,
         Habitat = factor(Habitat, levels = reefLevels)) %>%
  ggplot(aes(x = MDS1,
             y = MDS2)) +
  annotate(geom = "rect",
           ymin = 0.55, ymax = 1.75, xmin = 0.275, xmax = 0.775,
           col = "white",
           fill = "white") +
  geom_point(aes(col = Habitat),
             shape = 16,
             size = 2,
             show.legend = F) +
  scale_color_manual("",
                     values = setNames(
                       paste0(reefColors,"70"),
                       names(reefColors))) +
  ggnewscale::new_scale_color() +
  stat_chull(aes(col = Habitat,
                 fill = Habitat),
             alpha = 0.1,
             lwd = 0.05,
             show.legend = F) +
  geom_point(data = HabitatCentroids,
             aes(x = MDS1, y = MDS2),
             col = "black",
             size = 6,
             show.legend = F) +
  geom_point(data = HabitatCentroids,
             aes(x = MDS1, y = MDS2,
                 col = Habitat),
             size = 5) +
  geom_text(data = HabitatCentroids %>% left_join(mic_complabels),
            aes(x = MDS1, y = MDS2,
                label = Letters),
            col = "white",
            fontface = "bold",
            size = 4) +
  scale_color_manual("",
                     values = reefColors) +
  scale_fill_manual("",
                    values = reefColors,na.value = "transparent") +
  scale_y_continuous(name = paste0("MDS2 (", biplot_explainedVar[2],"%)"),
                     lim = c(min(c(ordsum[["sites"]][,"MDS2"],
                                   microbes_PCoA_envfit_df$MDSlab2),
                                 na.rm = T) - 0.05,
                             max(c(ordsum[["sites"]][,"MDS2"],
                                   microbes_PCoA_envfit_df$MDSlab2),
                                 na.rm = T) + 0.05)) +
  scale_x_continuous(name = paste0("MDS1 (", biplot_explainedVar[1],"%)"),
                     lim = c(min(c(ordsum[["sites"]][,"MDS1"],
                                   microbes_PCoA_envfit_df$MDSlab1),
                                 na.rm = T) - 0.05,
                             max(c(ordsum[["sites"]][,"MDS1"],
                                   microbes_PCoA_envfit_df$MDSlab1),
                                 na.rm = T) + 0.05)) +
  theme_bw() +
  theme(legend.position = c(0.81,0.85),
        legend.key.height = unit(0.001,"in"),
        legend.key.width = unit(0.001,"in"),
        legend.text = element_text(size = 10,
                                   hjust = 1),
        legend.text.position = "left",
        legend.key.spacing.y = unit(0.001,"in"),
        legend.spacing.x = unit(-0.1,"in"), 
        legend.background = element_blank(),
        legend.title = element_blank(),
        plot.title = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.border = element_rect(color="black", linewidth = 1)) +
  coord_flip()
pcoa_loadings_biplot_points_flipped

pcoa_loadings_biplot_vectors = microbes_PCoA_envfit_df %>%
  ggplot(aes(x = MDS1,
             y = MDS2)) +
  geom_segment(aes(x=xoffset,xend=MDSarrow1,
                   y=yoffset,yend=MDSarrow2),
               arrow = arrow(length = unit(0.3, "cm")),
               lwd = 0.75,
               col = "black") +
  ggrepel::geom_label_repel(data=microbes_PCoA_envfit_df,
                            aes(x=MDSarrow1,y=MDSarrow2,
                                label=species_print),
                            min.segment.length = 5,
                            size=4, max.overlaps = 10,
                            col = "black",label.size = 0,
                            fill = "transparent") +
  scale_y_continuous(breaks = 0,
                     lim = c(min(c(ordsum[["sites"]][,"MDS2"],
                                   microbes_PCoA_envfit_df$MDSlab2),
                                 na.rm = T) - 0.05,
                             max(c(ordsum[["sites"]][,"MDS2"],
                                   microbes_PCoA_envfit_df$MDSlab2),
                                 na.rm = T) + 0.05)) +
  scale_x_continuous(breaks = 0,
                     lim = c(min(c(ordsum[["sites"]][,"MDS1"],
                                   microbes_PCoA_envfit_df$MDSlab1),
                                 na.rm = T) - 0.05,
                             max(c(ordsum[["sites"]][,"MDS1"],
                                   microbes_PCoA_envfit_df$MDSlab1),
                                 na.rm = T) + 0.05)) +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.border = element_rect(color="black", linewidth = 1))
pcoa_loadings_biplot_vectors

pcoa_loadings_stackedBiPlot = cowplot::plot_grid(
  pcoa_loadings_biplot_points + coord_flip() + theme(legend.position = c(0.815,0.85), legend.text = element_text(size = 9)),
  pcoa_loadings_biplot_vectors + coord_flip(),
  nrow = 2,
  align = "hv",
  labels = c("a","b")
)
pcoa_loadings_stackedBiPlot




#### Figure 5b ----
## Directly relating MDS axes with environmental variables
siteScoresDf <- as.data.frame(
  vegan::scores(ord_get(microbes_PCoA), display = "sites")
)
directComparisons = samdat_tbl(microbes_PCoA) %>%
  mutate(MDS1 = siteScoresDf$MDS1,
         MDS2 = siteScoresDf$MDS2) 

## Together, distance from shore and residual Turbinaria %N explained 39.5% of the variability in microbial communities along MDS 1 (p < 0.01). 
summary(lm(MDS1 ~ Distance_from_shore + detrended, data = directComparisons)) # Explain over 1/3 the variability
# Call:
#   lm(formula = MDS1 ~ Distance_from_shore + detrended, data = directComparisons)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -1.15383 -0.27404 -0.01568  0.28776  1.07254 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)         -3.595e-01  5.002e-02  -7.187 2.30e-11 ***
#   Distance_from_shore  8.731e-04  9.441e-05   9.249  < 2e-16 ***
#   detrended           -2.984e+00  6.131e-01  -4.868 2.66e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.3993 on 162 degrees of freedom
# Multiple R-squared:  0.3952,	Adjusted R-squared:  0.3877 
# F-statistic: 52.92 on 2 and 162 DF,  p-value: < 2.2e-16


## Notably, sites that were close to shore had microbial communities more similar to offshore sites when nutrients were low, and more disparate microbial communities when nutrients were high (Fig. 5b). 
directComp_plot = directComparisons %>%
  ggplot(aes(x = Distance_from_shore,
             y = MDS1)) +
  stat_smooth(method = "lm",
              col = "black") +
  geom_point(size = 3.5,
             pch = 16,
             col = "black") +
  geom_point(aes(col = detrended),
             size = 3,
             pch = 16) +
  stat_smooth(method = "lm",
              col = "black",
              se = F) +
  scale_color_gradientn("Residual\nN (%)",
                        colors = rev(c(RColorBrewer::brewer.pal(11,"BrBG")[1],
                                       RColorBrewer::brewer.pal(11,"BrBG")[1],
                                       RColorBrewer::brewer.pal(11,"BrBG")[1],
                                       RColorBrewer::brewer.pal(11,"BrBG")[1],
                                       RColorBrewer::brewer.pal(11,"BrBG")))) +
  xlab("Distance from shore (m)") +
  ylab("Microbial MDS1") +
  CJsBasics::BasicTheme +
  theme(legend.position = c(0.875,0.16),
        legend.text = element_text(size = 10),
        legend.key.height = unit(0.1,"in"))
directComp_plot



fig5_cd = plot_grid(tukeyDispPlot,
                    directComp_plot + theme(legend.position = c(0.875,0.175)),
                    nrow = 2,
                    rel_heights = c(1,1),
                    labels = c("c","d"))
figure5_compiled = plot_grid(pcoa_loadings_stackedBiPlot,
                             fig5_cd,
                             nrow = 1,
                             rel_widths = c(1.2,1),
                             labels = c("a",""),
                             align = "hv",
                             axis = "tb")
figure5_compiled
ggsave("plots/figure5_4panel.jpg",
       figure5_compiled,
       width = 8, height = 8, units = "in", dpi = 300)

figure5_final = plot_grid(pcoa_loadings_biplot_points_flipped,
                          directComp_plot,
                          nrow = 1,
                          labels = c("a","b"))
figure5_final
ggsave("plots/figure5.jpg",
       figure5_final,
       width = 8, height = 4, units = "in", dpi = 300)
ggsave("plots/figure5.png",
       figure5_final,
       width = 8, height = 4, units = "in", dpi = 300)
ggsave("plots/figure5.tiff",
       figure5_final,
       width = 8, height = 4, units = "in", dpi = 300)


#### Compare distance among years ----
## ). Mantel tests on Bray-Curtis distance matrices revealed consistent spatial structure of lagoon microbial communities among years (p = 0.001 in all 3 comparisons; Mantel r = 0.45, 0.52, and 0.44 for 2021 vs 2022, 2021 vs 2023, and 2022 vs 2023, respectively; Supplementary materials S11).
siteCounter = c(unique(samdat_tbl(p2021)$Site),
                unique(samdat_tbl(p2022)$Site),
                unique(samdat_tbl(p2023)$Site)) %>% 
  table() %>% 
  as.data.frame() %>% 
  as_tibble() %>%
  rename("Site" = 1) %>%
  mutate(Site = as.character(Site)) %>%
  filter(Freq == 3)

j2021 = p2021 %>% 
  ps_filter(Site %in% siteCounter$Site) %>%
  phyloseq_validate(remove_undetected = TRUE) %>%
  tax_fix(unknowns = c("endosymbionts", "uncultured",
                       "Unknown_Family", "Unknown_Family Family")) %>%
  tax_transform(trans = "identity", rank = taxLevel)  %>%
  dist_calc("bray") 
j2022 = p2022 %>% 
  ps_filter(Site %in% siteCounter$Site) %>%
  phyloseq_validate(remove_undetected = TRUE) %>%
  tax_fix(unknowns = c("endosymbionts", "uncultured",
                       "Unknown_Family", "Unknown_Family Family")) %>%
  tax_transform(trans = "identity", rank = taxLevel)  %>%
  dist_calc("bray") 
j2023 = p2023 %>% 
  ps_filter(Site %in% siteCounter$Site) %>%
  phyloseq_validate(remove_undetected = TRUE) %>%
  tax_fix(unknowns = c("endosymbionts", "uncultured",
                       "Unknown_Family", "Unknown_Family Family")) %>%
  tax_transform(trans = "identity", rank = taxLevel)  %>%
  dist_calc("bray") 

b2021 = as.matrix(j2021@dist)
b2022 = as.matrix(j2022@dist)
b2023 = as.matrix(j2023@dist)

mantel(b2021, b2022)
# Mantel statistic based on Pearson's product-moment correlation 
# 
# Call:
# mantel(xdis = b2021, ydis = b2022) 
# 
# Mantel statistic r: 0.4486 
#       Significance: 0.001 
# 
# Upper quantiles of permutations (null model):
#    90%    95%  97.5%    99% 
# 0.0781 0.1032 0.1203 0.1456 
# Permutation: free
# Number of permutations: 999

mantel(b2021, b2023)
# Mantel statistic based on Pearson's product-moment correlation 
# 
# Call:
# mantel(xdis = b2021, ydis = b2023) 
# 
# Mantel statistic r: 0.5232 
#       Significance: 0.001 
# 
# Upper quantiles of permutations (null model):
#    90%    95%  97.5%    99% 
# 0.0650 0.0856 0.1026 0.1316 
# Permutation: free
# Number of permutations: 999

mantel(b2022, b2023)
# Mantel statistic based on Pearson's product-moment correlation 
# 
# Call:
# mantel(xdis = b2022, ydis = b2023) 
# 
# Mantel statistic r: 0.4393 
#       Significance: 0.001 
# 
# Upper quantiles of permutations (null model):
#    90%    95%  97.5%    99% 
# 0.0684 0.0922 0.1113 0.1301 
# Permutation: free
# Number of permutations: 999
