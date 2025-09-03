## This script imputes missing Turbinaria %N values from island-wide
## Turbinaria sampling 2016-2023

## Author: Christian John
## Date: 13 May 2025


#### Load packages ----
library(tidyverse) # Data wrangling


#### Import data ----
## ... Turbinaria data ----
## Import data
turb_all = read_csv("data/generated/turb_N_compiled.csv")


#### Data wrangling ----
## ... Filter to imputable sites ----
## i.e., only include sites with sufficient representation.
## If a site was sampled fewer than 3 times, remove it.
## Identify these cases: there are 12
site_filter <- turb_all %>% 
  group_by(Site, FieldSeason) %>% 
  tally() %>% 
  ungroup() %>% 
  tidyr::pivot_wider(names_from = FieldSeason, 
                     values_from=n) %>% 
  mutate(n_samps = rowSums(select(., 
                                  contains("-")),
                           na.rm=T)) %>% 
  filter(n_samps < 3)
site_filter_sites = site_filter %>% select(Site)
write_csv(site_filter_sites, file = "data/generated/turbSitesFilter")

## Generate final filtered dataset with poorly sampled sites removed.
## Pivot wider and then longer to add a NA when a site-date is missing
mydata = turb_all %>%
  filter(!(Site %in% site_filter$Site)) %>% 
  select(Site, Percent_N, FieldSeason) %>% 
  tidyr::pivot_wider(names_from = FieldSeason, 
                     values_from = Percent_N) %>% 
  tidyr::pivot_longer(cols=2:ncol(.),
                      names_to = 'FieldSeason',
                      values_to = 'Percent_N')

## 35 sites will have 1 Turbinaria %N value imputed:
mydata %>% 
  group_by(Site, FieldSeason) %>% 
  summarize(nacount = sum(is.na(Percent_N))) %>% 
  ungroup() %>% 
  filter(nacount > 0)

## Identify these sites and save a list of them
sitesWithMissingN = mydata %>%
  filter(is.na(Percent_N)) %>%
  pull(Site) %>%
  unique()
write_csv(data.frame("Site" = sitesWithMissingN),
          file = "data/generated/turbSitesWithMissingN.csv")
length(sitesWithMissingN)


## 192 total sites for the Turbinaria imputation:
## Identify unique sites and sampling dates
allSites = unique(mydata$Site)
allDates = unique(mydata$FieldSeason)
nSites = length(allSites)

## ... Preview time series ----
timeSeries = mydata %>%
  ggplot() +
  geom_point(aes(x = FieldSeason,
                 y = Percent_N,
                 col = Site)) +
  geom_line(aes(x = FieldSeason,
                y = Percent_N,
                group = Site,
                col = Site)) +
  scale_color_manual(values = rep(CJsBasics::KellyCols,40)) +
  xlab("Sampling date") +
  ylab("Turbinaria %N") +
  theme_bw() +
  theme(legend.position = "none") 
timeSeries


#### Simple lm imputation ----
## ... Summaries ----
## First some summaries
## How many sites in the full dataset, and how many will we work with here?
nSitesAllTurbs = length(unique(pull(turb_all, Site)))  
nSitesRetained = length(unique(pull(mydata, Site)))
nSitesAllTurbs-nSitesRetained
## 12. This matches the site_filter info above.

## Ask: How many points are missing, and what is their proportional contribution to the full dataset?
nDataPoints = mydata %>%
  nrow()
nDataPointsMissing = mydata %>%
  filter(is.na(Percent_N)) %>%
  nrow()
nDataPointsMissing/nDataPoints ## Answer: Less than 5%
(nDataPoints - nDataPointsMissing)/nDataPoints

## Most missing values occurred in 2016:
mydata %>%
  filter(is.na(Percent_N)) %>%
  pull(FieldSeason) %>%
  table()

## ... Identify field season N means ----
## Generate daily mean value across sampling dates
## Use mutate to keep all rows
mydata_meanImp = mydata %>%
  group_by(FieldSeason) %>%
  mutate(FieldSeason_N_mean = mean(Percent_N, na.rm=T)) %>%
  ungroup()

## ... Define imputation model ----
## Model N as a function of FieldSeason average and Site identity
mod_pars = lm(Percent_N ~ FieldSeason_N_mean + Site - 1,
              data = mydata_meanImp)
summary(mod_pars)
## Note, R2 in this model is 0.99. This is a good predictor of Turb N.
mod_pars_coeff = mod_pars$coefficients

## ... Preview model results ----
## Plot of site-wise coefficients
coeff_data = summary(mod_pars)$coefficients %>%
  as.data.frame() %>% 
  arrange(Estimate) %>%
  mutate(Site = rownames(.)) %>%
  filter(Site != "FieldSeason_N_mean")
coeffPlot = coeff_data %>%
  mutate(Sort = order(Estimate)) %>%
  ggplot() +
  geom_hline(yintercept = 0,
             lty = 2) +
  geom_errorbar(aes(x = Sort,
                    ymin = Estimate - `Std. Error`,
                    ymax = Estimate + `Std. Error`)) +
  geom_point(aes(x = Sort,
                 y = Estimate)) +
  xlab("") +
  ylab("Coefficient estimates +/- Std. Err.") +
  coord_flip() +
  CJsBasics::BasicTheme +
  theme(axis.text.y = element_blank())
coeffPlot


## ... Correct daily mean N based on missing values ----
## Generate a corrected mean estimate based on model results 
## and knowledge of which sites generated mean estimate
N_correctedSeasonMean = fastDummies::dummy_cols(mydata_meanImp,
                                                select_columns = "Site") %>%
  # Add in site coefficients from the lm based on site name
  mutate(SiteFullName = paste0("Site",Site)) %>%
  left_join(data.frame(SiteFullName = names(mod_pars_coeff),
                       N_coeff = mod_pars_coeff)) %>% 
  # Create a mean adjuster column; if site is not measured on a given day
  # we need to re-compute the date-wide mean accounting for that site's
  # coefficient.
  # If the site is MISSING it should be used in the correction
  # corrected mean = uC
  # estimated mean = uE
  # uC = uE + (nmissing/npresent)*(mean(missingCoeffs))
  # The mutate that follows allows us to include a site coefficient
  # if the site is being imputed
  mutate(SiteCoeff = N_coeff,
         inclCorr = ifelse(is.na(Percent_N),
                           1,
                           0),
         SiteCorr = SiteCoeff*inclCorr) %>%
  # Next, for each field season:
  # calculate the mean site correction
  # and weight it based on the representation of missing sites
  # This should be a small number
  group_by(FieldSeason) %>%
  mutate(nMissing = sum(is.na(Percent_N)),
         nPresent = sum(!is.na(Percent_N)),
         FieldSeasonCorrection = mean(SiteCorr,na.rm=T)*(nMissing/nPresent)) %>%
  ungroup() %>%
  # Then apply the mean field season correction to identify 
  # the corrected field season N mean
  mutate(N_mean_corrected = FieldSeason_N_mean + FieldSeasonCorrection) %>%
  rename(FieldSeason_N_mean_raw = FieldSeason_N_mean,
         FieldSeason_N_mean_corrected = N_mean_corrected) %>%
  select(Site, FieldSeason, SiteFullName, Percent_N, 
         FieldSeason_N_mean_raw, FieldSeason_N_mean_corrected)

## Check how much the time-wise mean N is adjusted after accounting 
## for missing values (Note that it is not very much)
N_correctedSeasonMean %>%
  group_by(FieldSeason) %>%
  summarize(N_mean_corrected = mean(FieldSeason_N_mean_corrected),
            N_mean_raw = mean(FieldSeason_N_mean_raw)) %>%
  ungroup() %>%
  ggplot(aes(x = N_mean_raw,
             y = N_mean_corrected)) +
  geom_point() +
  geom_abline(col = "grey80",
              linewidth = 0.5) +
  xlab("Raw mean %N") +
  ylab("Corrected mean %N") +
  theme_bw()

## View data
N_correctedSeasonMean %>%
  group_by(FieldSeason) %>%
  summarize(N_mean_corrected = mean(FieldSeason_N_mean_corrected),
            N_mean_raw = mean(FieldSeason_N_mean_raw)) %>%
  ungroup() %>%
  ggplot() +
  geom_histogram(aes(x = (N_mean_raw - N_mean_corrected),
                     fill = FieldSeason)) +
  theme_bw()
## Note that the changes are marginal 
## (2 orders of magnitude smaller than typical N variation)

## ... Predict N using corrected field season means ----
N_predictions = N_correctedSeasonMean %>%
  ## Use the corrected field season N mean as our updated mean
  mutate(FieldSeason_N_mean = FieldSeason_N_mean_corrected)

N_predictions$N_corrected_pred = predict(
  mod_pars, 
  newdata = N_predictions)

## ... Impute missing N values with predicted values ----
Nimput = N_predictions %>%
  mutate(wasImputed = ifelse(is.na(Percent_N),
                             "YES","NO"),
         Percent_N = case_when(is.na(Percent_N) ~ N_corrected_pred,
                               TRUE ~ Percent_N)) %>%
  select(Site, FieldSeason, Percent_N, N_corrected_pred, wasImputed)

## Identify spread of predicted N vs observed N in cases we observed N
## If the imputation model is effective, these points should
## fall along a 1:1 line (they do)
imputationFid_cor = cor.test(Nimput$Percent_N, Nimput$N_corrected_pred)
imputationRangeMin = min(c(Nimput$Percent_N,
                        Nimput$N_corrected_pred))
imputationRangeMax = max(c(Nimput$Percent_N,
                           Nimput$N_corrected_pred))
imputationFidelity = Nimput %>%
  filter(wasImputed == "NO") %>%
  ggplot(aes(x = Percent_N,
             y = N_corrected_pred)) +
  geom_abline() +
  geom_point(shape = 1, alpha = 0.5) +
  stat_smooth(method = "lm", se = F, col = "black") +
  xlim(imputationRangeMin,
       imputationRangeMax) +
  ylim(imputationRangeMin,
       imputationRangeMax) +
  coord_equal() + 
  annotate("text", x = 0.5, y = 0.90, 
           label = paste0("r = ", round(imputationFid_cor$estimate, 2))) +
  annotate("text", x = 0.5, y = 0.85,
           label = paste0("p = ", round(imputationFid_cor$p.value, 5))) +
  xlab("Observed N") +
  ylab("Predicted N using sampling date mean + site ID") +
  theme_bw()
imputationFidelity
ggsave(filename = "plots/supplement/s3_imputationFidelity.jpg",
       plot = imputationFidelity,
       width = 4, height = 4)

## ... Tidy up imputed dataset ----
## Keep Sample date, Site ID, Percent_N (potentially imputed), and an indicator of whether or not the N was imputed
Nimput_share = Nimput %>%
  select(Site, FieldSeason, Percent_N, wasImputed)
head(Nimput_share)
tail(Nimput_share)


#### Export results ----
## ... Write spreadsheets ----
## Imputed nitrogen data
write_csv(Nimput_share,
          "data/generated/turb_N_imputed.csv")
