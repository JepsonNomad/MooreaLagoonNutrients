## README

### Code to reproduce results in Lagoon Nutrients

This repository contains code to generate results and figures from John et al. (in review at *Limnology and Oceanography*) Terrigenous inputs link nutrient dynamics to microbial communities in a tropical lagoon. 

![Manuscript figure showing enriched N in bays, fringing reefs, and reef passes](plots/figure2/figure2.png)

### Repository structure

Please use the code in the `/scripts/` folder in order. 

#### Data wrangling 
* `000_reefCols.R`: Colors for plotting
* `001_compileRawTurbData.R`: Combining Turbinaria %N data analyzed at UGA + UCSB
* `002_compileRawWaterData.R`: Combining water column nutrients and fDOM data
* `003_imputeTurb.R`: Imputing missing Turbinaria %N values

#### Main results
* `010_figure1.R`: Produces the map and analyses of LTER **Turbinaria** %N data in Figure 1
* `011_figure2_abc_habitatN.R`: Generates the boxplots and relevant analyses in Figure 2 a, b, and c
* `012_figure2_de_map.R`: Generates the map and relevant analyses in Figure 2 d and e
* `013_figure3_waterChemistryOrdination.R`: Produces the water chemistry ordination
* `014_figure4_rainImpacts.R`: Analyses of temporal scale in nutrient enrichment
* `015_figure5_microbialCommunity.R`: Ordination and microbial community analyses

#### Supplementary materials
* `101_supplementS1.R`: Summarizes sampling effort across years
* `102_supplementS2.R`: Cross-calibrates **Turbinaria** %N across UGA and UCSB labs
* `103_supplementS5_figure2abc.R`: Reproduces Figure 2abc without imputed **Turbinaria** %N values
* `104_supplementS5_figure2de.R`: Reproduces Figure 2de without imputed **Turbinaria** %N values
* `105_supplementS9.R`: Boxplots of each water nutrient and fDOM variable by habitat
* `108_supplementS11_microbialcommunity.R`: Microbial community ordinations for individual years

