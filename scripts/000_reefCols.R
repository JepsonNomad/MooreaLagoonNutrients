## > Define parameters ----
# reefColors = c("Bay" = "#999999",
#                "Fringing reef" = "#911000",
#                "Mid lagoon" =  "#d7b814",
#                "Back reef" =  "#d7b814",
#                "Reef crest" = "#3b9012", # "darkgreen", # old one: #3b9ab2
#                "Reef pass" = "#1a1a1a")

reefColors = c("Bay" = "#999999",
               "Fringing reef" = "#910000",
               # "Mid lagoon" =  "#e48729",
               # "Mid lagoon" =  "#e9a919",
               # "Mid lagoon" =  "orange2",
               "Mid lagoon" =  "#ff9f3b",
               # "Back reef" =  "#578e3d",
               "Back reef" =  "#4876FF",
               "Reef crest" = "#5e1c9a", 
               "Reef pass" = "#1a1a1a")

reefLevels = c("Bay",
               "Fringing reef",
               "Mid lagoon",
               "Back reef",
               "Reef crest",
               "Reef pass")

scale_reef_fill = function(...){
  scale_fill_manual(values = reefColors) 
}
scale_reef_col = function(...){
  scale_color_manual(values = reefColors)
}

landColors = c("Agriculture" = "#b3fc5cff",
               "Buildings" = "black",  
               "Cleared" = "#b3fc5cff",
               "Dirt" = "orange4",
               "FallowBuffer" = "#b3fc5cff",
               "Forest" = "darkgreen",
               "Monoculture" = "orange2",
               "Pasture" = "#b3fc5cff",
               "Paved" = "black",       
               "Sand" = "beige",
               "Water" = "lightblue") 

## From Joshua Stevens' cartography website
bi_scale_cols = c("1-1" = "#e8e8e8",
                  "2-1" = "#dfb0d6",
                  "3-1" = "#be64ac",
                  "1-2" = "#ace4e4",
                  "2-2" = "#a5add3",
                  "3-2" = "#8c62aa",
                  "1-3" = "#5ac8c8",
                  "2-3" = "#5698b9",
                  "3-3" = "#3b4994")


## Create a legend for your color palette. 
generateReefColsLegend = function(myLabel = "Habitat",
                                  myCols = reefColors,
                                  outpath = "plots/legends/habitatPointCols.png"){
  library(tidyverse)
  myplot = data.frame("x" = 1:length(myCols), 
                      "y" = 1:length(myCols),
                      "a" = names(myCols),
                      "z" = myCols) %>%
    mutate(afac = factor(a, levels = names(myCols))) %>%
    ggplot() +
    geom_point(aes(x = x, y = y, col = afac)) +
    scale_color_manual(myLabel, values = myCols) +
    theme_void() +
    theme(legend.key.height = unit(0.001,"in"),
          legend.key.width = unit(0.001,"in"),
          legend.text.position = "left",
          legend.key.spacing.y = unit(0.001,"in"),
          legend.spacing.x = unit(-0.1,"in"), 
          legend.text = element_text(hjust = 1),
          panel.background = element_rect(color = "transparent",fill='transparent'), 
          plot.background = element_rect(color = "transparent",fill='transparent'), 
          legend.background = element_rect(color = "transparent",
                                           fill='transparent'), #transparent legend bg
          legend.box.background = element_rect(color = "transparent",
                                               fill='transparent')
    )
  myLegend = cowplot::get_legend(myplot)
  ggsave(file = outpath,
         myLegend,
         width = 1.5, height = 1, units = "in",
         dpi = 600, 
         bg='transparent')
  return(myLegend)
}
# generateReefColsLegend(myCols = reefColors[c(1,2,4,5,6)])

# generateReefColsLegend(myLabel = "Habitat",
#                        myCols = reefColors[c("Bay", "Fringing reef",
#                                              "Mid lagoon",
#                                              "Back reef",
#                                              "Reef pass")],
#                        outpath = "plots/legends/habitatCols.png")

# generateReefColsLegend(myLabel = "",
#                        myCols = reefColors[c("Bay", "Fringing reef",
#                                              "Mid lagoon",
#                                              "Back reef",
#                                              "Reef pass")],
#                        outpath = "plots/legends/habitatCols_noname.png")
