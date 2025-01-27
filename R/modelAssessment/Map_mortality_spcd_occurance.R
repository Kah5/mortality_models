### 
# Code to make species distribution maps of the eastern tree species
# ##
# Author: Kelly Heilman
# source of species distribution shapefiles: https://github.com/wpetry/USTreeAtlas
library(sf)
library(ggplot2)
library(tidyverse)
library(here)
library(FIESTA)
library(maps)
library(mapdata)

################################################################################
# Read in mortality data for 17 species
################################################################################
cleaned.data <- readRDS( "data/cleaned.data.mortality.TRplots.RDS")
unique(cleaned.data$SPCD)

# get the top species
nspp <- cleaned.data %>% group_by(SPCD) %>% summarise(n = n(), 
                                                      pct = n/nrow(cleaned.data)) %>% arrange (desc(`pct`))

nspp$cumulative.pct <- cumsum(nspp$pct)



# link up to the species table:
nspp$COMMON <- FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD),]$COMMON

nspp[1:17,]$COMMON


cleaned.data.17 <- cleaned.data %>% dplyr::select(-spp) %>% filter(SPCD %in% unique(nspp[1:17,]$SPCD))

################################################################################
# For each species get the combined species name that matches the little maps
################################################################################
plt.all <- left_join(cleaned.data.17, ref_species) 
plt.shp <- plt.all %>% mutate(GENUS = tolower(GENUS), 
                              shp.code = paste0(substr(GENUS, 1, 4), substr(SPECIES, 1,4)))

# fix some exceptions to the 4 letter rule:
# caryspp. is separated into multiple carya for the maps (maybe overlay each cayra species on the map later)
# Pinsus strobis is pinustrb instead of pinustrob (pinus strobiformis*)
# sugar maple is acersacr instead of acersacc (silver maple)


plt.shp$shp.code <- ifelse(plt.shp$SCIENTIFIC_NAME %in% "Pinus strobus", "pinustrb", plt.shp$shp.code)
plt.shp$shp.code <- ifelse(plt.shp$SCIENTIFIC_NAME %in% "Acer saccharum", "acersacr", plt.shp$shp.code)

################################################################################
# Function to create maps
################################################################################
# get the states that we need:
states <- map_data("state")
#9=CT, 25=MA, 33=NH, 23=ME, 50=VT, 44=RI, 42=PA, 39=OH, 54=WV
state_sub <- filter(states, region %in% c("connecticut","maine","new hampshire","vermont","new york", "new jersey",
                                          "rhode island","pennsylvania","ohio","west virginia", "massachusetts", "virginia", "delaware", 
                                          "north carolina", "kentucky", "tennessee", "michigan", "indiana", "district of columbia", "south carolina", 
                                          "georgia", "maryland"))

canada <- map_data("worldHires", "Canada")
spp <- "pinustrb"

# map out all the species on the plot:
ggplot() +
  geom_polygon(data = canada, 
               aes(x=long, y=lat, group = group), 
               color = "black", fill = "white") +
  geom_polygon(data = state_sub, 
               aes(x=long, y=lat, group = group), 
               color = "black", fill = "white") +
  #geom_sf(alpha = 0.75, aes(fill = as.character(Distribution)))+
  #scale_fill_manual(values = c("Species Distribution" = "forestgreen", "Outside Distribution" = "white"))+
  geom_jitter(data = plt.shp, aes(x = LONG_FIADB, y = LAT_FIADB, color = COMMON_NAME), size = 0.5)+theme_bw()+
  coord_sf(xlim = c(-85, -68), ylim = c(35, 48))+theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'lightblue'), 
                                                       #legend.position = c(0.87, 0.25),
                                                       legend.title = element_blank(),
                                                       legend.background = element_rect(fill = "white", color = "black"))+
  scale_color_manual(values = c(
    
    # maples
    "red maple" ="black", 
    "sugar maple" = "#878787", 
    # oaks
    "northern red oak" = "#a6cee3", 
    "white oak" = "#cab2d6", 
    "chestnut oak" = "#1f78b4", 
    "black oak" = "#8dd3c7",
    
    # other decidious
    "hickory spp."= "#6a3d9a",
    "American beech" = "#b15928", 
    "black cherry" = "#c51b7d", 
    "white ash" = "#f1b6da", 
    "yellow-poplar" = "#c2a5cf", 
    "yellow birch" = "#ffff99", 
    
    # conifers
    "red spruce"= "#e31a1c", 
    "balsam fir" = "#ff7f00", 
    "northern white-cedar" = "#fdbf6f", 
    "eastern white pine" = "#66c2a4", 
    "eastern hemlock" = "#00441b"
    
    
  ))

plt.pct.mortality <- plt.shp %>% group_by(state, county, pltnum, LAT_FIADB, LONG_FIADB, M) %>% 
  summarise(ntrees = n()) %>% spread(M, ntrees) %>%
  mutate(mort.ratio = ifelse(is.na(`1`),0,
                             ifelse(is.na(`0`), 1, `1`/(`0`+`1`)))) %>%
  mutate(pct.mortality = ifelse(mort.ratio > 1, 100, mort.ratio*100))

plt.pct.mortality %>% filter(pct.mortality > 100)
hist(plt.pct.mortality$pct.mortality)
ggplot(data = plt.pct.mortality, aes(x = pct.mortality))+geom_histogram()+theme_bw(base_size = 12)+xlab("% mortality on plot")
ggsave(height = 2, width = 3, units = "in", here("images/species_comp/all_observed_pct_mortality_hist.png"))

nonzero <- plt.pct.mortality %>% filter(pct.mortality> 0)
quantile(nonzero$pct.mortality)

plt.pct.mortality <- plt.pct.mortality %>% mutate(Mort.quantiles = cut(pct.mortality, 
                                                       breaks = c(0, 1, 2.5, 5, 10, 20, 30, 40, 50, 100), 
                                                       include.lowest=TRUE))

unique(plt.pct.mortality$Mort.quantiles)
large.size.font <- ggplot() +
  geom_polygon(data = canada, 
               aes(x=long, y=lat, group = group), 
               color = "black", fill = "white") +
  geom_polygon(data = state_sub, 
               aes(x=long, y=lat, group = group), 
               color = "black", fill = "white") +
  #geom_sf(alpha = 0.75, aes(fill = as.character(Distribution)))+
  #scale_fill_manual(values = c("Species Distribution" = "forestgreen", "Outside Distribution" = "white"))+
  geom_jitter(data = plt.pct.mortality, aes(x = LONG_FIADB, y = LAT_FIADB, color = Mort.quantiles), size = 5)+theme_bw()+
  coord_sf(xlim = c(-85, -68), ylim = c(35, 48))+theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'lightblue'), 
                                                       #legend.position = c(0.87, 0.25),
                                                       legend.title = element_blank(),
                                                       legend.background = element_rect(fill = "white", color = "black"))+
  scale_color_manual(values = c(
    "[0,1]" ="lightgrey",
    "(1,2.5]" = "#fcc5c0", 
    "(2.5,5]" = "#fa9fb5", 
    "(5,10]" = "#feb24c", 
    "(10,20]" = "#fd8d3c", 
    "(20,30]" = "#fc4e2a", 
    "(30,40]" = "#e31a1c", 
    "(40,50]" =  "#bd0026", 
    "(50,100]" = "#49006a" 
    
    
  ))


legend.large <- cowplot::get_legend(large.size.font)

mortality_percent_map <- ggplot() +
  geom_polygon(data = canada, 
               aes(x=long, y=lat, group = group), 
               color = "black", fill = "white") +
  geom_polygon(data = state_sub, 
               aes(x=long, y=lat, group = group), 
               color = "black", fill = "white") +
  #geom_sf(alpha = 0.75, aes(fill = as.character(Distribution)))+
  #scale_fill_manual(values = c("Species Distribution" = "forestgreen", "Outside Distribution" = "white"))+
  geom_jitter(data = plt.pct.mortality, aes(x = LONG_FIADB, y = LAT_FIADB, color = Mort.quantiles), size = 0.5)+theme_bw()+
  coord_sf(xlim = c(-84, -65), ylim = c(37, 47.5))+theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'lightblue'), 
                                                       #legend.position = c(0.87, 0.25),
                                                       legend.title = element_blank(),
                                                       legend.background = element_rect(fill = "white", color = "black"), 
                                                       legend.position = "none")+
  scale_color_manual(values = c(
    "[0,1]" ="lightgrey",
    "(1,2.5]" = "#fcc5c0", 
    "(2.5,5]" = "#fa9fb5", 
    "(5,10]" = "#feb24c", 
    "(10,20]" = "#fd8d3c", 
    "(20,30]" = "#fc4e2a", 
    "(30,40]" = "#e31a1c", 
    "(40,50]" =  "#bd0026", 
    "(50,100]" = "#49006a" 
   

  ))

png(height = 5, width = 8.5, res = 300, units = "in", here("images/species_comp/all_observed_pct_mortality_map.png") )
cowplot::plot_grid(legend.large, mortality_percent_map,  align = "v", rel_widths = c(0.15, 0.89))
dev.off()


# function to just map out the presence for the plots
plot.cored.spp.distn <- function(spp){
  dir.path = "C:/Users/KellyHeilman/OneDrive - USDA/Documents/USTreeAtlas/shp/" # where the distribution shape files are
  
  
  if(file.exists(paste0(dir.path, spp, "/"))){
    spp.distribution <- st_read(paste0(dir.path, spp, "/"))
    
    plt.spp <- plt.shp %>% filter(shp.code %in% spp)
    
    # plot distribution 
    spp.distribution %>% mutate(Distribution = ifelse(CODE == 1, "Species Distribution", "Outside Distribution"))%>%
      ggplot() +
      geom_polygon(data = canada, 
                   aes(x=long, y=lat, group = group), 
                   color = "black", fill = "white") +
      geom_polygon(data = state_sub, 
                   aes(x=long, y=lat, group = group), 
                   color = "black", fill = "white") +
      geom_sf(alpha = 0.75, aes(fill = as.character(Distribution)))+
      scale_fill_manual(values = c("Species Distribution" = "forestgreen", "Outside Distribution" = "white"))+
      geom_point(data = plt.spp, aes(x = LONG_FIADB, y = LAT_FIADB), size = 0.75)+theme_bw()+ggtitle(unique(plt.spp$COMMON_NAME))+
      coord_sf(xlim = c(-85, -68), ylim = c(35, 48))+theme(panel.grid = element_blank(), legend.position = "none", panel.background = element_rect(fill = 'lightblue'))
    
    ggsave(height = 6, width = 8, units = "in", here("images/species_comp/", paste0(spp, "_distribution_maps.png")))
  }
  else{
    cat(paste0("no range map exists for ", spp))
  }
}


# color by % live and dead for that taxa
plot.mort.spp.distn <- function(spp){
  dir.path = "C:/Users/KellyHeilman/OneDrive - USDA/Documents/USTreeAtlas/shp/" # where the distribution shape files are
  
  
  if(file.exists(paste0(dir.path, spp, "/"))){
    spp.distribution <- st_read(paste0(dir.path, spp, "/"))
    
    plt.spp <- plt.shp %>% filter(shp.code %in% spp)
    
    plt.spp <-  plt.spp %>% group_by(PLOT.ID, LONG_FIADB, LAT_FIADB, COMMON_NAME, M) %>% summarise(ntrees = n()) %>% ungroup()%>%
                group_by(PLOT.ID, LONG_FIADB, LAT_FIADB, COMMON_NAME) %>% spread(M, ntrees) %>% mutate(fraction.dead = ifelse(!is.na(`1`) & !is.na(`0`),`1`/(`0`+`1`), 
                                                                                                                    ifelse(is.na(`1`), 0, 
                                                                                                                                     ifelse(is.na(`0`), 1 ))))
    plt.spp$pct.mort <- ifelse(plt.spp$fraction.dead <= 0.1, "< 0.1",
                                      ifelse(plt.spp$fraction.dead >0.1 & plt.spp$fraction.dead <= 0.49, "0.1 - 0.49", 
                               ifelse(plt.spp$fraction.dead >0.49 & plt.spp$fraction.dead <0.9, "0.5 - 0.89", "> 0.9")))
    
    plt.spp$pct.mort <-factor(plt.spp$pct.mort, levels = c("< 0.1", "0.1 - 0.49", "0.5 - 0.89", "> 0.9"))
    # plot distribution 
   spp.map <-  spp.distribution %>% mutate(Distribution = ifelse(CODE == 1, "Species Distribution", "Outside Distribution"))%>%
      ggplot() +
      geom_polygon(data = canada, 
                   aes(x=long, y=lat, group = group), 
                   color = "black", fill = "white") +
      geom_polygon(data = state_sub, 
                   aes(x=long, y=lat, group = group), 
                   color = "black", fill = "white") +
      geom_sf(alpha = 0.75, aes(fill = Distribution))+
      scale_fill_manual(values = c("Species Distribution" = "forestgreen", "Outside Distribution" = "white"))+
      geom_point(data = plt.spp, aes(x = LONG_FIADB, y = LAT_FIADB, color = pct.mort), size = 0.5)+
      scale_color_manual(values = c("< 0.1" = "black", "0.1 - 0.49" = "#fdcc8a", "0.5 - 0.89" = "#fc8d59", "> 0.9" = "#d7301f"))+
      theme_bw()+ggtitle(unique(plt.spp$COMMON_NAME))+
      coord_sf(xlim = c(-85, -68), ylim = c(35, 48))+theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'lightblue'))
    spp.map 
    ggsave(height = 6, width = 8, units = "in", here("images/species_comp/", paste0(spp, "mortality_pct_distribution_maps.png")))
  spp.map
    }
  else{
    cat(paste0("no range map exists for ", spp))
  }
}

################################################################################
# Apply to all the species except hickory
################################################################################
# map out with presence only 
lapply(unique(plt.shp$shp.code), plot.cored.spp.distn)

# map out colored by mortality
species.maps <- lapply(unique(plt.shp$shp.code), plot.mort.spp.distn)
species.maps.legend <- cowplot::get_legend(species.maps[[1]])

png(height = 12, width = 14, units = "in", res = 150, "images/species_comp/all_spp_mortality.png")
cowplot::plot_grid(cowplot::plot_grid(species.maps[[1]] + theme(legend.position = "none"), 
                   species.maps[[2]]+ theme(legend.position = "none"), 
                   species.maps[[3]]+ theme(legend.position = "none"), 
                   species.maps[[4]]+ theme(legend.position = "none"), 
                   species.maps[[5]]+ theme(legend.position = "none"), 
                   species.maps[[7]]+ theme(legend.position = "none"), 
                   species.maps[[8]]+ theme(legend.position = "none"), 
                   species.maps[[9]]+ theme(legend.position = "none"), 
                   species.maps[[10]]+ theme(legend.position = "none"), 
                   species.maps[[11]]+ theme(legend.position = "none"), 
                   species.maps[[12]]+ theme(legend.position = "none"), 
                   species.maps[[13]]+ theme(legend.position = "none"), 
                   species.maps[[14]]+ theme(legend.position = "none"), 
                   species.maps[[15]]+ theme(legend.position = "none"), 
                   species.maps[[16]]+ theme(legend.position = "none"), 
                   species.maps[[17]]+ theme(legend.position = "none"), ncol = 4), 
                  species.maps.legend, rel_widths = c(0.85, 0.1), ncol = 2)
dev.off()
################################################################################
# Make a unique hickory map with multiple species distributions
################################################################################
# get all the hickory shapefiles:
dir.path = "C:/Users/KellyHeilman/OneDrive - USDA/Documents/USTreeAtlas/shp/" 
all.hickory <- list.files(path = dir.path, pattern = "cary")

# list of species in NE based on the TREEs of NA book
# "caryaqua" # not in range
# "carycord" # yes
# "caryflor" # not in range
# "caryglab" # yes 
# "caryilli" # borderline--maybe in ohio
# "carylaci" # yes
# "carymyri" # out of range
# "caryovat" # yes
# "carypall" # borderline -- maybe in WV
# "carytexa" # out of range
# "carytome" # yes

in.range.hickory <- all.hickory[all.hickory %in% c("carycord", "caryglab","carylaci", "caryovat", "carytome", "carypall", "caryilli")]

spp.distribution <- list()
for(i in 1:length(in.range.hickory)){
  spp.distribution[[i]] <- st_read(paste0(dir.path, in.range.hickory[i], "/"))
  spp.distribution[[i]]$`Hickory Species` <- ifelse(spp.distribution[[i]]$CODE == 1, in.range.hickory[i], "out of range")
}

# 1,3 is not in the NE
plt.spp <- plt.shp %>% filter(shp.code %in% "caryspp.")
# Carya in range

#spp.distribution[[4]] %>% 
 hickory.presence.map <-  ggplot() +
  geom_polygon(data = canada, 
               aes(x=long, y=lat, group = group), 
               color = "black", fill = "white") +
  geom_polygon(data = state_sub, 
               aes(x=long, y=lat, group = group), 
               color = "black", fill = "white") +
 
  geom_sf(data = spp.distribution[[2]], alpha = 0.5, aes(fill = `Hickory Species`))+
  geom_sf(data = spp.distribution[[3]], alpha = 0.5, aes(fill = `Hickory Species`))+
  geom_sf(data = spp.distribution[[4]], alpha = 0.5, aes(fill = `Hickory Species`))+
   geom_sf(data = spp.distribution[[5]], alpha = 0.5, aes(fill = `Hickory Species`))+
   geom_sf(data = spp.distribution[[6]], alpha = 0.5, aes(fill = `Hickory Species`))+
   geom_sf(data = spp.distribution[[7]], alpha = 0.5, aes(fill = `Hickory Species`)) +
  geom_sf(data = spp.distribution[[1]], alpha = 0.5, aes(fill = `Hickory Species`))+
  scale_fill_manual(values = c("carycord" = "forestgreen", 
                               "caryglab" = "#beaed4",
                               "caryilli" = "#fdc086",
                               "carylaci" = "#ffff99", 
                               "caryovat" = "#386cb0", 
                               "carypall" = "#f0027f",
                               "carytome" = "#bf5b17",
                               "out of range" = "white"))+
    geom_point(data = plt.spp, aes(x = LONG_FIADB, y = LAT_FIADB), size = 0.75)+
    #scale_color_manual(values = c("< 0.1" = "black", "0.1 - 0.49" = "#fdcc8a", "0.5 - 0.89" = "#fc8d59", "> 0.9" = "#d7301f"))+
    theme_bw()+ggtitle(unique(plt.spp$COMMON_NAME))+
    coord_sf(xlim = c(-85, -68), ylim = c(35, 48))+theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'lightblue'))
 hickory.presence.map
  ggsave(height = 6, width = 8, units = "in", here("images/species_comp/", paste0(spp, "_distribution_maps.png")))
  
# make the map with percent mortality
  plt.spp <-  plt.spp %>% group_by(PLOT.ID, LONG_FIADB, LAT_FIADB, COMMON_NAME, M) %>% summarise(ntrees = n()) %>% ungroup()%>%
    group_by(PLOT.ID, LONG_FIADB, LAT_FIADB, COMMON_NAME) %>% spread(M, ntrees) %>% mutate(fraction.dead = ifelse(!is.na(`1`) & !is.na(`0`),`1`/(`0`+`1`), 
                                                                                                                  ifelse(is.na(`1`), 0, 
                                                                                                                         ifelse(is.na(`0`), 1 ))))
  plt.spp$pct.mort <- ifelse(plt.spp$fraction.dead <= 0.1, "< 0.1",
                             ifelse(plt.spp$fraction.dead >0.1 & plt.spp$fraction.dead <= 0.49, "0.1 - 0.49", 
                                    ifelse(plt.spp$fraction.dead >0.49 & plt.spp$fraction.dead <0.9, "0.5 - 0.89", "> 0.9")))
  
  plt.spp$pct.mort <-factor(plt.spp$pct.mort, levels = c("< 0.1", "0.1 - 0.49", "0.5 - 0.89", "> 0.9"))
  
 all.hickory.spp.map <-  ggplot() +
    geom_polygon(data = canada, 
                 aes(x=long, y=lat, group = group), 
                 color = "black", fill = "white") +
    geom_polygon(data = state_sub, 
                 aes(x=long, y=lat, group = group), 
                 color = "black", fill = "white") +
    geom_sf(data = spp.distribution[[1]], alpha = 0.5, aes(fill = `Hickory Species`))+
    geom_sf(data = spp.distribution[[2]], alpha = 0.5, aes(fill = `Hickory Species`))+
    geom_sf(data = spp.distribution[[3]], alpha = 0.5, aes(fill = `Hickory Species`))+
    geom_sf(data = spp.distribution[[4]], alpha = 0.5, aes(fill = `Hickory Species`))+
    geom_sf(data = spp.distribution[[5]], alpha = 0.5, aes(fill = `Hickory Species`))+
    geom_sf(data = spp.distribution[[6]], alpha = 0.5, aes(fill = `Hickory Species`))+
    geom_sf(data = spp.distribution[[7]], alpha = 0.5, aes(fill = `Hickory Species`)) + 
    scale_fill_manual(values = c( 
                                 "caryglab" = "#beaed4",
                                 "caryilli" = "#fdc086",
                                 "carylaci" = "#ffff99",
                                 "carycord" = "forestgreen",
                                 "caryovat" = "#386cb0", 
                                 "carypall" = "#f0027f",
                                 "carytome" = "#bf5b17",
                                 "out of range" = "white"))+
    geom_point(data = plt.spp, aes(x = LONG_FIADB, y = LAT_FIADB, color = pct.mort), size = 0.75)+
  scale_color_manual(values = c("< 0.1" = "black", "0.1 - 0.49" = "#fdcc8a", "0.5 - 0.89" = "#fc8d59", "> 0.9" = "#d7301f"))+
    theme_bw()+ggtitle(unique(plt.spp$COMMON_NAME))+
    coord_sf(xlim = c(-85, -68), ylim = c(35, 48))+theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'lightblue'))
 all.hickory.spp.map 
  ggsave(height = 6, width = 8, units = "in", here("images/species_comp/", paste0(spp, "pct_mort_distribution_maps.png")))
  

# To do: Make one big or two big figures of species distributions with overlay
  hickory.legend <- cowplot::get_legend(hickory.presence.map)
  png(height = 15, width = 12, units = "in", res = 150, "images/species_comp/all_spp_mortality2.png")
  cowplot::plot_grid(cowplot::plot_grid(species.maps[[15]]+ theme(legend.position = "none")+ theme(axis.title = element_blank()),
                                        species.maps[[16]]+ theme(legend.position = "none")+ theme(axis.title = element_blank()),
                                        species.maps[[17]]+ theme(legend.position = "none")+ theme(axis.title = element_blank()),
                                        species.maps[[14]]+ theme(legend.position = "none")+ theme(axis.title = element_blank()),
                                        species.maps[[5]]+ theme(legend.position = "none")+ theme(axis.title = element_blank()),
                                        species.maps[[1]] + theme(legend.position = "none") + theme(axis.title  = element_blank()), 
                                        species.maps[[11]]+ theme(legend.position = "none")+ theme(axis.title = element_blank()),
                                        
                                        species.maps[[13]]+ theme(legend.position = "none")+ theme(axis.title = element_blank()),
                                        species.maps[[3]]+ theme(legend.position = "none") + theme(axis.title = element_blank()),
                                        species.maps[[7]]+ theme(legend.position = "none")+ theme(axis.title = element_blank()),
                                        species.maps[[8]]+ theme(legend.position = "none")+ theme(axis.title = element_blank()),
                                        species.maps[[10]]+ theme(legend.position = "none")+ theme(axis.title = element_blank()),
                                        
                                        species.maps[[2]]+ theme(legend.position = "none") + theme(axis.title = element_blank()),
                                         species.maps[[4]]+ theme(legend.position = "none")+ theme(axis.title = element_blank()),
                                         species.maps[[9]]+ theme(legend.position = "none")+ theme(axis.title = element_blank()),
                                        
                                        species.maps[[12]]+ theme(legend.position = "none")+ theme(axis.title = element_blank()),
                                        
                                        
                                      all.hickory.spp.map + theme(legend.position = "none")+ theme(axis.title = element_blank()),
                                      hickory.legend, ncol = 3), 
                     species.maps.legend, rel_widths = c(0.85, 0.15), ncol = 2)
  dev.off()

######################################################################################################################
# Map out species-level model predicted probability of mortality, by species, by plot
######################################################################################################################
  SPCD.df <- data.frame(SPCD = nspp[1:17, ]$SPCD, 
                        Species = nspp[1:17, ]$COMMON)
# loop to read in the predictions for pmort
  for(i in 1:17){# run for each of the 17 species
    common.name <- nspp[1:17, ] %>% filter(SPCD %in% SPCD.df[i,]$SPCD) %>% dplyr::select(COMMON)
     
    SPCD.id <- SPCD.df[i,]$SPCD
      model.number <- 6
      
      remper.correction <- 0.5
      model.name <- paste0("mort_model_",model.number,"_SPCD_", SPCD.id, "_remper_correction_", remper.correction)
      remp.cor <- 0.5
      
      SPCD.id <-  SPCD.df[i,]$SPCD
        
      cat(paste("mapping out put from stan mortality model ",model.number, " for SPCD", SPCD.df[i,]$SPCD, common.name$COMMON, " remper correction", remper.correction))
    
      fit.1 <- readRDS( paste0("SPCD_stanoutput_full/samples/model_",model.number,"_SPCD_",SPCD.id, "_remper_correction_", remper.correction, ".RDS"))
      
      # tload the data used to fit the model 
      load(paste0("SPCD_standata_general_full/SPCD_",SPCD.id, "remper_correction_", remper.correction,"model_",model.number, ".Rdata")) # load the species code data
      
      
      species.table <- unique(train.data[,c("SPCD","SPP")])
      species.table$COMMON <- FIESTA::ref_species[match(species.table$SPCD, FIESTA::ref_species$SPCD),]$COMMON_NAME
      species.table$SPP <- as.character(species.table$SPP)
      
      
      names(fit.1) <- c("alpha_SPP", colnames(mod.data$xM),
                        ## in sample predicted status
                        paste0("yrep[",1:mod.data$N, "]"),
                        paste0("psurv[",1:mod.data$N, "]"),
                        ## out of  sample predicted status
                        paste0("yhat[",1:mod.data$Nrep, "]"),
                        ## out of sample predicted prob mor
                        paste0("psurv.hat[",1:mod.data$Nrep, "]"),
                        ## in sample predicted status
                        paste0("log_lik[",1:mod.data$N, "]"),
                        "lp__")
      par.names = c("alpha_SPP", colnames(mod.data$xM)) #,
      
      nvariables <- length(par.names)
      
      species.table <- unique(train.data[,c("SPCD","SPP")])
      species.table$COMMON <- FIESTA::ref_species[match(species.table$SPCD, FIESTA::ref_species$SPCD),]$COMMON_NAME
      species.table$SPP <- as.character(species.table$SPP)
      
      fit_ssm_df <- as.data.frame(fit.1) 
     
     # get the probability of survival estimates for in-sample data  
      psurv.estimates <- fit_ssm_df %>% dplyr::select( paste0("psurv[",1:mod.data$N, "]")) 
      psurv.m <- reshape2::melt(psurv.estimates)
      
      psurv.quant <- psurv.m %>% group_by(variable) %>% summarise(median = quantile(value, 0.5, na.rm =TRUE),
                                                                  ci.lo = quantile(value, 0.005, na.rm =TRUE),
                                                                  ci.hi = quantile(value, 0.975, na.rm =TRUE))
      
      
      
      # add data from the training dataset
      psurv.quant$Mobs <- as.character(mod.data$y)
      psurv.quant$COMMON <- unique(species.table$COMMON)
      
      train.data.attributes <- train.data %>% select(date, state, county, pltnum, point, tree, PLOT.ID, spp, dbhcur, dbhold, Species, LAT_FIADB, LONG_FIADB)
      psurvival.train <- cbind(psurv.quant, train.data.attributes)
      
    
      # get the probability of survival estimates for out of sample data 
      psurv.hat.estimates <- fit_ssm_df %>% dplyr::select( paste0("psurv.hat[",1:mod.data$Nrep, "]")) 
      psurv.hat.m <- reshape2::melt(psurv.hat.estimates)
      
      psurv.hat.quant <- psurv.hat.m %>% group_by(variable) %>% summarise(median = quantile(value, 0.5, na.rm =TRUE),
                                                                  ci.lo = quantile(value, 0.005, na.rm =TRUE),
                                                                  ci.hi = quantile(value, 0.975, na.rm =TRUE))
      
      
      
      # add data from the training dataset
      psurv.hat.quant$Mobs <- as.character(mod.data$ytest)
      psurv.hat.quant$COMMON <- unique(species.table$COMMON)
      
      test.data.attributes <- test.data %>% select(date, state, county, pltnum, point, tree, PLOT.ID, spp, dbhcur, dbhold, Species, LAT_FIADB, LONG_FIADB)
      psurvival.test <- cbind(psurv.hat.quant, test.data.attributes)
      psurvival.test$dataset <- "held-out"
      psurvival.train$dataset <- "in-sample"
      psurvival.all <-  rbind(psurvival.train, psurvival.test)
      
      
      
      # calculate mortality
     pmort.summaries <-  psurvival.all %>% mutate(pmort.med = 1-median) %>% group_by(PLOT.ID, date, LAT_FIADB, LONG_FIADB, spp) %>% 
        summarise(mean.pmort = median(pmort.med), 
                  pct.mortality = median(pmort.med)*100,
                  sum.pmort = sum(pmort.med))
     # generate bins
     pmort.summaries <- pmort.summaries %>% mutate(Mort.quantiles = cut(pct.mortality, 
                                                                            breaks = c(0, 1, 2.5, 5, 10, 20, 30, 40, 50, 100), 
                                                                            include.lowest=TRUE))
    
     large.size.font <- ggplot() +
       #geom_sf(alpha = 0.75, aes(fill = as.character(Distribution)))+
       #scale_fill_manual(values = c("Species Distribution" = "forestgreen", "Outside Distribution" = "white"))+
       geom_jitter(data = pmort.summaries, aes(x = LONG_FIADB, y = LAT_FIADB, color = Mort.quantiles), size = 5)+theme_bw()+
       coord_sf(xlim = c(-85, -68), ylim = c(35, 48))+theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'lightblue'), 
                                                            #legend.position = c(0.87, 0.25),
                                                            legend.title = element_blank(),
                                                            legend.background = element_rect(fill = "white", color = "black"))+
       scale_color_manual(values = c(
         "[0,1]" ="lightgrey",
         "(1,2.5]" = "#fcc5c0", 
         "(2.5,5]" = "#fa9fb5", 
         "(5,10]" = "#feb24c", 
         "(10,20]" = "#fd8d3c", 
         "(20,30]" = "#fc4e2a", 
         "(30,40]" = "#e31a1c", 
         "(40,50]" =  "#bd0026", 
         "(50,100]" = "#49006a" 
         
         
       ))
     
     
     legend.large <- cowplot::get_legend(large.size.font)
  
     mortality_percent_map_spp <- ggplot() +
       geom_polygon(data = canada, 
                    aes(x=long, y=lat, group = group), 
                    color = "black", fill = "white") +
       geom_polygon(data = state_sub, 
                    aes(x=long, y=lat, group = group), 
                    color = "black", fill = "white") +
       #geom_sf(alpha = 0.75, aes(fill = as.character(Distribution)))+
       #scale_fill_manual(values = c("Species Distribution" = "forestgreen", "Outside Distribution" = "white"))+
       geom_jitter(data = pmort.summaries, aes(x = LONG_FIADB, y = LAT_FIADB, color = Mort.quantiles), size = 0.5)+theme_bw()+
       coord_sf(xlim = c(-84, -65), ylim = c(37, 47.5))+theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'lightblue'), 
                                                              #legend.position = c(0.87, 0.25),
                                                              legend.title = element_blank(),
                                                              legend.background = element_rect(fill = "white", color = "black"), 
                                                              legend.position = "none")+
       scale_color_manual(values = c(
         "[0,1]" ="lightgrey",
         "(1,2.5]" = "#fcc5c0", 
         "(2.5,5]" = "#fa9fb5", 
         "(5,10]" = "#feb24c", 
         "(10,20]" = "#fd8d3c", 
         "(20,30]" = "#fc4e2a", 
         "(30,40]" = "#e31a1c", 
         "(40,50]" =  "#bd0026", 
         "(50,100]" = "#49006a" 
         
         
       ))+ggtitle(paste0("Plot averaged predicted mortality probability for ", SPCD.df[i,]$Species))
     
    # png(height = 5, width = 8.5, res = 300, units = "in", here(paste0("images/predicted_mortality/average_predicted_pmort_SPCD_",SPCD.df[i,]$SPCD,"_map.png") ))
    cowplot::plot_grid(legend.large, mortality_percent_map_spp,  align = "v", rel_widths = c(0.15, 0.89))+
       theme(panel.background = element_rect(fill = "white", colour = NA),
             plot.background = element_rect(fill = "white", colour = NA))
     ggsave(height = 5, width = 8.5, units = "in", here(paste0("images/predicted_mortality/average_predicted_pmort_SPCD_",SPCD.df[i,]$SPCD,"_map.png") ))
     #dev.off()
     
     
    saveRDS(psurvival.all, paste0("SPCD_stanoutput_full/pmort_LL_model_",model.number,"_SPCD_",SPCD.id, "_remper_correction_", remper.correction, ".RDS"))
     
  }
  
### now read in all the species-level data frames and make one big plot level summary map:
  
pmort.files <- paste0("SPCD_stanoutput_full/", list.files("SPCD_stanoutput_full/", "pmort_LL"))
pmort.list <- lapply(pmort.files, readRDS) 
pmort.df <- do.call(rbind, pmort.list)


# calculate mortality
pmort.summaries <- pmort.df %>% mutate(pmort.med = 1-median) %>% group_by(PLOT.ID, date, LAT_FIADB, LONG_FIADB, spp) %>% 
  summarise(mean.pmort = median(pmort.med), 
            pct.mortality = median(pmort.med)*100,
            sum.pmort = sum(pmort.med))
# generate bins
pmort.summaries <- pmort.summaries %>% mutate(Mort.quantiles = cut(pct.mortality, 
                                                                   breaks = c(0, 1, 2.5, 5, 10, 20, 30, 40, 50, 100), 
                                                                   include.lowest=TRUE))

large.size.font <- ggplot() +
  #geom_sf(alpha = 0.75, aes(fill = as.character(Distribution)))+
  #scale_fill_manual(values = c("Species Distribution" = "forestgreen", "Outside Distribution" = "white"))+
  geom_jitter(data = pmort.summaries, aes(x = LONG_FIADB, y = LAT_FIADB, color = Mort.quantiles), size = 5)+theme_bw()+
  coord_sf(xlim = c(-85, -68), ylim = c(35, 48))+theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'lightblue'), 
                                                       #legend.position = c(0.87, 0.25),
                                                       
                                                       legend.background = element_rect(fill = "white", color = "black"))+
  scale_color_manual(values = c(
    "[0,1]" ="lightgrey",
    "(1,2.5]" = "#fcc5c0", 
    "(2.5,5]" = "#fa9fb5", 
    "(5,10]" = "#feb24c", 
    "(10,20]" = "#fd8d3c", 
    "(20,30]" = "#fc4e2a", 
    "(30,40]" = "#e31a1c", 
    "(40,50]" =  "#bd0026", 
    "(50,100]" = "#49006a" 
    
    
  ))+ guides(color=guide_legend(title="Predicted\n mortality\n probability"))


legend.large <- cowplot::get_legend(large.size.font)

mortality_percent_map_all <- ggplot() +
  geom_polygon(data = canada, 
               aes(x=long, y=lat, group = group), 
               color = "black", fill = "white") +
  geom_polygon(data = state_sub, 
               aes(x=long, y=lat, group = group), 
               color = "black", fill = "white") +
  #geom_sf(alpha = 0.75, aes(fill = as.character(Distribution)))+
  #scale_fill_manual(values = c("Species Distribution" = "forestgreen", "Outside Distribution" = "white"))+
  geom_jitter(data = pmort.summaries, aes(x = LONG_FIADB, y = LAT_FIADB, color = Mort.quantiles), size = 0.5)+theme_bw()+
  coord_sf(xlim = c(-84, -65), ylim = c(37, 47.5))+theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'lightblue'), 
                                                         #legend.position = c(0.87, 0.25),
                                                         legend.title = element_blank(),
                                                         legend.background = element_rect(fill = "white", color = "black"), 
                                                         legend.position = "none")+
  scale_color_manual(values = c(
    "[0,1]" ="lightgrey",
    "(1,2.5]" = "#fcc5c0", 
    "(2.5,5]" = "#fa9fb5", 
    "(5,10]" = "#feb24c", 
    "(10,20]" = "#fd8d3c", 
    "(20,30]" = "#fc4e2a", 
    "(30,40]" = "#e31a1c", 
    "(40,50]" =  "#bd0026", 
    "(50,100]" = "#49006a" 
    
    
  ))#+ggtitle(paste0("Plot averaged predicted mortality probability predicted from species-level models"))

# png(height = 5, width = 8.5, res = 300, units = "in", here(paste0("images/predicted_mortality/average_predicted_pmort_SPCD_",SPCD.df[i,]$SPCD,"_map.png") ))
cowplot::plot_grid(legend.large, mortality_percent_map_all,  align = "v", rel_widths = c(0.15, 0.89))+
  theme(panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA))
ggsave(height = 5, width = 8.5, units = "in", here(paste0("images/predicted_mortality/average_predicted_pmort_alltrees_map.png") ))
#dev.off()


# read in the psurvival predictions from the hierarchical models

psurvival.train <- readRDS("SPCD_stanoutput_joint/ll.train.pmort.RDS")

psurvival.test <- readRDS("SPCD_stanoutput_joint/ll.test.pmort.RDS")
psurvival.test$dataset <- "held-out"
psurvival.train$dataset <- "in-sample"
psurvival.all <-  rbind(psurvival.train, psurvival.test)



# calculate mortality
pmort.summaries.joint <-  psurvival.all %>% mutate(pmort.med = 1-median) %>% group_by(PLOT.ID, date, LAT_FIADB, LONG_FIADB, spp) %>% 
  summarise(mean.pmort = median(pmort.med), 
            pct.mortality = median(pmort.med)*100,
            sum.pmort = sum(pmort.med))
# generate bins
pmort.summaries.joint <- pmort.summaries.joint %>% mutate(Mort.quantiles = cut(pct.mortality, 
                                                                   breaks = c(0, 1, 2.5, 5, 10, 20, 30, 40, 50, 100), 
                                                                   include.lowest=TRUE))

large.size.font <- ggplot() +
  #geom_sf(alpha = 0.75, aes(fill = as.character(Distribution)))+
  #scale_fill_manual(values = c("Species Distribution" = "forestgreen", "Outside Distribution" = "white"))+
  geom_jitter(data = pmort.summaries, aes(x = LONG_FIADB, y = LAT_FIADB, color = Mort.quantiles), size = 5)+theme_bw()+
  coord_sf(xlim = c(-85, -68), ylim = c(35, 48))+theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'lightblue'), 
                                                       #legend.position = c(0.87, 0.25),
                                                       legend.title = element_blank(),
                                                       legend.background = element_rect(fill = "white", color = "black"))+
  scale_color_manual(values = c(
    "[0,1]" ="lightgrey",
    "(1,2.5]" = "#fcc5c0", 
    "(2.5,5]" = "#fa9fb5", 
    "(5,10]" = "#feb24c", 
    "(10,20]" = "#fd8d3c", 
    "(20,30]" = "#fc4e2a", 
    "(30,40]" = "#e31a1c", 
    "(40,50]" =  "#bd0026", 
    "(50,100]" = "#49006a" 
    
    
  ))


legend.large <- cowplot::get_legend(large.size.font)

mortality_percent_map_spp_joint <- ggplot() +
  geom_polygon(data = canada, 
               aes(x=long, y=lat, group = group), 
               color = "black", fill = "white") +
  geom_polygon(data = state_sub, 
               aes(x=long, y=lat, group = group), 
               color = "black", fill = "white") +
  #geom_sf(alpha = 0.75, aes(fill = as.character(Distribution)))+
  #scale_fill_manual(values = c("Species Distribution" = "forestgreen", "Outside Distribution" = "white"))+
  geom_jitter(data = pmort.summaries.joint, aes(x = LONG_FIADB, y = LAT_FIADB, color = Mort.quantiles), size = 0.5)+theme_bw()+
  coord_sf(xlim = c(-84, -65), ylim = c(37, 47.5))+theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'lightblue'), 
                                                         #legend.position = c(0.87, 0.25),
                                                         legend.title = element_blank(),
                                                         legend.background = element_rect(fill = "white", color = "black"), 
                                                         legend.position = "none")+
  scale_color_manual(values = c(
    "[0,1]" ="lightgrey",
    "(1,2.5]" = "#fcc5c0", 
    "(2.5,5]" = "#fa9fb5", 
    "(5,10]" = "#feb24c", 
    "(10,20]" = "#fd8d3c", 
    "(20,30]" = "#fc4e2a", 
    "(30,40]" = "#e31a1c", 
    "(40,50]" =  "#bd0026", 
    "(50,100]" = "#49006a" 
    
    
  ))+ggtitle(paste0("Plot averaged predicted mortality probability - hierarchical model"))

# png(height = 5, width = 8.5, res = 300, units = "in", here(paste0("images/predicted_mortality/average_predicted_pmort_SPCD_",SPCD.df[i,]$SPCD,"_map.png") ))
cowplot::plot_grid(legend.large, mortality_percent_map_spp_joint,  align = "v", rel_widths = c(0.15, 0.89))+
  theme(panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA))
ggsave(height = 5, width = 8.5, units = "in", here(paste0("images/predicted_mortality/average_predicted_pmort_joint_SPCD_map.png") ))
#dev.off()

cowplot::plot_grid(legend.large, cowplot::plot_grid(mortality_percent_map_all, mortality_percent_map_spp_joint + ggtitle(""), align = 
                                             "hv", 
                                           nrow = 2, 
                                           labels = "AUTO"),
                     align = "v", rel_widths = c(0.15, 0.89))+
  theme(panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA))
ggsave(height = 10, width = 8.5, units = "in", here(paste0("images/predicted_mortality/average_predicted_pmort_alltrees_map.png") ))

#-------------------------------------------------------------------------------------
# generate predictions for additional species

cleaned.data <- cleaned.data %>% filter(!is.na(ba) & !is.na(slope) & ! is.na(physio) & !is.na(aspect))%>% dplyr::select(state, county, pltnum, cndtn, point, tree, PLOT.ID, cycle, spp, dbhcur, dbhold, damage, Species, SPCD,
                                                                                                                        remper, LAT_FIADB, LONG_FIADB, elev, DIA_DIFF, annual.growth, M, relative.growth, si, physio:RD) %>% distinct()
plot.medians <- unique(cleaned.data %>% ungroup()%>% dplyr::select(PLOT.ID, si, ba, slope, aspect, MAP, MATmin, MATmax, damage.total, elev, Ndep.remper.avg, physio, RD)) %>% ungroup() %>% summarise(si.median = median(si, na.rm =TRUE), 
                                                                                                                                                                                                      RD.median = median(RD, na.rm = TRUE),
                                                                                                                                                                                                      ba.median = median(ba, na.rm =TRUE), 
                                                                                                                                                                                                      slope.median = median(slope, na.rm = TRUE), 
                                                                                                                                                                                                      aspect.median = median(aspect, na.rm = TRUE),
                                                                                                                                                                                                      damage.median = median(damage.total, na.rm =TRUE),
                                                                                                                                                                                                      elev.median = median(elev, na.rm =TRUE),
                                                                                                                                                                                                      Ndep.median = median(Ndep.remper.avg, na.rm =TRUE),
                                                                                                                                                                                                      physio.median = median(physio, na.rm = TRUE),
                                                                                                                                                                                                      
                                                                                                                                                                                                      MAP.median = median(MAP, na.rm =TRUE), 
                                                                                                                                                                                                      MATmin.median = median(MATmin, na.rm =TRUE), 
                                                                                                                                                                                                      MATmax.median = median(MATmax, na.rm =TRUE), 
                                                                                                                                                                                                      
                                                                                                                                                                                                      RD.sd = sd(RD, na.rm = TRUE),
                                                                                                                                                                                                      ba.sd = sd(ba, na.rm =TRUE),
                                                                                                                                                                                                      si.sd = sd(si, na.rm =TRUE), 
                                                                                                                                                                                                      slope.sd = sd(slope, na.rm =TRUE),
                                                                                                                                                                                                      aspect.sd = sd(aspect, na.rm = TRUE),
                                                                                                                                                                                                      damage.sd = sd(damage.total, na.rm =TRUE),
                                                                                                                                                                                                      elev.sd = sd(elev, na.rm =TRUE),
                                                                                                                                                                                                      Ndep.sd = sd(Ndep.remper.avg, na.rm =TRUE),
                                                                                                                                                                                                      physio.sd = sd(physio, na.rm = TRUE),
                                                                                                                                                                                                      
                                                                                                                                                                                                      MAP.sd = sd(MAP, na.rm =TRUE), 
                                                                                                                                                                                                      MATmin.sd = sd(MATmin, na.rm =TRUE), 
                                                                                                                                                                                                      MATmax.sd = sd(MATmax, na.rm =TRUE)
)

#View(cleaned.data %>% group_by(SPGRPCD, SPCD) %>% summarise(n()))

cleaned.data.full <- cleaned.data

# get all the observations from non-focal trees
SPCD.id <- nspp[18:126,]$SPCD

cleaned.data <- cleaned.data.full %>% filter(SPCD %in% SPCD.id)
remper.correction <- 0.5

# scale the cleaned data tree-level diameters by species?
cleaned.data <- cleaned.data %>% ungroup()  %>% group_by(SPCD) %>% 
  mutate(rempercur = ifelse(M ==1, remper*remper.correction, remper), 
         annual.growth = DIA_DIFF/rempercur) %>% mutate(DIA.median = median(dbhcur, na.rm =TRUE), 
                                                        DIA.sd = sd(dbhcur, na.rm = TRUE),  
                                                        BAL.median = median(BAL, na.rm=TRUE),
                                                        BAL.sd = sd(BAL, na.rm = TRUE),
                                                        RD.median = median(RD, na.rm=TRUE), 
                                                        RD.sd = sd(RD, na.rm =TRUE),
                                                        annual.growth.median = median(annual.growth, na.rm = TRUE), 
                                                        annual.growth.sd = sd(annual.growth, na.rm = TRUE)) %>% 
  ungroup() %>% mutate(DIA_scaled = (dbhcur - DIA.median)/DIA.sd, 
                       annual.growth.scaled = (annual.growth - annual.growth.median)/annual.growth.sd, 
                       RD.scaled = (RD-RD.median)/RD.sd,
                       BAL.scaled = (BAL-BAL.median)/BAL.sd,
                       si.scaled = (si - plot.medians$si.median)/plot.medians$si.sd, 
                       ba.scaled = (ba - plot.medians$ba.median)/plot.medians$ba.sd,
                       aspect.scaled = (aspect - plot.medians$aspect.median)/plot.medians$aspect.sd,
                       slope.scaled = (slope - plot.medians$slope.median)/plot.medians$slope.sd,
                       damage.scaled = (damage.total - plot.medians$damage.median)/plot.medians$damage.sd,
                       MAP.scaled = (MAP-plot.medians$MAP.median)/plot.medians$MAP.sd, 
                       elev.scaled = (elev-plot.medians$elev.median)/plot.medians$elev.sd,
                       Ndep.scaled = (Ndep.remper.avg- plot.medians$Ndep.median)/plot.medians$Ndep.sd,
                       physio.scaled = (physio-plot.medians$physio.median)/plot.medians$physio.sd,
                       MATmin.scaled = (MATmin- plot.medians$MATmin.median)/plot.medians$MATmin.sd,
                       MATmax.scaled = (MATmax - plot.medians$MATmax.median)/plot.medians$MATmax.sd)

SPP.df <- data.frame(SPCD = unique(cleaned.data$SPCD), 
                     SPP = 1:length(unique(cleaned.data$SPCD)))

cleaned.data <- left_join(cleaned.data, SPP.df) 
cleaned.data <- cleaned.data %>%  filter(!is.na(si) & !is.na(dbhcur)& !is.na(M) & !is.na(annual.growth.scaled) & !is.na(ppt.anom))
#summary(cleaned.data$BAL.scaled)
cleaned.data$S <- ifelse(cleaned.data$M == 1, 0, 1)


train.data <- cleaned.data

# model.6 data
# 6. All Fixed effects and all growth + Diameter interactions
mod.data.6 <- list(N = nrow(train.data), 
                   y = train.data$S, 
                   xM = as.matrix(train.data %>% dplyr::select(annual.growth.scaled, 
                                                               DIA_scaled, 
                                                               RD.scaled, 
                                                               ba.scaled, 
                                                               BAL.scaled, 
                                                               damage.scaled,
                                                               MATmax.scaled, 
                                                               MATmin.scaled, 
                                                               MAP.scaled,
                                                               ppt.anom, 
                                                               tmin.anom, 
                                                               tmax.anom, 
                                                               slope.scaled, 
                                                               aspect.scaled,
                                                               elev.scaled, 
                                                               Ndep.scaled,
                                                               physio.scaled) %>%
                                    # generate growth interactions
                                    mutate_at(.funs = list(growth.int = ~.*annual.growth.scaled), 
                                              .vars = vars(DIA_scaled:physio.scaled)) %>% 
                                    mutate_at(.funs = list(DIA.int = ~.*DIA_scaled), 
                                              .vars = vars(RD.scaled:physio.scaled)) ))


# generate predictions 
# read in the alpha population 
mualpha <- readRDS("SPCD_stanoutput_joint/alpha.p_model_6_1000samples.rds")
# read in the betas populations
mubetas <- readRDS("SPCD_stanoutput_joint/beta_model_6_1000samples.rds")

posterior.predict <- function(x){
 exp.x <- exp(-1*as.data.frame( mualpha) %>% dplyr::select(alpha) + rowSums(as.data.frame(mubetas) %>% dplyr::select(paste0("mu_beta[", 1:48,"]"))*mod.data.6$xM[x,]))
  p.mort <- exp.x / (1+ exp.x)
  p.mort$alpha
  }

pred.list <- lapply(1:mod.data.6$N, FUN = posterior.predict)
pmort.hat.oos <- do.call(cbind, pred.list)
pmort.hat.oos.m <- reshape2::melt(pmort.hat.oos)
colnames(pmort.hat.oos.m) <- c("sample", "treeid", "pmort")
pmort.hat.oos.summary <- pmort.hat.oos.m %>% group_by(treeid) %>% summarise(median = quantile(pmort, 0.5, na.rm =TRUE), 
                                                        ci.lo = quantile(pmort, 0.025, na.rm =TRUE), 
                                                        ci.hi = quantile(pmort, 0.975, na.rm =TRUE))

pmort.hat.oos.summary$S <- train.data$S
pmort.hat.oos.summary$SPCD <- train.data$SPCD
pmort.hat.oos.summary$COMMON <-  FIESTA::ref_species[match(pmort.hat.oos.summary$SPCD, FIESTA::ref_species$SPCD),]$COMMON
pmort.hat.oos.summary$GENUS <-  FIESTA::ref_species[match(pmort.hat.oos.summary$SPCD, FIESTA::ref_species$SPCD),]$GENUS
pmort.hat.oos.summary$SPECIES <-  FIESTA::ref_species[match(pmort.hat.oos.summary$SPCD, FIESTA::ref_species$SPCD),]$SPECIES


# AUC using mltools auc_roc function
# for in sample data
actuals = train.data$S
preds = as.vector(pmort.hat.oos.summary$median)
auc.oos.nonfocal.spp <-auc_roc(preds, actuals)


auc.oos.spp <- pmort.hat.oos.summary %>% group_by(SPCD, COMMON) %>% summarise(`AUC` = auc_roc(median, S), 
                                                                              `# observations` = n())

hist(auc.oos.spp$AUC)

auc.oos.spp %>% filter(!is.na(AUC)) %>% filter(AUC > 0.5)
auc.oos.spp %>% filter(!is.na(AUC)) %>% filter(AUC > 0.5) %>% ungroup()%>% summarise(total.trees = sum(`# observations`))
auc.oos.spp %>% filter(!is.na(AUC)) %>% filter(AUC < 0.5)
auc.oos.spp %>% filter(!is.na(AUC)) %>% filter(AUC < 0.5) %>% ungroup()%>% summarise(total.trees = sum(`# observations`))


auc.oos.spp$GENUS <-  FIESTA::ref_species[match(auc.oos.spp$SPCD, FIESTA::ref_species$SPCD),]$GENUS
auc.oos.spp$SPECIES <-  FIESTA::ref_species[match(auc.oos.spp$SPCD, FIESTA::ref_species$SPCD),]$SPECIES

nspp$GENUS <- FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD),]$GENUS
nSPP.GENUS <- nspp[1:17,]$GENUS

auc.oos.spp %>% filter(!is.na(AUC)) %>% arrange(desc(AUC)) %>% ungroup()|> gt() %>%
  data_color(
    columns = GENUS,
    rows =  GENUS %in% nSPP.GENUS,
    method = "factor",
    #palette = c("red", "green"),
    #domain = c(0, 50)
  )%>%
  gtsave("tab_1.html", path = getwd())

auc.oos.spp %>% filter(!is.na(AUC)) %>% arrange(desc(AUC)) %>% group_by(AUC < 0.5, GENUS %in% nSPP.GENUS) %>%
  summarise(`total trees` = sum(`# observations`)) %>% ungroup()|> gt() %>%
  # data_color(
  #   columns = GENUS,
  #   rows =  GENUS %in% nSPP.GENUS,
  #   method = "factor",
  #   #palette = c("red", "green"),
  #   #domain = c(0, 50)
  # )%>%
  gtsave("summary_oos_species.html", path = getwd())

