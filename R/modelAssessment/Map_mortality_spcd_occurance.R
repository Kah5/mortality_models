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
  
#########################################################################################
# Generate covariate data figures for Figure 1 ----
#########################################################################################

# plot PRISM climate data for each location----
#all.monthly <- readRDS("C:/Users/KellyHeilman/Box/01. kelly.heilman Workspace/mortality/Eastern-Mortality/mortality_models/data/PRISM_monthly_NE_periodic.RDS" )

all.seas <- readRDS("C:/Users/KellyHeilman/Box/01. kelly.heilman Workspace/mortality/Eastern-Mortality/mortality_models/data/seasonal_climate_periodic_NE_plots.RDS" )
# 
# all.seas
# join up to get the statecode information
trimmed.data <- cleaned.data %>% dplyr::select(PLOT.ID, state, remper, date) %>% distinct()
unique(trimmed.data$state)
# state.id <- 9
# state.trim <- trimmed.data %>% filter(state %in% state.id )

# just get the plots
climate.plots <- all.seas %>% filter( PLOT.ID %in% trimmed.data$PLOT.ID )#

state.summary <- climate.plots %>% group_by(lat, lon, PLOT.ID, year) %>% 
  select(-ID, -cn, -plot, -invyr, -plot_status_cd, -cycle, -elev) %>% distinct() %>%
  group_by(statecd, year) %>%
  summarise(mean_PPT = mean(yr_PPT, na.rm = TRUE), 
            mean_Tmax = mean(yr_MeanTmax, na.rm = TRUE))
state.summary$State <- ref_statecd[match(state.summary$statecd, ref_statecd$VALUE),]$MEANING

region.summary <- climate.plots %>% group_by(lat, lon, PLOT.ID, year) %>% 
  select(-ID, -cn, -plot, -invyr, -plot_status_cd, -cycle, -elev) %>% distinct() %>%
  group_by(year) %>%
  summarise(mean_PPT_all = mean(yr_PPT, na.rm = TRUE), 
            mean_Tmax_all = mean(yr_MeanTmax, na.rm = TRUE))

ggplot()+geom_line(data = state.summary, aes(x = year, y = mean_PPT, group = State), color = "lightblue")+
  geom_line(data = region.summary, aes(x = year, y = mean_PPT_all), color = "blue", linewidth = 1.1)+
  theme_minimal(base_size = 14)+ylab("Mean Precipitation (mm)")+xlab("Year")

ggplot()+geom_line(data = state.summary, aes(x = year, y = mean_Tmax, group = State, color = State), alpha = 0.5)+
  geom_line(data = region.summary, aes(x = year, y = mean_Tmax_all), color = "red", linewidth = 1.1)+
  theme_minimal(base_size = 14)+ylab("Average Annual Tmax (C)")+xlab("Year")+
  geom_label(data = state.summary %>% filter(year == last(year)), aes(label = State, 
                                                               x = year + 6, 
                                                               y = mean_Tmax, 
                                                               color = State))+
  theme(legend.position = "none")


plot.rempers <- cleaned.data %>% select(state, remper, date, PLOT.ID) %>% distinct() %>%
  group_by(state, remper, date) %>%
  summarise(n()) %>%
  mutate(T1 = date-remper, 
         T2 = date)
plotcommon.remper <- plot.rempers %>% group_by(state) %>% mutate(n.max = max(`n()`)) %>% 
  filter(`n()` == n.max) %>%
  rename("statecd" = "state")

state.summary.remper <- left_join(state.summary, plotcommon.remper) %>%
  mutate(`State Remeasurement` = ifelse(year >=T1 & year <= T2, 0.2, 0.1))
remper.vals <- state.summary.remper %>% filter(year %in% T1:T2) %>% group_by(statecd, State, T1, T2) %>% 
  summarise(PPT.val = mean(mean_PPT), 
            Tmax.val = mean(mean_Tmax))
state.scales = c(
  "Maine" = "#00441b",
  "Vermont" = "#004529",
  "New Hampshire" = "#006837",
  "New York" = "#238443",
  
 
  "Connecticut" = "#fd8d3c",
  "New Jersey" = "#fc4e2a",
  "Maryland" = "#a50026",
  
  "Ohio" = "#d73027",         
  "Pennsylvania" = "#f46d43",
  "West Virginia" = "#bd0026"
)







state.locations <- state.summary.remper  %>% filter(year == last(year)) %>% arrange(mean_Tmax)
state.locations$Yloc <- c(11.9, 12.4, 13, 13.5, 15.0, 15.9, 16.5, 17.0, 17.5, 18.1)

ggplot()+geom_line(data = state.summary.remper , aes(x = year, y = mean_Tmax, group = State, color = State, linewidth =   `State Remeasurement`))+
  scale_linewidth(range = c(0.1, 1))+
  geom_line(data = region.summary, aes(x = year, y = mean_Tmax_all), color = "black", linewidth = 1.1)+
  theme_bw()+ylab("Average Annual Tmax (C)")+xlab("Year")+
  geom_label(data = state.locations, aes(label = State, 
                                                                      x = year + 6, 
                                                                      y = Yloc, 
                                                                      color = State))+
  theme(legend.position = "none")+scale_color_manual(values = state.scales)+xlim(1900, 2038)
ggsave(height = 5, width = 8, units = "in", dpi = 300, "images/all_state_Tmax_time_series.png")


# do the same with precipitation

state.scales = c(
  "Maine" = "#00441b",
  "Vermont" = "#004529",
  "New Hampshire" = "#006837",
  "New York" = "#238443",
  
  
  "Connecticut" = "#fd8d3c",
  "New Jersey" = "#fc4e2a",
  "Maryland" = "#a50026",
  
  "Ohio" = "#d73027",         
  "Pennsylvania" = "#f46d43",
  "West Virginia" = "#bd0026"
)

state.scales.ppt <- c("Ohio"="#543005",
"Vermont"="#8c510a",
"Maryland"="#bf812d",
"Pennsylvania"="#dfc27d",
"New York"="#fec44f",
"New Jersey"="#66c2a4",
"New Hampshire"="#80cdc1",
"Connecticut" = "#35978f",
"Maine"="#01665e",
"West Virginia"="#003c30")

state.locations.ppt <- state.summary.remper  %>%  filter(year == last(year)) %>% arrange(mean_PPT)
state.locations.ppt$yloc <- c(900, 950, 1000, 1050, 1100, 1150, 1200, 1260, 1350, 1400)

ggplot()+geom_line(data = state.summary.remper, aes(x = year, y = mean_PPT, group = State, color = State), alpha = 0.9)+
  scale_linewidth(range = c(0.1, 1))+
  geom_line(data = region.summary, aes(x = year, y = mean_PPT_all), color = "black", linewidth = 1.1)+
  theme_bw(base_size = 12)+ylab("Mean Precipitation (mm)")+xlab("Year")+
xlim(1900, 2038)+
  geom_label(data = state.locations.ppt, aes(label = State, 
                                         x = year + 8, 
                                         y = yloc, 
                                         color = State))+
  scale_color_manual(values = state.scales.ppt)+
  theme(legend.position = "none")
ggsave(height = 5, width = 8, units = "in", dpi = 300, "images/all_state_PPT_time_series.png")

# get the N deposition data
boxdir <- "C:/Users/KellyHeilman/Box/01. kelly.heilman Workspace/mortality/Eastern-Mortality/mortality_models/"
N.oxidized <- read_delim(paste0(boxdir,"data/Ndep/Atmospheric_Oxidized.txt"))
N.oxidized$county <- as.numeric(N.oxidized$CountyFIPS)
N.oxidized$state <- as.numeric(N.oxidized$StateFIPS)

N.reduced <- read_delim(paste0(boxdir, "data/Ndep/Atmospheric_Reduced.txt"))
N.reduced$county <- as.numeric(N.reduced$CountyFIPS)
N.reduced$statecd <- as.numeric(N.reduced$StateFIPS)

Ntotal <- N.reduced
Ntotal[,5:92] <- N.reduced[,5:92]+N.oxidized[,5:92] # add together to get the total:
N.total.dep <- Ntotal %>% dplyr::select(statecd, county,`AREA (HA)`, y1930:y2017)
N.total.dep


Ndep.melt <- reshape2::melt(N.total.dep, id.vars = c( "county", "statecd", "AREA (HA)") )
Ndep.melt$year <- substring(Ndep.melt$variable, 2)

Ndep.totals <- Ndep.melt %>% filter(statecd %in% unique(plot.rempers$state)) %>% group_by(statecd, year) %>%
  mutate(areal.value = value/`AREA (HA)`) %>% 
  summarise(Avg.Ndep = mean(value), 
            Total.State.Ndep = sum(value)) %>%
  mutate(year = as.numeric(year))

Ndep.totals$State <- ref_statecd[match(Ndep.totals$statecd, ref_statecd$VALUE),]$MEANING

Ndep.total.remper <- left_join(Ndep.totals, plotcommon.remper) %>%
  mutate(`State Remeasurement` = ifelse(year >=T1 & year <= T2, 0.2, 0.1))


state.locations.ndep <- Ndep.total.remper %>% filter(year == T2) %>% arrange(Avg.Ndep)
state.locations.ndep$Ndep.loc <- c(9, 10, 11, 13, 14, 15, 16, 17, 19, 20)

state.scales.Ndep <- c("Maine"="#543005",
                      "Vermont"="#8c510a",
                      "New Hampshire"="#bf812d",
                      "New York" ="#b2abd2",
                      "Ohio"="#b2abd2",
                      "Maryland"="#8073ac",
                      "West Virginia"="#807dba",
                      "Connecticut" = "#6a51a3",
                      "Pennsylvania"="#542788",
                      "New Jersey"="#2d004b")



ggplot()+geom_line(data = Ndep.total.remper, aes(x = year, y = Avg.Ndep, group = State, color = State), alpha = 0.9)+
  scale_linewidth(range = c(0.1, 1))+
  #geom_line(data = region.summary, aes(x = year, y = mean_PPT_all), color = "black", linewidth = 1.1)+
  theme_bw(base_size = 12)+ylab("Reduced N deposition (kg/ha/year)")+xlab("Year")+
  xlim(1900, 2038)+
  geom_label(data = state.locations.ndep, aes(label = State,
                                             x = 1920,
                                             y = Ndep.loc,
                                             color = State))+
  scale_color_manual(values = state.scales.Ndep)+
  theme(legend.position = "none")
ggsave(height = 5, width = 8, units = "in", dpi = 300, "images/all_state_Ndep_time_series.png")

#### Plot pests and disturbances over time ----
# read in the pest disturbance datafiles:
# spongy moth --defoliation records by state
# spruce budworm data --defoliation recoreds by state
# beech bark disease: make table from: https://borealisdata.ca/dataset.xhtml?persistentId=doi:10.7939/DVN/10835
# additional beech bark information: https://www.nrs.fs.usda.gov/pubs/jrnl/2007/nrs_2007_morin_001.pdf
# Citation: By downloading these data the user agrees to cite the following references in any and all publications in which these data are used:
# Cale, Jonathan A; Morin, Randall S. "County-level detection dates for beech scale in Canada and the United States", http://dx.doi.org/10.7939/DVN/10835 V1 [Version]
# 
# Cale, Jonathan A; Garrison-Johnston; Mariann T; Teale, Stephen A; Castello, John D. Beech bark disease: Over a century of research revisited. 2017. Forest Ecology and Management. DOI: 10.1016/j.foreco.2017.03.031.
# 
# Morin RS, Liebhold AM, Tobin PC, Gottschalk KW, Luzader E. Spread of beech bark disease in the eastern United States and its relationship to regional forest composition. 2007. Canadian Journal of Forest Research. 37: 726-736. DOI: 10.1139/X06-281

# hemlock wooley adelgid infestations:
# hard to track down: using this data here: https://hemlock-woolly-adelgid-national-initiative-gmsts.hub.arcgis.com/pages/distribution-maps
# 

# spongy month defoliation records by state:
spongy <- read.csv("data/NE_spongy_moth_outbreaks.csv") %>% 
  rename( "statecd"="STATECD") %>% left_join(., plotcommon.remper) %>%
  group_by(statecd)%>%
  mutate(max.def = max(Acres.Defoliated))
state.locations.spongy <- spongy %>% filter(Acres.Defoliated == max.def ) %>% arrange(Acres.Defoliated)
state.locations.spongy$acre.loc <- c(10000,120000, 130000,  600000, 650000, 700000, 1500000, 1900000, 2400000, 30000000, 44000000)

spongy$State <- factor(spongy$State, levels = rev(c(state.locations.spongy$State)))

state.scales.spongy <- c( "Ohio" = "#fdae61", 
                          "West Virginia" = "#d6604d", 
                          "Maryland" = "#f4a582", 
                          "Vermont" = "#7fcdbb", 
                          "New Hampshire" = "#41b6c4", 
                          "New York" = "#1d91c0",
                          "Connecticut" = "#225ea8" ,
                           "Maine" = "#253494",
                          "Pennsylvania"="#081d58")


ggplot()+geom_area(data = spongy %>% filter(!is.na(State)), aes(x = year, y = Acres.Defoliated, group = State, fill = State), alpha = 0.9, position = "stack")+
  #scale_linewidth(range = c(0.1, 1))+
  #geom_line(data = region.summary, aes(x = year, y = mean_PPT_all), color = "black", linewidth = 1.1)+
  theme_bw(base_size = 12)+ylab("Acres Defoliated (Lymantria Dispar)")+xlab("Year")+
  xlim(1900, 2038)+
  geom_label_repel(data = state.locations.spongy, aes(label = State,
                                              x = year,
                                              y = Acres.Defoliated,
                                              color = State, 
                                              nudge_x = ifelse(year < 1990, year - 25, year + 20), 
                                              nudge_y = Acres.Defoliated+950000), 
                   box.padding = 1, max.overlaps = Inf, min.segment.length = 0, segment.size = 0.5, 
                  
                 
                  direction = "y",
                  
                   segment.color = "black")+
  scale_fill_manual(values = state.scales.spongy)+
  scale_color_manual(values = state.scales.spongy)+
  theme(legend.position = "none")
ggsave(height = 5, width = 8, units = "in", dpi = 300, "images/all_state_Spongy_time_series.png")

# spruce budworm defoliation records by state:
budworm <- read.csv("data/NE_spruce_budworm_outbreaks.csv") %>% 
  rename( "statecd"="STATECD") %>% left_join(., plotcommon.remper) %>%
  mutate(K.Acres.Defoliated = Spruce.Budworm.Acres.Defoliated/1000)%>%
  group_by(statecd)%>%
  mutate(max.def = max(Spruce.Budworm.Acres.Defoliated))
state.locations.budworm <- budworm %>% filter(Year == T1) %>% arrange(Spruce.Budworm.Acres.Defoliated)
#state.locations.budworm$acre.loc <- c(10000,120000, 130000,  600000, 650000, 700000, 1500000, 1900000, 2400000, 30000000, 44000000)

budworm$State <- factor(budworm$State, levels = rev(c(state.locations.budworm$State)))

state.scales.budworm <- c( "Ohio" = "#fdae61", 
                          "West Virginia" = "#d6604d", 
                          "Maryland" = "#f4a582", 
                          "Vermont" = "#7fcdbb", 
                          "New Hampshire" = "#41b6c4", 
                          "New York" = "#1d91c0",
                          "Connecticut" = "#225ea8" ,
                          "Maine" = "#253494",
                          "Pennsylvania"="#081d58")


ggplot()+geom_area(data = budworm %>% filter(!is.na(State) & ! State %in% "Total"), aes(x = Year, y = K.Acres.Defoliated, group = State, fill = State), alpha = 0.9, position = "stack")+
  #scale_linewidth(range = c(0.1, 1))+
  #geom_line(data = region.summary, aes(x = year, y = mean_PPT_all), color = "black", linewidth = 1.1)+
  theme_bw(base_size = 12)+ylab("Spruce Budworm Defoliation (1000 acres)")+xlab("Year")+
  xlim(1900, 2038)+
  geom_label_repel(data = state.locations.budworm, aes(label = State,
                                                      x = Year,
                                                      y = K.Acres.Defoliated,
                                                      color = State, 
                                                      nudge_x = ifelse(State == "Vermont", 1990, Year - 25), 
                                                      nudge_y = Spruce.Budworm.Acres.Defoliated+95),
                   box.padding = 5, max.overlaps = Inf, min.segment.length = 0, segment.size = 0.5, 
                  
                   segment.color = "black")+
  scale_fill_manual(values = state.scales.budworm)+
  scale_color_manual(values = state.scales.budworm)+
  theme(legend.position = "none")
ggsave(height = 5, width = 8, units = "in", dpi = 300, "images/all_state_budworm_time_series.png")

### Beech Bark disease spread ---
beech.bark <- read.csv("data/beech_bark_spread/Cale_Morin-BeechScaleDatesCanadaUS.csv") %>%
  filter(COUNTRY %in% "USA") %>%
  rename("State" = "PRVSTTNAME") %>% 
  filter(State %in% unique(state.summary.remper$State)) # just gets the remper states

# since this is just the time it was first observed, we want county area or something:
# Get county data for all the states
counties <- counties(state = c("DE","RI","ME", "MA","MD", "NH", "NJ", "NY", "OH", "PA", "VT", "WV"), cb = TRUE) %>% 
  data.frame()%>% 
  rename("State" = "STATE_NAME", 
         "CONAME" = "NAME") %>% select(State, STUSPS, CONAME, ALAND, AWATER)
# separate handling for CT, which changed thier county names around 2017
CT.counties <- counties(state = "CT", cb = FALSE, year = 2015) %>% 
  data.frame()%>% 
  rename( "CONAME" = "NAME") %>%
  mutate(State = "Connecticut", 
         STUSPS = "CT") %>% select(State, STUSPS, CONAME, ALAND, AWATER)

counties <- rbind(counties, CT.counties)

# join to beech bark data
beech.bark.cos <- counties %>% left_join(.,beech.bark %>% select(-REFERENCE)) %>% select(-STUSPS) 
beech.bark.cos %>% filter(is.na(ALAND)) # make sure all have matches

# now need to fill out timeseries:
beech.scale <- beech.bark.cos  %>% distinct() %>%
  expand_grid(year = 1900:2020) %>% # expand the DF to include years 
  mutate(BB.county = ifelse(SCALEYR <= year, "Beech Bark Present", "No observed Beech Bark"))%>% # categorical variable that says if there is BB present in this year
  # calculate rolling total of land where BB is present in each year
  group_by(State, year) %>%
           group_by(State, year, BB.county) %>%
           summarise(Total.Land = sum(ALAND)) %>% spread(BB.county, Total.Land, fill = 0) %>%
  mutate(Total.State = `Beech Bark Present` + `No observed Beech Bark`)%>%
  mutate(Percent.infested = (`Beech Bark Present`/Total.State)*100)

# get the labels of the first year:

state.locations.beech.scale <- beech.scale %>% left_join(., beech.bark.cos %>% group_by(State) %>% summarise(first.yr = min(SCALEYR))) %>% 
  filter(year == first.yr) %>% arrange(first.yr)

ggplot()+geom_line(data = beech.scale, aes(x = year, y = Percent.infested, group = State, color = State), linewidth = 2)+
  #scale_linewidth(range = c(0.1, 1))+
  #geom_line(data = region.summary, aes(x = year, y = mean_PPT_all), color = "black", linewidth = 1.1)+
  theme_bw(base_size = 12)+ylab("Beech Scale Presence (% land area)")+xlab("Year")+
  xlim(1900, 2020)+
  geom_label_repel(data = state.locations.beech.scale, aes(label = State,
                                                       x = year,
                                                       y = Percent.infested,
                                                       color = State,
                                                      # nudge_x = year,#ifelse(State == "Vermont", 1990, year - 25),
                                                       nudge_y = Percent.infested+10
                                                       ),
                   box.padding = 1, max.overlaps = Inf, min.segment.length = 0, segment.size = 0.5,
                   direction = "y", segment.color = "black")+
  #scale_fill_manual(values = state.scales.budworm)+
  scale_color_manual(values = state.scales.budworm)+
  theme(legend.position = "none")
ggsave(height = 5, width = 8, units = "in", dpi = 300, "images/all_state_Beech_scale_time_series.png")

### Beech Bark disease spread ---
HWA <- read.csv("data/HWA_invested_counties_2023_3_6_24.csv")  # just gets the remper states

# since this is just the time it was first observed, we want county area or something:
# Get county data for all the states
HWA.cos <- counties %>% mutate(CNTYNAME = toupper(CONAME)) %>% 
  rename("ST" = "STUSPS")%>% left_join(.,HWA %>% select(-X) %>% distinct() ) # some duplicates in the data

HWA.cos

HWA.infest <- HWA.cos  %>% distinct() %>%
  expand_grid(year = 1900:2020) %>% # expand the DF to include years 
  mutate(HWA.county = ifelse(year_txt <= year, "HWA Present", "No observed HWA"))%>% # categorical variable that says if there is BB present in this year
  # calculate rolling total of land where BB is present in each year
  group_by(State, year) %>%
  group_by(State, year, HWA.county) %>%
  summarise(Total.Land = sum(ALAND)) %>% 
  spread(HWA.county, Total.Land, fill = 0) %>%
  mutate(Total.State = `HWA Present` + `No observed HWA`)%>%
  mutate(Percent.infested = (`HWA Present`/Total.State)*100)

# get the labels of the first year:

state.locations.HWA.infest <- HWA.infest %>% 
  left_join(., HWA.cos %>% group_by(State)%>% summarise(first.yr = min(year_txt, na.rm =TRUE))) %>% 
  filter(year == first.yr) %>% arrange(first.yr)

# set color scale:
state.scales.budworm <- c("#b11e02",
                          "#ff355e", 
                          "#fd5b78",
                          "#ff6037", 
                          "#ffcc33",
                          "goldenrod",
                          "#4aa682", 
                          "darkgreen", 
                          "#16d0cb", 
                          "#50bfe6",
                          "#4067a8",
                          "#9c3ed6", 
                          "#ee34d2") 
names(state.scales.budworm) <- state.locations.HWA.infest$State

ggplot()+geom_line(data = HWA.infest, aes(x = year, y = Percent.infested, group = State, color = State), linewidth = 2)+
  #scale_linewidth(range = c(0.1, 1))+
  #geom_line(data = region.summary, aes(x = year, y = mean_PPT_all), color = "black", linewidth = 1.1)+
  theme_bw(base_size = 12)+ylab("Hemlock Wooley Adelgid Presence (% land area)")+xlab("Year")+
  xlim(1900, 2020)+
  geom_label_repel(data = state.locations.HWA.infest, aes(label = State,
                                                           x = year,
                                                           y = Percent.infested,
                                                           color = State,
                                                           nudge_x = ifelse(State == "Vermont", 1970, 
                                                                            ifelse(State %in% c("Ohio", "Maine", "West Virginia", "New Hampshire"), year+10, year - 25)),
                                                           nudge_y = ifelse(State == "New Hampshire",-25, Percent.infested)
  ),
  box.padding =2, max.overlaps = Inf, min.segment.length = 0, segment.size = 0.5,
  direction = "y", 
  segment.color = "black")+
  #scale_fill_manual(values = state.scales.budworm)+
  scale_color_manual(values = state.scales.budworm)+
  theme(legend.position = "none")
ggsave(height = 5, width = 8, units = "in", dpi = 300, "images/all_state_HWA_time_series.png")


# alternate labelling:
state.locations.HWA.infest <- state.locations.HWA.infest %>%
  mutate(State.year = paste(State, year))

ggplot()+geom_line(data = HWA.infest, aes(x = year, y = Percent.infested, group = State, color = State), linewidth = 2)+
  #scale_linewidth(range = c(0.1, 1))+
  #geom_line(data = region.summary, aes(x = year, y = mean_PPT_all), color = "black", linewidth = 1.1)+
  theme_bw(base_size = 12)+ylab("Hemlock Wooley Adelgid Presence (% land area)")+xlab("Year")+
  xlim(1900, 2020)+
  geom_label_repel(data = state.locations.HWA.infest, aes(label = State.year,
                                                          x = year,
                                                          y = Percent.infested,
                                                          color = State,
                                                          nudge_x = ifelse(State %in% c("Vermont", "Ohio"),2035, 
                                                                           ifelse(State %in% c("Maine", "New Hampshire"), year+10, year - 35)),
                                                          nudge_y = ifelse(State == "Ohio",12, Percent.infested)
  ),
  box.padding =2, max.overlaps = Inf, min.segment.length = 0, segment.size = 0.5,
  direction = "y", 
  segment.color = "black")+
  #scale_fill_manual(values = state.scales.budworm)+
  scale_color_manual(values = state.scales.budworm)+
  theme(legend.position = "none")+xlim(1900, 2030)
ggsave(height = 5, width = 8, units = "in", dpi = 300, "images/all_state_HWA_time_series_alternate.png")



# plot up the other drivers of veg change
# make up state-level summarise of variables
cleaned.data %>% group_by(SPCD) %>% 
  summarise(elevation = mean(elev), 
            damage = mean(damage), 
            ba.mean = median(ba),
            RD.mean = median(RD), 
            BAL.mean = median(BAL),
            n.dead = sum(M)) 

# read in historic damage data:
damage.table <- readRDS("data/N.DAMAGE.table.RDS")
ggplot(data = damage.table %>% filter(! damage_agent %in% "None"), aes(x = damage_agent, y = n.by.damage))+geom_point()

#####################################################################################################state######################################################################################################################
#Map out species-level model predicted probability of mortality for Figure 1 ----
#by species, by plot----
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

