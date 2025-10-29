library(tidyverse)
library(ggplot2)
library(cowplot)
library(posterior)
library(FIESTA)
library(sf)
output.dir <- output.folder <- "C:/Users/KellyHeilman/Box/01. kelly.heilman Workspace/mortality/Eastern-Mortality/mortality_models/"

# # get the complete species list
nspp <- data.frame(SPCD = c(316, 318, 833, 832, 261, 531, 802, 129, 762,  12, 541,  97, 621, 400, 371, 241, 375))
nspp$Species <- paste(FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD),]$GENUS, FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD),]$SPECIES)
# link up to the species table:
nspp$COMMON <- FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD),]$COMMON


# read in the test data for all the species
spp.table <- data.frame(SPCD.id = nspp[1:17,]$SPCD, 
                        spp = 1:17, 
                        COMMON = nspp[1:17,]$COMMON)
# set the species order using the factors:
SP.TRAITS <- read.csv("data/NinemetsSpeciesTraits.csv") %>% filter(COMMON_NAME %in% unique(nspp[1:17,]$COMMON))
# order the trait db by softwood-hardwood, then shade tolerance, then name (this puts all the oaks together b/c hickory and red oak have the same tolerance values)
SP.TRAITS <- SP.TRAITS %>% group_by(SFTWD_HRDWD) %>% arrange(desc(SFTWD_HRDWD), desc(ShadeTol), desc(COMMON_NAME))

SP.TRAITS$Color <- c(# softwoods
  "#b2df8a",
  "#003c30", 
  "#b2182b", 
  "#fee090", 
  "#33a02c",
  
  
  # sugar  maples
  "#a6cee3",
  "#1f78b4",
  
  # red maple
  "#e31a1c",
  # yellow birch
  "#fdbf6f",
  # oaks
  "#cab2d6",
  "#8073ac",
  "#6a3d9a",
  
  # hickory
  "#7f3b08",
  # white ash
  "#bababa",
  # black cherry
  "#4d4d4d",
  # yellow poplar
  "#ff7f00",
  "#fccde5" # paper birch
  
  
)

SP.TRAITS$`Shade Tolerance`  <- ifelse(SP.TRAITS$ShadeTol >=4, "High", 
                                       ifelse(SP.TRAITS$ShadeTol <=2.5, "Low", "Moderate"))

# set up custom colors for species
sppColors <- SP.TRAITS$Color 
names(sppColors) <- unique(SP.TRAITS$COMMON_NAME) 

species_fill <- scale_fill_manual(values = sppColors)
species_color <- scale_color_manual(values = sppColors)

# common species ordering scheme:
disturb.species.order <- c(
  # spruce fir forests
  "balsam fir", 
  "red spruce", 
  "northern white-cedar", 
  # HWA & beech
  "eastern hemlock", 
  "American beech", 
  
  # spongy moth susceptible
 # "black oak", 
  "chestnut oak", 
  "northern red oak",
  "white oak", 
  "yellow birch", 
  "paper birch",

  
  # spongy moth resistant
  "hickory spp.", 
  "eastern white pine", 
  "red maple", 
  "sugar maple", 
  
  # spongy moth immune
  "black cherry", 
  "white ash", 
  "yellow-poplar"
)

# summaries the mortality rates (observed) by species:-----

# unfiltered tree remeasurement data
TREE.remeas <- readRDS( "data/unfiltered_TREE.remeas.rds")

# BA.issues <- TREE.remeas %>% 
# mutate(ba_sq_ft_cur = ((dbhcur^2))*0.005454, 
#        ba_sq_ft_old = ((dbhold^2))*0.005454) %>% 
#   group_by(PLOT.ID, point, date, cndtn, SPCD) %>%
#   
#   mutate(SPCD_BA = sum(ba_sq_ft_cur*volfac, na.rm =TRUE), 
#          SPCD_BA_old = sum(ba_sq_ft_old*volfac, na.rm = TRUE)) %>% 
#   ungroup() %>%
#   mutate(M = ifelse(status %in% c(1,3), 0, 1)) %>% 
#   group_by(PLOT.ID, point, date, cndtn) %>%
#  summarise(BA_total = sum(SPCD_BA, na.rm =TRUE), 
#            BA_total_wrong = sum(SPCD, na.rm =TRUE), 
#          BA_total_old = sum(SPCD_BA_old, na.rm = TRUE), 
#          total_mort = sum(M), 
#          ba_tre_sum = sum(ba_sq_ft_cur, na.rm =TRUE), 
#          ba_tre_sum_old = sum(ba_sq_ft_old, na.rm =TRUE))
# 
# ggplot(data = BA.issues)+geom_point(aes(x = ba_tre_sum, y = total_mort))
# ggplot(data = BA.issues)+geom_point(aes(x = ba_tre_sum_old, y = total_mort))
# ggplot(data = BA.issues)+geom_point(aes(x = ba_tre_sum, y = BA_total))
# ggplot(data = BA.issues)+geom_point(aes(x = ba_tre_sum_old, y = ba_tre_sum))
# 
# ggplot(data = BA.issues)+geom_point(aes(x = BA_total_wrong, y = total_mort))
# ggplot(data = BA.issues)+geom_point(aes(x = BA_total_wrong, y = BA_total))
# ggplot(data = BA.issues)+geom_point(aes(x = BA_total, y = total_mort))
# ggplot(data = BA.issues)+geom_point(aes(x = BA_total_wrong, y = total_mort))
#   


all.remeas <- TREE.remeas %>%
  # do the filtering section
  filter( exprem > 0 & # if exprem == 0, these could be modeled plots?
            dbhold >= 5 & # need an initial dbh greater than 5
            ! remper == 0 & # if remper is listed as zero, filter out
            DIA_DIFF >= 0 & # filter diameter differences >= 0
            # !status == 3 & # keep the cut trees for this
            SPCD %in% nspp[1:17,]$SPCD & # filter species in the top 17 of all species
            !is.na(status) & # filter out trees with no status
            !is.na(elev)) #&# filter out plots with no FIADB lat long for elevation
#!is.na(relative.growth)) # remove any trees with NA for growth

# all.remeas %>% left_join(., BA.issues) %>% 
#   ggplot()+geom_boxplot(aes(x = status, y = BA_total, group = status))+facet_wrap(~SPCD)
# 
# all.remeas %>% left_join(., BA.issues) %>% 
#   ggplot()+geom_boxplot(aes(x = status, y = BA_total_old, group = status))+facet_wrap(~SPCD, scales = "free_y")
# 
# all.remeas %>% left_join(., BA.issues) %>% 
#   ggplot()+geom_boxplot(aes(x = status, y = BA_total_wrong, group = status))+facet_wrap(~SPCD, scales = "free_y")
# 
# all.remeas %>% left_join(., BA.issues) %>% 
#   ggplot()+geom_boxplot(aes(x = status, y = ba_tre_sum_old, group = status))+facet_wrap(~SPCD, scales = "free_y")
# 
# all.remeas %>% left_join(., BA.issues) %>% 
#   ggplot()+geom_boxplot(aes(x = status, y = ba_tre_sum, group = status))+facet_wrap(~SPCD, scales = "free_y")
# 

# mortality rates across the whole species distributions:----
total.tree.sums <- all.remeas %>% group_by(SPCD, Species, remper)%>%
  summarise(total_trees_volfac = sum(volfac, na.rm = TRUE), 
            total_n = n())


# logging only
all.remeas %>% group_by(SPCD, Species, remper)%>% filter(status == 3) %>%
  summarise(logged_trees_volfac = sum(volfac, na.rm = TRUE),
            n_cut = n())%>%
  left_join(., total.tree.sums) %>%# join to total trees
  
  # for each remper period, get the mortality rate per year 
  mutate(all_cut_rate_volfac = ((logged_trees_volfac/total_trees_volfac)*100)/remper, 
         all_cut_rate = ((n_cut/total_n)*100)/remper)%>%
  ungroup()%>%
  
  # average all the remper mortality rates together to get a single species value
  group_by(SPCD, Species)%>%
  summarise(cut_volfac_mort = mean(all_cut_rate_volfac, na.rm =TRUE), 
            cut_n_mort = mean(all_cut_rate, na.rm =TRUE))

# mortality rates by species: non logging dead trees
mort.rate.species <- all.remeas %>% group_by(SPCD, Species, remper)%>% filter(status %in% c(2, 4, 5, 6)) %>%
  summarise(nonlog_dead_trees_volfac = sum(volfac, na.rm = TRUE), 
            n_nonlog_dead = n()) %>%
  left_join(., total.tree.sums) %>%# join to total trees
  
  # for each remper period, get the mortality rate per year 
  mutate(nonlog_mort_rate_volfac = ((nonlog_dead_trees_volfac/total_trees_volfac)*100)/remper, 
         nonlog_mort_rate = ((n_nonlog_dead/total_n)*100)/remper)%>%
  ungroup()%>%
  
  # average all the remper mortality rates together to get a single species value
  group_by(SPCD, Species)%>%
  summarise(species_volfac_mort = mean(nonlog_mort_rate_volfac, na.rm =TRUE), 
            species_n_mort = mean(nonlog_mort_rate, na.rm =TRUE))

mort.rate.species$Species <- factor(mort.rate.species$Species, levels = disturb.species.order)

ggplot(data = mort.rate.species)+
  geom_bar(aes(x = Species, y = species_volfac_mort, fill = Species), stat = "identity")+
  theme_bw(base_size = 16)+
  species_fill+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        #panel.grid = element_blank(), 
        legend.position = "none", 
        axis.title.x = element_blank())+
  ylab("Average Species Mortality Rate (% per year)")



# mortality rates by species and states distributions:----
st.total.tree.sums <- all.remeas %>% group_by(SPCD, Species,state, remper)%>%
  summarise(total_trees_volfac = sum(volfac, na.rm = TRUE), 
            total_n = n())

# get the number of plots in each state:
st.total.plt.spcd <- all.remeas %>% 
  select(SPCD, Species, state,pltnum, PLOT.ID) %>% distinct()%>%
  group_by(SPCD, Species, state)%>%
  summarise(n_plots_SPCD = n())


# mortality rates by species: non logging dead trees
st.mort.rate.species <- all.remeas %>% group_by(SPCD, Species, state, remper)%>% filter(status %in% c(2, 4, 5, 6)) %>%
  summarise(nonlog_dead_trees_volfac = sum(volfac, na.rm = TRUE), 
            n_nonlog_dead = n()) %>%
  left_join(., st.total.tree.sums) %>%# join to total trees
  
  # for each remper period, get the mortality rate per year 
  mutate(nonlog_mort_rate_volfac = ((nonlog_dead_trees_volfac/total_trees_volfac)*100)/remper, 
         nonlog_mort_rate = ((n_nonlog_dead/total_n)*100)/remper)%>%
  ungroup()%>%
  
  # average all the remper mortality rates together to get a single species value
  group_by(SPCD, Species, state)%>%
  summarise(species_volfac_mort = mean(nonlog_mort_rate_volfac, na.rm =TRUE), 
            species_n_mort = mean(nonlog_mort_rate, na.rm =TRUE))


state.df <- ref_statecd %>% 
  rename("state" = "VALUE", 
         "region" = "MEANING") %>%
  mutate(region = tolower(region))

st.mort.df <- st.mort.rate.species %>% left_join(., state.df) %>%
  left_join(., st.total.plt.spcd)
st.mort.df$Species <- factor(st.mort.df$Species, levels = disturb.species.order)

summary(st.mort.df$species_volfac_mort)
quantile(st.mort.df$species_volfac_mort, c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99))
hist(st.mort.df$species_volfac_mort)

st.mort.df <- st.mort.df %>% mutate(
  cut.mort.rate = cut(species_volfac_mort, breaks = c(0, 0.25, 0.5, 0.75, 1, 1.5, 2, 4, 8))
)
cut.mort.values <- data.frame(
  cut.mort.rate = c("(0,0.25]", "(0.25,0.5]", "(0.5,0.75]", "(0.75,1]", "(1,1.5]", "(1.5,2]", "(2,4]", "(4,8]"), 
  mortality.rate = c("<0.25", "0.25 - 0.5", "0.5 - 0.75", "0.75 - 1", "1 - 1.5", "1.5 - 2", "2 - 4", "4-6"), 
  hex.colors = c("#fff7f3",
                 "#fde0dd",
                 "#fcc5c0",
                 "#fa9fb5",
                 "#f768a1",
                 "#dd3497",
                 "#ae017e",
                 "#7a0177"))

st.mort.df <- st.mort.df %>% left_join(., cut.mort.values)
st.mort.df$mortality.rate <- factor(st.mort.df$mortality.rate, levels = c("<0.25", "0.25 - 0.5", "0.5 - 0.75", "0.75 - 1", "1 - 1.5", "1.5 - 2", "2 - 4", "4-6"))

val.vec <- as.vector(cut.mort.values$hex.colors)
names(val.vec) <- as.vector(cut.mort.values$mortality.rate)
fill_mort_rate <- scale_fill_manual(values = val.vec, name = "Mortality\nRate\n(%/year)", drop = FALSE)

st.mort.df <- st.mort.df %>% left_join(.,
                                       data.frame(region = unique(st.mort.df$region),
                                                  Abbrev = c("ME", "NH", "NY", "VT", "WV", "CT", "MD", "OH", "PA", "NJ")))


mort.rate.state.species <- ggplot(data = st.mort.df %>% filter(n_plots_SPCD > 50))+
  geom_bar(aes(x = Species, y = species_volfac_mort, fill = Species), stat = "summary", fun = median,  color = "black")+
  theme_bw(base_size = 16)+
  species_fill+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        axis.title.x = element_blank(), 
        legend.position = "none", 
        panel.grid = element_blank())+
  geom_point(aes(x = Species, y = species_volfac_mort, group = region), alpha = 0.9) +
  #geom_text(aes(x = Species, y = species_volfac_mort, label = region),position = position_jitter(width = 0.1), alpha = 0.6)
  ggrepel::geom_text_repel(aes(x = Species, y = species_volfac_mort, label = Abbrev), color = "black", size = 2.5, segment.color = "grey", max.overlaps = 12) +  
  ylab("Mortality Rate (% per year)")


ggsave(paste0(output.dir, "images/barplot_mortality_rate_state.png"), plot = mort.rate.state.species, 
       width = 10, height = 6, units = "in")

ggsave(paste0(output.dir, "images/barplot_mortality_rate_state.svg"), plot = mort.rate.state.species, 
       width = 10, height = 6, units = "in")


ggplot(data = st.mort.df %>% filter(n_plots_SPCD > 50))+
  geom_bar(aes(x = region, y = species_volfac_mort), stat = "summary", fun = median,  color = "black")+
  theme_bw(base_size = 16)+
  species_fill+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        axis.title.x = element_blank(), 
        legend.position = "none", 
        panel.grid = element_blank())+
  geom_point(aes(x = region, y = species_volfac_mort, group = Species, color = Species), alpha = 0.9, position =  position_jitter(width = 0.25)) +
  species_color+
  #geom_text(aes(x = Species, y = species_volfac_mort, label = region),position = position_jitter(width = 0.1), alpha = 0.6)
  #ggrepel::geom_text_repel(aes(x = Species, y = species_volfac_mort, label = Abbrev), color = "black", size = 2.5, segment.color = "grey", max.overlaps = 12) +  
  ylab("Mortality Rate (% per year)")

st.mort.df$region <- factor(st.mort.df$region, levels = c("maine", "new hampshire", "vermont", "new york", "connecticut", 
                                                          "new jersey", "pennsylvania", "ohio", "maryland", "west virginia"))
barplot.state.spp.colors <- ggplot(data = st.mort.df %>% filter(n_plots_SPCD > 50))+
  geom_bar(aes(x = region, y = species_volfac_mort, fill = Species), stat = "identity", position = "dodge", fun = median,  color = "black", width = 0.75)+
  theme_bw(base_size = 16)+
  species_fill+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        axis.title.x = element_blank(), 
        legend.position = "none", 
        panel.grid = element_blank())+
  #geom_point(aes(x = region, y = species_volfac_mort, group = Species, color = Species), alpha = 0.9, position =  position_jitter(width = 0.25)) +
  species_color+
  #geom_text(aes(x = Species, y = species_volfac_mort, label = region),position = position_jitter(width = 0.1), alpha = 0.6)
  #ggrepel::geom_text_repel(aes(x = Species, y = species_volfac_mort, label = Abbrev), color = "black", size = 2.5, segment.color = "grey", max.overlaps = 12) +  
  ylab("Mortality Rate (% per year)")

ggsave(paste0(output.dir, "images/barplot_mortality_rate_state_species.png"), plot = barplot.state.spp.colors, 
       width = 10, height = 6, units = "in")

ggsave(paste0(output.dir, "images/barplot_mortality_rate_state_species.svg"), plot = barplot.state.spp.colors, 
       width = 10, height = 6, units = "in")


# plot up all the mortality rates
barplot.mort <- ggplot(data = st.mort.df )+
  geom_bar(aes(x = Species, y = species_volfac_mort, fill = mortality.rate), stat = "identity", color = "black")+
  theme_bw(base_size = 16)+
  facet_wrap(~region)+
  fill_mort_rate+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        #panel.grid = element_blank(), 
        #legend.position = "none", 
        axis.title.x = element_blank())+
  ylab("Average Species Mortality Rate (% per year)")
mort.rate.legend <- get_legend(barplot.mort)

# get state-level information and plot up observed mortality rates:

library(maps)
library(mapdata)
states <- map_data("state")
#9=CT, 25=MA, 33=NH, 23=ME, 50=VT, 44=RI, 42=PA, 39=OH, 54=WV
state_sub <- filter(states, region %in% c("connecticut","maine","new hampshire","vermont","new york", "new jersey",
                                          "rhode island","pennsylvania","ohio","west virginia", "massachusetts", "virginia", "delaware", 
                                          "north carolina", "kentucky", "tennessee", "michigan", "indiana", "district of columbia", "south carolina", 
                                          "georgia", "maryland"))


canada <- map_data("worldHires", "Canada")

# join up mortality data to state_sub
state_sub <- state_sub %>% left_join(., st.mort.df)
state_sub$mortality.rate <- factor(state_sub$mortality.rate, levels = c("> 0.25", "0.25 - 0.5", "0.5 - 0.75", "0.75 - 1", "1 - 1.5", "1.5 - 2", "2 - 4", "4-6"))

# pull the terrain background from stamenmaps
height <- max(state_sub$lat) - min(state_sub$lat)
width <- max(state_sub$long) - min(state_sub$long)
ne_borders <- c(bottom  = min(state_sub$lat)  - 0.1 * height, 
                top     = max(state_sub$lat)  + 0.1 * height,
                left    = min(state_sub$long) - 0.1 * width,
                right   = max(state_sub$long) + 0.1 * width)
#map <- get_stadiamap(ne_borders, zoom = 6, maptype = "stamen_terrain")


library(rnaturalearth)
library(terra)
output.dir <- "C:/Users/KellyHeilman/Box/01. kelly.heilman Workspace/mortality/Eastern-Mortality/mortality_models/"
nat.earth <- terra::rast(paste0(output.dir,"data/NaturalEarth/NE1_50M_SR_W/NE1_50M_SR_W/NE1_50M_SR_W.tif"))
domain <- ext(-87, -65, 36, 48)

nat.crop <- crop(nat.earth, y=domain)
rast.table <- as.data.frame(nat.crop, xy = TRUE)


raster_data <- as.data.frame(nat.crop, xy = TRUE) 

# Combine relevant layers using rgb (you may need to adjust band indices)
raster_data$rgb <- rgb(raster_data[, 3], raster_data[, 4], raster_data[, 5], maxColorValue = 255) # Assuming RGB bands


lakes = ne_download(scale = 50, type = 'lakes', category = 'physical', returnclass = "sf")

# Convert sf object to a data frame with geometry as coordinates
lakes_df <- lakes %>%
  st_as_sf() %>%
  st_cast("POLYGON") %>% # Ensure geometry is in polygon format
  st_coordinates() %>%
  as.data.frame()

# Add grouping variables for polygons
lakes_df <- lakes_df %>%
  rename(long = X, lat = Y) %>%
  mutate(group = interaction(L1, L2)) # Grouping for polygons

# create a basemap from natural earth raster:
NE_basemap <- ggplot(data = raster_data, aes(x = x, y = y)) +
  geom_tile(fill = raster_data$rgb, alpha = 0.75) +
  #geom_sf(data = lakes, fill = "lightblue") + theme_bw()+
  coord_sf(xlim = c(-85.5, -67.5), ylim = c(37, 47.5))+
  theme_void() 


mortality_rate_maps <- list()

for(i in 1:length(nspp$SPCD)){
  
  state_sub_species <- state_sub %>% filter(Species %in% nspp[i,]$COMMON &  n_plots_SPCD >= 25)
  
  
  
  mortality_rate_maps[[i]] <- NE_basemap + 
    geom_polygon(data = canada, 
                 aes(x=long, y=lat, group = group), 
                 color = "black", fill = "white") +
    geom_polygon(data = lakes_df , 
                 aes(x = long, y = lat, group = group), 
                 color = "black", fill = "lightblue") +
    geom_polygon(data = state_sub, 
                 aes(x=long, y=lat, group = group), 
                 color = "black", fill = "white") +
    
    geom_polygon(data = state_sub_species, 
                 aes(x=long, y=lat, group = group, fill = mortality.rate), 
                 color = "black", alpha = 0.9) +
    coord_sf(xlim = c(-85, -67.5), ylim = c(37, 47.5))+theme(panel.grid = element_blank(), #panel.background = element_rect(fill = 'lightblue'), 
                                                             legend.position = "none",
                                                             axis.title  = element_blank(),
                                                             #legend.title = element_blank(),
                                                             legend.background = element_rect(fill = "white", color = "black") 
    )+
    ggtitle(paste0(disturb.species.order[i]))+
    fill_mort_rate
  
}


mort.maps.species <- plot_grid(
  plot_grid(plotlist = mortality_rate_maps, ncol = 3, align = "hv"), 
  mort.rate.legend, rel_widths = c(1, 0.1), ncol = 2
)

ggsave(paste0(output.dir, "images/map_mortality_rate_state.png"), plot = mort.maps.species, 
       width = 12, height = 13, units = "in")

ggsave(paste0(output.dir, "images/map_mortality_rate_state.svg"), plot = mort.maps.species, 
       width = 12, height = 13, units = "in")

# observed vs predicted mortality rates -----
# predictions compiled in R/hierarchicalModel/posterior_prediction_joint_longchains.R
# read in predicted mortality by species by tree:
spp.tree.pred.mort <- do.call(rbind, lapply(list.files(path = paste0(
  output.folder,
  "SPCD_stanoutput_joint_v3/predicted_mort/Mort_10yr/"), 
  pattern = "tree_10yr_mortality_", 
  full.names = TRUE), 
  readRDS))

# scale tree-level predictions of mortality by volfac to estimate mortality on trees per acre scale:

# combine with volfac estimates by tree:
# volfacs are listed in all.remeas, so we need to join spp.tree.pred.mort

spp.tree.pred.mort.volfac <- left_join(spp.tree.pred.mort, 
                                       all.remeas %>% select(state, county, pltnum, cndtn, point, tree, SPCD,status, dbhcur, dbhold, volfac, LAT_FIADB, LONG_FIADB)%>% distinct())

# map up probability of survival for each species
ggplot(data = spp.tree.pred.mort.volfac %>% filter(COMMON %in% c("balsam fir", "red spruce", "northern white-cedar")), aes(x =LONG_FIADB,  y = LAT_FIADB, color = 1-p10year.med), size = 0.1)+
  geom_point()+facet_wrap(~COMMON)+ scale_color_viridis_c(option = "magma")

# get the median 1 year mortality rate for each tree:
spp.tree.mort.probs <- spp.tree.pred.mort.volfac %>% mutate(p1year.survival = p10year.med^(1/10))%>%
  mutate(p1year.mortality = 1- p1year.survival) %>%
  rename("Species"="COMMON")


ggplot(data = spp.tree.mort.probs)+
  geom_density(aes(y = p1year.mortality*100, fill = Species))+
  facet_wrap(~Species, scales = "free")+species_fill

# get plot-level averages of species rate estimates for the modelled rate:--
avg.plt.mort.annual <- spp.tree.mort.probs %>% group_by(state, PLOT.ID, Species, SPCD, LONG_FIADB, LAT_FIADB)%>%
  summarise(pMort_annual_plot = median(p1year.mortality, na.rm =TRUE), 
            pMort_annual_plot.sd = sd(p1year.mortality, na.rm =TRUE), 
            nTrees = n())



ggplot(data = avg.plt.mort.annual)+
  geom_density(aes(y = pMort_annual_plot*100, fill = Species))+
  facet_wrap(~Species, scales = "free")+species_fill

ggplot(data = avg.plt.mort.annual)+
  geom_boxplot(aes(x = Species, y = pMort_annual_plot*100, fill = Species))+
  species_fill

ggplot(data = avg.plt.mort.annual)+
  geom_density(aes(y = pMort_annual_plot*100, fill = Species))+
  facet_wrap(~Species, scales = "free")+species_fill

ggplot(data = avg.plt.mort.annual)+
  geom_point(aes(x = LONG_FIADB, y = LAT_FIADB, color = pMort_annual_plot*100), size = 0.5)+
  facet_wrap(~Species)+
  scale_color_viridis_c()
hist(avg.plt.mort.annual$pMort_annual_plot*100)

# Create discrete categories from continuous_var
avg.plt.mort.annual$pMort_annual_plot_discrete <- cut(avg.plt.mort.annual$pMort_annual_plot*100,
                                breaks = c(0, 0.5, 1, 1.5, 2, 5, 7.5, 10, 20  ),
                                labels = c("< 0.5", "0.5 - 1", "1 - 1.5", "1.5 - 2", "2 - 5", "5 - 7.5", "7.5 - 10", "> 10"))
nspp$COMMON

disturb.species.order <- c("balsam fir" ,"red spruce","northern white-cedar","eastern hemlock",     
                          "American beech" ,"chestnut oak","northern red oak",    
                           "white oak", "yellow birch" ,"paper birch" ,"hickory spp.",        
                            "eastern white pine","red maple","sugar maple" ,"black cherry",      
                            "white ash" ,"yellow-poplar"   )

species.plot.mortality.list <- list()
for(i in 1:length(disturb.species.order)){
  
species.plot.mortality.list[[i]] <- ggplot() + 
  geom_polygon(data = canada, 
               aes(x=long, y=lat, group = group), 
               color = "black", fill = "white") +
  geom_polygon(data = lakes_df , 
               aes(x = long, y = lat, group = group), 
               color = "black", fill = "lightblue") +
  geom_polygon(data = state_sub, 
               aes(x=long, y=lat, group = group), 
               color = "black", fill = "white") +
 
  geom_point(data = avg.plt.mort.annual %>% filter(Species %in% disturb.species.order[i]) , aes(x = LONG_FIADB, y = LAT_FIADB, color = pMort_annual_plot_discrete), size = 0.1)+
  scale_color_manual(values =  c(
    "#74add1",
    "#abd9e9",
   
    "#fcbba1",
    "#fc9272",
    "#fb6a4a",
    "#ef3b2c",
    "#cb181d",
    "#99000d"), name = "Mortality Probability\n(%/year)")+
  

  coord_sf(xlim = c(-85, -67.5), ylim = c(37, 47.5))+theme(panel.grid = element_blank(), #panel.background = element_rect(fill = 'lightblue'), 
                                                          legend.position = "none",
                                                           axis.title  = element_blank(),
                                                           #legend.title = element_blank(),
                                                           legend.background = element_rect(fill = "white", color = "black") 
  ) + ggtitle(paste0(disturb.species.order[i]))
}

mort.pred.legend <- get_legend(ggplot() + 
                                 geom_point(data = avg.plt.mort.annual  , aes(x = LONG_FIADB, y = LAT_FIADB, color = pMort_annual_plot_discrete), size = 3)+
                                 scale_color_manual(values =  c(
                                   "#74add1",
                                   "#abd9e9",
                                   
                                   "#fcbba1",
                                   "#fc9272",
                                   "#fb6a4a",
                                   "#ef3b2c",
                                   "#cb181d",
                                   "#99000d"), name = "Mortality Probability\n(%/year)")+
                                 
                                 
                                 coord_sf(xlim = c(-85, -67.5), ylim = c(37, 47.5))+theme(panel.grid = element_blank(), #panel.background = element_rect(fill = 'lightblue'), 
                                                                                          #legend.position = "none",
                                                                                          axis.title  = element_blank(),
                                                                                          #legend.title = element_blank(),
                                                                                          legend.background = element_rect(fill = "white", color = "black") 
                                 ) )
mort.pred.maps.species <- plot_grid(
  plot_grid(plotlist = species.plot.mortality.list, ncol = 3, align = "hv"), 
  mort.pred.legend, rel_widths = c(1, 0.2), ncol = 2
)

ggsave(paste0(output.dir, "images/map_annual_mortality_prob_plot.png"), plot = mort.pred.maps.species, 
       width = 12, height = 13, units = "in")

ggsave(paste0(output.dir, "images/map_annual_mortality_prob_plot.svg"), plot = mort.pred.maps.species, 
       width = 12, height = 13, units = "in")


# ideas for improving this figure:
# summarise by county (need to scale by volfac)
# alternative colors?

# Looking at four attibutes of tree mortality-----
# after Davis et al. hemlock pollen decline:
# Quaternary history and the stability of forest communities. Pages 132–153 in D. C. West, H. H. Shugart, and D. B. Botkin, eds. Forest succession: concepts and applications. Springer, New York.

# Specificity: Is it only one species or a group?--
# Synchrony: Is it happening at the same time across the range?--
# Rapidity: ??--
  # dia_diff as metric? lower diameter differences mean less growth prior to mortality?
  # may not be able to answer
# Uniformity: Is it happening across the range at similar rates?--

# Specificity: Are Species proability of mortality correlated with one another?
nmort.df <- spp.tree.mort.probs %>% group_by(state, PLOT.ID, Species, SPCD, LONG_FIADB, LAT_FIADB)%>%
  summarise(pMort_annual_plot = median(p1year.mortality, na.rm =TRUE), 
            pMort_annual_plot.sd = sd(p1year.mortality, na.rm =TRUE), 
            nTrees = n(), 
            Nmort = sum(mortality.draw))


# set up a dataframe with plots as rows and columns the pmort for different species
pmort.spread <- avg.plt.mort.annual %>% ungroup()%>% 
  mutate(pMort.pct  = pMort_annual_plot*100)%>%
  select(LAT_FIADB, LONG_FIADB, PLOT.ID,Species,pMort.pct)%>%
  group_by(LAT_FIADB, LONG_FIADB, PLOT.ID) %>% 
  spread(Species,pMort.pct,  fill = NA)
#p_correlations <- cor(pmort.spread[,4:ncol(pmort.spread)], use = "pairwise.complete")
# library(ggcorrplot)
# library(Hmisc)
# library(GGally)
# 
# 
# ggpairs(pmort.spread, columns = 4:ncol(pmort.spread))

library(ggcorrplot)
library(Hmisc)
p_correlations <- rcorr(as.matrix(pmort.spread[,4:ncol(pmort.spread)]), type = "spearman" ) # or "spearman"

# lets look at correlations for species that co-occur in at least 10 plots
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
melted.correlation <- get_lower_tri(p_correlations$r[disturb.species.order,disturb.species.order]) %>% reshape2::melt()%>%
  rename("r"="value") %>% left_join(., 
                                    get_lower_tri(p_correlations$n[disturb.species.order,disturb.species.order]) %>% reshape2::melt())%>%
  rename("n" = "value")%>%
  left_join(., 
            get_lower_tri(p_correlations$P[disturb.species.order,disturb.species.order]) %>% reshape2::melt())%>%
  rename("pval" = "value")%>%
  mutate(R_revised = ifelse(n>=50, round(r, digits = 1), NA), 
         P_revised = ifelse(n >=50, pval, NA))%>%
  mutate(R_revised = ifelse(R_revised == 1, NA, R_revised))%>%
  mutate(sig.label = ifelse(P_revised <= 0.05, R_revised, NA))%>%
  filter(!is.na(R_revised))# if the species are colocated in less than 10 plots omit the correlation


pmort_correlation <-  ggplot(data = melted.correlation, aes(Var2, Var1, fill = R_revised))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-0.6,0.6), space = "Lab", 
                       name="Spearman\nCorrelation") +
  theme_bw(base_size = 16)+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                    hjust = 0))+
  coord_fixed()+
  scale_x_discrete(position = "top",
                   limits = levels(melted.correlation$Var2))+
  scale_y_discrete(position = "left",
                   limits = rev(levels(melted.correlation$Var1)))+
 
  geom_text(aes(Var2, Var1, label = sig.label), color = "black", size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    #panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.9, 0.65),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))#+coord_flip()
  

ggsave(filename = paste0(output.dir, "images/posterior_plot_pMort_species_correlations.png"), 
       plot = pmort_correlation, 
       height = 6, width = 6, units = "in", 
       dpi = 350)

# fourteen species pairs have R_revised >= 0.3
melted.correlation %>% filter(R_revised >=0.3)
# northern region:
# balsam fir - Red spruce
# balsam fir - eastern hemlock

# eastern hemlock - yellow birch
# American beech - red maple

# American beech - red maple

# oaks:
# chestnut oak - N red oak
# chestnut oak - white oak

# chestnut oak - E. white pine
# chestnut oak - red maple

# N red oak - paper birch
# hickor spp. - northern red oak

# white oak - yellow birch 
# white oak - red maple

# red maple - paper birch
# sugar maple - red maple

# four species pairs have correlations >= 0.4
melted.correlation %>% filter(R_revised >=0.4)


# Synchrony and uniformity across space:
# moran's global 
# moran's local 
# Ripley's K-function:



data_built <- ggplot_build(p)

# The data for the filled contours is in the second layer
data_2d <- data_built$data[[2]]

# The head() function shows the structure
head(data_2d)

library(spdep)
polygon.contour.list <- polygon.buffer.list <- buffer.map.list <- moran.local.list <- polygon_hull.list <- list()
for(i in 1:17){
  
plt_mort_sf <- st_as_sf(avg.plt.mort.annual %>% filter(Species %in% nspp[i,]$COMMON), coords = c("LONG_FIADB", "LAT_FIADB"))
# get the spatial weights matrix 

plot_neighbors <- knn2nb(knearneigh(plt_mort_sf, k = 10)) # 10 nearest neighbors

# get spatial weights 
plot_weights <- nb2listw(plot_neighbors, style = "W")

# global Moran's I
moran_global <- moran.test(plt_mort_sf$pMort_annual_plot, plot_weights)
print(moran_global)

# Moran plot 
moran.plot(plt_mort_sf$pMort_annual_plot, plot_weights, #labels = as.character(tree_data$id),
           main = "Moran Plot for Tree Mortality Risk")


#Calculate Local Moran's I
local_moran_result <- localmoran(plt_mort_sf$pMort_annual_plot, plot_weights)

# Add the local Moran statistics to the sf object
plt_mort_sf$local_moran_I <- local_moran_result[, 1]
plt_mort_sf$p_value <- local_moran_result[, 5]

# get the type of local cluster
# standardization and  filter by significance
plt_mort_sf$cluster_type <- "Not Significant"
plt_mort_sf$standardized_risk <- scale(plt_mort_sf$pMort_annual_plot)
plt_mort_sf$quadrant <- local_moran_result[, 3] # Quadrant for the moran plot

# assign labels
significant_clusters <- which(plt_mort_sf$p_value <= 0.05)
plt_mort_sf$cluster_type[significant_clusters] <- ifelse(
  plt_mort_sf$standardized_risk[significant_clusters] > 0 & plt_mort_sf$quadrant[significant_clusters] > 0, "High-High",
  ifelse(plt_mort_sf$standardized_risk[significant_clusters] < 0 & plt_mort_sf$quadrant[significant_clusters] < 0, "Low-Low",
         ifelse(plt_mort_sf$standardized_risk[significant_clusters] > 0 & plt_mort_sf$quadrant[significant_clusters] < 0, "High-Low",
                ifelse(plt_mort_sf$standardized_risk[significant_clusters] < 0 & plt_mort_sf$quadrant[significant_clusters] > 0, "Low-High", "Not Significant"))))

# plot the  clusters
# ggplot(plt_mort_sf) +
#   geom_sf(aes(color = cluster_type), size = 3) +
#   scale_color_manual(values = c("High-High" = "red", "Low-Low" = "blue", "High-Low" = "pink", "Low-High" = "lightblue", "Not Significant" = "grey")) +
#   theme_minimal()
# 

# identify hotspots of mortality risk and create a polygon
high_high_clusters <- plt_mort_sf %>%
  filter(p_value < 0.05 & local_moran_I > 0)# & pMort_annual_plot >= 0.01)



# Convex hull around all HH cluster points
if (nrow(high_high_clusters) > 0) {
  hh_polygon_chull <- st_convex_hull(st_union(high_high_clusters))
} else {
  message("No HH clusters found to create a polygon.")
  hh_polygon_chull <- NULL
}
if (!is.null(hh_polygon_chull)) {
  ggplot() +
    geom_sf(data = plt_mort_sf, aes(color = pMort_annual_plot), size = 2) +
    geom_sf(data = high_high_clusters, color = "red", size = 3) +
    geom_sf(data = hh_polygon_chull, fill = "red", alpha = 0.3, color = NA) +
    labs(title = "HH Clusters and Convex Hull Polygon") +
    theme_minimal()
}

# Create a buffer around each High-High cluster point, then union
# You might want to adjust the buffer distance (dist)
if (nrow(high_high_clusters) > 0) {
  hh_polygon_buffer <- st_union(st_buffer(high_high_clusters, dist = 0.5))
} else {
  message("No HH clusters found to create a polygon.")
  hh_polygon_buffer <- NULL
}

gg.buffer.map <- ggplot() + 
  geom_polygon(data = canada, 
               aes(x=long, y=lat, group = group), 
               color = "black", fill = "white") +
  geom_polygon(data = lakes_df , 
               aes(x = long, y = lat, group = group), 
               color = "black", fill = "lightblue") +
  geom_polygon(data = state_sub, 
               aes(x=long, y=lat, group = group), 
               color = "black", fill = "white")+
    geom_sf(data = plt_mort_sf, aes(color = pMort_annual_plot_discrete), size = 2) +
    #geom_sf(data = high_high_clusters, color = "red", size = 3) +
    geom_sf(data = hh_polygon_buffer, fill = "red", alpha = 0.3, color = NA) +
    #labs(title = "HH Clusters and Convex Hull Polygon") +
    theme_minimal()+
  scale_color_manual(values =  c(
    "#74add1",
    "#abd9e9",
    
    "#fcbba1",
    "#fc9272",
    "#fb6a4a",
    "#ef3b2c",
    "#cb181d",
    "#99000d"), name = "Mortality Probability\n(%/year)")+
  
  
  coord_sf(xlim = c(-85, -67.5), ylim = c(37, 47.5))+theme(panel.grid = element_blank(), #panel.background = element_rect(fill = 'lightblue'), 
                                                           #legend.position = "none",
                                                           axis.title  = element_blank(),
                                                           #legend.title = element_blank(),
                                                           legend.background = element_rect(fill = "white", color = "black") 
  ) +ggtitle(paste0("mortality risk hotspots for ", nspp[i,]$COMMON))


# # use spatstat package to get KDE
library(spatstat)
spp.mort.rate <- avg.plt.mort.annual %>% filter(Species %in% nspp[i,]$COMMON)
# spp_mort_high_sf <- st_as_sf(spp.mort.rate %>% filter(pMort_annual_plot >= mean(spp.mort.rate$pMort_annual_plot)), coords = c("LONG_FIADB", "LAT_FIADB"), crs = 4269)%>%
#   st_transform(., crs = 5070)
# 
# # use state spatial polygons as the windown for the density object using as.owin 
# us_polygon_df <-map("usa", fill = TRUE, plot = FALSE) %>% st_as_sf() %>%
#   st_transform(., crs = 5070)
# us_state_owin <- as.owin(us_polygon_df)
# 
# ppp_obj_us <- as.ppp(st_coordinates(spp_mort_high_sf), W = us_state_owin)
# kde_result_us <- density(ppp_obj_us, sigma = 12)


# Convert the im object to a data frame for ggplot2
# kde_df <- as.data.frame(kde_result_us)

# Plot with filled contours
# ggplot(kde_df %>% filter(value >=0), aes(x = x, y = y, z = value)) +
#   geom_contour_filled() +
#   scale_fill_viridis_d() # Or other color scales
# 
# plot(kde_result_us)

# Define contour levels (e.g., at 25%, 50%, 75% of the maximum density)
# levels <- quantile(kde_result_us$v, probs = c(0.25, 0.50, 0.75), na.rm =TRUE)
# 
# # Extract contours
# contours_list <- contour(kde_result_us, levels = levels, draw = FALSE)


#mean(spp.mort.rate$pMort_annual_plot)
# if there are any points with > 1% mortality per year:
avg.plt.mort.annual %>% filter(pMort_annual_plot >= 0.005) %>% ungroup()%>% dplyr::select(Species) %>% distinct()

# get 2D KDE of above average mortality rates form stat_density2d
p <- ggplot(data = spp.mort.rate %>% filter(pMort_annual_plot >= quantile(spp.mort.rate$pMort_annual_plot, 0.75)), aes(x = LONG_FIADB, y = LAT_FIADB)) +
  geom_point(alpha = 0.5) + # Optional: show the original points
  coord_sf(xlim = c(-85, -60), ylim = c(30, 50))+
  
  stat_density2d(aes(fill = ..level..,), alpha = 0.5, geom = "polygon", contour_var = "count",#, contour_var = "count",
                 breaks =  c(0,0.5, 1, 2, 5, 10, 15, 20, 25, 30, 60)) +
  scale_fill_viridis_c() + # A color scale for the density levels
  #coord_equal() + # Ensures correct aspect ratio for spatial data
  labs(title = "Density of higher than average mortality",
       x = "Longitude",
       y = "Latitude") +
  theme_minimal()
p


p2 <- ggplot(data = spp.mort.rate %>% filter(pMort_annual_plot >= quantile(spp.mort.rate$pMort_annual_plot, 0.5)), aes(x = LONG_FIADB, y = LAT_FIADB)) +
  geom_point(alpha = 0.5) + # Optional: show the original points
  coord_sf(xlim = c(-85, -60), ylim = c(30, 50))+
  #geom_hdr(probs = c(0.25, 0.5, 0.75, 0.95))
  
   geom_density2d(aes(color = ..level..,), alpha = 0.5)+ 
  #                breaks = density_levels) +
  scale_fill_viridis_c() + # A color scale for the density levels
  #coord_equal() + # Ensures correct aspect ratio for spatial data
  labs(title = "Density of higher than average mortality",
       x = "Longitude",
       y = "Latitude") +
  theme_minimal()







data_built <- ggplot_build(p2)

# The data for the filled contours is in the second layer
data_2d <- data_built$data[[2]]

# set a polygon IDs for each contour piece
data_2d$pol_id <- paste0(data_2d$group)

# Split and create sf polygons based on the IDs
polygons_list <- lapply(unique(data_2d$pol_id), function(x) {
  # Subset data for the current polygon
  polygon_data <- data_2d[data_2d$pol_id == x, ]
  
  # close polygon
  closed_polygon <- rbind(polygon_data, polygon_data[1, ])
  
  # sf polygon 
  polygon_sf <- st_polygon(list(as.matrix(closed_polygon[, c("x", "y")])))
  
  # density level 
  level_data <- unique(polygon_data[, grepl("level", names(polygon_data))])
  
  # sf object with the level data
  st_as_sf(level_data, geometry = st_sfc(polygon_sf))
})

# combine species contours into a df
final_polygons_sf <- do.call(rbind, polygons_list) %>% 
  mutate(Species = nspp[i,]$COMMON, 
         mean.plt.mort.spp = mean(spp.mort.rate$pMort_annual_plot))

final_polygons_sf %>% filter(nlevel >=0.5) %>% ggplot()+geom_sf(aes(fill = nlevel))+scale_fill_viridis_c()+
  geom_point(data = spp.mort.rate %>% filter(pMort_annual_plot >= 0.01), aes(x = LONG_FIADB, y = LAT_FIADB), size = 0.1)




# save each species buffers:
polygon.contour.list[[i]] <- final_polygons_sf
polygon.buffer.list[[i]] <- hh_polygon_buffer
polygon_hull.list[[i]] <- hh_polygon_chull
buffer.map.list[[i]] <- gg.buffer.map
moran.local.list[[i]] <- plt_mort_sf

rm(hh_polygon_buffer, hh_polygon_chull, gg.buffer.map, final_polygons_sf)
}

names(polygon.contour.list) <- nspp$COMMON
names(polygon.buffer.list) <- nspp$COMMON
names(polygon_hull.list) <- nspp$COMMON
#spatial_buffer_points <- st_sfc(polygon.buffer.list[[1]], polygon.buffer.list[[2]], polygon.buffer.list[[3]])


poly.contours <- do.call(rbind, polygon.contour.list)%>% #data.frame()%>%
  group_by(Species) %>% mutate(max.level = max(level))

poly.contours %>% as.data.frame ()%>% group_by(Species)%>% summarise(max.level = max(level), min (level))%>% 
  arrange(max.level)


polygon.contour.list[[1]] %>% filter(nlevel > 0.5) %>% ggplot()+geom_sf(aes(fill = level))+scale_fill_viridis_c()
polygon.contour.list[[2]] %>% filter(nlevel > 0.5) %>% ggplot()+geom_sf(aes(fill = level))+scale_fill_viridis_c()
polygon.contour.list[[3]]  %>% filter(level >= 15) %>% ggplot()+geom_sf(aes(fill = level))+scale_fill_viridis_c()
polygon.contour.list[[4]]   %>% filter(level >= 15)%>% ggplot()+geom_sf(aes(fill = level))+scale_fill_viridis_c()
polygon.contour.list[[5]]  %>% filter(level >= 15)%>% ggplot()+geom_sf(aes(fill = level))+scale_fill_viridis_c()
polygon.contour.list[[6]]  %>% filter(level >= 15)%>% ggplot()+geom_sf(aes(fill = level))+scale_fill_viridis_c()
polygon.contour.list[[7]]   %>% filter(level >= 15)%>% ggplot()+geom_sf(aes(fill = level))+scale_fill_viridis_c()
polygon.contour.list[[8]]   %>% filter(level >= 15)%>% ggplot()+geom_sf(aes(fill = level))+scale_fill_viridis_c()
polygon.contour.list[[9]]   %>% filter(level >= 15)%>% ggplot()+geom_sf(aes(fill = level))+scale_fill_viridis_c()
polygon.contour.list[[10]] %>% filter(level >= 15)%>% ggplot()+geom_sf(aes(fill = level))+scale_fill_viridis_c()
polygon.contour.list[[11]]  %>% filter(level >= 15)%>% ggplot()+geom_sf(aes(fill = level))+scale_fill_viridis_c()
polygon.contour.list[[12]]   %>% filter(level >= 15)%>% ggplot()+geom_sf(aes(fill = level))+scale_fill_viridis_c()
polygon.contour.list[[13]]   %>% filter(level >= 15)%>% ggplot()+geom_sf(aes(fill = level))+scale_fill_viridis_c()
polygon.contour.list[[14]]  %>% filter(level >= 15) %>% ggplot()+geom_sf(aes(fill = level))+scale_fill_viridis_c()
polygon.contour.list[[15]]   %>% filter(level >= 15)%>% ggplot()+geom_sf(aes(fill = level))+scale_fill_viridis_c()
polygon.contour.list[[16]]   %>% filter(level >= 15)%>% ggplot()+geom_sf(aes(fill = level))+scale_fill_viridis_c()
polygon.contour.list[[17]]   %>% filter(level >= 15)%>% ggplot()+geom_sf(aes(alpha = level), fill = "forestgreen")#+scale_fill_viridis_c()


# plot the point-based buffers all together:
base.map <- ggplot() + 
  geom_polygon(data = canada, 
               aes(x=long, y=lat, group = group), 
               color = "black", fill = "white") +
  geom_polygon(data = lakes_df , 
               aes(x = long, y = lat, group = group), 
               color = "black", fill = "lightblue") +
  geom_polygon(data = state_sub, 
               aes(x=long, y=lat, group = group), 
               color = "black", fill = "white")+
  theme_minimal()

# plot kde based polygons for each species:
all.kde <- base.map +
  geom_sf(data = poly.contours %>% filter(nlevel >= 0.75) , aes(fill = Species, group = Species), alpha = 0.7) +
  species_fill+coord_sf(xlim = c(-85, -67.5), ylim = c(37, 47.5))+theme(panel.grid = element_blank(), #panel.background = element_rect(fill = 'lightblue'), 
                                                                       legend.position = "none",
                                                                       axis.title  = element_blank(),
                                                                       #legend.title = element_blank(),
                                                                       legend.background = element_rect(fill = "white", color = "black") )


TREE.remeas %>% dplyr::select(state, stname, date, remper) %>% distinct() %>%
  filter(!remper == 0) %>%
  group_by(state, stname, date)%>%
  summarise(avg.remper = median(remper))%>%
  mutate(date2 = date - avg.remper)

# most mortality here is 1981-1995 (1997)
northern.kde <- base.map +
  geom_sf(data = poly.contours %>% 
            filter(Species %in% c("balsam fir", "northern white-cedar",  "red spruce"))%>%
                     filter(nlevel > 0.75), aes(fill = Species, group = Species), alpha = 0.75, color = NA) +
  species_fill+coord_sf(xlim = c(-85, -67.5), ylim = c(37, 47.5))+theme(panel.grid = element_blank(), #panel.background = element_rect(fill = 'lightblue'), 
                                                                        legend.position = "none",
                                                                        axis.title  = element_blank(),
                                                                        #legend.title = element_blank(),
                                                                        legend.background = element_rect(fill = "white", color = "black") )



  


# peaks of mortality in PA for Oaks, in NY, NH, and VT, and ME for birches
# most mortality here is 1976-1989 for PA and WV and 1980-1997ish for birches
# Oaks:
# PA: remper 1976-1989 # peak SM defoliation = 1990
# NJ: remper 1972-1987 # peak SM defoliation = 1981
# Birches:
# NY: remper = 1980-1993 # peak SM defoliation = 1980
# VT: remper = 1982-1997 # peak SM defoliation = 1953, but some later on
# NH: remper = 1983 - 1997 # peak SM defoliation = 1981, and large peaks in 1980s


oak.kde <- base.map +
  geom_sf(data = poly.contours %>% 
            filter(Species %in% c( "chestnut oak","northern red oak", "white oak", "paper birch", "yellow birch"))%>%
            filter(nlevel >= 0.70), aes(fill = Species, group = Species), alpha = 0.75, color = NA) +
  species_fill+coord_sf(xlim = c(-85, -67.5), ylim = c(37, 47.5))+theme(panel.grid = element_blank(), #panel.background = element_rect(fill = 'lightblue'), 
                                                                        legend.position = "none",
                                                                        axis.title  = element_blank(),
                                                                        #legend.title = element_blank(),
                                                                        legend.background = element_rect(fill = "white", color = "black") )


  
  

# most mortality is in PA, NY, VT, NH, an Maine for both
# could be HWA related for hemlock in PA and NY
# PA = 1976 - 1989 (HWA first detected 1979)
# NY = 1980 - 1993 (HWA first detected 1984)
# VT = 1982 - 1997 (HWA first detected 2008)-Hemlock looper?
# NH = 1983 - 1997 (HWA first detected 2004)-Hemlock looper?


beech.hemlock<- base.map +
  geom_sf(data = poly.contours %>% 
            filter(Species %in% c("eastern hemlock", "American beech"))%>%
            filter(nlevel > 0.70), aes(fill = Species, group = Species), alpha = 0.75, color = NA) +
  species_fill+coord_sf(xlim = c(-85, -67.5), ylim = c(37, 47.5))+theme(panel.grid = element_blank(), #panel.background = element_rect(fill = 'lightblue'), 
                                                                        legend.position = "none",
                                                                        axis.title  = element_blank(),
                                                                        #legend.title = element_blank(),
                                                                        legend.background = element_rect(fill = "white", color = "black") )





# hickory mortality is WV, OH = 
spongy.resistant.kde <- base.map +
  geom_sf(data = poly.contours %>% 
            filter(Species %in% c("sugar maple", "red maple", "eastern white pine", "hickory spp."))%>%
            filter(nlevel >= 0.7), aes(fill = Species, group = Species), alpha = 0.75, color = NA) +
  species_fill+coord_sf(xlim = c(-85, -67.5), ylim = c(37, 47.5))+theme(panel.grid = element_blank(), #panel.background = element_rect(fill = 'lightblue'), 
                                                                        legend.position = "none",
                                                                        axis.title  = element_blank(),
                                                                        #legend.title = element_blank(),
                                                                        legend.background = element_rect(fill = "white", color = "black") )






spongy.immune.kde <- base.map +
  geom_sf(data = poly.contours %>% 
            filter(Species %in% c("black cherry", "yellow-poplar", "white ash"))%>%
            filter(nlevel >=0.7), aes(fill = Species, group = Species), alpha = 0.75, color = NA) +
  species_fill+coord_sf(xlim = c(-85, -67.5), ylim = c(37, 47.5))+theme(panel.grid = element_blank(), #panel.background = element_rect(fill = 'lightblue'), 
                                                                        legend.position = "none",
                                                                        axis.title  = element_blank(),
                                                                        #legend.title = element_blank(),
                                                                        legend.background = element_rect(fill = "white", color = "black") )








# all species point-based buffers together:
base.map +
  geom_sf(data = polygon.buffer.list[["balsam fir"]], aes(fill = "balsam fir"), alpha = 0.5, color = NA) +
  geom_sf(data = polygon.buffer.list[["red spruce"]], aes(fill = "red spruce"), alpha = 0.5, color = NA) +
  geom_sf(data = polygon.buffer.list[["northern white-cedar"]], aes(fill = "northern white-cedar"), alpha = 0.5, color = NA) +
  
  geom_sf(data = polygon.buffer.list[["eastern hemlock"]], aes(fill = "eastern hemlock"), alpha = 0.5, color = NA) +
  geom_sf(data = polygon.buffer.list[["American beech"]], aes(fill = "American beech"), alpha = 0.5,color = NA) +
  
  
  geom_sf(data = polygon.buffer.list[["northern red oak"]], aes(fill = "northern red oak"), alpha = 0.5, color = NA) +
  geom_sf(data = polygon.buffer.list[["chestnut oak"]], aes(fill = "chestnut oak"), alpha = 0.5,color = NA) +
  geom_sf(data = polygon.buffer.list[["white oak"]], aes(fill = "white oak"), alpha = 0.5, color = NA) +
  geom_sf(data = polygon.buffer.list[["yellow birch"]], aes(fill = "yellow birch"), alpha = 0.5,color = NA) +
  geom_sf(data = polygon.buffer.list[["paper birch"]], aes(fill = "paper birch"), alpha = 0.5,color = NA) +
  geom_sf(data = polygon.buffer.list[["hickory spp."]], aes(fill = "hickory spp."), alpha = 0.5,color = NA) +
  geom_sf(data = polygon.buffer.list[["eastern white pine"]], aes(fill = "eastern white pine"), alpha = 0.5,color = NA) +
  geom_sf(data = polygon.buffer.list[["sugar maple"]], aes(fill = "sugar maple"), alpha = 0.5,color = NA) +
  geom_sf(data = polygon.buffer.list[["red maple"]], aes(fill = "red maple"), alpha = 0.5,color = NA) +
  geom_sf(data = polygon.buffer.list[["black cherry"]], aes(fill = "black cherry"), alpha = 0.5,color = NA) +
  geom_sf(data = polygon.buffer.list[["yellow-poplar"]], aes(fill = "yellow-poplar"), alpha = 0.5,color = NA) +
  geom_sf(data = polygon.buffer.list[["white ash"]], aes(fill = "white ash"), alpha = 0.5,color = NA) +
  species_fill+coord_sf(xlim = c(-85, -67.5), ylim = c(37, 47.5))+theme(panel.grid = element_blank(), #panel.background = element_rect(fill = 'lightblue'), 
                                                                        #legend.position = "none",
                                                                        axis.title  = element_blank(),
                                                                        #legend.title = element_blank(),
                                                                        legend.background = element_rect(fill = "white", color = "black") 
  )


spruce.fir.hotspots <- base.map +
  geom_sf(data = polygon.buffer.list[["balsam fir"]], aes(fill = "balsam fir"), alpha = 0.9, color = NA) +
  geom_sf(data = polygon.buffer.list[["red spruce"]], aes(fill = "red spruce"), alpha = 0.9, color = NA) +
  geom_sf(data = polygon.buffer.list[["northern white-cedar"]], aes(fill = "northern white-cedar"), alpha = 0.5, color = NA) +
  
  #geom_sf(data = polygon.buffer.list[["eastern hemlock"]], aes(fill = "eastern hemlock"), alpha = 0.5, color = NA) +
  #geom_sf(data = polygon.buffer.list[["American beech"]], aes(fill = "American beech"), alpha = 0.5,color = NA) +
  
  
 
  species_fill+coord_sf(xlim = c(-85, -67.5), ylim = c(37, 47.5))+theme(panel.grid = element_blank(), #panel.background = element_rect(fill = 'lightblue'), 
                                                                        #legend.position = "none",
                                                                        axis.title  = element_blank(),
                                                                        #legend.title = element_blank(),
                                                                        legend.background = element_rect(fill = "white", color = "black") 
  )

oak.hotspots <- base.map +
  geom_sf(data = polygon.buffer.list[["northern red oak"]], aes(fill = "northern red oak"), alpha = 0.5, color = NA) +
  geom_sf(data = polygon.buffer.list[["chestnut oak"]], aes(fill = "chestnut oak"), alpha = 0.5,color = NA) +
  geom_sf(data = polygon.buffer.list[["white oak"]], aes(fill = "white oak"), alpha = 0.5, color = NA) +
  geom_sf(data = polygon.buffer.list[["yellow birch"]], aes(fill = "yellow birch"), alpha = 0.5,color = NA) +
  geom_sf(data = polygon.buffer.list[["paper birch"]], aes(fill = "paper birch"), alpha = 0.5,color = NA) +
  
  species_fill+coord_sf(xlim = c(-85, -67.5), ylim = c(37, 47.5))+theme(panel.grid = element_blank(), #panel.background = element_rect(fill = 'lightblue'), 
                                                                        #legend.position = "none",
                                                                        axis.title  = element_blank(),
                                                                        #legend.title = element_blank(),
                                                                        legend.background = element_rect(fill = "white", color = "black") 
  )



beech.hemlock.hotspots <- base.map +
  geom_sf(data = polygon.buffer.list[["eastern hemlock"]], aes(fill = "eastern hemlock"), alpha = 0.75, color = "black") +
  geom_sf(data = polygon.buffer.list[["American beech"]], aes(fill = "American beech"), alpha = 0.75,color = "black") +
  geom_sf(data = polygon.buffer.list[["eastern white pine"]], aes(fill = "eastern white pine"), alpha = 0.75,color = "black") +
  species_fill+coord_sf(xlim = c(-85, -67.5), ylim = c(37, 47.5))+theme(panel.grid = element_blank(), #panel.background = element_rect(fill = 'lightblue'), 
                                                                        #legend.position = "none",
                                                                        axis.title  = element_blank(),
                                                                        #legend.title = element_blank(),
                                                                        legend.background = element_rect(fill = "white", color = "black") 
  )

beech.hemlock.hotspots



maple.hotspots <- base.map +
  geom_sf(data = polygon.buffer.list[["sugar maple"]], aes(fill = "sugar maple"), alpha = 0.5, color = "black") +
  geom_sf(data = polygon.buffer.list[["red maple"]], aes(fill = "red maple"), alpha = 0.5,color = "black") +
 # geom_sf(data = polygon.buffer.list[["eastern white pine"]], aes(fill = "eastern white pine"), alpha = 0.5,color = "black") +
  species_fill+coord_sf(xlim = c(-85, -67.5), ylim = c(37, 47.5))+theme(panel.grid = element_blank(), #panel.background = element_rect(fill = 'lightblue'), 
                                                                        #legend.position = "none",
                                                                        axis.title  = element_blank(),
                                                                        #legend.title = element_blank(),
                                                                        legend.background = element_rect(fill = "white", color = "black") 
  )

maple.hotspots

other.hotspots <- base.map +
  geom_sf(data = polygon.buffer.list[["black cherry"]], aes(fill = "black cherry"), alpha = 0.7,color = "black") +
  geom_sf(data = polygon.buffer.list[["yellow-poplar"]], aes(fill = "yellow-poplar"), alpha = 0.7,color = "black") +
  geom_sf(data = polygon.buffer.list[["white ash"]], aes(fill = "white ash"), alpha = 0.7,color = "black") +
  geom_sf(data = polygon.buffer.list[["hickory spp."]], aes(fill = "hickory spp."), alpha = 0.7,color = "black") +
  # geom_sf(data = polygon.buffer.list[["eastern white pine"]], aes(fill = "eastern white pine"), alpha = 0.5,color = "black") +
  species_fill+coord_sf(xlim = c(-85, -67.5), ylim = c(37, 47.5))+theme(panel.grid = element_blank(), #panel.background = element_rect(fill = 'lightblue'), 
                                                                        #legend.position = "none",
                                                                        axis.title  = element_blank(),
                                                                        #legend.title = element_blank(),
                                                                        legend.background = element_rect(fill = "white", color = "black") 
  )

other.hotspots

spruce.fir.hothulls <- base.map +
  geom_sf(data = polygon_hull.list[["balsam fir"]], aes(fill = "balsam fir"), alpha = 0.9, color = NA) +
  geom_sf(data = polygon_hull.list[["red spruce"]], aes(fill = "red spruce"), alpha = 0.9, color = NA) +
  geom_sf(data = polygon_hull.list[["northern white-cedar"]], aes(fill = "northern white-cedar"), alpha = 0.5, color = NA) +
  
  #geom_sf(data = polygon.buffer.list[["eastern hemlock"]], aes(fill = "eastern hemlock"), alpha = 0.5, color = NA) +
  #geom_sf(data = polygon.buffer.list[["American beech"]], aes(fill = "American beech"), alpha = 0.5,color = NA) +
  
  
  
  species_fill+coord_sf(xlim = c(-85, -67.5), ylim = c(37, 47.5))+theme(panel.grid = element_blank(), #panel.background = element_rect(fill = 'lightblue'), 
                                                                        #legend.position = "none",
                                                                        axis.title  = element_blank(),
                                                                        #legend.title = element_blank(),
                                                                        legend.background = element_rect(fill = "white", color = "black") 
  )



# geom_boxplot()# to do this, we need to handle some of the predicted trees have a volfac == 0
all.remeas %>% select(state, volfac)%>% distinct()%>% group_by(state) %>% summarise(nunique_volfac = n())

# states 9, 23, 24, 33, and 50 may have fewer plot designs and volfacs are either zero, sawtimber, poletimber, etc
all.remeas %>% select(state, volfac, dbhcur) %>% filter(state %in% c(9, 23, 24, 33, 50))%>% distinct()%>% 
  group_by(state, volfac) %>% 
  summarise(min_dbh = min(dbhcur))

all.remeas %>% filter(dbhcur>5)%>% select(state, PLOT.ID, point, date, cndtn) %>% 
  distinct()%>% group_by(state, PLOT.ID, date, cndtn)%>% summarise(npoints = n())%>% 
  ungroup()%>% group_by(state, date, cndtn) %>% summarise(max_npoints = max(npoints), min_npoints = min(npoints))

all.remeas %>% select(state, volfac, dbhcur) %>% filter(!state %in% c(9, 23, 24, 33, 50))%>% distinct()%>% 
  group_by(state, volfac) %>% 
  summarise(min_dbh = min(dbhcur), 
            dbhcur/volfac)
# formula for TPA
#TPA = (BAF / 0.005454*DIA^2)/Npoints
#volfac*Npoints = BAF/0.005454*DIA^2
(volfac*Npoints)*(0.005454*DIA^2)

# the trees with volfac == 0 are in mostly fixed radius plot states (9, 23, 33, 50), except for 50
all.remeas %>% filter(volfac == 0 & status %in% c(2,4,5,6))%>% group_by(state, mortfac > 0) %>% summarise(n())

# if status == 1 and volfac == 0, assume it is not counted

BAFestimated <- left_join(all.remeas, all.remeas %>% select(state, PLOT.ID, point) %>% distinct()%>%
  group_by(state, PLOT.ID) %>% summarise(Npoints = n())) %>%
  ungroup()%>%
  #filter(Npoints >1) %>% 
  mutate(BAFest = (volfac*Npoints)*(0.005454*dbhcur^2))                                           

# state 36, new york == 10 point design, 37.5 BAF?
# state 39, ohio == 5 point design, 18.7 BAF?
# state 42, Pennsylvannia == 5 point design, 18.7 BAF?
# state 54, West Virginia == 10 point design, 37.5BAF or 5 point design, 18.7 BAF?

ggplot(BAFestimated, aes( y = BAFest, fill = state, group = state))+geom_histogram()+facet_wrap(~state)+ylim(0,50)
ggplot(BAFestimated, aes( y = Npoints, fill = state, group = state))+geom_histogram()+facet_wrap(~state)

ggplot(BAFestimated %>% filter(state == 54), aes( y = BAFest, fill = state, group = state))+geom_histogram()+facet_wrap(~state)
ggplot(BAFestimated %>% filter(state == 54), aes( y = BAFest, fill = state, group = state))+geom_histogram()+facet_wrap(~Npoints >5)

# state 34, new jersey is a mix of fixed radius and variable radius designs
ggplot(BAFestimated %>% filter(state == 34), aes( y = BAFest, fill = state, group = state))+geom_histogram()+facet_wrap(~state)
ggplot(BAFestimated %>% filter(state == 34), aes( y = BAFest, fill = state, group = state))+geom_histogram()+facet_wrap(~Npoints >5)
ggplot(BAFestimated %>% filter(state == 34), aes( y = volfac, fill = state, group = state))+geom_histogram()+facet_wrap(~Npoints >1)

BAFestimated %>% filter(volfac == 0)%>% group_by(state) %>% summarise(n())
# fixed radius volfacs--
# volfac == 10.0196, for trees < 11.0 in
# volfac == 4.9928, for trees >= 11 inches
# all of these are single point plots:
BAFestimated %>% filter(state == 34) %>% filter(volfac == 10.0196 | volfac == 4.9928) %>% 
  group_by(volfac, point)%>% summarise(maxdbh = max(dbhcur), 
                                mindbh = min(dbhcur))

BAFestimated %>% filter(state == 34) %>% filter(!volfac == 10.0196 & !volfac == 4.9928) #%>% 
  
ggplot(BAFestimated %>% filter(state == 34 & !volfac == 10.0196 & !volfac == 4.9928), aes( y = BAFest, fill = state, group = state))+geom_histogram()+facet_wrap(~Npoints >5)
ggplot(BAFestimated %>% filter(state == 34 & !volfac == 10.0196 & !volfac == 4.9928), aes( y = Npoints, fill = state, group = state))+geom_histogram()+facet_wrap(~Npoints >5)


BAFestimated %>% filter(Npoints > 10)%>% View()

# states 34, 36, 39, 42, and 54 likely have variable radius plot designs
all.remeas %>% group_by(state,volfac == 0) %>% summarise(n = n())

# total trees in each species on the stand and in the plots

  
group_by(SPCD, Species, state, remper, PLOT.ID, LONG_FIADB, LAT_FIADB) %>%
  
  summarise(n_trees_volfac = sum(volfac, na.rm = TRUE), 
            n_count = n())




total.spp.predicted.plt <- spp.tree.pred.mort.volfac %>% rename("Species" = "COMMON")%>% 
  group_by(SPCD, Species, state, PLOT.ID, LAT_FIADB, LONG_FIADB)%>%
  summarise(n_predicted_plot = n(), 
            volfac_predicted_plot = sum(volfac, na.rm =TRUE))

spp.predicted.dead.plot <- spp.tree.pred.mort.volfac %>% rename("Species" = "COMMON")%>% 
  mutate(predicted.status = ifelse(survival.draw == 1, "1", "2"))%>%
  filter(predicted.status == "2")%>%
  group_by(SPCD, Species, state, PLOT.ID, LAT_FIADB, LONG_FIADB)%>%
  
  summarise(predicted_dead_trees_volfac = sum(volfac, na.rm = TRUE), 
            n_predicted_dead = n())%>%
  left_join(., total.spp.predicted.plt )



# get the plot-level predicted mortality rates by each species
predicted.mort.rates.spp.plot <- spp.predicted.dead.plot %>% 
  group_by(SPCD, Species, state, PLOT.ID, LAT_FIADB, LONG_FIADB)%>%
  # for each remper period, get the mortality rate per year 
  summarise(predicted_10yrmort_rate_volfac = ((predicted_dead_trees_volfac/volfac_predicted_plot)), 
            predicted_10yrmort_rate = ((n_predicted_dead/n_predicted_plot)))%>%
  ungroup()%>%
  # 10 year survival rates
  mutate(predicted_10yr_surv_rate_volfac = 1-predicted_10yrmort_rate_volfac, 
         predicted_10yr_surv_rate = 1-predicted_10yrmort_rate)%>%
  # get annualized predicted survival and mortality rates
  mutate(predicted_1year_surv_rate_volfac = predicted_10yr_surv_rate_volfac^(1/10),
         predicted_1year_surv_rate = predicted_10yr_surv_rate^(1/10))%>%
  mutate(predicted_1year_mort_rate_volfac = 1 - predicted_1year_surv_rate_volfac, 
         predicted_1year_mort_rate = 1 - predicted_1year_surv_rate)

hist(predicted.mort.rates.spp.plot$predicted_1year_mort_rate_volfac)

  
# plot up 10 year survival rates by plot, by species:
ggplot(data = predicted.mort.rates.spp.plot %>% filter(Species %in% c("balsam fir", "red spruce", "northern white-cedar")))+
  geom_point(aes(x =LONG_FIADB,  y = LAT_FIADB, color = predicted_10yrmort_rate*100))+
  facet_wrap(~Species)+ scale_color_viridis_c(option = "magma")


ggplot(data = predicted.mort.rates.spp.plot %>% filter(Species %in% c("balsam fir", "red spruce", "northern white-cedar")))+
  geom_density(aes(y = predicted_1year_surv_rate*100, fill = Species))+
  facet_wrap(~Species)#+ylab("predicted_1year_mort_rate_volfac")

ggplot(data = predicted.mort.rates.spp.plot %>% filter(Species %in% c("balsam fir", "red spruce", "northern white-cedar")))+
  geom_density(aes(y = predicted_10yrmort_rate*100, fill = Species))+
  facet_wrap(~Species)+ylab("Predicted 10-year mortality rates")


ggplot(data = predicted.mort.rates.spp.plot)+
  geom_density(aes(y = predicted_10yrmort_rate*100, fill = Species))+
  facet_wrap(~Species)+ylab("Predicted 10-year mortality rates")+
  species_fill




ggplot(data = predicted.mort.rates.spp.plot)+
  geom_density(aes(y = predicted_10yrmort_rate_volfac*100, fill = Species))+
  facet_wrap(~Species)+ylab("Predicted 10-year mortality rates (% of stems/species/plot)")+
  species_fill

predicted.mort.rates.spp.plot %>% filter(Species %in% c("balsam fir", "red spruce", "northern white-cedar"))%>%
  ggplot(aes(x = Species, y = predicted_10yrmort_rate_volfac, fill = Species)) +
  geom_dots(layout = "weave") +
  stat_slabinterval()



ggplot(data = predicted.mort.rates.spp.plot %>% filter(Species %in% c("balsam fir", "red spruce", "northern white-cedar")))+
  geom_point(aes(x =LONG_FIADB,  y = LAT_FIADB, color = predicted_10yrmort_rate*100))+
  facet_wrap(~Species)+ scale_color_viridis_c(option = "magma")

# get state-level mortality rate estimates for the modelled estimates---
# sum up volfac by state
# sum up number of plots by state

st.mort.rate.species <- all.remeas %>% group_by(SPCD, Species, state, remper)%>% filter(status %in% c(2, 4, 5, 6)) %>%
  summarise(nonlog_dead_trees_volfac = sum(volfac, na.rm = TRUE), 
            n_nonlog_dead = n()) %>%
  left_join(., st.total.tree.sums) %>%# join to total trees
  
  # for each remper period, get the mortality rate per year 
  mutate(nonlog_mort_rate_volfac = ((nonlog_dead_trees_volfac/total_trees_volfac)*100)/remper, 
         nonlog_mort_rate = ((n_nonlog_dead/total_n)*100)/remper)%>%
  ungroup()%>%
  
  # average all the remper mortality rates together to get a single species value
  group_by(SPCD, Species, state)%>%
  summarise(species_volfac_mort = mean(nonlog_mort_rate_volfac, na.rm =TRUE), 
            species_n_mort = mean(nonlog_mort_rate, na.rm =TRUE))


# get the total trees predicted for each species: 
# predictions were made on standard 10 year intervals
total.spp.predicted <- spp.tree.pred.mort.volfac %>% rename("Species" = "COMMON")%>% 
  group_by(SPCD, Species, state)%>%
  summarise(n_predicted_total = n(), 
            volfac_predicted_total = sum(volfac, na.rm =TRUE))

spp.predicted.dead <- spp.tree.pred.mort.volfac %>% rename("Species" = "COMMON")%>% 
  mutate(predicted.status = ifelse(survival.draw == 1, "1", "2"))%>%
  filter(predicted.status == "2")%>%
  group_by(SPCD, Species, state)%>%
  
  summarise(predicted_dead_trees_volfac = sum(volfac, na.rm = TRUE), 
            n_predicted_dead = n())%>%
  left_join(., total.spp.predicted )
  

predicted.mort.rates.spp.state <- spp.predicted.dead %>% group_by(SPCD, Species, state)%>%
  # for each remper period, get the mortality rate per year 
  summarise(predicted_mort_rate_volfac = ((predicted_dead_trees_volfac/volfac_predicted_total)*100)/10, 
         predicted_mort_rate = ((n_predicted_dead/n_predicted_total)*100)/10)%>%
  ungroup()


# combine with the observed mortality rates by species:
p.o.mort.rates <- left_join(st.mort.rate.species, predicted.mort.rates.spp.state)

ggplot(data = p.o.mort.rates, aes(x = species_volfac_mort, y = predicted_mort_rate_volfac, color = Species))+
  geom_point()+
  species_color+geom_abline(aes(slope = 1, intercept = 0), linetype = "dashed")+
  xlim(0, 6)+ylim(0,6)

# looking at distibution of mortality rates within the context of disturbances----

PLOT <- read_delim(paste0(output.dir,"data/formatted_older_matching_plts_PLOT.txt"))
colnames(PLOT)
unique(PLOT$typcur)
TREE.typcd <- TREE.remeas %>% left_join(., PLOT %>% select(PLOT.ID, cndtn, typcur, typold, date) %>% distinct())

# northern white cedar mortality is in N white cedar forest type (14), and spruce fir types
# red spruce and balsam fir (13), balsam fir(11), and red spruce (19)

TREE.typcd %>% filter(Species %in% "northern white-cedar" & dbhold > 5)%>%
  filter(status == c(1, 2, 4, 5)) %>% 
  mutate(M = ifelse(status == 1, 0, 1))%>%
  group_by( typcur, M)%>%
  summarise(ntree = n()) %>% 
  spread(M, ntree)%>% arrange(desc(`1`)) 

# get plots with at least 1 dead nwc tree:
Trees.maine <- TREE.typcd %>% filter(dbhold > 5)%>%
  filter(status %in% c(1, 2,4, 5)) %>% 
  mutate(M = ifelse(status == 1, 0, 1))%>%
  group_by( PLOT.ID, stname, Species, typcur, M)%>%
  summarise(ntree = n()) %>% 
  ungroup() %>% filter(stname %in% c("ME", "VT", "NH"))

Dead.NWC.plots <- Trees.maine %>% filter(Species %in% "northern white-cedar" ) %>% 
  group_by(PLOT.ID, M)%>%
  summarise(totals = sum(ntree)) %>% 
  spread(M, totals) %>% filter(!is.na(`1`))

noDead.NWC.plots <- Trees.maine %>% filter(Species %in% "northern white-cedar"  )%>% 
  group_by(PLOT.ID, M)%>%
  summarise(totals = sum(ntree)) %>% 
  spread(M, totals) %>% filter(is.na(`1`) | `1`==0)

# get the mortality fraction for other species for plots with dead N. white cedar
Trees.maine %>% filter(PLOT.ID %in% Dead.NWC.plots$PLOT.ID) %>% 
  mutate(plot.type = "Dead white-cedar") %>% rbind(., 
                                                   Trees.maine %>% filter(PLOT.ID %in% noDead.NWC.plots$PLOT.ID) %>% 
                                                     mutate(plot.type = "No Dead white-cedar")) %>% 
  filter(M == 1)%>%
  ggplot()+geom_bar(aes(x = plot.type, y = ntree, fill = Species), position = "stack", stat = "identity")+facet_wrap(~typcur)

TREE.typcd %>% filter(Species %in% "northern white-cedar" & dbhold > 5 & !status == 3) %>% 
  mutate(M = ifelse(status == 1, "live", "dead"))%>%
  ggplot()+geom_boxplot(aes(x = as.factor(typcur), y = dbhcur,  fill = M), alpha = 0.75, position = "dodge")

TREE.typcd %>% filter(Species %in% "northern white-cedar" & dbhold > 5 & !status == 3) %>% 
  mutate(M = ifelse(status == 1, "live", "dead"))%>%
  group_by(M, damage)%>%
  summarise(n())

TREE.typcd %>% filter(Species %in% "northern white-cedar" & dbhold > 5 & !status == 3)%>%
  #filter(damage == 50)%>%
  group_by(PLOT.ID) %>% 
  ggplot()+geom_jitter(aes(x = LONG_FIADB,y = LAT_FIADB, color = as.factor(damage)))


TREE.typcd %>% filter(Species %in% "balsam fir")%>%
  filter(status %in% c(1, 2, 4, 5)) %>% 
  mutate(M = ifelse(status == 1, 0, 1))%>%
  group_by( typcur, M)%>%
  summarise(ntree = n()) %>% 
  spread(M, ntree)%>% arrange(desc(`1`)) 



TREE.typcd %>%filter(Species %in% "balsam fir") %>%
  filter(status %in% c(1, 2, 4, 5)) %>% 
  mutate(M = ifelse(status == 1, 0, 1))%>%
  group_by( typcur, stname, M)%>%
  summarise(ntree = n()) %>% 
  spread(M, ntree)%>% arrange(desc(`1`)) 


# most of the Northern white cedar tree mortality is in Maine
TREE.typcd %>% filter(Species %in% "northern white-cedar" & dbhold > 5)%>%
  filter(status == c(1, 2,4, 5)) %>% 
  mutate(M = ifelse(status == 1, 0, 1))%>%
  group_by( stname, M)%>%
  summarise(ntree = n()) %>% 
  spread(M, ntree)%>% arrange(desc(`1`)) 


TREE.typcd %>% filter(Species %in% "northern white-cedar"& dbhold > 5)%>%
  filter(status == c(1, 2, 4, 5)) %>% 
  mutate(M = ifelse(status == 1, 0, 1))%>%
  group_by(PLOT.ID, M)%>%
  summarise(ntree = n()) %>% 
  spread(M, ntree) %>%
  mutate(total.NWC = sum(c(`0`,`1`), na.rm =TRUE)) %>% 
  arrange(desc(`1`))




# get summary of damages for later use:
N.DAMAGE <- cleaned.data %>% group_by(SPCD, damage) %>% summarise(n.by.damage = n())
N.DAMAGE$SPECIES <- ref_species[match(N.DAMAGE$SPCD, ref_species$SPCD),]$COMMON
ref_damage<- ref_codes %>% filter(VARIABLE %in% "AGENTCD")
N.DAMAGE$damage_agent <- ref_damage[match(N.DAMAGE$damage, ref_damage$VALUE),]$MEANING
N.DAMAGE$damage_agent <- ifelse(N.DAMAGE$damage == 0, "None", N.DAMAGE$damage_agent)





maine_var <- readRDS(paste0(output.folder, "/SPCD_stanoutput_joint_v3/predicted_mort/","variance_explained_state_23.RDS"))






color.pred.class.2 <- c(
  
  "Size" = "darkgreen",
  "Change in Size" = "#66a61e",
  "Growth x Size" = "#a1d99b" ,
  
  "Climate" = "darkblue",
  "Climate x G & S" = "#1d91c0",
  
  "N deposition" = "#67000d" , 
  "Ndep x G & S"= "#bd0026",
  
  "% Damage" = "#762a83", 
  "Damage x G & S" = "#9970ab",
  
  
  "Competition" = "#8c510a",
  "Competition x G & S" = "#d95f02" ,
  
  "Site Conditions" = "#e6ab02",
  "Site x G & S" = "yellow3")

var.part.fill <- scale_fill_manual(values = color.pred.class.2, name = "",
                                   labels = c("Change in Size" = expression(Delta ~ "D" ), 
                                              "Growth x Size" = expression(Delta ~ "D x D"), 
                                              "Size" = "Diameter (D)", 
                                              "Climate x G & S" = expression("Climate x " ~Delta ~ "D or D"), 
                                              "Ndep x G & S" = expression("N dep. x " ~Delta ~ "D or D"),
                                              "Damage x G & S" = expression("Damage x " ~Delta ~ "D or D"), 
                                              "Competition x G & S" = expression("Compeition x " ~Delta ~ "D or D"),
                                              "Site x G & S" = expression("Site x " ~Delta ~ "D or D")
                                   ))

var.part.color <- scale_color_manual(values = color.pred.class.2, name = "",
                                     labels = c("Change in Size" = expression(Delta ~ "D" ), 
                                                "Growth x Size" = expression(Delta ~ "D x D"), 
                                                "Size" = "Diameter (D)", 
                                                "Climate x G & S" = expression("Climate x " ~Delta ~ "D or D"), 
                                                "Ndep x G & S" = expression("N dep. x " ~Delta ~ "D or D"),
                                                "Damage x G & S" = expression("Damage x " ~Delta ~ "D or D"), 
                                                "Competition x G & S" = expression("Compeition x " ~Delta ~ "D or D"),
                                                "Site x G & S" = expression("Site x " ~Delta ~ "D or D")
                                     ))

# centralized code for the manuscript figures:
input.folder <- "C:/Users/KellyHeilman/Box/01. kelly.heilman Workspace/mortality/Eastern-Mortality/mortality_models/SPCD_stanoutput_joint_v3/"
output.folder <- "C:/Users/KellyHeilman/Box/01. kelly.heilman Workspace/mortality/Eastern-Mortality/mortality_models/SPCD_stanoutput_joint_v3/"

# re-invisioning the betas effects plots ----
# get all the covariates using posterior package
model.no <- 6
mod.data <- readRDS(paste0(input.folder, "all_SPCD_model_",model.no,".RDS"))
ncovar <- length(colnames(mod.data$xM))

full.model <- data.frame(Covariates = colnames(mod.data$xM), 
                         id = 1:length(colnames(mod.data$xM)))

# read in the betas, marginal responses, and variance parsing summaries ---
betas.df <- readRDS(paste0(input.folder, "SPCD_stanoutput_joint_v3/samples/u_betas_model_", model.no, "_5000samples.rds"))
marginal_response_df <- readRDS(paste0(input.folder, "SPCD_stanoutput_joint_v3/all_marginal_responses.RDS"))
interaction_response_df <- readRDS(paste0(input.folder, "SPCD_stanoutput_joint_v3/interaction_responses.RDS"))
var_summary <- read.csv(paste0(output.folder, "SPCD_stanoutput_joint_v3/predicted_mort/regional_var/variance_partitioning_summary_by_predictor_all.csv"))


betas.quant <- betas.df %>% summarise_draws(median, ~quantile(., probs = c(0.025, 0.975))) %>%
  rename(`ci.lo` = "2.5%", `ci.hi` = "97.5%") %>%
  mutate(remper.cor = 0.5)
# relabel u_betas to meaningful species ids names
betas.quant$spp <- rep(1:17, ncovar)
betas.quant$cov <- rep(1:ncovar, each = 17)


covariate_names <- c(colnames(mod.data$xM))  # Replace with your covariate names
betas.quant$Covariate <- rep(covariate_names, each = 17)
betas.quant$Species <- rep(nspp[1:17,]$COMMON, ncovar)


# reorder by the value of the covariate
#betas.quant <- betas.quant %>% arrange(by = median) 
# betas.quant$parameter <- factor(betas.quant$parameter, levels = betas.quant$parameter)

# get overlapping zero to color the error bars
betas.quant$`significance` <- ifelse(betas.quant$ci.lo < 0 & betas.quant$ci.hi < 0, "significant", 
                                     ifelse(betas.quant$ci.lo > 0 & betas.quant$ci.hi > 0, "significant", "not overlapping zero"))



betas.quant$Covariate <- factor(betas.quant$Covariate, levels = unique(betas.quant$Covariate))
# order species by hardwood softwood, then shade tolence
betas.quant$Species <- factor(betas.quant$Species, levels = SP.TRAITS$COMMON_NAME)


ggplot(data = na.omit(betas.quant), aes(x = Species, y = median, color = significance))+geom_point()+
  geom_errorbar(data = na.omit(betas.quant), aes(x = Species , ymin = ci.lo, ymax = ci.hi, color = significance), width = 0.1)+
  geom_abline(aes(slope = 0, intercept = 0), color = "grey", linetype = "dashed")+
  facet_wrap(~Covariate, scales= "free_y")+
  theme_bw(base_size = 14)+
  theme( axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.grid  = element_blank(), legend.position = "none")+
  ylab("Effect on Survival")+xlab("Parameter")+
  scale_color_manual(values = c("not overlapping zero"="darkgrey", "significant"="black"))
#ggsave(height = 10, width = 15,dpi = 350, units = "in",paste0(output.folder,"SPCD_stanoutput_joint_v3/images/Estimated_effects_on_mortality_model_model_",model.no,"_all_species_betas.png"))



betas.quant$Species <- factor(betas.quant$Species, levels = rev(disturb.species.order))

betas.sig.df <- betas.quant %>% filter(significance %in% "significant") %>%
  left_join(., unique(var_summary[,c("Covariate", "Predictor","predictor.class2")]))

plot_beta_effects <- function(pred.group, betas.sig){
  
  df.covar <- betas.sig %>% filter(predictor.class2 %in%  pred.group)
  strip.fill <- as.character(color.pred.class.2[unique(df.covar$predictor.class2)])
  
  ggplot()+
    geom_errorbar(data = na.omit(df.covar), aes(y = Predictor , xmin = ci.lo, xmax = ci.hi, color = Species, linetype = significance), position = position_dodge(width = 0.9), width = 0)+
    geom_point(data = na.omit(df.covar), aes(y = Predictor, x = median, color = Species, shape = significance), position = position_dodge(width = 0.9), size = 1.5)+
    
    
    facet_wrap(~predictor.class2, scales= "free_y", ncol = 1)+
    geom_vline(xintercept = 0, color = "grey", linetype = "dashed")+
    theme_bw(base_size = 14)+
    theme(axis.text.y = element_text(angle = 90, hjust = 0.5),
          panel.grid  = element_blank(), 
          strip.background = element_rect(fill = c(strip.fill)), 
          strip.text = element_text(color = "white"), 
          axis.title.y = element_blank()
    )+
    xlab("Effect on Survival")+ylab("Parameter")+
    species_color+
    scale_shape_manual(values = c("significant" = 19,  "not overlapping zero" = 1), drop = FALSE, name = "Significance")+
    scale_linetype_manual(values = c("significant" = "solid",  "not overlapping zero" = "dotted"), 
                          drop = FALSE,  name = "Significance")+
    guides(color = guide_legend(reverse = TRUE))
  
  
}


plot_beta_effects(pred.group = "Change in Size", betas.sig = betas.sig.df)
plot_beta_effects(pred.group = "Competition", betas.sig = betas.sig.df)

plot_beta_effects(pred.group = "Competition", betas.sig = betas.sig.df)

beta.sig.plots <- lapply(unique(betas.sig.df$predictor.class2), function(x){
  plot_beta_effects(pred.group = x, betas.sig = betas.sig.df)
})

beta.sig.plots


# plot for all the betas, not just significant ones:
betas.all <- betas.quant %>% 
  left_join(., unique(var_summary[,c("Covariate", "Predictor","predictor.class2")]))

beta.all.plots <- lapply(unique(betas.all$predictor.class2), function(x){
  plot_beta_effects(pred.group = x, betas.sig = betas.all)
})


beta.all.plots[[2]]

# plots for the betas that have at least one significant species effect:
betas.sig.all <- betas.quant %>% filter(Covariate %in% unique(betas.sig.df$Covariate))%>% 
  left_join(., unique(var_summary[,c("Covariate", "Predictor","predictor.class2")]))

# group by variables
# list.vars <- list(c("Change in Size", "Size", "Growth x Size"),
# c("% Damage", "Competition", "Competition x G & S"),
# c("Climate", "Climate x G & S"),
# c("N deposition", "Ndep x G & S", "Site Conditions"))

beta.sig.all.plots <- lapply(unique(betas.sig.all$predictor.class2), function(x){
  plot_beta_effects(pred.group = x, betas.sig = betas.sig.all)
})
library(patchwork)
# set up the layouts using patchwork
growth.plt <- (beta.sig.all.plots[[1]] / beta.sig.all.plots[[2]]/ beta.sig.all.plots[[8]] & coord_cartesian(xlim=c(-1, 3.567)))+ plot_layout(ncol = 1, guides = "collect", axes = "collect_x") 

compete.plt <- ( (beta.sig.all.plots[[3]]) / beta.sig.all.plots[[12]]/beta.sig.all.plots[[4]] /beta.sig.all.plots[[9]]  & coord_cartesian(xlim=c(-0.7, 0.5))) + 
  plot_layout(ncol = 1, guides = "collect", axes = "collect_x", heights = c(2, 1,1,1)) 

climate.plt <- ( beta.sig.all.plots[[5]] / beta.sig.all.plots[[10]] & coord_cartesian(xlim=c(-1, 1))) + plot_layout(ncol = 1, guides = "collect", axes = "collect_x") 

site.ndep.plt <- (  beta.sig.all.plots[[7]]/ beta.sig.all.plots[[11]]/beta.sig.all.plots[[6]]   & coord_cartesian(xlim=c(-0.5, 0.7))) + plot_layout(ncol = 1, guides = "collect", axes = "collect_x") 

# save all the plots
save_plot(paste0(output.folder,"images/significant_betas_growth_dia.png"), 
          growth.plt, base_width = 5, base_height = 8) 

save_plot(paste0(output.folder,"images/significant_betas_growth_dia.svg"),
          growth.plt, base_width = 5, base_height = 8)

save_plot(paste0(output.folder,"images/significant_betas_competition.png"), 
          compete.plt, base_width = 5, base_height = 10) 

save_plot(paste0(output.folder,"images/significant_betas_competition.svg"),
          compete.plt, base_width = 5, base_height = 10)

save_plot(paste0(output.folder,"images/significant_betas_climate.png"), 
          climate.plt, base_width = 5, base_height = 10) 

save_plot(paste0(output.folder,"images/significant_betas_climate.svg"),
          climate.plt, base_width = 5, base_height = 10)


save_plot(paste0(output.folder,"images/significant_betas_siteNdep.png"), 
          site.ndep.plt, base_width = 5, base_height = 8) 

save_plot(paste0(output.folder,"images/significant_betas_siteNdep.svg"),
          site.ndep.plt, base_width = 5, base_height = 8)


#########################################################################
# SVG of variance partitions

# read in regional variance partitioning:
var_summary <- readRDS(paste0(output.dir, "SPCD_stanoutput_joint_v3/predicted_mort/var_summary_regional.rds"))

color.pred.class <- c(
  "Size" = "darkgreen",
  "Change in Size" = "#66a61e",
  "Site x G & S"= "#fdbf6f",
  "Competition x G & S"="#d95f02" ,
  "Climate x G & S"="#7570b3" ,
  "Climate"= "#e7298a" ,
  "Growth x Size" =  "#1b9e77",
  "Site Conditions" = "#e6ab02",
  "Competition" = "#a6761d"
)

color.pred.class.2 <- c(
  
  "Size" = "darkgreen",
  "Change in Size" = "#66a61e",
  "Growth x Size" = "#a1d99b" ,
  
  "Climate" = "#081d58",
  "Climate x G & S" = "#1d91c0",
  
  "N deposition" = "#67000d" , 
  "Ndep x G & S"= "#bd0026",
  
  "% Damage" = "#762a83", 
  "Damage x G & S" = "#9970ab",
  
  
  "Competition" = "#8c510a",
  "Competition x G & S" = "#d95f02" ,
  
  "Site Conditions" = "#e6ab02",
  "Site x G & S" = "yellow3")



var.summary.region <- var_summary
var.summary.region$COMMON <- factor(var.summary.region$COMMON, 
                                    levels = rev(c("balsam fir", "red spruce", "northern white-cedar", 
                                                   "eastern hemlock", "American beech", 
                                                   "black oak", "chestnut oak", "northern red oak", "white oak", "yellow birch", "paper birch", 
                                                   "hickory spp.", "eastern white pine", "red maple", "sugar maple", 
                                                   "black cherry", "white ash", "yellow-poplar")))
var.summary.region$predictor.class <- factor(var.summary.region$predictor.class, 
                                             levels = c(
                                               
                                               "Size" ,
                                               "Change in Size" ,
                                               "Climate" ,
                                               "Site Conditions",
                                               "Competition",
                                               "N deposition", 
                                               "% Damage", 
                                               "Growth x Size",
                                               
                                               "Site x G & S",
                                               "Competition x G & S" ,
                                               
                                               "Ndep x G & S",
                                               "Climate x G & S",
                                               "Damage x G & S"))


var.summary.region$predictor.class2 <- factor(var.summary.region$predictor.class2, 
                                              levels = c(
                                                
                                                "Size",
                                                "Change in Size",
                                                "Growth x Size" ,
                                                
                                                "Climate",
                                                "Climate x G & S",
                                                
                                                "N deposition", 
                                                "Ndep x G & S" ,
                                                
                                                "% Damage", 
                                                "Damage x G & S",
                                                
                                                
                                                "Competition",
                                                "Competition x G & S",
                                                
                                                "Site Conditions",
                                                "Site x G & S"
                                                
                                                
                                              ))




ggplot(data = var.summary.region)+
  geom_bar(aes(x = COMMON, y = mean, fill = predictor.class2), stat = "identity", position = "stack")+
  theme_bw(base_size = 16)+ theme(axis.text.x = element_text(angle = 60, hjust = 1), 
                                  panel.background = element_blank())+
  ylab("Across-tree variance in  \n p(survival) explained")+
  scale_fill_manual(values = color.pred.class.2, name = "",
                    labels = c("Change in Size" = expression(Delta ~ "D" ), 
                               "Growth x Size" = expression(Delta ~ "D x D"), 
                               "Size" = "Diameter (D)", 
                               "Climate x G & S" = expression("Climate x " ~Delta ~ "D or D"), 
                               "Ndep x G & S" = expression("N dep. x " ~Delta ~ "D or D"),
                               "Damage x G & S" = expression("Damage x " ~Delta ~ "D or D"), 
                               "Competition x G & S" = expression("Compeition x " ~Delta ~ "D or D"),
                               "Site x G & S" = expression("Site x " ~Delta ~ "D or D")
                    ),
                    guide = guide_legend(direction = "horizontal",
                                         ncol = 14,nrow = 1,reverse = TRUE,
                                         label.position="top", label.hjust = 0,
                                         label.vjust = 0.5,
                                         label.theme = element_text(angle = 90)))+
  coord_flip()+
  xlab("")+
  theme(panel.grid = element_blank(), 
        legend.position = "top", 
        panel.border = element_blank(), 
        axis.ticks.y = element_blank())


plt.regional.plt.all <- ggplot(data = var.summary.region, aes(x = COMMON, y = mean_logit, fill = predictor.class2))+
  #geom_bar(aes(x = COMMON, y = mean_logit, fill = predictor.class2), stat = "identity", position = "stack", width = 1)+
  stat_summary(geom = "col", fun = sum, width = 0.7,
               position = "stack") +
  theme_bw(base_size = 16)+ theme(axis.text.x = element_text(angle = 60, hjust = 1), 
                                  panel.background = element_blank())+
  ylab("Across-tree variance in  \n p(survival) explained")+
  scale_fill_manual(values = color.pred.class.2, name = "",
                    labels = c("Change in Size" = expression(Delta ~ "D" ), 
                               "Growth x Size" = expression(Delta ~ "D x D"), 
                               "Size" = "Diameter (D)", 
                               "Climate x G & S" = expression("Climate x " ~Delta ~ "D or D"), 
                               "Ndep x G & S" = expression("N dep. x " ~Delta ~ "D or D"),
                               "Damage x G & S" = expression("Damage x " ~Delta ~ "D or D"), 
                               "Competition x G & S" = expression("Compeition x " ~Delta ~ "D or D"),
                               "Site x G & S" = expression("Site x " ~Delta ~ "D or D")
                    ),
                    guide = guide_legend(direction = "horizontal",
                                         ncol = 14,nrow = 1,reverse = TRUE,
                                         label.position="top", label.hjust = 0,
                                         label.vjust = 0.5,
                                         label.theme = element_text(angle = 90)))+
  coord_flip()+
  xlab("")+
  theme(panel.grid = element_blank(), 
        legend.position = "top", 
        panel.border = element_blank(), 
        axis.ticks.y = element_blank())

ggsave(plot = plt.regional.plt.all, height = 8, width = 7.5, units = "in", dpi = 500, 
       paste0(output.dir, "images/regional_species_var_partitioning.png"))

ggsave(plot = plt.regional.plt.all, height = 8, width = 7.5, units = "in", dpi = 500, 
       paste0(output.dir, "images/regional_species_var_partitioning.svg"))


plt.regional.plt.all2 <- ggplot(data = var.summary.region, aes(x = COMMON, y = mean, fill = predictor.class2))+
  #geom_bar(aes(x = COMMON, y = mean_logit, fill = predictor.class2), stat = "identity", position = "stack", width = 1)+
  stat_summary(geom = "col", fun = sum, width = 0.7,
               position = "stack") +
  theme_bw(base_size = 16)+ theme(axis.text.x = element_text(angle = 60, hjust = 1), 
                                  panel.background = element_blank())+
  ylab("Across-tree variance in  \n p(survival) explained")+
  scale_fill_manual(values = color.pred.class.2, name = "",
                    labels = c("Change in Size" = expression(Delta ~ "D" ), 
                               "Growth x Size" = expression(Delta ~ "D x D"), 
                               "Size" = "Diameter (D)", 
                               "Climate x G & S" = expression("Climate x " ~Delta ~ "D or D"), 
                               "Ndep x G & S" = expression("N dep. x " ~Delta ~ "D or D"),
                               "Damage x G & S" = expression("Damage x " ~Delta ~ "D or D"), 
                               "Competition x G & S" = expression("Compeition x " ~Delta ~ "D or D"),
                               "Site x G & S" = expression("Site x " ~Delta ~ "D or D")
                    ),
                    guide = guide_legend(direction = "horizontal",
                                         ncol = 14,nrow = 1,reverse = TRUE,
                                         label.position="top", label.hjust = 0,
                                         label.vjust = 0.5,
                                         label.theme = element_text(angle = 90)))+
  coord_flip()+
  xlab("")+
  theme(panel.grid = element_blank(), 
        legend.position = "top", 
        panel.border = element_blank(), 
        axis.ticks.y = element_blank())

ggsave(plot = plt.regional.plt.all2, height = 8, width = 7.5, units = "in", dpi = 500, 
       paste0(output.dir, "images/regional_species_var_partitioning_meanpsurv.png"))

ggsave(plot = plt.regional.plt.all2, height = 8, width = 7.5, units = "in", dpi = 500, 
       paste0(output.dir, "images/regional_species_var_partitioning_meanpsurv.svg"))


# some summaries for putting numbers in the paper
var.summary.region %>% group_by(COMMON) %>% arrange(COMMON, desc(mean_logit)) %>% 
  filter(predictor.class2 %in% c("Size", "Change in Size", "Growth x Size"))%>%
  group_by(COMMON)%>%
  summarise(sum.logit = sum(mean_logit))

var.summary.region %>% group_by(COMMON) %>% arrange(COMMON, desc(mean_logit)) %>% 
  filter(!predictor.class2 %in% c("Size", "Change in Size", "Growth x Size"))%>%
  group_by(COMMON) %>% 
  group_by(predictor.class2, COMMON)%>%
  summarise(total_cat = sum(mean_logit)*100)%>%
  select(total_cat,  predictor.class2, COMMON) %>%
  
  filter(total_cat > 2.5) %>% 
  arrange( desc(total_cat)) %>%
  group_by(COMMON)%>% View()


var.summary.region %>% group_by(COMMON) %>% arrange(COMMON, desc(mean_logit)) %>% 
  filter(!predictor.class2 %in% c("Size", "Change in Size", "Growth x Size"))%>%
  mutate(logit_percent = round(mean_logit*100, digits = 2))%>%
  group_by(predictor.class2, COMMON)%>%
  mutate(total_cat = sum(mean_logit)*100)%>%
  select(Covariate, logit_percent, total_cat,  predictor.class2, COMMON) %>%
  filter(predictor.class2 %in%  "Competition x G & S")%>%
  #filter(total_cat > 5) %>% 
  arrange(desc(total_cat), desc(logit_percent)) %>% 
  filter(total_cat >1)%>% View()

#########################################################################
# Marginal effects, but plotting effects for species in affected by different pests



marginal_response_df$predictor.class2
#interaction_response_df<- left_join(interaction_response_df, Covariate.types.df)

spruce.fir <- c("balsam fir", "red spruce",  "northern white-cedar")
mixed <- c("American beech", "eastern hemlock")
spongy.susceptible<- c("chestnut oak", "white oak", "northern red oak", "paper birch", "yellow birch")
spongy.immune <- c("white ash", "yellow-poplar", "black cherry")
spongy.resist <- c("sugar maple", "red maple", "eastern white pine", "hickory spp.")





# set up groups of species to look at:
spruce.fir <- c("balsam fir", "red spruce",  "northern white-cedar")
mixed <- c("American beech", "eastern hemlock")
spongy.susceptible<- c("chestnut oak", "white oak",  "northern red oak", "paper birch", "yellow birch")
spongy.immune <- c("white ash", "yellow-poplar", "black cherry")
spongy.resist <- c("sugar maple", "red maple", "eastern white pine", "hickory spp.")

interaction.terms <- unique(interaction_response_df$Covariate)
all_marginal_response_df <- marginal_response_df

# set up species colors:
SP.TRAITS <- read.csv("data/NinemetsSpeciesTraits.csv") %>% filter(COMMON_NAME %in% unique(nspp[1:17,]$COMMON))
# order the trait db by softwood-hardwood, then shade tolerance, then name (this puts all the oaks together b/c hickory and red oak have the same tolerance values)
SP.TRAITS <- SP.TRAITS %>% group_by(SFTWD_HRDWD) %>% arrange(desc(SFTWD_HRDWD), desc(ShadeTol), desc(COMMON_NAME))

SP.TRAITS$Color <- c(# softwoods
  "#b2df8a",
  "#003c30", 
  "#b2182b", 
  "#fee090", 
  "#33a02c",
  
  
  # sugar  maples
  "#a6cee3",
  "#1f78b4",
  
  # red maple
  "#e31a1c",
  # yellow birch
  "#fdbf6f",
  # oaks
  "#cab2d6",
  "#8073ac",
  "#6a3d9a",
  
  # hickory
  "#7f3b08",
  # white ash
  "#bababa",
  # black cherry
  "#4d4d4d",
  # yellow poplar
  "#ff7f00",
  "#fccde5" # paper birch
  
  
)

SP.TRAITS$`Shade Tolerance`  <- ifelse(SP.TRAITS$ShadeTol >=4, "High", 
                                       ifelse(SP.TRAITS$ShadeTol <=2.5, "Low", "Moderate"))

# set up custom colors for species
sppColors <- SP.TRAITS$Color 
names(sppColors) <- unique(SP.TRAITS$COMMON_NAME) 

species_fill <- scale_fill_manual(values = sppColors)
species_color <- scale_color_manual(values = sppColors)



marginal_response_df

# covariates that are scaled by species:
#Ndep_Diff_per_yr
species.scaling <- train.data %>% select(Species, SPCD, Ndep_Diff_per_yr.median, Ndep_Diff_per_yr.sd) %>% distinct()%>%
  mutate(Covariate = "Ndep.scaled", 
         Clean_Name = "delta N deposition",
         Units = "kg/m^2/year") %>% 
  rename("Val.mean" = "Ndep_Diff_per_yr.median",
         "Val.sd" = "Ndep_Diff_per_yr.sd")%>%
  #DIA_DIFF_scaled
  rbind(., train.data %>% select(Species, SPCD, DIA.DIFF.median, DIA.DIFF.sd) %>% distinct()%>%
          mutate(Covariate = "DIA_DIFF_scaled", 
                 Clean_Name = "Diameter Difference",
                 Units = "Inches")%>% 
          rename("Val.mean" = "DIA.DIFF.median",
                 "Val.sd" = "DIA.DIFF.sd"))%>%
  #BAL
  rbind(.,train.data %>% select(Species, SPCD, BAL.median, BAL.sd) %>% distinct()%>%
          mutate(Covariate = "BAL.scaled",
                 Clean_Name = "Basal Area Larger Than",
                 Units = "ft^2")%>% 
          rename("Val.mean" = "BAL.median",
                 "Val.sd" = "BAL.sd")
  )%>%
  #DIA
  rbind(., train.data %>% select(Species, SPCD, DIA.median, DIA.sd) %>% distinct()%>%
          mutate(Covariate = "DIA_scaled", 
                 Clean_Name = "Diameter",
                 Units = "Inches")%>% 
          rename("Val.mean" = "DIA.median",
                 "Val.sd" = "DIA.sd") 
  )%>%
  #plt_ba_sq_ft_old
  rbind(., train.data %>% select(Species, SPCD, plt_ba_sq_ft_old.median, plt_ba_sq_ft_old.sd) %>% distinct()%>%
          mutate(Covariate = "ba.scaled", 
                 Clean_Name = "Plot Basal Area",
                 Units = "ft^2")%>% 
          rename("Val.mean" = "plt_ba_sq_ft_old.median",
                 "Val.sd" = "plt_ba_sq_ft_old.sd")
  )%>%
  select(Species, SPCD, Covariate, Val.mean, Val.sd, Clean_Name, Units)

# covariates scaled across space:
plot.medians <- readRDS("data/plot.medians_SPCD_all.rds")
plot.scaled <- plot.medians %>% select(damage.median, damage.sd) %>% 
  mutate(Covariate = "damage.scaled",
         Clean_Name = "Damage",
         Units = "% of stems")%>%
  rename("Val.mean" = "damage.median",
         "Val.sd" = "damage.sd")%>%
  
  #damage.scaled = (damage.total - plot.medians$damage.median)/plot.medians$damage.sd,
  #MAP.scaled = (MAP-plot.medians$MAP.median)/plot.medians$MAP.sd,
  rbind(., plot.medians %>% select(MAP.median, MAP.sd) %>% 
          mutate(Covariate = "MAP.scaled", 
                 Clean_Name = "Mean Annual Precipitation",
                 Units = "mm")%>%
          rename("Val.mean" = "MAP.median",
                 "Val.sd" = "MAP.sd"))%>%
  #MATmax.scaled = (MATmax - plot.medians$MATmax.median)/plot.medians$MATmax.sd)
  rbind(.,plot.medians %>% select(MATmax.median, MATmax.sd) %>% 
          mutate(Covariate = "MATmax.scaled", 
                 Clean_Name = "Average Annual Maximum Temperature",
                 Units = "Degrees C")%>%
          rename("Val.mean" = "MATmax.median",
                 "Val.sd" = "MATmax.sd"))%>%
  
  
  rbind(.,plot.medians %>% select(slope.median, slope.sd) %>% 
          mutate(Covariate = "slope.scaled", 
                 Clean_Name = "Slope",
                 Units = "%")%>%
          rename("Val.mean" = "slope.median",
                 "Val.sd" = "slope.sd"))%>%
  expand_grid(., Species = unique(species.scaling$Species))%>%
  left_join(., unique(species.scaling[,c("Species", "SPCD")]))%>%
  select(Species, SPCD, Covariate, Val.mean, Val.sd, Clean_Name, Units)

# data.frame with the median and sd for each species to rescale the variable by:
combined.scaled.main <- rbind(plot.scaled, species.scaling)

marginal_response_df_unscaled <- marginal_response_df %>% left_join(., combined.scaled.main)%>%
  mutate(Raw.value = ifelse(!is.na(Val.mean), (Value*Val.sd) + Val.mean, Value))


#Ndep X dia_diff
#p1.value   p2.value p2.rank             covariate     Pred.1          Pred.2
#1  -1.12091712  0.1224049 

interaction_response_df_unscaled <- interaction_response_df %>% 
  # join the first precdictor with the raw values of the covariate
  left_join(., combined.scaled.main %>% rename("Pred.1" = "Covariate"))%>%
  mutate(Raw.value = ifelse(!is.na(Val.mean), (p1.value*Val.sd) + Val.mean, p1.value)) %>%
  left_join(., combined.scaled.main %>% rename("Pred.2" = "Covariate", 
                                               "Val.mean.2" = "Val.mean", 
                                               "Val.sd.2" = "Val.sd", 
                                               "Clean_Name.2" = "Clean_Name", 
                                               "Units.2" = "Units"))%>%
  mutate(Raw.value.p2 = ifelse(!is.na(Val.mean.2), (p2.value*Val.sd.2) + Val.mean.2, p2.value))

interaction_response_df_unscaled %>% select(Species, p2.rank, Pred.2, Raw.value.p2) %>% distinct()
# need to join the interaction plots by the median and sd values and convert p1 and p2 values
# function to plot up the effects on 10 year mortality probabilities
# using the unscaled "raw" values for covariates
plot_main_region_effects <- function(species.group, covar, ymax.spp){
  
  
  
  if(! covar %in% interaction.terms){  
    
    df.species <- marginal_response_df_unscaled %>% filter(Species %in%  species.group & Covariate %in% covar)
    strip.fill <- as.character(color.pred.class.2[unique(df.species$predictor.class2)])
    
    #Ndepx growth= -2.23104221
    # Ndep = 
    
    p1 <-  ggplot(df.species , aes(x = Raw.value, y = 1-(mean)^10, color = Species)) +
      geom_line(size = 1) +
      geom_ribbon(aes(ymin = 1-(ci.lo)^10, ymax = 1-(ci.hi)^10, fill = Species), alpha = 0.2, color = NA) +
      facet_wrap(~Predictor) +
      labs(
        x = paste(df.species$predictor.class2, "(", unique(df.species$Units),")"),
        y = "10-year Mortality Probability",
        #title = "Effect of Predictors on Probability of Mortality"
      ) +
      species_color + species_fill+
      #var.part.fill + var.part.color+
      theme_bw(base_size = 14)+theme(panel.grid = element_blank(), 
                                     legend.key.size = unit(1, "cm"), 
                                     legend.position = "none", 
                                     strip.background = element_rect(fill = strip.fill), 
                                     strip.text = element_text(color = "white")
      )+coord_cartesian(ylim = c(0, ymax.spp))
    
    
  }else{
    
    df.species <- interaction_response_df_unscaled %>% filter(Species %in%  species.group & Covariate %in% covar)
    
    
    df.species$p2.rank <- factor(df.species$p2.rank, levels = c("high", "median", "low"))
    strip.fill <- as.character(color.pred.class.2[unique(df.species$predictor.class2)])
    pred.name <-  marginal_response_df %>% filter(Covariate %in% unique(df.species$Pred.1))%>% select(Predictor)%>% distinct()
    
    if(unique(df.species$Pred.2) %in% "DIA_DIFF_scaled"){
      
      p1 <-  ggplot(data = df.species) +
        geom_line(aes(x = Raw.value, y = 1-(mean)^10, color = Species, group = interaction(Species, p2.rank), linetype = p2.rank), size = 1) +
        geom_ribbon(aes(x = Raw.value, ymin = 1-(ci.lo.10)^10, ymax = 1-(ci.hi.90)^10, fill = Species, group = interaction(Species, p2.rank)), color = NA, alpha = 0.25) +
        #facet_grid(cols = vars(predictor.class2), rows = vars(Predictor)) +
        facet_wrap(~Predictor)+
        labs(
          x = paste(pred.name$Predictor, "(", unique(df.species$Units),")"),
          y = "10-year Mortality Probability",
          #title = "Effect of Predictors on Probability of Mortality"
        ) + scale_linetype_manual(values = c("low"= "dotted" ,
                                             "median" = "dashed",
                                             "high"= "solid" ), 
                                  name = expression(Delta ~ "Diameter"))+
        species_color + species_fill+
        # scale_fill_manual(values = c("low"="#a1dab4" ,
        #                              "median" = "#41b6c4",
        #                              "high"= "#225ea8" ), 
        #                   name = expression(Delta ~ "Diameter"))+
        
        #species_fill + species_color + 
        #named_species_linetype +
        theme_bw(base_size = 14)+theme(panel.grid = element_blank(), 
                                       legend.key.size = unit(1, "cm"), 
                                       strip.background = element_rect(fill = strip.fill)
        )+coord_cartesian(ylim = c(0, ymax.spp))
      
      
    }else{
      
      p1 <-  ggplot(data = df.species) +
        geom_line(aes(x = Raw.value, y = 1-(mean)^10, color = Species, group = interaction(Species, p2.rank), linetype = p2.rank, size = p2.rank)) +
        geom_ribbon(aes(x = Raw.value, ymin = 1-(ci.lo.10)^10, ymax = 1-(ci.hi.90)^10, fill = Species, group = interaction(Species, p2.rank)), color = NA, alpha = 0.25) +
        #facet_grid(cols = vars(predictor.class2), rows = vars(Predictor)) +
        facet_wrap(~Predictor)+
        labs(
          x = paste(pred.name$Predictor, "(", unique(df.species$Units),")"),
          y = "10-year Mortality Probability",
          #title = "Effect of Predictors on Probability of Mortality"
        ) + 
        species_color + species_fill+
        scale_size_manual(values = c("low"=0.25 ,
                                     "median" = 1.2,
                                     "high"= 2.2 ), 
                          name = "Diameter", 
                          labels = c("low" = "small", 
                                     "median" = "medium", 
                                     "high" = "large"))+
        scale_linetype_manual(values = c("low"="solid" ,
                                         "median" = "dashed",
                                         "high"= "longdash" ), 
                              name = "Diameter", 
                              labels = c("low" = "small", 
                                         "median" = "medium", 
                                         "high" = "large"))+
        
        
        theme_bw(base_size = 14)+theme(panel.grid = element_blank(), 
                                       legend.key.size = unit(1, "cm"), 
                                       strip.background = element_rect(fill = strip.fill)
        )+coord_cartesian(ylim = c(0, ymax.spp))
      
      
      
      
      
      
    }
  }
  return(p1)
}


# Plot up the important marginal effects for different groups of species---
# select main effects and interactions to plot for each species group:
# read in the variance partitioning summary by species
#var_summary <- readRDS(paste0(output.folder, "variance_partitioning_summary_region_species.RDS"))




plot_main_region_effects(species.group = c("northern red oak"), 
                         covar = "MAP.scaled_DIA.int", 
                         ymax.spp = 0.01)

plot_main_region_effects(species.group = c("northern red oak", "hickory spp."), 
                         covar = c("tmax.anom_DIA.int"),
                         ymax.spp = 0.02)

plot_main_region_effects(species.group = c("chestnut oak"), 
                         covar = c("Ndep.scaled"),
                         ymax.spp = 0.1)

plot_main_region_effects(species.group = c("red maple", "yellow birch"), 
                         covar = c("MATmax.scaled_growth.int"),
                         ymax.spp = 0.25)

plot_main_region_effects(species.group = c("balsam fir"), 
                         covar = c("BAL.scaled_DIA.int"),
                         ymax.spp = 0.25)

plot_main_region_effects(species.group = c("balsam fir", "yellow-poplar", "yellow birch"), 
                         covar = c("ba.scaled_growth.int"),
                         ymax.spp = 0.5)


plot_main_region_effects(species.group = c( "paper birch"), 
                         covar = c("damage.scaled_DIA.int"),
                         ymax.spp = 0.05)

plot_main_region_effects(species.group = c("paper birch"), 
                         covar = c("damage.scaled_growth.int"),
                         ymax.spp = 0.25)

plot_main_region_effects(species.group = c("eastern white pine"), 
                         covar = c("Ndep.scaled_growth.int"),
                         ymax.spp = 0.09)
# spruce fir:----

# get all of the covariates that contribute to 1% or greater of the species
sprucefir.top5 <- var_summary %>% filter(Species %in% spruce.fir) %>% ungroup()%>%
  group_by( Covariate, predictor.class2, Species)%>%
  summarise(pct.region = median(mean_logit)*100)%>% arrange(Species, desc(pct.region)) %>% 
  group_by(Species)%>% #View()
  filter(pct.region >= 1 ) %>% ungroup() #%>% select(Covariate) %>% distinct()

# D X deltaD
# deltaD
# D
# Ndep x deltaD
# damage x deltaD
# damage
# matmax x D
# ba x D

region.plt.list <-list()
for(h in 1:length(unique(sprucefir.top5$Covariate))){
  region.plt.list[[h]] <- plot_main_region_effects(species.group = spruce.fir, 
                                                   covar = unique(sprucefir.top5$Covariate)[h], 
                                                   ymax.spp = 0.5)
}



# make separate legends
dia.diff.legend <- get_legend(region.plt.list[[3]]+scale_color_discrete(guide = "none")+scale_fill_discrete(guide = "none")+theme_bw(base_size = 18)+ theme(legend.key.width = unit(2,"cm")))
#diameter.legend <- get_legend(region.plt.list[[4]]+scale_color_discrete(guide = "none")+scale_fill_discrete(guide = "none")+theme_bw(base_size = 18)+ theme(legend.key.width = unit(2,"cm")))
species.legend <- get_legend(region.plt.list[[1]]+scale_linetype_discrete(guide = "none")+scale_size_discrete(guide = "none")+theme_bw(base_size = 18)+ theme(legend.key.width = unit(2,"cm")))

spruce.fir.effects <-plot_grid(
  plot_grid(region.plt.list[[1]]+theme(legend.position = "none"),
            region.plt.list[[2]]+theme(legend.position = "none"),
            region.plt.list[[5]]+theme(legend.position = "none"),
            
            
            ncol = 3, align = "hv"),
  
  plot_grid(region.plt.list[[3]]+theme(legend.position = "none"),
            region.plt.list[[4]]+theme(legend.position = "none"),
            region.plt.list[[6]]+theme(legend.position = "none"),
            ncol = 3, align = "hv"),
  plot_grid(species.legend, dia.diff.legend, ncol = 2),
  rel_heights = c(1,1,0.3),
  nrow = 4)

save_plot(paste0(output.folder,"images/spruce.fir_regional_marginal_effects_1pct.png"), 
          spruce.fir.effects, base_width = 16, base_height = 15) 

save_plot(paste0(output.folder,"images/spruce.fir_regional_marginal_effects_1pct.svg"),
          spruce.fir.effects, base_width = 16, base_height = 15)




# spongy moth susceptible:----
spongy.susceptible.top5 <- var_summary %>% filter(Species %in% spongy.susceptible) %>% ungroup()%>%
  group_by( Covariate, predictor.class2, Species)%>%
  summarise(pct.region = median(mean_logit)*100)%>% arrange(Species, desc(pct.region)) %>% 
  group_by(Species)%>% #View()
  filter(pct.region >= 1 ) %>% ungroup() #%>% select(Covariate) %>% distinct()

unique(spongy.susceptible.top5$Covariate)

region.plt.list <-list()
for(h in 1:length(unique(spongy.susceptible.top5$Covariate))){
  region.plt.list[[h]] <- plot_main_region_effects(species.group = spongy.susceptible, 
                                                   covar = unique(spongy.susceptible.top5$Covariate)[h], 
                                                   ymax.spp = 0.05)
}



# make separate legends
dia.diff.legend <- get_legend(region.plt.list[[6]]+scale_color_discrete(guide = "none")+scale_fill_discrete(guide = "none")+theme_bw(base_size = 18)+ theme(legend.key.width = unit(2,"cm")))
diameter.legend <- get_legend(region.plt.list[[8]]+scale_color_discrete(guide = "none")+scale_fill_discrete(guide = "none")+theme_bw(base_size = 18)+ theme(legend.key.width = unit(2,"cm")))
species.legend <- get_legend(region.plt.list[[1]]+scale_linetype_discrete(guide = "none")+scale_size_discrete(guide = "none")+theme_bw(base_size = 18)+ theme(legend.key.width = unit(2,"cm")))

spongy.susceptible.effects <-plot_grid(
  plot_grid(region.plt.list[[2]]+theme(legend.position = "none"),
            region.plt.list[[1]]+theme(legend.position = "none"),
            region.plt.list[[6]]+theme(legend.position = "none"),
            region.plt.list[[5]]+theme(legend.position = "none"),
            region.plt.list[[9]]+theme(legend.position = "none"),
            ncol = 6, align = "hv"),
  
  plot_grid(
    region.plt.list[[7]]+theme(legend.position = "none"),
    region.plt.list[[3]]+theme(legend.position = "none"),
    region.plt.list[[4]]+theme(legend.position = "none"),
    region.plt.list[[10]]+theme(legend.position = "none"),
    region.plt.list[[11]]+theme(legend.position = "none"),
    region.plt.list[[8]]+theme(legend.position = "none"),
    
    ncol = 6, align = "hv"),
  plot_grid(species.legend, diameter.legend,  dia.diff.legend, ncol = 2),
  rel_heights = c(1,1,0.3),
  nrow = 4)

save_plot(paste0(output.folder,"images/spongy.susceptible_dominant_marginal_variance_effects.png"), 
          spongy.susceptible.effects, base_width = 16, base_height = 15) 

save_plot(paste0(output.folder,"images/spongy.susceptible_dominant_marginal_variance_effects.svg"),
          spongy.susceptible.effects, base_width = 16, base_height = 15)




# spongy moth resistant:----

spongy.resist.top5 <- var_summary %>% filter(Species %in% spongy.resist) %>% ungroup()%>%
  group_by( Covariate, predictor.class2, Species)%>%
  summarise(pct.region = median(mean_logit)*100)%>% arrange(Species, desc(pct.region)) %>% 
  group_by(Species)%>% #View()
  filter(pct.region >= 1 ) %>% ungroup()

unique(spongy.resist.top5$Covariate)

region.plt.list <-list()
for(h in 1:length(unique(spongy.resist.top5$Covariate))){
  region.plt.list[[h]] <- plot_main_region_effects(species.group = spongy.resist, 
                                                   covar = unique(spongy.resist.top5$Covariate)[h], 
                                                   ymax.spp = 0.05)
}



# make separate legends
dia.diff.legend <- get_legend(region.plt.list[[4]]+scale_color_discrete(guide = "none")+scale_fill_discrete(guide = "none")+theme_bw(base_size = 18)+ theme(legend.key.width = unit(2,"cm")))
diameter.legend <- get_legend(region.plt.list[[5]]+scale_color_discrete(guide = "none")+scale_fill_discrete(guide = "none")+theme_bw(base_size = 18)+ theme(legend.key.width = unit(2,"cm")))
species.legend <- get_legend(region.plt.list[[1]]+scale_linetype_discrete(guide = "none")+scale_size_discrete(guide = "none")+theme_bw(base_size = 18)+ theme(legend.key.width = unit(2,"cm")))

spongy.resist.effects <-plot_grid(
  plot_grid(region.plt.list[[2]]+theme(legend.position = "none"),
            region.plt.list[[1]]+theme(legend.position = "none"),
            region.plt.list[[6]]+theme(legend.position = "none"),
            region.plt.list[[4]]+theme(legend.position = "none"),
            
            
            ncol = 4, align = "hv"),
  
  plot_grid(
    
    region.plt.list[[3]]+theme(legend.position = "none"),
    region.plt.list[[7]]+theme(legend.position = "none"),
    region.plt.list[[5]]+theme(legend.position = "none"),
    
    
    ncol = 4, align = "hv"),
  plot_grid(species.legend, diameter.legend, dia.diff.legend, ncol = 3),
  rel_heights = c(1,1,0.3),
  nrow = 4)

save_plot(paste0(output.folder,"images/spongy.resist_dominant_marginal_variance_effects.png"), 
          spongy.resist.effects, base_width = 16, base_height = 15) 

save_plot(paste0(output.folder,"images/spongy.resist_dominant_marginal_variance_effects.svg"),
          spongy.resist.effects, base_width = 16, base_height = 15)

# spongy moth immune:----

spongy.immune.top5 <- var_summary %>% filter(Species %in% spongy.immune) %>% ungroup()%>%
  group_by( Covariate, predictor.class2, Species)%>%
  summarise(pct.region = median(mean_logit)*100)%>% arrange(Species, desc(pct.region)) %>% 
  group_by(Species)%>% #View()
  filter(pct.region >= 1 ) %>% ungroup()

unique(spongy.immune.top5$Covariate)

region.plt.list <-list()
for(h in 1:length(unique(spongy.immune.top5$Covariate))){
  region.plt.list[[h]] <- plot_main_region_effects(species.group = spongy.immune, 
                                                   covar = unique(spongy.immune.top5$Covariate)[h], 
                                                   ymax.spp = 0.05)
}



# make separate legends
dia.diff.legend <- get_legend(region.plt.list[[2]]+scale_color_discrete(guide = "none")+scale_fill_discrete(guide = "none")+theme_bw(base_size = 18)+ theme(legend.key.width = unit(2,"cm")))
#diameter.legend <- get_legend(region.plt.list[[10]]+scale_color_discrete(guide = "none")+scale_fill_discrete(guide = "none")+theme_bw(base_size = 18)+ theme(legend.key.width = unit(2,"cm")))
species.legend <- get_legend(region.plt.list[[1]]+scale_linetype_discrete(guide = "none")+scale_size_discrete(guide = "none")+theme_bw(base_size = 18)+ theme(legend.key.width = unit(2,"cm")))

spongy.immune.effects <-plot_grid(
  plot_grid(region.plt.list[[3]]+theme(legend.position = "none"),
            region.plt.list[[1]]+theme(legend.position = "none"),
            region.plt.list[[2]]+theme(legend.position = "none"),
            region.plt.list[[4]]+theme(legend.position = "none"),
            
            
            ncol = 4, align = "hv"),
  
  
  
  
  plot_grid(species.legend,  dia.diff.legend, ncol = 3),
  rel_heights = c(1,0.3),
  nrow = 2)

save_plot(paste0(output.folder,"images/spongy.immune_dominant_marginal_variance_effects.png"), 
          spongy.immune.effects, base_width = 16, base_height = 9) 

save_plot(paste0(output.folder,"images/spongy.immune_dominant_marginal_variance_effects.svg"),
          spongy.immune.effects, base_width = 16, base_height = 9)

# hemlock, beech, and :----
mixed.top5 <- var_summary %>% filter(Species %in% mixed) %>% ungroup()%>%
  group_by( Covariate, predictor.class2, Species)%>%
  summarise(pct.region = median(mean_logit)*100)%>% arrange(Species, desc(pct.region)) %>% 
  group_by(Species)%>% #View()
  filter(pct.region >= 1 ) %>% ungroup() #%>% select(Covariate) %>% distinct()


region.plt.list <-list()
for(h in 1:length(unique(mixed.top5$Covariate))){
  region.plt.list[[h]] <- plot_main_region_effects(species.group = mixed, 
                                                   covar = unique(mixed.top5$Covariate)[h], 
                                                   ymax.spp = 0.15)
}


region.plt.list[[3]]

# make separate legends
#dia.diff.legend <- get_legend(region.plt.list[[8]]+scale_color_discrete(guide = "none")+scale_fill_discrete(guide = "none")+theme_bw(base_size = 18)+ theme(legend.key.width = unit(2,"cm")))
#diameter.legend <- get_legend(region.plt.list[[9]]+scale_color_discrete(guide = "none")+scale_fill_discrete(guide = "none")+theme_bw(base_size = 18)+ theme(legend.key.width = unit(2,"cm")))
species.legend <- get_legend(region.plt.list[[1]]+scale_linetype_discrete(guide = "none")+scale_size_discrete(guide = "none")+theme_bw(base_size = 18)+ theme(legend.key.width = unit(2,"cm")))

mixed.effects <- plot_grid(
  plot_grid(region.plt.list[[3]]+theme(legend.position = "none"),
            region.plt.list[[1]]+theme(legend.position = "none"),
            region.plt.list[[2]]+theme(legend.position = "none"),
            region.plt.list[[4]]+theme(legend.position = "none"),
            
            region.plt.list[[5]]+theme(legend.position = "none"),
            
            
            ncol = 5, align = "hv"),
  
  plot_grid(species.legend, ncol = 1),
  rel_heights = c(1,0.3),
  nrow = 2)

save_plot(paste0(output.folder,"images/mixed_dominant_marginal_variance_effects.png"), 
          mixed.effects, base_width = 16, base_height = 9) 

save_plot(paste0(output.folder,"images/mixed_dominant_marginal_variance_effects.svg"),
          mixed.effects, base_width = 16, base_height = 9)





## Select variance partitioning across the region:---------------------
# read in the state-level variance partitioning:
st.var.files <- list.files(paste0(output.folder,"SPCD_stanoutput_joint_v3/predicted_mort/"), pattern = "variance_partitioning_summary_by_predictor_state_", full.names = TRUE)
last_three_str_sub <- str_sub(st.var.files, -3, -1)
Covariate_names <- read.csv(paste0(output.dir, "data/model_covariate_types_v2.csv"))

st.var.all.df <- do.call(rbind, lapply(st.var.files, function(x){
  read.csv(x)%>% mutate(state = str_sub(x, -6, -5)) %>% rename("Covariate" = "predictor")%>% left_join(., Covariate_names)
}))

state.variance.df <- st.var.all.df %>% left_join(.,state.df %>% mutate(state = as.character(state)))
state.variance.df
#var.summary.region <- var_summary
state.variance.df$COMMON <- factor(state.variance.df$Species, 
                                   levels = rev(c("balsam fir", "red spruce", "northern white-cedar", 
                                                  "eastern hemlock", "American beech", 
                                                  "black oak", "chestnut oak", "northern red oak", "white oak", "yellow birch", "paper birch", 
                                                  "hickory spp.", "eastern white pine", "red maple", "sugar maple", 
                                                  "black cherry", "white ash", "yellow-poplar")))
state.variance.df$predictor.class <- factor(state.variance.df$predictor.class, 
                                            levels = c(
                                              
                                              "Size" ,
                                              "Change in Size" ,
                                              "Climate" ,
                                              "Site Conditions",
                                              "Competition",
                                              "N deposition", 
                                              "% Damage", 
                                              "Growth x Size",
                                              
                                              "Site x G & S",
                                              "Competition x G & S" ,
                                              
                                              "Ndep x G & S",
                                              "Climate x G & S",
                                              "Damage x G & S"))


state.variance.df$predictor.class2 <- factor(state.variance.df$predictor.class2, 
                                             levels = c(
                                               
                                               "Size",
                                               "Change in Size",
                                               "Growth x Size" ,
                                               
                                               "Climate",
                                               "Climate x G & S",
                                               
                                               "N deposition", 
                                               "Ndep x G & S" ,
                                               
                                               "% Damage", 
                                               "Damage x G & S",
                                               
                                               
                                               "Competition",
                                               "Competition x G & S",
                                               
                                               "Site Conditions",
                                               "Site x G & S"
                                               
                                               
                                             ))



spruce.fir.state <- ggplot(data = state.variance.df %>% filter(COMMON %in% spruce.fir))+
  geom_bar(aes(x = region, y = mean, fill = predictor.class2), stat = "identity", position = "stack")+
  theme_bw(base_size = 16)+ theme(axis.text.x = element_text(angle = 60, hjust = 1), 
                                  panel.background = element_blank())+
  ylab("Across-tree variance in  \n p(survival) explained")+
  scale_fill_manual(values = color.pred.class.2, name = "",
                    labels = c("Change in Size" = expression(Delta ~ "D" ), 
                               "Growth x Size" = expression(Delta ~ "D x D"), 
                               "Size" = "Diameter (D)", 
                               "Climate x G & S" = expression("Climate x " ~Delta ~ "D or D"), 
                               "Ndep x G & S" = expression("N dep. x " ~Delta ~ "D or D"),
                               "Damage x G & S" = expression("Damage x " ~Delta ~ "D or D"), 
                               "Competition x G & S" = expression("Compeition x " ~Delta ~ "D or D"),
                               "Site x G & S" = expression("Site x " ~Delta ~ "D or D")
                    ),
                    guide = guide_legend(direction = "horizontal",
                                         ncol = 14,nrow = 1,reverse = TRUE,
                                         label.position="top", label.hjust = 0,
                                         label.vjust = 0.5,
                                         label.theme = element_text(angle = 90)))+
  coord_flip()+
  xlab("")+
  theme(panel.grid = element_blank(), 
        legend.position = "top", 
        panel.border = element_blank(), 
        axis.ticks.y = element_blank())+
  facet_wrap(~COMMON, scales = "free_y", ncol = 1)

save_plot(paste0(output.folder,"images/spruce_fir_state_variance_effects.png"), 
          spruce.fir.state, base_width = 16, base_height = 15) 

save_plot(paste0(output.folder,"images/spruce_fir_state_variance_effects.svg"),
          spruce.fir.state, base_width = 16, base_height = 15)


spongy.susceptible.state <- ggplot(data = state.variance.df %>% filter(COMMON %in% spongy.susceptible))+
  geom_bar(aes(x = region, y = mean, fill = predictor.class2), stat = "identity", position = "stack")+
  theme_bw(base_size = 16)+ theme(axis.text.x = element_text(angle = 60, hjust = 1), 
                                  panel.background = element_blank())+
  ylab("Across-tree variance in  \n p(survival) explained")+
  scale_fill_manual(values = color.pred.class.2, name = "",
                    labels = c("Change in Size" = expression(Delta ~ "D" ), 
                               "Growth x Size" = expression(Delta ~ "D x D"), 
                               "Size" = "Diameter (D)", 
                               "Climate x G & S" = expression("Climate x " ~Delta ~ "D or D"), 
                               "Ndep x G & S" = expression("N dep. x " ~Delta ~ "D or D"),
                               "Damage x G & S" = expression("Damage x " ~Delta ~ "D or D"), 
                               "Competition x G & S" = expression("Compeition x " ~Delta ~ "D or D"),
                               "Site x G & S" = expression("Site x " ~Delta ~ "D or D")
                    ),
                    guide = guide_legend(direction = "horizontal",
                                         ncol = 14,nrow = 1,reverse = TRUE,
                                         label.position="top", label.hjust = 0,
                                         label.vjust = 0.5,
                                         label.theme = element_text(angle = 90)))+
  coord_flip()+
  xlab("")+
  theme(panel.grid = element_blank(), 
        legend.position = "right", 
        panel.border = element_blank(), 
        axis.ticks.y = element_blank())+
  facet_wrap(~COMMON, scales = "free_y", ncol = 1)

save_plot(paste0(output.folder,"images/spongy_susceptible_state_variance_effects.png"), 
          spongy.susceptible.state , base_width = 16, base_height = 15) 

save_plot(paste0(output.folder,"images/spongy_susceptible_state_variance_effects.svg"),
          spongy.susceptible.state , base_width = 16, base_height = 15)



beech.hemlock.pine.pop.st <-  ggplot(data = state.variance.df %>% filter(COMMON %in% c("American beech", "eastern hemlock", "yellow-poplar", "eastern white pine")))+
  geom_bar(aes(x = region, y = mean, fill = predictor.class2), stat = "identity", position = "stack")+
  theme_bw(base_size = 16)+ theme(axis.text.x = element_text(angle = 60, hjust = 1), 
                                  panel.background = element_blank())+
  ylab("Across-tree variance in  \n p(survival) explained")+
  scale_fill_manual(values = color.pred.class.2, name = "",
                    labels = c("Change in Size" = expression(Delta ~ "D" ), 
                               "Growth x Size" = expression(Delta ~ "D x D"), 
                               "Size" = "Diameter (D)", 
                               "Climate x G & S" = expression("Climate x " ~Delta ~ "D or D"), 
                               "Ndep x G & S" = expression("N dep. x " ~Delta ~ "D or D"),
                               "Damage x G & S" = expression("Damage x " ~Delta ~ "D or D"), 
                               "Competition x G & S" = expression("Compeition x " ~Delta ~ "D or D"),
                               "Site x G & S" = expression("Site x " ~Delta ~ "D or D")
                    ),
                    guide = guide_legend(direction = "horizontal",
                                         ncol = 14,nrow = 1,reverse = TRUE,
                                         label.position="top", label.hjust = 0,
                                         label.vjust = 0.5,
                                         label.theme = element_text(angle = 90)))+
  coord_flip()+
  xlab("")+
  theme(panel.grid = element_blank(), 
        legend.position = "right", 
        panel.border = element_blank(), 
        axis.ticks.y = element_blank())+
  facet_wrap(~COMMON, scales = "free_y", ncol = 1)

save_plot(paste0(output.folder,"images/mixed_other_state_variance_effects.png"), 
          beech.hemlock.pine.pop.st , base_width = 16, base_height = 15) 

save_plot(paste0(output.folder,"images/mixed_other_state_variance_effects.svg"),
          beech.hemlock.pine.pop.st , base_width = 16, base_height = 15)

state.variance.df %>% filter(COMMON %in% spruce.fir)%>%
  mutate(logit_percent = round(mean_logit*100, digits = 2))%>%
  group_by(predictor.class2, COMMON, region, Covariate)%>%
  summarise(total_cat = sum(mean_logit)*100)%>%
  #select(Covariate, logit_percent, total_cat,  predictor.class2, COMMON, region) %>%
  filter(COMMON %in% "red spruce" )%>%
  filter(total_cat > 1) %>% 
  arrange(region, desc(total_cat))%>% View()

# plot up regional values against the mortality rates
train.data %>% filter(Species %in%  "red spruce") %>% 
  group_by(state) %>% 
  summarise(median(MATmax))

train.data.values <- train.data %>% left_join(., state.df)
ggplot(train.data.values %>% filter(Species %in% "red spruce"), aes(x = region, y = MATmax))+geom_boxplot()
ggplot(train.data.values %>% filter(Species %in% "red spruce"), aes(x = region, y = Difference_per_yr))+geom_boxplot()

ggplot(train.data.values %>% filter(Species %in% spruce.fir), aes(x = region, y = damage.total, fill = Species))+geom_boxplot()
ggplot(train.data.values %>% filter(Species %in% spruce.fir), 
       aes(x = region, y = Difference_per_yr, fill = Species))+geom_boxplot()

ggplot(train.data.values %>% filter(Species %in% spruce.fir), 
       aes(x = region, y = MATmax, fill = Species))+geom_boxplot()
ggplot(train.data.values %>% filter(Species %in% spruce.fir), 
       aes(x = region, y = plt_ba_sq_ft_old, fill = Species))+geom_boxplot()

ggplot(train.data.values %>% filter(Species %in% spongy.susceptible), aes(x = region, y = damage.total, fill = Species))+geom_boxplot()


ggplot(train.data.values %>% filter(Species %in% spruce.fir), aes(x = region, y = damage.total, fill = Species))+geom_boxplot()
ggplot(train.data.values %>% filter(Species %in% spongy.susceptible), aes(x = region, y = damage.total, fill = Species))+geom_boxplot()
