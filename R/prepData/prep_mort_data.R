library(rstan)
library(MASS)
library(here)
library(tidyverse)
library(gt)
library(terra)
library(tidyverse)
library(odbc)


boxdir <- "C:/Users/KellyHeilman/Box/01. kelly.heilman Workspace/mortality/Eastern-Mortality/mortality_models/"
# lets start with just remeasured trees from plots where we have tree cores (and therefore lat-longs)
load(paste0(boxdir, "data/Tree Records Master.RData"))
# loads all.trees and radial.inc but we just need all.trees for the PLOT.ID 

# loads internal CN connection that will allow us to get plot lat and lons
CN.connect <- read.csv("C:/Users/KellyHeilman/Box/tree core project/Tree Core plots needing matching to their old Periodic data.csv")

# load the plot-level information and climate data, called plots
load(paste0(boxdir, "data/Tree core plots plus 800m PRISM.RData"))

colnames(plots)

TREE <- read_delim(paste0(boxdir, "data/formatted_older_matching_plts_TREE.txt"))
TREE$dbhcur <- as.numeric(TREE$dbhcur)
TREE$dbhold <- as.numeric(TREE$dbhold)

#TREE <- TREE %>% filter(PLOT.ID %in% unique(plots$PLOT.ID))
TREE$SPCD <- TREE$spp
TREE$Species <- FIESTA::ref_species[match(TREE$SPCD, FIESTA::ref_species$SPCD),]$COMMON

# TREE %>% filter(SPCD %in% nspp[1:15,]$spp & status %in% c(1,2)) %>% group_by(SPCD, status, Species) %>% summarise(n()) %>% 
#  ungroup()%>% mutate(STATUS = ifelse(status == 2, "dead", "live"))%>% dplyr::select(Species, SPCD, `n()`, STATUS)%>%
#   spread(STATUS, `n()`) %>% ungroup() |> gt()


summary(TREE$dbhcur)
summary(TREE$crcls)
unique(TREE$crcls)
unique(TREE$state)
TREE %>% group_by(status) %>% summarise(n())

PLOT <- read_delim(paste0(boxdir,"data/formatted_older_matching_plts_PLOT.txt"))
colnames(PLOT)
PLOT$long
PLOT$lat
unique(PLOT$state)

PLOT %>% filter( state == 10) %>% summarise(n())
# including DE adds 250
#------------------------------------------------------
# get the periodic plot lat longs that match:
#------------------------------------------------------


# get periodic plot information from all the plots:
# 09 CT
# 10 DE
# 23 ME
# 24 MD
# 25 MA
# 33 NH
# 34 NJ
# 36 NY
# 39 OH
# 42 PA
# 44 RI
# 50 VT
# 54 WV

ocon <- dbConnect(odbc(), "fiadb01p")


NE_plot <- dbGetQuery(ocon, "SELECT cn, statecd, unitcd, countycd, plot, invyr, plot_status_cd, cycle, lat, lon, elev, designcd
                      FROM fs_fiadb.plot
                      WHERE statecd =  ANY(09, 10, 23, 24, 25, 33, 34, 36, 39, 42, 44, 50, 54) 
                      and invyr < 2000") %>%
  as_tibble()%>%
  rename_with(tolower)
unique(NE_plot$statecd)
length(unique(NE_plot$cn))
# kindCD or design cd

# RI_tree <- dbGetQuery(ocon, "SELECT cn, statecd, unitcd, countycd, plot, subp, tree, invyr, dia, statuscd, spcd
#                       FROM fs_fiadb.tree
#                       WHERE statecd =  ANY(44) 
#                       and invyr < 2000") %>%
#   as_tibble()%>%
#   rename_with(tolower)
# View(RI_tree %>% filter( countycd == 5, plot == 1))
# 
# RI_tree %>% mutate(st_ut_ct_plt_subp_tree = paste0(statecd,"_", unitcd, "_",countycd, "_", plot, "_", subp, "_", tree)) %>%
#   select(st_ut_ct_plt_subp_tree, invyr, dia, spcd) %>% 
#   # then spread by diameter

unique(NE_plot$statecd)
length(unique(NE_plot$cn))
#length(unique(NE_plot$))
#NE_plot %>% filter(statecd == 25 | statecd == 44) %>% group_by(statecd, invyr, designcd, kindcd) %>% summarise(n())
# kindCD or design cd


# when done, disconnect from ORACLE FIADB
dbDisconnect(ocon)
rm(ocon)

summary(NE_plot$lat)
summary(NE_plot$elev)

# get matching plot ids
NE_plot$PLOT.ID <- as.numeric(paste0(NE_plot$statecd, 
                                       sprintf("%03d", NE_plot$countycd),
                                       sprintf("%04d", NE_plot$plot)))


NE_plot %>% filter(PLOT.ID %in% unique(PLOT$PLOT.ID))

# filter only plots with tree cores
#PLOT <- PLOT %>% filter(PLOT.ID %in% unique(plots$PLOT.ID))

colnames(NE_plot) <- tolower(colnames(NE_plot))
colnames(NE_plot)[9:10] <- c("LAT_FIADB", "LONG_FIADB")
colnames(NE_plot)[12] <- "PLOT.ID"

NE_plot_ll <- left_join(PLOT, NE_plot)
length(NE_plot_ll$LAT_FIADB)
summary(NE_plot_ll$LAT_FIADB)

# 7418 plots out of 39940that are missing LLS


# JOIN up elevation data to this dataset
PRISM <- rast(paste0(boxdir,"data/PRISM_us_dem_800m_asc.asc"))
plot(PRISM)

# make a spatial vector of the lat long points to extract the prism elevation from
spplots <- vect(cbind( c(NE_plot_ll$LONG_FIADB), c(NE_plot_ll$LAT_FIADB)), crs="+proj=longlat +datum=WGS84")

#PRISM 
# now need to extract using LL of the plots
# elevation is in meters
# fast extraction here!
spplots.elev <- extract(PRISM, spplots)
NE_plot_ll$elev <- spplots.elev$PRISM_us_dem_800m_asc

summary(NE_plot_ll$elev)

#unique(PLOT %>% filter(remper >= 20) %>% dplyr::select(state, date))
summary(NE_plot_ll$exprem)

NE_plot_ll %>% group_by(state) %>% summarise(mean (elev, na.rm  = TRUE))


NE_plot_ll %>% group_by(state, exprem > 0) %>% summarise(n())

PLOT.remper <- unique( NE_plot_ll %>% dplyr::select(PLOT.ID, mdate, date, remper, exprem, LAT_FIADB, LONG_FIADB, elev))




#-----------------------------------------------------------------------------
# Join up with tree level data
TREE.remeas <- left_join(TREE, PLOT.remper)

TREE.remeas$DIA_DIFF <- TREE.remeas$dbhcur - TREE.remeas$dbhold
TREE.remeas$remper <- as.numeric(TREE.remeas$remper)
#hist(TREE.remeas$DIA_DIFF)
unique(TREE.remeas$status)
colnames(TREE.remeas)
#--------------------------------------------------------------------------------------
# generate map with 15 species on it
#--------------------------------------------------------------------------------------
library(ggmap)
library(maps)
library(mapdata)

nspp <- TREE.remeas %>% group_by(spp) %>% summarise(n()) %>% arrange(desc(`n()`))
nspp$Species <- FIESTA::ref_species[match(nspp$spp, FIESTA::ref_species$SPCD),]$COMMON_NAME
nspp[1:15,]

#TREE.remeas$LAT <- CN.connect[match(TREE.remeas$PLOT.ID, CN.connect$PLOT.ID),]$LAT_ACTUAL_NAD83
#TREE.remeas$LON <- CN.connect[match(TREE.remeas$PLOT.ID, CN.connect$PLOT.ID),]$LON_ACTUAL_NAD83

unique(TREE.remeas$state)

TREEs.15.spp <- TREE.remeas %>% filter(spp %in% nspp[1:17,]$spp)

occurance.number <- TREEs.15.spp %>% group_by(LAT_FIADB, LONG_FIADB, PLOT.ID, spp)%>% summarise(n=n()) 
occurance.number$Species <- FIESTA::ref_species[match(occurance.number$spp, FIESTA::ref_species$SPCD),]$COMMON_NAME
#ggplot()+geom_point(data = TREEs.15.spp, aes(LON, LAT, color = spp))

ggplot()+geom_jitter(data = occurance.number, aes(LONG_FIADB, LAT_FIADB, color = Species))
# get state information
states <- map_data("state")
#9=CT, 25=MA, 33=NH, 23=ME, 50=VT, 44=RI, 42=PA, 39=OH, 54=WV
state_sub <- filter(states, region %in% c("connecticut","maine","new hampshire","vermont","new york", "new jersey",
                                          "rhode island","pennsylvania","ohio","west virginia", "massachusetts"))


# pull the terrain background from stamenmaps
height <- max(state_sub$lat) - min(state_sub$lat)
width <- max(state_sub$long) - min(state_sub$long)
ne_borders <- c(bottom  = min(state_sub$lat)  - 0.1 * height, 
                top     = max(state_sub$lat)  + 0.1 * height,
                left    = min(state_sub$long) - 0.1 * width,
                right   = max(state_sub$long) + 0.1 * width)
#map <- get_stadiamap(ne_borders, zoom = 6, maptype = "stamen_terrain")


library(maps)
library(mapdata)
states <- map_data("state")
#9=CT, 25=MA, 33=NH, 23=ME, 50=VT, 44=RI, 42=PA, 39=OH, 54=WV
state_sub <- filter(states, region %in% c("connecticut","maine","new hampshire","vermont","new york", "new jersey",
                                          "rhode island","pennsylvania","ohio","west virginia", "massachusetts", "virginia", "delaware", 
                                          "north carolina", "kentucky", "tennessee", "michigan", "indiana", "district of columbia", "south carolina", 
                                          "georgia", "maryland"))

canada <- map_data("worldHires", "Canada")

ggplot() +
  geom_polygon(data = canada, 
               aes(x=long, y=lat, group = group), 
               color = "black", fill = "white") +
  geom_polygon(data = state_sub, 
               aes(x=long, y=lat, group = group), 
               color = "black", fill = "white")+
  geom_jitter(data = occurance.number, aes(LONG_FIADB, LAT_FIADB, color = Species))+
  coord_sf(xlim = c(-85, -68), ylim = c(35, 48))+theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'lightblue'))
ggsave(height = 8, width = 10, "images/cored_plots_17_species_map.png")




# # Groups:   exprem > 0 [2]
# `exprem > 0` `mortfac > 0`  `n()`
# <lgl>        <lgl>          <int>
#   1 FALSE        FALSE         185222
# 2 TRUE         FALSE         241877
# 3 TRUE         TRUE           10834

# do all the trees have a pmort in them?
# 
TREE.remeas %>% group_by(PLOT.ID, point, tree) %>% mutate(mortfac.tree = ifelse(mortfac > 0, 0, 1))%>% 
  ungroup() %>% group_by(PLOT.ID) %>% mutate(mortfac.total = sum(mortfac.tree))%>% 
  filter(mortfac.total == 0 )

saveRDS(TREE.remeas, "data/unfiltered_TREE.remeas.rds")

####################################################################################################
# explore the consequences of our filtering
####################################################################################################
# save a table of diameter differences by species and status
TREE.remeas %>% filter( exprem > 0 & dbhold > 0 & ! remper == 0) %>% 
  #dplyr::select(PLOT.ID, state, spp, remper, status, DIA_DIFF, dbhold, dbhcur, crcls, point, state, county, pltnum, tree, date) %>%
  group_by(PLOT.ID, point, state, county, pltnum, tree, date) %>% 
  mutate(annual.growth = ifelse(status %in% c(2, 3, 4) & ! is.na(dbhold) & ! remper == 0, 
                                DIA_DIFF/(remper/2), DIA_DIFF/remper),
         M = ifelse(status %in%  c(2, 3, 4), 1, 0), 
         relative.growth = (annual.growth/dbhold)*100) %>% ungroup() %>% 
  filter(SPCD %in% nspp[1:17,]$spp)%>%
  mutate(Tree.status = ifelse(M == 1, "dead", "live"), 
         DIAMETER_diff = ifelse(DIA_DIFF > 0, "positive", 
                                ifelse(DIA_DIFF == 0,"zero", "negative")))%>%
  
  group_by(Species, SPCD, Tree.status, DIAMETER_diff) %>% summarise(n()) %>% ungroup() %>% 
  group_by(Species, SPCD, Tree.status) %>% spread(DIAMETER_diff, `n()`) %>% ungroup() |> gt() |> 
  gtsave("images/filtering_exploration/DIA_DIFF_by_spcd_live_dead_status.png")


# save a table of diameter differences by type of mortality
TREE.remeas %>% filter( exprem > 0 & dbhold > 0 & ! remper == 0) %>% 
  #dplyr::select(PLOT.ID, state, spp, remper, status, DIA_DIFF, dbhold, dbhcur, crcls, point, state, county, pltnum, tree, date) %>%
  group_by(PLOT.ID, point, state, county, pltnum, tree, date) %>% 
  mutate(annual.growth = ifelse(status %in% c(2, 3, 4) & ! is.na(dbhold) & ! remper == 0, 
                                DIA_DIFF/(remper/2), DIA_DIFF/remper),
         M = ifelse(status %in%  c(2, 3, 4), 1, 0), 
         relative.growth = (annual.growth/dbhold)*100) %>% ungroup() %>% 
  filter(SPCD %in% nspp[1:17,]$spp)%>%
  mutate(Tree.status = ifelse(M == 1, "dead", "live"), 
         DIAMETER_diff = ifelse(DIA_DIFF > 0, "positive", 
                                ifelse(DIA_DIFF == 0,"zero", "negative")), 
         mort.status = ifelse(status == 1, "live", 
                              ifelse(status == 2, "dead (not salvagable)", 
                                     ifelse(status == 3, "cut", 
                                            ifelse(status == 4, "dead (salvagable)", "snag")))))%>%
  
  group_by(Species, SPCD, mort.status, DIAMETER_diff) %>% summarise(n()) %>% ungroup() %>% 
  group_by(Species, SPCD) %>% spread(mort.status, `n()`) %>% select(Species, DIAMETER_diff, live, `dead (salvagable)`, `dead (not salvagable)`, cut, snag )%>% ungroup() |> gt() |> 
  gtsave("images/filtering_exploration/DIA_DIFF_by_spcd_status.png")

# save a table of diameter differences by type of mortality
TREE.remeas %>% filter( exprem > 0 & dbhold > 0 & ! remper == 0) %>% 
  #dplyr::select(PLOT.ID, state, spp, remper, status, DIA_DIFF, dbhold, dbhcur, crcls, point, state, county, pltnum, tree, date) %>%
  group_by(PLOT.ID, point, state, county, pltnum, tree, date) %>% 
  mutate(annual.growth = ifelse(status %in% c(2, 3, 4) & ! is.na(dbhold) & ! remper == 0, 
                                DIA_DIFF/(remper/2), DIA_DIFF/remper),
         M = ifelse(status %in%  c(2, 3, 4), 1, 0), 
         relative.growth = (annual.growth/dbhold)*100) %>% ungroup() %>% 
  filter(SPCD %in% nspp[1:17,]$spp)%>%
  mutate(Tree.status = ifelse(M == 1, "dead", "live"), 
         DIAMETER_diff = ifelse(DIA_DIFF > 0, "positive", 
                                ifelse(DIA_DIFF == 0,"zero", "negative")), 
         mort.status = ifelse(status == 1, "live", 
                              ifelse(status == 2, "dead (not salvagable)", 
                                     ifelse(status == 3, "cut", 
                                            ifelse(status == 4, "dead (salvagable)", "snag")))))%>%
  
  group_by( mort.status, DIAMETER_diff) %>% summarise(n()) %>% ungroup() %>% 
  spread(mort.status, `n()`) %>% select(DIAMETER_diff, live, `dead (salvagable)`, `dead (not salvagable)`, cut, snag )%>% ungroup() |> gt() |> 
  gtsave("images/filtering_exploration/DIA_DIFF_by_status.png")


#########################################################################################
# consequences of filtering for trees > 5 in
# save a table of diameter differences by species and status
TREE.remeas %>% filter( exprem > 0 & dbhold > 0 & ! remper == 0 & !status == 5) %>% 
  #dplyr::select(PLOT.ID, state, spp, remper, status, DIA_DIFF, dbhold, dbhcur, crcls, point, state, county, pltnum, tree, date) %>%
  group_by(PLOT.ID, point, state, county, pltnum, tree, date) %>% 
  mutate(annual.growth = ifelse(status %in% c(2, 3, 4) & ! is.na(dbhold) & ! remper == 0, 
                                DIA_DIFF/(remper/2), DIA_DIFF/remper),
         M = ifelse(status %in%  c(2, 3, 4), 1, 0), 
         relative.growth = (annual.growth/dbhold)*100) %>% ungroup() %>% 
  filter(SPCD %in% nspp[1:17,]$spp)%>%
  mutate(Tree.status = ifelse(M == 1, "dead", "live"), 
         DIAMETER_diff = ifelse(DIA_DIFF > 0, "positive", 
                                ifelse(DIA_DIFF == 0,"zero", "negative")))%>%
  
  group_by(Species, SPCD, dbhold > 5) %>% summarise(n()) %>% ungroup() %>% 
  group_by(Species, SPCD) %>% spread(`dbhold > 5`, `n()`) %>%
  rename("dhbold < 5"=`FALSE`, 
         "dbhold > 5"=`TRUE`) %>% ungroup() |> gt() |> 
  grand_summary_rows(
    columns = c("dhbold < 5", "dbhold > 5"),
    fns = list(
      total ~ sum(.)
    ),
    fmt = ~ fmt_number(., use_seps = FALSE, decimals = 0)
  )|>
  gtsave("images/filtering_exploration/TREE_size_dead_status.png")


#View(TREE.remeas2)
#unique(TREE.remeas %>% filter(dbhcur > 100)%>% select(state, dat
# filter out the trees that are not modeled exprem > 0 and 
TREE_growth.mort <- TREE.remeas %>% filter(DIA_DIFF >= 0 & exprem > 0 & dbhold > 0 & ! remper == 0) %>% 
  #dplyr::select(PLOT.ID, state, spp, remper, status, DIA_DIFF, dbhold, dbhcur, crcls, point, state, county, pltnum, tree, date) %>%
  group_by(PLOT.ID, point, state, county, pltnum, tree, date) %>% 
  mutate(annual.growth = ifelse(status %in% c(2, 3, 4) & ! is.na(dbhold) & ! remper == 0, 
                                DIA_DIFF/(remper/2), DIA_DIFF/remper),
         M = ifelse(status %in%  c(2, 3, 4), 1, 0), 
         relative.growth = (annual.growth/dbhold)*100)
#hist(TREE_growth.mort$annual.growth)

# # revised model
# TREE_remeas <- readRDS("data/periodic.nonharvest.data.RDS")
# summary(TREE_remeas$ba)
# # get the annualized growth rate:
# 
# TREE_remeas$INVYR <- substr(TREE_remeas$mdate, 1,2)
# TREE_remeas$DIA_DIFF <- TREE_remeas$dbhcur - TREE_remeas$dbhold
# hist(TREE_remeas$DIA_DIFF)


hist(TREE_growth.mort$annual.growth)
head(TREE_growth.mort$spp)

# join up with si information from PLOT


cleaned.data <- TREE_growth.mort %>% filter(!is.na(annual.growth) & !is.na(status) & DIA_DIFF >=0 & !is.na(elev) )
# View(cleaned.data %>% group_by(M, spp) %>% summarise(n()))
# View(cleaned.data %>% group_by(M, state) %>% summarise(n()))
# View(cleaned.data %>% group_by(M) %>% summarise(n()))
# 
# # only 2400 tree mortality events detected in the valid plots of this dataset:
# A tibble: 2 × 2
# M  `n()`
# <dbl>  <int>
#   1     0 136674
# 2     1   2400
# select top species in the TREE ring datasets:
#stan.tr <- readRDS("data/formatted/STANdata.RDS")
#unique(stan.tr$dbh.mat$SPCD)

#high.SPP <- stan.tr$dbh.mat %>% group_by(SPCD) %>% summarise(n = n()) %>% filter(n>1000) %>% dplyr::select(SPCD)

#cleaned.data <- cleaned.data %>% filter(SPCD %in% unique(high.SPP$SPCD))

# --------------------------Merge with plot covariate data ------------------
# get plot-level information:


unique.PLOT.cov <- unique(PLOT %>% dplyr::select(PLOT.ID, cndtn, physio, ba, slope, aspect, si, siage, mdate, cycle))
TREE_growth.cov <- left_join(TREE_growth.mort, unique.PLOT.cov)

# calculate the % damaged trees on the plot
pct.damages <- TREE.remeas %>% group_by(PLOT.ID) %>% mutate(n.tree = n()) %>% ungroup() %>% group_by(PLOT.ID, damage)%>%
  mutate(n.damage = n(), 
         pct.damage.by.cd = (n.damage/n.tree)*100) %>% ungroup()%>% dplyr::select(PLOT.ID,damage, pct.damage.by.cd) %>% 
  group_by(PLOT.ID) %>% distinct() %>% spread(damage, pct.damage.by.cd, fill = 0)

summary(pct.damages)

colnames(pct.damages)[2:10] <- paste0("damage.",colnames(pct.damages)[2:10])

TREE_growth.cov <- left_join(TREE_growth.cov, pct.damages)

ggplot(TREE_growth.cov, aes(M, 100-damage.0, group = M))+geom_violin()
TREE_growth.cov <- TREE_growth.cov %>% mutate(damage.total = 100-damage.0)


# get climate data from PRISM:
# note that we need to reextract it for the plots without tree cores:
# need to use adapt PRISM extraction code with terra packages
# code to extract PRISM climate is in get_PRISM_periodic_timeseries.R
length(unique(TREE_growth.cov$PLOT.ID)) # 11,128 unique plots

# read in RDS:
plots.clim <- readRDS(paste0(boxdir, "data/seasonal_climate_periodic_NE_plots.RDS"))


# colnames(plots)
# plots.clim.m <- reshape2::melt(plots, id.vars = c("PLOT.ID","STATECD" , "COUNTYCD" ,     "PLOT" ,        
#                                                   "MAX_INVYR" ,    "COUNTYCD_NIMS", "PLOT_NIMS"  ))
# plots.clim.m <- plots.clim.m %>% filter(!variable %in% c("tmax.2015.", "tmin.2015.", "ppt.2015."))
# 
# # split up the long variable name into climate, month, and year
# clim.month.split <- strsplit(as.character(plots.clim.m$variable),".",fixed = TRUE)
# clim.month <- do.call(rbind, clim.month.split)
# 
# clim.month <- data.frame(clim.month)
# colnames(clim.month) <- c("climate", "year","month")
# plots.clim.m <- cbind(plots.clim.m %>% dplyr::select(-variable), clim.month)#separate(variable, c("climate", "month")) #%>% rename(CLIMNA = value)
# 
# 
# # get the months and climate and years
# plots.clim.m <- plots.clim.m %>% 
#   #separate(variable, into =c("climate", "year", "month")) %>%
#   rename(PRISM = value, 
#          climate.variable = climate, 
#          Year = year)
# 
# plots.clim.m$Year <- as.numeric(plots.clim.m$Year)
# 
# #plots.clim.m %>% group_by(PLOT.ID, climate.variable, month)
# # prefilter each climate variable so we can match it up the CLIMNA more efficiently
# tmax.prism <- plots.clim.m %>% filter(climate.variable %in% "tmax")
# tmin.prism <- plots.clim.m %>% filter(climate.variable %in% "tmin")
# tave.prism <- plots.clim.m %>% filter(climate.variable %in% "tave")
# ppt.prism <- plots.clim.m %>% filter(climate.variable %in% "ppt")

# generate climate summaries for each variable, 
# get MAP and MAT
# also get annual values to be aggregated


MAP <- plots.clim %>% group_by(PLOT.ID) %>% filter(!year %in% 1895) %>% summarise(MAP = mean(wateryr_PPT, na.rm =TRUE))
MATmin <- plots.clim %>% group_by(PLOT.ID) %>% filter(!year %in% 1895) %>% summarise(MATmin = mean(yr_MeanTmin, na.rm =TRUE))
MATmax <- plots.clim %>% group_by(PLOT.ID) %>% filter(!year %in% 1895) %>% summarise(MATmax = mean(yr_MeanTmax, na.rm =TRUE))

PLOT.ID.remper <- unique(TREE_growth.cov %>% ungroup()%>% dplyr::select(PLOT.ID, date, remper))

annual.ppt.remp <- left_join(plots.clim, PLOT.ID.remper)

annual.ppt.remp <- annual.ppt.remp %>% rename(Year = year)
# get the remper anomoly from the site mean
annual.ppt.anom <- annual.ppt.remp %>% ungroup()%>% group_by(PLOT.ID, remper, date)%>%
  mutate(MAP = mean(wateryr_PPT, na.rm=TRUE),
         MAP.sd = sd(wateryr_PPT, na.rm =TRUE), 
         MATmin = mean(yr_MeanTmin, na.rm =TRUE),
         MATmin.sd = sd(yr_MeanTmin, na.rm =TRUE),
         MATmax = mean(yr_MeanTmax, na.rm =TRUE), 
         MATmax.sd = sd(yr_MeanTmax, na.rm =TRUE)) %>% 
  
  ungroup()%>% group_by(PLOT.ID, remper, date) %>%
  mutate(ppt.anomoly = (wateryr_PPT - MAP)/MAP.sd, 
         tmax.anomoly = (yr_MeanTmax- MATmax)/MATmax.sd, 
         tmin.anomoly = (yr_MeanTmin - MATmin)/MATmin.sd,
         date2 = as.numeric(date) - as.numeric(remper)) %>% ungroup()%>%
  group_by(PLOT.ID) %>% filter(Year >=date2 & Year <= date) %>%
  summarise(ppt.anom = mean(ppt.anomoly), 
            tmax.anom = mean(tmax.anomoly), 
            tmin.anom = mean(tmin.anomoly))

hist(annual.ppt.anom$ppt.anom)
hist(annual.ppt.anom$tmax.anom)
hist(annual.ppt.anom$tmin.anom)

annual.ppt.anom
# join up to the tree database
TREE_growth.cov <- left_join(TREE_growth.cov, MATmin)
TREE_growth.cov <- left_join(TREE_growth.cov, MATmax)
TREE_growth.cov <- left_join(TREE_growth.cov, MAP)
TREE_growth.cov <- left_join(TREE_growth.cov, annual.ppt.anom)

rm(annual.ppt, annual.tmin, annual.tmax, ppt.prism, tmin.prism, tmax.prism, plots.clim.m, clim.month.split, clim.month, radial.inc, spplots.elev)
rm(MAP, MATmax, MATmin, annual.ppt.anom, annual.ppt.remp, occurance.number, PRISM, spplots)

# do some renaming to make things easier
#TREE_growth.cov <- TREE_growth.cov %>% rename(SPCD = spp)

# calculate BAL basal area larger than
TREE.remeas.BAL <- TREE.remeas %>%
  mutate(ba_sq_ft = ((dbhcur^2))*0.005454) %>% 
  group_by(PLOT.ID, point, date, cndtn) %>% 
  mutate(BAL = sitreeE::PBAL(ba_sq_ft*volfac)) %>% select(state, county, pltnum, cndtn, point, tree, SPCD, status, 
                                                   dbhcur, dbhold, PLOT.ID, BAL, ba_sq_ft, volfac, cycle)

TREE_growth.cov <- left_join(TREE_growth.cov, TREE.remeas.BAL) %>% distinct()

# calculate species-level BA for like species and non-like species for each tree
PLOT.remeas.BAL <- TREE.remeas %>%
  mutate(ba_sq_ft = ((dbhcur^2))*0.005454) %>% 
  group_by(PLOT.ID, point, date, cndtn, SPCD) %>%
  mutate(SPCD_BA = sum(ba_sq_ft*volfac, na.rm =TRUE), 
         SPCD_density = sum(volfac, na.rm = TRUE)) %>% ungroup() %>%
  group_by(PLOT.ID, point, date, cndtn) %>%
  mutate(BA_total = sum(SPCD, na.rm =TRUE), 
         density_total = sum(volfac, na.rm = TRUE)) %>%
  mutate(non_SPCD_BA = BA_total - SPCD_BA, 
         non_SPCD_density = density_total - SPCD_density) %>%  
  select(PLOT.ID, point, date, cndtn, SPCD, SPCD_BA, BA_total, non_SPCD_BA, SPCD_density, density_total, non_SPCD_density)

TREE_growth.cov <- left_join(TREE_growth.cov, PLOT.remeas.BAL) %>% distinct()

cleaned.data <- TREE_growth.cov %>% filter(!is.na(annual.growth) & !is.na(status) & DIA_DIFF >=0 & !is.na(elev))%>% distinct()

# ggplot(cleaned.data, aes(non_SPCD_BAL, M))+geom_point()+facet_wrap(~SPCD)
# ggplot(cleaned.data, aes(elev, M))+geom_point()+facet_wrap(~SPCD)
# ggplot(cleaned.data, aes(BAL, M))+geom_point()+facet_wrap(~SPCD)
# ggplot(cleaned.data, aes(DIA_DIFF, M))+geom_point()+facet_wrap(~SPCD)
# ggplot(cleaned.data, aes(si, M))+geom_point()+facet_wrap(~SPCD)
# ggplot(cleaned.data, aes(dbhcur, M))+geom_point()+facet_wrap(~SPCD)

# View(cleaned.data %>% group_by(M, SPCD) %>% summarise(n()))
# 
# ggplot(cleaned.data, aes(LONG_FIADB, LAT_FIADB, color = M))+geom_jitter(size = 0.5)


summary(cleaned.data$elev)
saveRDS(cleaned.data, "data/cleaned.data.mortality.TRplots.RDS")

#----------------------------------------------------------------------
# Get the N deposition data for all the plots
# N deposition timeseries/spatial variation from...
# the time series go from 1930-2017 by state and county...so it should be simple to match up 
N.oxidized <- read_delim(paste0(boxdir,"data/Ndep/Atmospheric_Oxidized.txt"))
N.oxidized$county <- as.numeric(N.oxidized$CountyFIPS)
N.oxidized$state <- as.numeric(N.oxidized$StateFIPS)

N.reduced <- read_delim(paste0(boxdir, "data/Ndep/Atmospheric_Reduced.txt"))
N.reduced$county <- as.numeric(N.reduced$CountyFIPS)
N.reduced$state <- as.numeric(N.reduced$StateFIPS)

Ntotal <- N.reduced
Ntotal[,5:92] <- N.reduced[,5:92]+N.oxidized[,5:92] # add together to get the total:

# get the county and state for each plot and then put the timeseries together by plotid for all trees:
N.total.dep <- Ntotal %>% dplyr::select(state, county, y1930:y2017)
N.total.dep

Ndep.plot.df <- left_join(cleaned.data %>% ungroup() %>% dplyr::select(PLOT.ID, cycle, date, remper, county, state) %>% distinct(), N.total.dep)
Ndep.melt <- reshape2::melt(Ndep.plot.df,id.vars = c( "date", "PLOT.ID", "cycle", "remper", "county", "state") )
Ndep.melt$year <- substring(Ndep.melt$variable, 2)

Ndep.remper <- Ndep.melt %>% group_by(PLOT.ID, cycle) %>% mutate(prev.date = date - remper) %>% group_by(PLOT.ID, cycle)%>%
                              filter(year >= prev.date & year <= date) %>% group_by(PLOT.ID, cycle, county, state) %>% 
                              mutate(Ndep.remper.avg = mean(value, na.rm =TRUE), 
                                     Nyears = n())
#View(Ndep.remper)

rm(test, plots.clim)
# save Ndep average for the remper
saveRDS(Ndep.remper, "data/Ndep.average.remper.NE.periodic.RDS")

Ndep.remper <- Ndep.remper %>% dplyr::select(date, PLOT.ID, cycle, remper, county, state, Ndep.remper.avg) %>% distinct()
# join up to cleaned data
cleaned.data <- left_join(cleaned.data, Ndep.remper)
colnames(cleaned.data)

# joining tables leads to alot of NA values (missing counties??)
cleaned.data <- cleaned.data %>% filter(!is.na(Ndep.remper.avg))

# lets check the plot anyways
map.ndep.trees <- ggplot(cleaned.data, aes(x = LONG_FIADB, y = LAT_FIADB, color = Ndep.remper.avg))+geom_point()

png(height = 6, width = 10, units = "in", res = 150, "images/Ndep_map_trees_NE_FIA.png")
map.ndep.trees
dev.off()

saveRDS(cleaned.data, "data/cleaned.data.mortality.TRplots.RDS")


#----------------------------------------------------------------
# get all the plot and tree information for the whole NE:

TREE <- read_delim("data/formatted_older_matching_plts_TREE.txt")
TREE$dbhcur <- as.numeric(TREE$dbhcur)
TREE$dbhold <- as.numeric(TREE$dbhold)

#TREE <- TREE %>% filter(PLOT.ID %in% unique(plots$PLOT.ID))
TREE$SPCD <- TREE$spp
TREE$Species <- FIESTA::ref_species[match(TREE$SPCD, FIESTA::ref_species$SPCD),]$COMMON


summary(TREE$dbhcur)
summary(TREE$crcls)
unique(TREE$crcls)


PLOT <- read_delim("data/formatted_older_matching_plts_PLOT.txt")
colnames(PLOT)
summary(PLOT$ba)
PLOT$typcur

TPA <- TREE %>% group_by(state, county, pltnum, SPCD) %>% summarise(TPA = sum(volfac, na.rm =TRUE))

# Get TPA and BA for each plot and put on the log scale

PLOT.tpa <- left_join(PLOT, TPA)
PLOT.tpa.all <- PLOT.tpa %>% dplyr::select(state, county, pltnum,SPCD, typcur, ba, TPA) %>% filter(!is.na(TPA)) %>%
  filter(TPA >0 & ba > 0) %>% distinct()

unique(PLOT.tpa.all$typcur)

cleaned.data <- readRDS( "data/cleaned.data.mortality.TRplots.RDS")
test.m <- left_join(cleaned.data, PLOT[, c("state", "county", "pltnum", "typcur")])


PLOT.tpa.all.5 <- PLOT.tpa.all %>% group_by( SPCD ) %>% summarise(n()) %>% filter(`n()` >5 )
# create data for stan model

library(quantreg)

 # fit a species-level quantile regression
fit.quantile.tmax.regression <- function(type = 12){
  prepped <- PLOT.tpa.all %>% filter(SPCD == type) %>% mutate(baa = ba/TPA, 
                                                                ba.log = log(baa), 
                                                                TPA.log = log(TPA))
  
  model <- rq(TPA.log ~ ba.log, data = prepped, tau = 0.99)
  model.coefs <-  data.frame(SPCD = type, 
                             alpha = coef(model)[1],
                             beta = coef(model)[2])
  model.coefs
}

fit.quantile.tmax.regression(type = unique(PLOT.tpa.all$SPCD[2]))
all.model.coeffs <- list()

for(i in 1:length(unique(PLOT.tpa.all.5$SPCD))){
  all.model.coeffs[[i]] <- fit.quantile.tmax.regression(type = unique(PLOT.tpa.all.5$SPCD)[i])
}

# use the coefficients to get the RD for each plot:
all.model.coeffs.df <- do.call(rbind, all.model.coeffs)
all.model.coeffs.df

# now do the estimates on all the trees
TREE.RD <- left_join(TREE.remeas.BAL %>% filter(!is.na(volfac)), all.model.coeffs.df) 
Nmax <- TREE.RD %>% mutate(TPA.log = log(volfac), 
                                BA.log = log(ba_sq_ft)) %>% 
  mutate(Nmax = exp(alpha + BA.log*beta))%>% mutate(N.Nmax = volfac/Nmax)
RD.by.plot <- Nmax %>% group_by(state, county, pltnum, PLOT.ID, cycle)%>%
  summarise(RD = sum(N.Nmax, na.rm = TRUE))
hist(RD.by.plot$RD)
summary(RD.by.plot$RD)
#RD.by.plot %>% filter(RD > 1)
# why are there some very large values??


saveRDS (RD.by.plot, "data/Relative_periodic_plot_densities.rds")
cleaned.data2 <- left_join(cleaned.data, RD.by.plot)
cleaned.data2 <- cleaned.data2  %>% ungroup() %>% 
  mutate(BAL.ratio = SPCD_BA/BA_total, 
         density.ratio = SPCD_density/total_density)

saveRDS(cleaned.data2, "data/cleaned.data.mortality.TRplots.RDS")
ggplot(cleaned.data2, aes(x = LONG_FIADB, y = LAT_FIADB))+geom_point()

unique(cleaned.data2$SPCD)
colnames(cleaned.data2)

rm(RD.by.plot, annual.ppt.anom, annual)
# ------------------------------------------------------------------------------
# make some intital exploratory data plots
# ------------------------------------------------------------------------------
mort.17 <- cleaned.data2 %>% filter(SPCD %in% c(12, 97, 129, 241, 261, 316, 318, 371, 400, 531, 541, 621, 762, 802, 832, 833, 837))


ggplot(mort.17, aes(x= as.character(M), y = non_SPCD_BA))+geom_violin()+facet_wrap(~SPCD)
ggplot(mort.17, aes(x= as.character(M), y = SPCD_BA))+geom_violin()+facet_wrap(~SPCD)
ggplot(mort.17, aes(x= as.character(M), y = BAL.ratio))+geom_violin()+facet_wrap(~SPCD, scales = "free_y")
ggplot(mort.17, aes(x= as.character(M), y = density.ratio))+geom_violin()+facet_wrap(~SPCD, scales = "free_y")
ggplot(mort.17, aes(x= as.character(M), y = RD))+geom_boxplot()+facet_wrap(~SPCD)
ggplot(mort.17, aes(x= as.character(M), y = BAL))+geom_boxplot()+facet_wrap(~SPCD)
ggplot(mort.17, aes(x= as.character(M), y = damage.total))+geom_boxplot()+facet_wrap(~SPCD)
ggplot(mort.17, aes(x= as.character(M), y = ba))+geom_boxplot()+facet_wrap(~SPCD)
ggplot(mort.17, aes(x= as.character(M), y = elev))+geom_boxplot()+facet_wrap(~SPCD)
ggplot(mort.17, aes(x= as.character(M), y = annual.growth))+geom_boxplot()+facet_wrap(~SPCD)
ggplot(mort.17, aes(x= as.character(M), y = dbhcur))+geom_boxplot()+facet_wrap(~SPCD)

library(ggbeeswarm)
cleaned.data2 <- readRDS("data/cleaned.data.mortality.TRplots.RDS")
#g + 
 # ggdist::stat_dots(side = "both")

g + 
  ggdist::stat_dots(position = position_nudge(x = -.25))

nspp[1:20,]$spp
# filter to just get the target species
cleaned.data.spp <- cleaned.data2 %>% filter(SPCD %in% nspp[1:17,]$SPCD)

# get common names
cleaned.data.spp$COMMON <- FIESTA::ref_species[match(cleaned.data.spp$SPCD, FIESTA::ref_species$SPCD),]$COMMON
cleaned.data.spp <- cleaned.data.spp %>% mutate(`Tree Status` = ifelse(M == 1, "Dead", "Live"))
# ggplot(data = cleaned.data.spp, aes(x = factor(Species), y = ppt.anom, fill = as.character(M)))+ 
#   ggbeeswarm::geom_quasirandom(size = 1, width = .33, alpha = .3) 

# ggplot(data = cleaned.data.spp, aes(x = factor(Species), y = ppt.anom, fill = as.character(M)))+ 
#   ggbeeswarm::geom_quasirandom(size = 1, width = .33, alpha = .3) 
# 
# cleaned.data.spp 

# color scheme
my_pal <- rcartocolor::carto_pal(n = 8, name = "Bold")[c(1, 3, 7, 2)]

## theme for horizontal charts
theme_flip <-
  theme(
    axis.text.x = element_text(face = "plain", family = "Ariel", size = 16),
    axis.text.y = element_text(face = "bold", family = "Ariel", size = 16),
    panel.grid.major.x = element_line(color = "grey90", size = .6),
    panel.grid.major.y = element_blank(),
    legend.position = "top", 
    legend.text = element_text(family = "Ariel", size = 18),
    legend.title = element_text(face = "bold", size = 18, margin = margin(b = 25))
  )

g_ridges <- 
  ggplot(cleaned.data.spp , aes(ppt.anom, fct_rev(`Tree Status`), color =`Tree Status`, fill = `Tree Status`)) + 
  coord_cartesian(clip = "off") +
  scale_y_discrete(expand = c(.07, .07)) +
  scale_color_manual(values = my_pal, guide = "none") +
  scale_fill_manual(values = my_pal, guide = "none") + facet_wrap(~COMMON)+
  theme_flip

g_ridges +
  ggridges::geom_density_ridges(
    alpha = .7, size = 1.5
  )

g_ridges + 
  ggridges::stat_density_ridges(
    quantile_lines = TRUE, quantiles = 2, 
    color = "black", alpha = .8, size = 1.5
  ) 


#######################################################################################
histogram.ppt.plt <- ggplot(cleaned.data.spp, aes(x = forcats::fct_rev(`Tree Status`), y = ppt.anom, 
                 color = `Tree Status`, fill = `Tree Status`)) +
  geom_boxplot(
    width = .2, fill = "white",
    size = 1, outlier.shape = NA
  ) +
  ggdist::stat_halfeye(
   # adjust = .33,
    #width = .67, 
    color = NA,
    position = position_nudge(x = .15)
  ) +
  geom_point(
    position = position_nudge(x = -.22),
    shape = 95, size = 24, alpha = .25
  )+
  coord_flip() +
  scale_x_discrete(expand = c(.07, .07)) +
  scale_y_continuous(breaks = 1:9) +
  scale_color_manual(values = my_pal, guide = "none") +
  scale_fill_manual(values = my_pal, guide = "none") +
  facet_wrap(~COMMON)+
  theme_flip

histogram.ppt.plt
ggsave(height = 12, width = 12, units = "in", "images/PPT_Mortality_plot.png")


histogram.growth.plt <- ggplot(cleaned.data.spp, aes(x = forcats::fct_rev(`Tree Status`), y = annual.growth, 
                                                  color = `Tree Status`, fill = `Tree Status`)) +
  geom_boxplot(
    width = .2, fill = "white",
    size = 1, outlier.shape = NA
  ) +
  ggdist::stat_halfeye(
    # adjust = .33,
    #width = .67, 
    color = NA,
    position = position_nudge(x = .15)
  ) +
  geom_point(
    position = position_nudge(x = -.22),
    shape = 95, size = 24, alpha = .25
  )+
  coord_flip() +
  scale_x_discrete(expand = c(.07, .07)) +
  scale_y_continuous(breaks = 1:9) +
  scale_color_manual(values = my_pal, guide = "none") +
  scale_fill_manual(values = my_pal, guide = "none") +
  facet_wrap(~COMMON)+
  theme_flip

histogram.growth.plt 
ggsave(height = 12, width = 12, units = "in", "images/Growth_Mortality_plot.png")



ggplot()+geom_violin(data = cleaned.data.spp, aes(x = factor(Species), y = ppt.anom, fill = as.character(M)))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot()+geom_violin(data = cleaned.data.spp, aes(x = factor(Species), y = tmax.anom, fill = as.character(M)))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot()+geom_violin(data = cleaned.data.spp, aes(x = factor(Species), y = tmin.anom, fill = as.character(M)))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggplot()+geom_violin(data = cleaned.data.spp, aes(x = factor(Species), y = damage.total, fill = as.character(M)))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot()+geom_violin(data = cleaned.data.spp, aes(x = factor(Species), y = MAP, fill = as.character(M)))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot()+geom_violin(data = cleaned.data.spp, aes(x = factor(Species), y = MATmin, fill = as.character(M)))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot()+geom_violin(data = cleaned.data.spp, aes(x = factor(Species), y = MATmax, fill = as.character(M)))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot()+geom_violin(data = cleaned.data.spp, aes(x = factor(Species), y = RD, fill = as.character(M)))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot()+geom_violin(data = cleaned.data.spp, aes(x = factor(Species), y = BAL, fill = as.character(M)))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot()+geom_violin(data = cleaned.data.spp, aes(x = factor(Species), y = aspect, fill = as.character(M)))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot()+geom_violin(data = cleaned.data.spp, aes(x = factor(Species), y = slope, fill = as.character(M)))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot()+geom_violin(data = cleaned.data.spp %>% filter(relative.growth < 10), aes(x = factor(Species), y = relative.growth, fill = as.character(M)))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot()+geom_violin(data = cleaned.data.spp %>% filter(relative.growth < 10), aes(x = factor(Species), y = si, fill = as.character(M)))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot()+geom_violin(data = cleaned.data.spp %>% filter(relative.growth < 10), aes(x = factor(Species), y = elev, fill = as.character(M)))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
