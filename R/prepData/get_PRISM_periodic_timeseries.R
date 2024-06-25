#------------------------------------------------------
# get the periodic plot lat longs that match:
#------------------------------------------------------
library(tidyverse)
library(odbc)
library(terra)

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

NE_plot <- dbGetQuery(ocon, "SELECT cn, statecd, unitcd, countycd, plot, invyr, plot_status_cd, cycle, lat, lon, elev
                      FROM fs_fiadb.plot
                      WHERE statecd =  ANY(09, 10, 23, 24, 25, 33, 34, 36, 39, 42, 44, 50, 54) 
                      and invyr < 1999") %>%
  as_tibble()%>%
  rename_with(tolower)
unique(NE_plot$statecd)
length(unique(NE_plot$cn))
NE_plot

# when done, disconnect from ORACLE FIADB
dbDisconnect(ocon)
rm(ocon)

summary(NE_plot$lat)
#summary(NE_plot$elev)

# get matching plot ids
NE_plot$PLOT.ID <- as.numeric(paste0(NE_plot$statecd, 
                                     sprintf("%03d", NE_plot$countycd),
                                     sprintf("%04d", NE_plot$plot)))


length(unique(NE_plot$PLOT.ID))
length(unique(NE_plot$cn))
nrow(unique(NE_plot[, c("lat", "lon")]))
# 

#-------------------------------------------------------------------------------
# Get prism time series for each unique lat long
#-------------------------------------------------------------------------------

### PRISM download Feb 13, 2024



# Search for PRISM PPT files
PRISM.path <-  "./PRISM_data/"
ppt.path <-  "./PRISM_data/PRISM_ppt_stable_4kmM2_189501_198012_bil/"
ppt.path.new <-  "./PRISM_data/PRISM_ppt_stable_4kmM3_198101_202307_bil/"


pptFiles.old <- list.files(path = ppt.path, pattern = glob2rx("*ppt*.bil"), full.names = TRUE)
pptFiles.new <- list.files(path = ppt.path.new, pattern = glob2rx("*ppt*.bil"), full.names = TRUE)
pptFiles <- c(pptFiles.old, pptFiles.new)

# do the same for Tmax
PRISM.path <-  "./PRISM_data/"
tmax.path <-  "./PRISM_data/PRISM_tmax_stable_4kmM3_189501_198012_bil/"
tmax.path.new <-  "./PRISM_data/PRISM_tmax_stable_4kmM3_198101_202307_bil/"


tmaxFiles.old <- list.files(path = tmax.path, pattern = glob2rx("*tmax*.bil"), full.names = TRUE)
tmaxFiles.new <- list.files(path = tmax.path.new, pattern = glob2rx("*tmax*.bil"), full.names = TRUE)
tmaxFiles <- c(tmaxFiles.old, tmaxFiles.new)


# do the same for Tmin
PRISM.path <-  "./PRISM_data/"
tmin.path <-  "./PRISM_data/PRISM_tmin_stable_4kmM3_189501_198012_bil/"
tmin.path.new <-  "./PRISM_data/PRISM_tmin_stable_4kmM3_198101_202307_bil/"


tminFiles.old <- list.files(path = tmin.path, pattern = glob2rx("*tmin*.bil"), full.names = TRUE)
tminFiles.new <- list.files(path = tmin.path.new, pattern = glob2rx("*tmin*.bil"), full.names = TRUE)
tminFiles <- c(tminFiles.old, tminFiles.new)

# make one big raster out of the precip monthly files (make take a little while)
ppt.all.rast <- rast(pptFiles)
ppt.extract <- extract(ppt.all.rast, NE_plot[, c("lon", "lat")], cells=FALSE, method="simple")

ppt.extract$lon <- NE_plot$lon
ppt.extract$lat <- NE_plot$lat

saveRDS(ppt.extract, "data/PPT_extracted.ppt.data_v1.rds")


# make one big raster out of the tmax monthly files (make take a little while)
tmax.all.rast <- rast(tmaxFiles)
tmax.extract <- extract(tmax.all.rast, NE_plot[, c("lon", "lat")], cells=FALSE, method="simple")

tmax.extract$lon <- NE_plot$lon
tmax.extract$lat <- NE_plot$lat

saveRDS(tmax.extract, "data/tmax_extracted.tmax.data_v1.rds")

# make one big raster out of the tmin monthly files (make take a little while)
tmin.all.rast <- rast(tminFiles)
tmin.extract <- extract(tmin.all.rast, NE_plot[, c("lon", "lat")], cells=FALSE, method="simple")

tmin.extract$lon <- NE_plot$lon
tmin.extract$lat <- NE_plot$lat

saveRDS(tmin.extract, "data/tmin_extracted.tmin.data_v1.rds")
ppt.extract <- readRDS("data/PPT_extracted.ppt.data_v1.rds")
tmin.extract <- readRDS("data/tmin_extracted.tmin.data_v1.rds")
tmax.extract <- readRDS("data/tmax_extracted.tmax.data_v1.rds")

# get the seasonal varaibles
colnames(ppt.extract) 
tail(colnames(tmin.extract)) 


tail(colnames(tmax.extract))
years <- rep(1895:2022, each = 12)
months <- rep(1:12, 128)

colnames.clim <- c(paste0(years, "_", months), "2023_1", "2023_2", "2023_3", "2023_4", "2023_5", "2023_6", "2023_7")
length(ppt.extract)
length(colnames.clim)
colnames(ppt.extract)[2:(1 + length(colnames.clim))] <- colnames.clim
colnames(tmin.extract)[2:(1 + length(colnames.clim))] <- colnames.clim
colnames(tmax.extract)[2:(1 + length(colnames.clim))] <- colnames.clim


 


# reorganize the data frames

yr.month <- data.frame(do.call(rbind, str_split(colnames.clim, "_",  n = 2)))
yr.month$yr.month <- colnames.clim
colnames(yr.month)<- c("year", "month", "variable")


ppt.m <- melt(ppt.extract, id.vars = c("lon", "lat", "ID"))
ppt.long <- left_join(ppt.m, yr.month)
colnames(ppt.long) <- c("lon", "lat", "ID", "yr.month", "ppt",
                        "year", "month")

tmin.m <- melt(tmin.extract, id.vars = c("lon", "lat", "ID"))
tmin.long <- left_join(tmin.m, yr.month)
colnames(tmin.long) <- c("lon", "lat", "ID", "yr.month", "tmin",
                        "year", "month")

tmax.m <- melt(tmax.extract, id.vars = c("lon", "lat", "ID"))
tmax.long <- left_join(tmax.m, yr.month)
colnames(tmax.long) <- c("lon", "lat", "ID", "yr.month", "tmax",
                        "year", "month")
temp.long <- left_join(tmax.long, tmin.long)

temp.ppt <- left_join(temp.long, ppt.long)


# assign water year
wtr_yr <- function(df, start_month=9) {
  # Year offset
  offset = ifelse(as.numeric(df$month) >= start_month - 1, 1, 0)
  # Water year
  adj.year = as.numeric(df$year) + offset
  # Return the water year
  adj.year
}

temp.ppt$water_year <- wtr_yr(temp.ppt)
# get total water year
water_year <- temp.ppt %>% filter(!water_year == 2023) %>% group_by(lon, lat, ID, water_year)%>% summarise(wateryr_PPT = sum(ppt, na.rm =TRUE), 
                                                                                                           wateryr_MeanTmax = mean(tmax, na.rm =TRUE), 
                                                                                                           wateryr_MeanTmin = mean(tmin, na.rm =TRUE) ) 

colnames(water_year)[4] <- "year"
# get annual means
year_means <- temp.ppt %>% filter(!year == 2023) %>% group_by(lon, lat, ID, year)%>% summarise(yr_PPT = sum(ppt, na.rm =TRUE), 
                                                                                               yr_MeanTmax = mean(tmax, na.rm =TRUE), 
                                                                                               yr_MeanTmin = mean(tmin, na.rm =TRUE) ) 



# FAll September october and november # dont use water year use calender year for fall
SeptOctNov <- temp.ppt %>% filter(month %in% 9:11 & !year == 2023)%>% group_by(lon, lat, ID, year)%>% summarise(
                                                                                                  Precip_SeptOctNov = sum(ppt), 
                                                                                                  Tmax_SeptOctNov = mean(tmax),
                                                                                                  Tmin_SeptOctNov = mean(tmin, na.rm=TRUE)) 


# Summer
JunJulyAug <- temp.ppt %>% filter(month %in% 6:8 & !water_year == 2023)%>% group_by(lon, lat, ID, water_year)%>% summarise( Precip_JunJulyAug = sum(ppt, na.rm = TRUE), 
                                                                                                             Tmax_JunJulyAug = mean(tmax, na.rm = TRUE),
                                                                                                             Tmin_JunJulAug = mean(tmin, na.rm=TRUE)) 

colnames(JunJulyAug)[4] <- "year"
# Spring 
Spring <- temp.ppt %>% filter(month %in% 9:11 & !water_year == 2023)%>% group_by(lon, lat, ID, water_year)%>% summarise( Precip_MarAprMay = sum(ppt, na.rm = TRUE), 
                                                                                                          Tmax_MarAprMay = mean(tmax, na.rm = TRUE),
                                                                                                          Tmin_MarAprMay = mean(tmin, na.rm=TRUE))
colnames(Spring)[4] <- "year"
# december january february - winter:

winter <- temp.ppt %>% filter(month %in% c(1,2,12) & !water_year == 2023) %>% group_by(lon, lat, ID, water_year) %>%summarise( Precip_DecJanFeb = sum(ppt, na.rm = TRUE), 
                                                                                                                         Tmax_DecJanFeb = mean(tmax, na.rm = TRUE),
                                                                                                                         Tmin_DecJanFeb = mean(tmin, na.rm=TRUE))
colnames(winter)[4] <- "year"

winter$year <- as.numeric(winter$year)
SeptOctNov$year <- as.numeric(SeptOctNov$year)
year_means$year <- as.numeric(year_means$year)
#water_year$year <- as.character(water_year$year)
#year$year <- as.character(winter$year)
# join all of them up into a seasonal timeseries dataset
seasonal.climate <- left_join(water_year, year_means, by= c("lon", "lat", "ID", "year")) %>%
                                left_join(., winter, , by= c("lon", "lat", "ID", "year")) %>%
                                left_join(., Spring, , by= c("lon", "lat", "ID", "year")) %>%
                                left_join(., JunJulyAug, , by= c("lon", "lat", "ID", "year")) %>%
                                left_join(., SeptOctNov, , by= c("lon", "lat", "ID", "year"))

# save
ggplot(seasonal.climate, aes(x = year, y = wateryr_PPT, group = ID))+geom_line()
NE_plot$ID <- 1:length(NE_plot$cn)
seasonal.climate.IDS <- left_join(seasonal.climate, NE_plot)

saveRDS(seasonal.climate.IDS, "data/seasonal_climate_periodic_NE_plots.RDS")

# then join up all the monthly data
tmax.monthly <- temp.ppt %>% select(lon, lat, ID, tmax, year, month) %>% group_by(lon, lat, ID, year) %>%
                spread(month, tmax)
colnames(tmax.monthly)[5:16] <- paste0("tmax_", colnames(tmax.monthly)[5:16])
tmax.monthly <- tmax.monthly %>% select(lat, lon, ID, year, tmax_1, tmax_2, tmax_3, tmax_4, tmax_5, tmax_6, tmax_7, 
                       tmax_8, tmax_9, tmax_10, tmax_11, tmax_12)

tmin.monthly <- temp.ppt %>% select(lon, lat, ID, tmin, year, month) %>% group_by(lon, lat, ID, year) %>%
  spread(month, tmin)
colnames(tmin.monthly)[5:16] <- paste0("tmin_", colnames(tmin.monthly)[5:16])
tmin.monthly <- tmin.monthly %>% select(lat, lon, ID,year,  tmin_1, tmin_2, tmin_3, tmin_4, tmin_5, tmin_6, tmin_7, 
                                        tmin_8, tmin_9, tmin_10, tmin_11, tmin_12)

ppt.monthly <- temp.ppt %>% select(lon, lat, ID, ppt, year, month) %>% group_by(lon, lat, ID, year) %>%
  spread(month, ppt)
colnames(ppt.monthly)[5:16] <- paste0("ppt_", colnames(ppt.monthly)[5:16])
ppt.monthly <- ppt.monthly %>% select(lat, lon, ID,year,  ppt_1, ppt_2, ppt_3, ppt_4, ppt_5, ppt_6, ppt_7, 
                                        ppt_8, ppt_9, ppt_10, ppt_11, ppt_12)


all.monthly <- left_join(tmax.monthly, tmin.monthly ) %>% left_join(., ppt.monthly ) %>% 
  left_join(., NE_plot)

# save
saveRDS(all.monthly, "data/PRISM_monthly_NE_periodic.RDS" )
head(all.monthly)

