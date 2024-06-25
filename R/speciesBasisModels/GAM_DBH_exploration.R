library(MASS)
library(here)
library(tidyverse)
library(gt)
library(FIESTA)
library(dplyr)
library(mltools)

cleaned.data <- readRDS( "C:/Users/KellyHeilman/Box/01. kelly.heilman Workspace/mortality/Eastern-Mortality/mortality_models/data/cleaned.data.mortality.TRplots.RDS")
cleaned.data <- cleaned.data %>% dplyr::select(state, county, pltnum, cndtn, point, tree, PLOT.ID, cycle, spp, dbhcur, dbhold, damage, Species, SPCD,
                                               remper, LAT_FIADB, LONG_FIADB, elev, DIA_DIFF, annual.growth, M, relative.growth, si, physio:RD) %>% distinct()
# get summary of damages for later use:
N.DAMAGE <- cleaned.data %>% group_by(SPCD, damage) %>% summarise(n.by.damage = n())
N.DAMAGE$SPECIES <- ref_species[match(N.DAMAGE$SPCD, ref_species$SPCD),]$COMMON
ref_damage<- ref_codes %>% filter(VARIABLE %in% "AGENTCD")
N.DAMAGE$damage_agent <- ref_damage[match(N.DAMAGE$damage, ref_damage$VALUE),]$MEANING
N.DAMAGE$damage_agent <- ifelse(N.DAMAGE$damage == 0, "None", N.DAMAGE$damage_agent)
#saveRDS(N.DAMAGE, "data/N.DAMAGE.table.RDS")


nspp <- cleaned.data %>% group_by(SPCD) %>% summarise(n = n(), 
                                                      pct = n/nrow(cleaned.data)) %>% arrange (desc(`pct`))

nspp$cumulative.pct <- cumsum(nspp$pct)



# link up to the species table:
nspp$COMMON <- FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD),]$COMMON
nspp$Species <- paste(FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD),]$GENUS, FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD),]$SPECIES)

#View(nspp)

nspp[1:17,]$COMMON

library(gt)
nspp[1:17,] %>% mutate(pct = round(pct, 3), 
                       cumulative.pct = round(cumulative.pct, 3)) %>% rename(`# of trees` = "n", 
                                                                             `% of trees` = "pct",
                                                                             `cumulative %` = "cumulative.pct", 
                                                                             `Common name` = "COMMON") %>%
  dplyr::select(Species, `Common name`, SPCD, `# of trees`, `% of trees`, `cumulative %`)|> gt()
# 15 species make up >75% of the total trees in the cored plots, so lets focus on those

# only 210 tree mortality events detected in this dataset:
# M  `n()`
# <dbl>  <int>
#   1     0 315132
# 2     1    210
# Split into training and testing


cleaned.data$SPGRPCD <- FIESTA::ref_species[match(cleaned.data$SPCD, FIESTA::ref_species$SPCD),]$E_SPGRPCD

SPGRP.df <- FIESTA::ref_codes %>% filter(VARIABLE %in% "SPGRPCD") %>% filter(VALUE %in% unique(cleaned.data$SPGRPCD))
cleaned.data$SPGRPNAME <- SPGRP.df[match(cleaned.data$SPGRPCD, SPGRP.df$VALUE),]$MEANING


#View(cleaned.data %>% filter(SPCD %in% unique(nspp[1:17,]$SPCD))%>% group_by( SPGRPNAME, SPCD) %>% summarise(n()))
# next to 97 (red spruce), 241 (white ceder), 531 (fagus grandifolia), 
# select species 318--red maple

# center and scale the covariate data
# for covariates at the plot level, scale by the unique plots so the # of trees on the plot doesnt affect the mean and sd values:
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


# for each species group, fit a model, plot the outputs, and save the results
# we source a function from another script
#source("R/SPCD_run_stan.R")

# this runs a stan model and saves the outputs
# SPGRPCD 2 throws uncerialize socklist error--too big of a diataset to parallelisze?
cleaned.data.full %>% group_by(SPCD) %>% summarise(n())
SPCD.df <- data.frame(SPCD = nspp[1:17, ]$SPCD, 
                      spcd.id = 1:17)
remper.cor.vector <- c(  0.5)
model.number <- 1
model.list <- 6
m <- 1

for(i in 1:17){# run for each of the 17 species
  common.name <- nspp[1:17, ] %>% filter(SPCD %in% SPCD.df[i,]$SPCD) %>% dplyr::select(COMMON)
  SPCD.id <- nspp[i,]$SPCD
  #for(m in 1:length(model.list)){  # run each of the 9 models
  model.number <- model.list[m]
  for (j in 1:length(remper.cor.vector)){ # for the gowoth only model explore the consequences of other assumptions about remeasurement period
    cat(paste("running gam mortality models for SPCD", nspp[i,]$SPCD, common.name$COMMON, " remper correction", remper.cor.vector[j]))
    model.no  <- model.number
    remper.correction <- remper.cor.vector[j]
    load(paste0("C:/Users/KellyHeilman/Box/01. kelly.heilman Workspace/mortality/Eastern-Mortality/mortality_models/SPCD_standata_general_full/SPCD_",SPCD.id, "remper_correction_", remper.correction,"model_",model.no, ".Rdata")) # load the species code data
    
    model.df <- data.frame(S = mod.data$y, 
                           mod.data$xM 
    )
    m.1 <- gam(formula = S ~ annual.growth.scaled , data =  model.df, family = "binomial")
    
    m.dbh.1 <- gam(formula = S ~ annual.growth.scaled + s(DIA_scaled, bs = "bs"), data =  model.df, family = "binomial")
    m.dbh.2 <- gam(formula = S ~ annual.growth.scaled + s(DIA_scaled, bs = "bs") + RD.scaled + ba.scaled + 
                     BAL.scaled + damage.scaled, data =  model.df, family = "binomial")
    m.dbh.3 <- gam(formula = S ~  annual.growth.scaled + s(DIA_scaled, bs = "bs") + RD.scaled + ba.scaled + 
                     BAL.scaled + damage.scaled + 
                     MATmax.scaled + MATmin.scaled + MAP.scaled + 
                     ppt.anom + tmin.anom + tmax.anom, data =  model.df, family = "binomial")
    m.dbh.4 <- gam(formula = S ~  annual.growth.scaled + s(DIA_scaled, bs = "bs") + RD.scaled + ba.scaled + 
                     BAL.scaled + damage.scaled + 
                     MATmax.scaled + MATmin.scaled + MAP.scaled + 
                     ppt.anom + tmin.anom + tmax.anom +
                     slope.scaled + aspect.scaled + elev.scaled + Ndep.scaled + physio.scaled, data =  model.df, family = "binomial")
    m.dbh.5 <- gam(formula = S ~  annual.growth.scaled + s(DIA_scaled, bs = "bs") + RD.scaled + ba.scaled + 
                     BAL.scaled + damage.scaled + 
                     MATmax.scaled + MATmin.scaled + MAP.scaled + 
                     ppt.anom + tmin.anom + tmax.anom + 
                     slope.scaled + aspect.scaled + elev.scaled + Ndep.scaled + physio.scaled + 
                     DIA_scaled_growth.int + 
                     RD.scaled_growth.int     +  ba.scaled_growth.int    +   BAL.scaled_growth.int    +
                     damage.scaled_growth.int +  MATmax.scaled_growth.int +  MATmin.scaled_growth.int +
                     MAP.scaled_growth.int    +  ppt.anom_growth.int      +  tmin.anom_growth.int     +
                     tmax.anom_growth.int     +  slope.scaled_growth.int  +  aspect.scaled_growth.int +
                     elev.scaled_growth.int   + Ndep.scaled_growth.int    + physio.scaled_growth.int +
                     RD.scaled_DIA.int        +  ba.scaled_DIA.int        +  BAL.scaled_DIA.int     +
                     damage.scaled_DIA.int  +   MATmax.scaled_DIA.int     + MATmin.scaled_DIA.int +   
                     MAP.scaled_DIA.int       +  ppt.anom_DIA.int         +  tmin.anom_DIA.int  +      
                     tmax.anom_DIA.int        +  slope.scaled_DIA.int     +  aspect.scaled_DIA.int+    
                     elev.scaled_DIA.int      +  Ndep.scaled_DIA.int      +  physio.scaled_DIA.int   , data =  model.df, family =  binomial )
    
    png(height = 5, width = 9, units = "in", res = 100, paste0("C:/Users/KellyHeilman/Box/01. kelly.heilman Workspace/mortality/Eastern-Mortality/mortality_models/images/Sample_diameter_survival_gams_cr_", SPCD.id, ".png"))
    par(mfrow =c(2,3))
    plot(m.dbh.1, rug = TRUE, main = paste("model 2, df =", df = m.dbh.1$smooth[[1]]$bs.dim - 1))
    plot(m.dbh.2, rug = TRUE, main = paste("model 3, df =", df = m.dbh.2$smooth[[1]]$bs.dim - 1))
    plot(m.dbh.3, rug = TRUE, main = paste("model 4, df =", df = m.dbh.3$smooth[[1]]$bs.dim - 1))
    plot(m.dbh.4, rug = TRUE, main = paste("model 5, df =", df = m.dbh.4$smooth[[1]]$bs.dim - 1))
    plot(m.dbh.5, rug = TRUE, main = paste("model 6, df =", df = m.dbh.5$smooth[[1]]$bs.dim - 1))
    dev.off()    
    
    
    base::save(m.1, m.dbh.1, m.dbh.2, m.dbh.3, m.dbh.4, 
               m.dbh.5, file = paste0("C:/Users/KellyHeilman/Box/01. kelly.heilman Workspace/mortality/Eastern-Mortality/mortality_models/dbh_smooth_gams/BDBH_GAM_models_SPCD_cr",SPCD.id,".Rdata"))
    
  }
}

