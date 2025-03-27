#library(rstan)
library(MASS)
library(here)
library(tidyverse)
library(gt)
library(FIESTA)
library(dplyr)
library(ggcorrplot)
library(car)

options(mc.cores = parallel::detectCores())
#saveRDS(cleaned.data2, "data/cleaned.data.mortality.TRplots.RDS")

cleaned.data <- readRDS( "data/cleaned.data.mortality.TRplots.RDS")

# get summary of damages for later use:
N.DAMAGE <- cleaned.data %>% group_by(SPCD, damage) %>% summarise(n.by.damage = n())
N.DAMAGE$SPECIES <- ref_species[match(N.DAMAGE$SPCD, ref_species$SPCD),]$COMMON
ref_damage<- ref_codes %>% filter(VARIABLE %in% "AGENTCD")
N.DAMAGE$damage_agent <- ref_damage[match(N.DAMAGE$damage, ref_damage$VALUE),]$MEANING
N.DAMAGE$damage_agent <- ifelse(N.DAMAGE$damage == 0, "None", N.DAMAGE$damage_agent)
saveRDS(N.DAMAGE, "data/N.DAMAGE.table.RDS")


nspp <- cleaned.data %>% group_by(SPCD) %>% summarise(n = n(), 
                                                      pct = n/nrow(cleaned.data)) %>% arrange (desc(`pct`))

nspp$cumulative.pct <- cumsum(nspp$pct)



# link up to the species table:
nspp$COMMON <- FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD),]$COMMON

View(nspp)

nspp[1:17,]$COMMON

library(gt)
nspp[1:17,] |> gt()
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


View(cleaned.data %>% filter(SPCD %in% unique(nspp[1:17,]$SPCD))%>% group_by( SPGRPNAME, SPCD) %>% summarise(n()))
# next to 97 (red spruce), 241 (white ceder), 531 (fagus grandifolia), 
# select species 318--red maple

# center and scale the covariate data
# for covariates at the plot level, scale by the unique plots so the # of trees on the plot doesnt affect the mean and sd values:
plot.means <- unique(cleaned.data %>% ungroup()%>% 
                       dplyr::select(PLOT.ID, si, ba, slope, aspect, MAP, MATmin, MATmax, damage.total, elev, Ndep.remper.avg)) %>% 
  ungroup() %>% summarise(si.mean = mean(si, na.rm =TRUE), 
                          ba.mean = mean(ba, na.rm =TRUE), 
                          slope.mean = mean(slope, na.rm = TRUE), 
                          aspect.mean = mean(aspect, na.rm = TRUE),
                          damage.mean = mean(damage.total, na.rm =TRUE),
                          elev.mean = mean(elev, na.rm =TRUE),
                          Ndep.mean = mean(Ndep.remper.avg, na.rm =TRUE),
                          
                          MAP.mean = mean(MAP, na.rm =TRUE), 
                          MATmin.mean = mean(MATmin, na.rm =TRUE), 
                          MATmax.mean = mean(MATmax, na.rm =TRUE), 
                          
                          ba.sd = sd(ba, na.rm =TRUE),
                          si.sd = sd(si, na.rm =TRUE), 
                          slope.sd = sd(slope, na.rm =TRUE),
                          aspect.sd = sd(aspect, na.rm = TRUE),
                          damage.sd = sd(damage.total, na.rm =TRUE),
                          elev.sd = sd(elev, na.rm =TRUE),
                          Ndep.sd = sd(Ndep.remper.avg, na.rm =TRUE),
                          
                          MAP.sd = sd(MAP, na.rm =TRUE), 
                          MATmin.sd = sd(MATmin, na.rm =TRUE), 
                          MATmax.sd = sd(MATmax, na.rm =TRUE)
  )

View(cleaned.data %>% group_by(SPGRPCD, SPCD) %>% summarise(n()))

cleaned.data.full <- cleaned.data

#-----------------------------------------------------------------------------------------
# Make the species level datasets for the top 17 species to run the model
#-----------------------------------------------------------------------------------------
# note: this now takes in the data built in prep_species_stan_data.R and remakes it for the glms
length(unique(cleaned.data$SPCD))

nspp[1:17,]$COMMON
# save these as .RDA files so we can just load, run the model, and 
SPCD.id <- 316#unique(cleaned.data$SPCD)[25]
set.seed(22)
remper.correction <- 0.5
stan2glm.data <- function(SPCD.id, remper.correction){
  
  load( paste0("SPCD_standata_general_full_standardized_v2/SPCD_",SPCD.id,"remper_correction_",remper.correction,"model_9.Rdata"))
  
  mod.data <- list(N = nrow(train.data), 
                   Nspp = length(unique(train.data$SPCD)),
                   SPP = train.data$SPP,
                   y = train.data$M,
                   si = train.data$si.scaled, 
                   slope = train.data$slope.scaled, 
                   aspect = train.data$aspect.scaled,
                   elev = train.data$elev.scaled, 
                   Ndep = train.data$Ndep.scaled,
                   DIA = train.data$DIA_scaled, 
                   MATmax = train.data$MATmax.scaled, 
                   MATmin= train.data$MATmin.scaled, 
                   MAP = train.data$MAP.scaled, 
                   annual_growth = train.data$annual.growth.scaled,
                   BAL = train.data$BAL.scaled, 
                   BA = train.data$ba.scaled, 
                   damage = train.data$damage.total,
                   MAPanom = train.data$ppt.anom, 
                   MATminanom = train.data$tmin.anom, 
                   MATmaxanom = train.data$tmax.anom,
                   PHYSIO = train.data$physio, 
                   RD = train.data$RD, 
                   SPCD.BA = train.data$SPCD.BA.scaled,
                   non_SPCD.BA.scaled = train.data$non_SPCD.BA.scaled,
                   prop.focal.ba.scaled = train.data$prop.focal.ba.scaled, 
                   DIA.diff  = train.data$DIA_DIFF,
                   
                   
                   N_test = nrow(test.data),
                   si_test = test.data$si.scaled, 
                   DIA_test = test.data$DIA_scaled, 
                   aspect_test = test.data$aspect.scaled, 
                   elev_test = test.data$elev.scaled, 
                   Ndep_test = test.data$Ndep.scaled,
                   slope_test = test.data$slope.scaled,
                   BAL_test = test.data$BAL.scaled,
                   BA_test = test.data$ba.scaled,
                   damage_test = test.data$damage.scaled, 
                   
                   MATmax_test = test.data$MATmax.scaled, 
                   MATmin_test = test.data$MATmin.scaled, 
                   MAP_test = test.data$MAP.scaled,
                   annual_growth_test = test.data$annual.growth.scaled, 
                   MAPanom_test = test.data$ppt.anom, 
                   MATminanom_test = test.data$tmin.anom, 
                   MATmaxanom_test = test.data$tmax.anom, 
                   PHYSIO_test = test.data$physio, 
                   RD_test = test.data$RD, 
                   SPCD.BA_test = test.data$SPCD.BA.scaled,
                   non_SPCD.BA.scaled_test = test.data$non_SPCD.BA.scaled,
                   prop.focal.ba.scaled_test = test.data$prop.focal.ba.scaled, 
                   DIA.diff_test  = test.data$DIA_DIFF)
  
  
  model.name <- paste0("simple_logistic_SPCD_", SPCD.id, "remper_",remper.correction)
  
  save(train.data, test.data, mod.data, model.name, file = paste0("SPCD_GLM_standata/SPCD_",SPCD.id,"remper_correction_",remper.correction, ".Rdata"))
}

# write the data for all 17 different species:
for(i in 1:length(unique(nspp[1:17,]$SPCD))){
  cat(i)
  stan2glm.data(SPCD.id = nspp[i,]$SPCD, remper.correction = 0.5)
}

#---------------------------------------------------------------------
# variable selection with glm models
#---------------------------------------------------------------------
#for(i in 1:15){# unique(cleaned.data.full$SPGRPCD)){
SPCD.df <- data.frame(SPCD = nspp[1:17,]$SPCD, 
                      spcd.id = 1:17)
remper.cor.vector <- c(0.5)
j <- 1
i <- 1

for(i in 1:length(SPCD.df$SPCD)){
  SPCD.id <- SPCD.df[i,]$SPCD
  common.name <- nspp[1:17,] %>% filter(SPCD %in% SPCD.id) %>% dplyr::select(COMMON)
  
  #if(SPCD.id == 621){
  # cat("Not running for yellow poplar, not enough data")
  #}else{
  
  for (j in 1:length(remper.cor.vector)){
    cat(paste("running glm mortality model for SPCD", SPCD.df[i,]$SPCD, common.name$COMMON, " remper correction", remper.cor.vector[j]))
    
    remper.correction <- remper.cor.vector[j]
    load(paste0("SPCD_GLM_standata/SPCD_", SPCD.id, "remper_correction_", remper.correction, ".Rdata")) # load the species code data
    covariate.data <- data.frame(M = mod.data$y, 
                                 annual.growth = mod.data$annual_growth, 
                                 dia.diff = mod.data$DIA.diff,
                                 DIA = mod.data$DIA, 
                                 si = mod.data$si, 
                                 slope = mod.data$slope, 
                                 aspect = mod.data$aspect, 
                                 MATmax = mod.data$MATmax, 
                                 MATmin= mod.data$MATmin, 
                                 MAP = mod.data$MAP, 
                                 BAL = mod.data$BAL, 
                                 damage = mod.data$damage, 
                                 MAPanom = mod.data$MAPanom, 
                                 MATmaxanom = mod.data$MATmaxanom, 
                                 MATminanom = mod.data$MATminanom, 
                                 PHYSIO = mod.data$PHYSIO, 
                                 BA = mod.data$BA, 
                                 RD = mod.data$RD, 
                                 elev = mod.data$elev, 
                                 Ndep = mod.data$Ndep, 
                                 SPCD.BA = mod.data$SPCD.BA,
                                 non_SPCD.BA = mod.data$non_SPCD.BA.scaled,
                                 prop.focal.ba = mod.data$prop.focal.ba.scaled 
    )
    
    test.covariate.data <- data.frame(M = test.data$M,
                                      annual.growth = mod.data$annual_growth_test, 
                                      dia.diff = mod.data$DIA.diff_test,
                                      DIA = mod.data$DIA_test, 
                                      si = mod.data$si_test, 
                                      slope = mod.data$slope_test, 
                                      aspect = mod.data$aspect_test, 
                                      MATmax = mod.data$MATmax_test, 
                                      MATmin= mod.data$MATmin_test, 
                                      MAP = mod.data$MAP_test, 
                                      BAL = mod.data$BAL_test, 
                                      damage = mod.data$damage_test, 
                                      MAPanom = mod.data$MAPanom_test, 
                                      MATmaxanom = mod.data$MATmaxanom_test, 
                                      MATminanom = mod.data$MATminanom_test, 
                                      PHYSIO = mod.data$PHYSIO_test, 
                                      BA = mod.data$PHYSIO_test, 
                                      RD = mod.data$RD_test, 
                                      elev = mod.data$elev_test, 
                                      Ndep = mod.data$Ndep_test, 
                                      SPCD.BA = mod.data$SPCD.BA_test,
                                      non_SPCD.BA = mod.data$non_SPCD.BA.scaled_test,
                                      prop.focal.ba = mod.data$prop.focal.ba.scaled_test )
    
    
    
    glm.A <-  glm(M ~ MAP , data = covariate.data, family = "binomial")
    glm.B <-  glm(M ~ MATmaxanom, data = covariate.data, family = "binomial")
    glm.C <-  glm(M ~ MATminanom, data = covariate.data, family = "binomial")
    glm.D <-  glm(M ~ MAPanom, data = covariate.data, family = "binomial")
    glm.E <-  glm(M ~ BAL, data = covariate.data, family = "binomial")
    glm.F <-  glm(M ~ damage , data = covariate.data, family = "binomial")
    glm.G <-  glm(M ~ si, data = covariate.data, family = "binomial")
    glm.H <-  glm(M ~ PHYSIO, data = covariate.data, family = "binomial")
    glm.I <-  glm(M ~ BA, data = covariate.data, family = "binomial")
    glm.J <-  glm(M ~ RD , data = covariate.data, family = "binomial")
    glm.K <-  glm(M ~ elev, data = covariate.data, family = "binomial")
    glm.L <-  glm(M ~ Ndep, data = covariate.data, family = "binomial")
    glm.M <-  glm(M ~ SPCD.BA, data = covariate.data, family = "binomial")
    glm.N <-  glm(M ~ non_SPCD.BA, data = covariate.data, family = "binomial")
    glm.O <-  glm(M ~ prop.focal.ba, data = covariate.data, family = "binomial")
    
    
    
    glm.1 <-  glm(M ~ dia.diff, data = covariate.data, family = "binomial")
    glm.2 <-  glm(M ~ dia.diff + DIA, data = covariate.data, family = "binomial")
    glm.3 <-  glm(M ~ dia.diff + DIA + exp(DIA), data = covariate.data, family = "binomial")
    
    glm.4 <-  glm(M ~ dia.diff + DIA + exp(DIA)+ aspect, data = covariate.data, family = "binomial")
    glm.5 <-  glm(M ~ dia.diff + DIA + exp(DIA)+ aspect + slope , data = covariate.data, family = "binomial")
    glm.6 <-  glm(M ~ dia.diff + DIA + exp(DIA) + aspect + slope + MATmax , data = covariate.data, family = "binomial")
    
    glm.7 <-  glm(M ~ dia.diff + DIA + exp(DIA)+ aspect + slope + MATmax 
                  + MATmin , data = covariate.data, family = "binomial")
    glm.8 <-  glm(M ~ dia.diff + DIA + exp(DIA)+ aspect + slope + MATmax 
                  + MATmin + MATmax , data = covariate.data, family = "binomial")
    glm.9 <-  glm(M ~ dia.diff + DIA + exp(DIA) + aspect + slope + MATmax 
                  + MATmin + MAP , data = covariate.data, family = "binomial")
    
    glm.10 <-  glm(M ~ dia.diff + DIA + exp(DIA)+ aspect + slope + MATmax 
                   + MATmin + MAP + MATmaxanom , data = covariate.data, family = "binomial")
    
    glm.11 <-  glm(M ~ dia.diff + DIA + exp(DIA)+ aspect + slope + MATmax 
                   + MATmin+ MAP + MATmaxanom 
                   + MATminanom , data = covariate.data, family = "binomial")
    
    
    glm.12 <-  glm(M ~ dia.diff + DIA + exp(DIA)+ aspect + slope + MATmax 
                   + MATmin+ MAP + MATmaxanom 
                   + MATminanom + MAPanom , data = covariate.data, family = "binomial")
    
    glm.13 <-  glm(M ~ dia.diff + DIA + exp(DIA)+ aspect + slope + MATmax 
                   + MATmin+ MAP + MATmaxanom 
                   + MATminanom + MAPanom + BAL , data = covariate.data, family = "binomial")
    
    glm.14<-  glm(M ~ dia.diff + DIA + exp(DIA)+ aspect + slope + MATmax 
                  + MATmin+ MAP + MATmaxanom 
                  + MATminanom + MAPanom + BAL + damage, data = covariate.data, family = "binomial")
    
    glm.15 <-  glm(M ~ dia.diff + DIA + exp(DIA)+ aspect + slope + MATmax 
                   + MATmin+ MAP + MATmaxanom 
                   + MATminanom + MAPanom + BAL + damage + si, data = covariate.data, family = "binomial")
    
    glm.16 <-  glm(M ~ dia.diff + DIA + exp(DIA)+ aspect + slope + MATmax 
                   + MATmin+ MAP + MATmaxanom 
                   + MATminanom + MAPanom + BAL + damage + si + PHYSIO, data = covariate.data, family = "binomial")
    
    glm.17 <-  glm(M ~ dia.diff + DIA + exp(DIA)+ aspect + slope + MATmax 
                   + MATmin+ MAP + MATmaxanom 
                   + MATminanom + MAPanom + BAL + damage + si +  PHYSIO + BA , data = covariate.data, family = "binomial")
    
    glm.18 <-  glm(M ~ dia.diff + DIA + exp(DIA)+ aspect + slope + MATmax 
                   + MATmin+ MAP + MATmaxanom 
                   + MATminanom + MAPanom + BAL + damage + si +  PHYSIO + BA + RD , data = covariate.data, family = "binomial")
    
    glm.19 <-  glm(M ~ dia.diff + DIA + exp(DIA)+ aspect + slope + MATmax 
                   + MATmin+ MAP + MATmaxanom 
                   + MATminanom + MAPanom + BAL + damage + si +  PHYSIO + BA + RD +
                     elev, data = covariate.data, family = "binomial")
    
    glm.20 <-  glm(M ~ dia.diff + DIA + exp(DIA)+ aspect + slope + MATmax 
                   + MATmin+ MAP + MATmaxanom 
                   + MATminanom + MAPanom + BAL + damage + si +  PHYSIO + BA + RD +
                     elev + Ndep, data = covariate.data, family = "binomial")
    
    # add in the 3 additional variables
    glm.20.a <-  glm(M ~ dia.diff + DIA + exp(DIA)+ aspect + slope + MATmax 
                     + MATmin+ MAP + MATmaxanom 
                     + MATminanom + MAPanom + BAL + damage + si +  PHYSIO + BA + RD +
                       elev + Ndep + SPCD.BA, data = covariate.data, family = "binomial")
    
    glm.20.b <-  glm(M ~ dia.diff + DIA + exp(DIA)+ aspect + slope + MATmax 
                     + MATmin+ MAP + MATmaxanom 
                     + MATminanom + MAPanom + BAL + damage + si +  PHYSIO + BA + RD +
                       elev + Ndep + SPCD.BA + non_SPCD.BA, data = covariate.data, family = "binomial")
    
    glm.20.c <-  glm(M ~ dia.diff + DIA + exp(DIA)+ aspect + slope + MATmax 
                     + MATmin+ MAP + MATmaxanom 
                     + MATminanom + MAPanom + BAL + damage + si +  PHYSIO + BA + RD +
                       elev + Ndep + SPCD.BA + non_SPCD.BA + prop.focal.ba, data = covariate.data, family = "binomial")
    
    
    glm.21 <-  glm(M ~ dia.diff + DIA + exp(DIA)+ aspect + slope + MATmax 
                   + MATmin+ MAP + MATmaxanom 
                   + MATminanom + MAPanom + BAL + damage + si +  PHYSIO + BA + RD + elev + Ndep + 
                     SPCD.BA + non_SPCD.BA + prop.focal.ba +
                     # annual growth interactcovarions
                     DIA*dia.diff + aspect*dia.diff  + slope*dia.diff  + MATmax*dia.diff  
                   + MATmin*dia.diff  + MAP*dia.diff  + MATmaxanom*dia.diff  
                   + MATminanom*dia.diff  + MAPanom*dia.diff  + BAL*dia.diff  + damage*dia.diff  + si*dia.diff +
                     +  PHYSIO*dia.diff + BA*dia.diff + RD*dia.diff + elev*dia.diff + Ndep*dia.diff +
                     SPCD.BA*dia.diff + non_SPCD.BA*dia.diff + prop.focal.ba*dia.diff , data = covariate.data, family = "binomial")
    
    
    
    # all growth + diameter interaction terms
    glm.22 <-  glm(M ~ dia.diff + DIA + exp(DIA)+ aspect + slope + MATmax 
                   + MATmin+ MAP + MATmaxanom 
                   + MATminanom + MAPanom + BAL + damage + si +  PHYSIO + BA + RD + elev + Ndep + 
                     SPCD.BA + non_SPCD.BA + prop.focal.ba +
                     # annual growth interactcovarions
                     DIA*dia.diff + aspect*dia.diff  + slope*dia.diff  + MATmax*dia.diff  
                   + MATmin*dia.diff  + MAP*dia.diff  + MATmaxanom*dia.diff  
                   + MATminanom*dia.diff  + MAPanom*dia.diff  + BAL*dia.diff  + damage*dia.diff  + si*dia.diff +
                     +  PHYSIO*dia.diff + BA*dia.diff + RD*dia.diff + elev*dia.diff + Ndep*dia.diff +
                     SPCD.BA*dia.diff + non_SPCD.BA*dia.diff + prop.focal.ba*dia.diff +
                     # all diameter interactions
                     aspect*DIA  + slope*DIA  + MATmax*DIA  
                   + MATmin*DIA  + MAP*DIA  + MATmaxanom*DIA  
                   + MATminanom*DIA  + MAPanom*DIA  + BAL*DIA  + damage*DIA  + si*DIA + 
                     PHYSIO*DIA + BA*DIA + RD*DIA + elev*DIA + Ndep*DIA +
                     SPCD.BA*DIA + non_SPCD.BA*DIA + prop.focal.ba*DIA , data = covariate.data, family = "binomial")
    
    glm.23 <-  glm(M ~ dia.diff + DIA + exp(DIA)+ aspect + slope + MATmax 
                   + MATmin+ MAP + MATmaxanom 
                   + MATminanom + MAPanom + BAL + damage + si +  PHYSIO + BA + RD + elev + Ndep + 
                     SPCD.BA + non_SPCD.BA + prop.focal.ba +
                     # annual growth interactcovarions
                     DIA*dia.diff + aspect*dia.diff  + slope*dia.diff  + MATmax*dia.diff  
                   + MATmin*dia.diff  + MAP*dia.diff  + MATmaxanom*dia.diff  
                   + MATminanom*dia.diff  + MAPanom*dia.diff  + BAL*dia.diff  + damage*dia.diff  + si*dia.diff +
                     +  PHYSIO*dia.diff + BA*dia.diff + RD*dia.diff + elev*dia.diff + Ndep*dia.diff +
                     SPCD.BA*dia.diff + non_SPCD.BA*dia.diff + prop.focal.ba*dia.diff +
                     # all diameter interactions
                     aspect*DIA  + slope*DIA  + MATmax*DIA  
                   + MATmin*DIA  + MAP*DIA  + MATmaxanom*DIA  
                   + MATminanom*DIA  + MAPanom*DIA  + BAL*DIA  + damage*DIA  + si*DIA + 
                     PHYSIO*DIA + BA*DIA + RD*DIA + elev*DIA + Ndep*DIA +
                     SPCD.BA*DIA + non_SPCD.BA*DIA + prop.focal.ba*DIA + 
                     # all aspect interactions
                     slope*aspect  + MATmax*aspect  
                   + MATmin*aspect  + MAP*aspect  + MATmaxanom*aspect  
                   + MATminanom*aspect  + MAPanom*aspect  + BAL*aspect  + damage*aspect  + si*aspect +
                     PHYSIO*aspect + BA*aspect + RD*aspect + elev*aspect + Ndep*aspect+
                     SPCD.BA*aspect + non_SPCD.BA*aspect + prop.focal.ba*aspect, data = covariate.data, family = "binomial")
    
    
    
    glm.24 <-  glm(M ~ dia.diff + DIA + exp(DIA)+ aspect + slope + MATmax 
                   + MATmin+ MAP + MATmaxanom 
                   + MATminanom + MAPanom + BAL + damage + si +  PHYSIO + BA + RD + elev + Ndep + 
                     SPCD.BA + non_SPCD.BA + prop.focal.ba +
                     # annual growth interactcovarions
                     DIA*dia.diff + aspect*dia.diff  + slope*dia.diff  + MATmax*dia.diff  
                   + MATmin*dia.diff  + MAP*dia.diff  + MATmaxanom*dia.diff  
                   + MATminanom*dia.diff  + MAPanom*dia.diff  + BAL*dia.diff  + damage*dia.diff  + si*dia.diff +
                     +  PHYSIO*dia.diff + BA*dia.diff + RD*dia.diff + elev*dia.diff + Ndep*dia.diff +
                     SPCD.BA*dia.diff + non_SPCD.BA*dia.diff + prop.focal.ba*dia.diff +
                     # all diameter interactions
                     aspect*DIA  + slope*DIA  + MATmax*DIA  
                   + MATmin*DIA  + MAP*DIA  + MATmaxanom*DIA  
                   + MATminanom*DIA  + MAPanom*DIA  + BAL*DIA  + damage*DIA  + si*DIA + 
                     PHYSIO*DIA + BA*DIA + RD*DIA + elev*DIA + Ndep*DIA +
                     SPCD.BA*DIA + non_SPCD.BA*DIA + prop.focal.ba*DIA + 
                     # all aspect interactions
                     slope*aspect  + MATmax*aspect  
                   + MATmin*aspect  + MAP*aspect  + MATmaxanom*aspect  
                   + MATminanom*aspect  + MAPanom*aspect  + BAL*aspect  + damage*aspect  + si*aspect +
                     PHYSIO*aspect + BA*aspect + RD*aspect + elev*aspect + Ndep*aspect+
                     SPCD.BA*aspect + non_SPCD.BA*aspect + prop.focal.ba*aspect +
                     # all slope interactions
                     MATmax*slope  
                   + MATmin*slope  + MAP*slope  + MATmaxanom*slope  
                   + MATminanom*slope  + MAPanom*slope  + BAL*slope  + damage*slope  + si*slope + 
                     PHYSIO*slope + BA*slope + RD*slope + elev*slope + Ndep*slope +
                     SPCD.BA*slope + non_SPCD.BA*slope + prop.focal.ba*slope, data = covariate.data, family = "binomial")
    
    glm.25 <-  glm(M ~ dia.diff + DIA + exp(DIA)+ aspect + slope + MATmax 
                   + MATmin+ MAP + MATmaxanom 
                   + MATminanom + MAPanom + BAL + damage + si +  PHYSIO + BA + RD + elev + Ndep + 
                     SPCD.BA + non_SPCD.BA + prop.focal.ba +
                     # annual growth interactcovarions
                     DIA*dia.diff + aspect*dia.diff  + slope*dia.diff  + MATmax*dia.diff  
                   + MATmin*dia.diff  + MAP*dia.diff  + MATmaxanom*dia.diff  
                   + MATminanom*dia.diff  + MAPanom*dia.diff  + BAL*dia.diff  + damage*dia.diff  + si*dia.diff +
                     +  PHYSIO*dia.diff + BA*dia.diff + RD*dia.diff + elev*dia.diff + Ndep*dia.diff +
                     SPCD.BA*dia.diff + non_SPCD.BA*dia.diff + prop.focal.ba*dia.diff +
                     # all diameter interactions
                     aspect*DIA  + slope*DIA  + MATmax*DIA  
                   + MATmin*DIA  + MAP*DIA  + MATmaxanom*DIA  
                   + MATminanom*DIA  + MAPanom*DIA  + BAL*DIA  + damage*DIA  + si*DIA + 
                     PHYSIO*DIA + BA*DIA + RD*DIA + elev*DIA + Ndep*DIA +
                     SPCD.BA*DIA + non_SPCD.BA*DIA + prop.focal.ba*DIA + 
                     # all aspect interactions
                     slope*aspect  + MATmax*aspect  
                   + MATmin*aspect  + MAP*aspect  + MATmaxanom*aspect  
                   + MATminanom*aspect  + MAPanom*aspect  + BAL*aspect  + damage*aspect  + si*aspect +
                     PHYSIO*aspect + BA*aspect + RD*aspect + elev*aspect + Ndep*aspect+
                     SPCD.BA*aspect + non_SPCD.BA*aspect + prop.focal.ba*aspect +
                     # all slope interactions
                     MATmax*slope  
                   + MATmin*slope  + MAP*slope  + MATmaxanom*slope  
                   + MATminanom*slope  + MAPanom*slope  + BAL*slope  + damage*slope  + si*slope + 
                     PHYSIO*slope + BA*slope + RD*slope + elev*slope + Ndep*slope +
                     SPCD.BA*slope + non_SPCD.BA*slope + prop.focal.ba*slope+ 
                     # all MATmax interactions
                     
                     MATmin*MATmax  + MAP*MATmax  + MATmaxanom*MATmax  
                   + MATminanom*MATmax  + MAPanom*MATmax  + BAL*MATmax  + damage*MATmax  + si*MATmax +
                     PHYSIO*MATmax + BA*MATmax + RD*MATmax + elev*MATmax +  Ndep*MATmax +
                     SPCD.BA*MATmax + non_SPCD.BA*MATmax + prop.focal.ba*MATmax, data = covariate.data, family = "binomial")
    
    
    glm.26 <-  glm(M ~ dia.diff + DIA + exp(DIA)+ aspect + slope + MATmax 
                   + MATmin+ MAP + MATmaxanom 
                   + MATminanom + MAPanom + BAL + damage + si +  PHYSIO + BA + RD + elev + Ndep + 
                     SPCD.BA + non_SPCD.BA + prop.focal.ba +
                     # annual growth interactcovarions
                     DIA*dia.diff + aspect*dia.diff  + slope*dia.diff  + MATmax*dia.diff  
                   + MATmin*dia.diff  + MAP*dia.diff  + MATmaxanom*dia.diff  
                   + MATminanom*dia.diff  + MAPanom*dia.diff  + BAL*dia.diff  + damage*dia.diff  + si*dia.diff +
                     +  PHYSIO*dia.diff + BA*dia.diff + RD*dia.diff + elev*dia.diff + Ndep*dia.diff +
                     SPCD.BA*dia.diff + non_SPCD.BA*dia.diff + prop.focal.ba*dia.diff +
                     # all diameter interactions
                     aspect*DIA  + slope*DIA  + MATmax*DIA  
                   + MATmin*DIA  + MAP*DIA  + MATmaxanom*DIA  
                   + MATminanom*DIA  + MAPanom*DIA  + BAL*DIA  + damage*DIA  + si*DIA + 
                     PHYSIO*DIA + BA*DIA + RD*DIA + elev*DIA + Ndep*DIA +
                     SPCD.BA*DIA + non_SPCD.BA*DIA + prop.focal.ba*DIA + 
                     # all aspect interactions
                     slope*aspect  + MATmax*aspect  
                   + MATmin*aspect  + MAP*aspect  + MATmaxanom*aspect  
                   + MATminanom*aspect  + MAPanom*aspect  + BAL*aspect  + damage*aspect  + si*aspect +
                     PHYSIO*aspect + BA*aspect + RD*aspect + elev*aspect + Ndep*aspect+
                     SPCD.BA*aspect + non_SPCD.BA*aspect + prop.focal.ba*aspect +
                     # all slope interactions
                     MATmax*slope  
                   + MATmin*slope  + MAP*slope  + MATmaxanom*slope  
                   + MATminanom*slope  + MAPanom*slope  + BAL*slope  + damage*slope  + si*slope + 
                     PHYSIO*slope + BA*slope + RD*slope + elev*slope + Ndep*slope +
                     SPCD.BA*slope + non_SPCD.BA*slope + prop.focal.ba*slope+ 
                     # all MATmax interactions
                     
                     MATmin*MATmax  + MAP*MATmax  + MATmaxanom*MATmax  
                   + MATminanom*MATmax  + MAPanom*MATmax  + BAL*MATmax  + damage*MATmax  + si*MATmax +
                     PHYSIO*MATmax + BA*MATmax + RD*MATmax + elev*MATmax +  Ndep*MATmax +
                     SPCD.BA*MATmax + non_SPCD.BA*MATmax + prop.focal.ba*MATmax+ 
                     # all MATmin interactions
                     
                     MAP*MATmin + MATminanom*MATmin
                   + MATminanom*MATmin + MAPanom*MATmin + BAL*MATmin + damage*MATmin + si*MATmin +
                     PHYSIO*MATmin + BA*MATmin + RD*MATmin + elev*MATmin + Ndep*MATmin +
                     SPCD.BA*MATmin + non_SPCD.BA*MATmin + prop.focal.ba*MATmin, data = covariate.data, family = "binomial")
    
    
    
    glm.27 <-  glm(M ~ dia.diff + DIA + exp(DIA)+ aspect + slope + MATmax 
                   + MATmin+ MAP + MATmaxanom 
                   + MATminanom + MAPanom + BAL + damage + si +  PHYSIO + BA + RD + elev + Ndep + 
                     SPCD.BA + non_SPCD.BA + prop.focal.ba +
                     # annual growth interactcovarions
                     DIA*dia.diff + aspect*dia.diff  + slope*dia.diff  + MATmax*dia.diff  
                   + MATmin*dia.diff  + MAP*dia.diff  + MATmaxanom*dia.diff  
                   + MATminanom*dia.diff  + MAPanom*dia.diff  + BAL*dia.diff  + damage*dia.diff  + si*dia.diff +
                     +  PHYSIO*dia.diff + BA*dia.diff + RD*dia.diff + elev*dia.diff + Ndep*dia.diff +
                     SPCD.BA*dia.diff + non_SPCD.BA*dia.diff + prop.focal.ba*dia.diff +
                     # all diameter interactions
                     aspect*DIA  + slope*DIA  + MATmax*DIA  
                   + MATmin*DIA  + MAP*DIA  + MATmaxanom*DIA  
                   + MATminanom*DIA  + MAPanom*DIA  + BAL*DIA  + damage*DIA  + si*DIA + 
                     PHYSIO*DIA + BA*DIA + RD*DIA + elev*DIA + Ndep*DIA +
                     SPCD.BA*DIA + non_SPCD.BA*DIA + prop.focal.ba*DIA + 
                     # all aspect interactions
                     slope*aspect  + MATmax*aspect  
                   + MATmin*aspect  + MAP*aspect  + MATmaxanom*aspect  
                   + MATminanom*aspect  + MAPanom*aspect  + BAL*aspect  + damage*aspect  + si*aspect +
                     PHYSIO*aspect + BA*aspect + RD*aspect + elev*aspect + Ndep*aspect+
                     SPCD.BA*aspect + non_SPCD.BA*aspect + prop.focal.ba*aspect +
                     # all slope interactions
                     MATmax*slope  
                   + MATmin*slope  + MAP*slope  + MATmaxanom*slope  
                   + MATminanom*slope  + MAPanom*slope  + BAL*slope  + damage*slope  + si*slope + 
                     PHYSIO*slope + BA*slope + RD*slope + elev*slope + Ndep*slope +
                     SPCD.BA*slope + non_SPCD.BA*slope + prop.focal.ba*slope+ 
                     # all MATmax interactions
                     
                     MATmin*MATmax  + MAP*MATmax  + MATmaxanom*MATmax  
                   + MATminanom*MATmax  + MAPanom*MATmax  + BAL*MATmax  + damage*MATmax  + si*MATmax +
                     PHYSIO*MATmax + BA*MATmax + RD*MATmax + elev*MATmax +  Ndep*MATmax +
                     SPCD.BA*MATmax + non_SPCD.BA*MATmax + prop.focal.ba*MATmax+ 
                     # all MATmin interactions
                     
                     MAP*MATmin + MATminanom*MATmin
                   + MATminanom*MATmin + MAPanom*MATmin + BAL*MATmin + damage*MATmin + si*MATmin +
                     PHYSIO*MATmin + BA*MATmin + RD*MATmin + elev*MATmin + Ndep*MATmin +
                     SPCD.BA*MATmin + non_SPCD.BA*MATmin + prop.focal.ba*MATmin+ 
                     # all MAP interatction
                     MATminanom*MAP +
                     MATmaxanom*MAP + MAPanom*MAP + BAL*MAP + damage*MAP + si*MAP +
                     PHYSIO*MAP + BA*MAP + RD*MAP + elev*MAP + Ndep*MAP +
                     SPCD.BA*MAP + non_SPCD.BA*MAP + prop.focal.ba*MAP, data = covariate.data, family = "binomial")
    
    
    glm.28 <-  glm(M ~ dia.diff + DIA + exp(DIA)+ aspect + slope + MATmax 
                   + MATmin+ MAP + MATmaxanom 
                   + MATminanom + MAPanom + BAL + damage + si +  PHYSIO + BA + RD + elev + Ndep + 
                     SPCD.BA + non_SPCD.BA + prop.focal.ba +
                     # annual growth interactcovarions
                     DIA*dia.diff + aspect*dia.diff  + slope*dia.diff  + MATmax*dia.diff  
                   + MATmin*dia.diff  + MAP*dia.diff  + MATmaxanom*dia.diff  
                   + MATminanom*dia.diff  + MAPanom*dia.diff  + BAL*dia.diff  + damage*dia.diff  + si*dia.diff +
                     +  PHYSIO*dia.diff + BA*dia.diff + RD*dia.diff + elev*dia.diff + Ndep*dia.diff +
                     SPCD.BA*dia.diff + non_SPCD.BA*dia.diff + prop.focal.ba*dia.diff +
                     # all diameter interactions
                     aspect*DIA  + slope*DIA  + MATmax*DIA  
                   + MATmin*DIA  + MAP*DIA  + MATmaxanom*DIA  
                   + MATminanom*DIA  + MAPanom*DIA  + BAL*DIA  + damage*DIA  + si*DIA + 
                     PHYSIO*DIA + BA*DIA + RD*DIA + elev*DIA + Ndep*DIA +
                     SPCD.BA*DIA + non_SPCD.BA*DIA + prop.focal.ba*DIA + 
                     # all aspect interactions
                     slope*aspect  + MATmax*aspect  
                   + MATmin*aspect  + MAP*aspect  + MATmaxanom*aspect  
                   + MATminanom*aspect  + MAPanom*aspect  + BAL*aspect  + damage*aspect  + si*aspect +
                     PHYSIO*aspect + BA*aspect + RD*aspect + elev*aspect + Ndep*aspect+
                     SPCD.BA*aspect + non_SPCD.BA*aspect + prop.focal.ba*aspect +
                     # all slope interactions
                     MATmax*slope  
                   + MATmin*slope  + MAP*slope  + MATmaxanom*slope  
                   + MATminanom*slope  + MAPanom*slope  + BAL*slope  + damage*slope  + si*slope + 
                     PHYSIO*slope + BA*slope + RD*slope + elev*slope + Ndep*slope +
                     SPCD.BA*slope + non_SPCD.BA*slope + prop.focal.ba*slope+ 
                     # all MATmax interactions
                     
                     MATmin*MATmax  + MAP*MATmax  + MATmaxanom*MATmax  
                   + MATminanom*MATmax  + MAPanom*MATmax  + BAL*MATmax  + damage*MATmax  + si*MATmax +
                     PHYSIO*MATmax + BA*MATmax + RD*MATmax + elev*MATmax +  Ndep*MATmax +
                     SPCD.BA*MATmax + non_SPCD.BA*MATmax + prop.focal.ba*MATmax+ 
                     # all MATmin interactions
                     
                     MAP*MATmin + MATminanom*MATmin
                   + MATminanom*MATmin + MAPanom*MATmin + BAL*MATmin + damage*MATmin + si*MATmin +
                     PHYSIO*MATmin + BA*MATmin + RD*MATmin + elev*MATmin + Ndep*MATmin +
                     SPCD.BA*MATmin + non_SPCD.BA*MATmin + prop.focal.ba*MATmin+ 
                     # all MAP interatction
                     MATminanom*MAP +
                     MATmaxanom*MAP + MAPanom*MAP + BAL*MAP + damage*MAP + si*MAP +
                     PHYSIO*MAP + BA*MAP + RD*MAP + elev*MAP + Ndep*MAP +
                     SPCD.BA*MAP + non_SPCD.BA*MAP + prop.focal.ba*MAP + 
                     
                     
                     # all MATminanom interatction
                     
                     MATmaxanom*MATminanom + MAPanom*MATminanom + BAL*MATminanom + damage*MATminanom + si*MATminanom +
                     PHYSIO*MATminanom + BA*MATminanom + RD*MATminanom + elev*MATminanom + Ndep*MATminanom +
                     SPCD.BA*MATminanom + non_SPCD.BA*MATminanom + prop.focal.ba*MATminanom, data = covariate.data, family = "binomial")
    
    glm.29 <-  glm(M ~ dia.diff + DIA + exp(DIA)+ aspect + slope + MATmax 
                   + MATmin+ MAP + MATmaxanom 
                   + MATminanom + MAPanom + BAL + damage + si +  PHYSIO + BA + RD + elev + Ndep + 
                     SPCD.BA + non_SPCD.BA + prop.focal.ba +
                     # annual growth interactcovarions
                     DIA*dia.diff + aspect*dia.diff  + slope*dia.diff  + MATmax*dia.diff  
                   + MATmin*dia.diff  + MAP*dia.diff  + MATmaxanom*dia.diff  
                   + MATminanom*dia.diff  + MAPanom*dia.diff  + BAL*dia.diff  + damage*dia.diff  + si*dia.diff +
                     +  PHYSIO*dia.diff + BA*dia.diff + RD*dia.diff + elev*dia.diff + Ndep*dia.diff +
                     SPCD.BA*dia.diff + non_SPCD.BA*dia.diff + prop.focal.ba*dia.diff +
                     # all diameter interactions
                     aspect*DIA  + slope*DIA  + MATmax*DIA  
                   + MATmin*DIA  + MAP*DIA  + MATmaxanom*DIA  
                   + MATminanom*DIA  + MAPanom*DIA  + BAL*DIA  + damage*DIA  + si*DIA + 
                     PHYSIO*DIA + BA*DIA + RD*DIA + elev*DIA + Ndep*DIA +
                     SPCD.BA*DIA + non_SPCD.BA*DIA + prop.focal.ba*DIA + 
                     # all aspect interactions
                     slope*aspect  + MATmax*aspect  
                   + MATmin*aspect  + MAP*aspect  + MATmaxanom*aspect  
                   + MATminanom*aspect  + MAPanom*aspect  + BAL*aspect  + damage*aspect  + si*aspect +
                     PHYSIO*aspect + BA*aspect + RD*aspect + elev*aspect + Ndep*aspect+
                     SPCD.BA*aspect + non_SPCD.BA*aspect + prop.focal.ba*aspect +
                     # all slope interactions
                     MATmax*slope  
                   + MATmin*slope  + MAP*slope  + MATmaxanom*slope  
                   + MATminanom*slope  + MAPanom*slope  + BAL*slope  + damage*slope  + si*slope + 
                     PHYSIO*slope + BA*slope + RD*slope + elev*slope + Ndep*slope +
                     SPCD.BA*slope + non_SPCD.BA*slope + prop.focal.ba*slope+ 
                     # all MATmax interactions
                     
                     MATmin*MATmax  + MAP*MATmax  + MATmaxanom*MATmax  
                   + MATminanom*MATmax  + MAPanom*MATmax  + BAL*MATmax  + damage*MATmax  + si*MATmax +
                     PHYSIO*MATmax + BA*MATmax + RD*MATmax + elev*MATmax +  Ndep*MATmax +
                     SPCD.BA*MATmax + non_SPCD.BA*MATmax + prop.focal.ba*MATmax+ 
                     # all MATmin interactions
                     
                     MAP*MATmin + MATminanom*MATmin
                   + MATminanom*MATmin + MAPanom*MATmin + BAL*MATmin + damage*MATmin + si*MATmin +
                     PHYSIO*MATmin + BA*MATmin + RD*MATmin + elev*MATmin + Ndep*MATmin +
                     SPCD.BA*MATmin + non_SPCD.BA*MATmin + prop.focal.ba*MATmin+ 
                     # all MAP interatction
                     MATminanom*MAP +
                     MATmaxanom*MAP + MAPanom*MAP + BAL*MAP + damage*MAP + si*MAP +
                     PHYSIO*MAP + BA*MAP + RD*MAP + elev*MAP + Ndep*MAP +
                     SPCD.BA*MAP + non_SPCD.BA*MAP + prop.focal.ba*MAP + 
                     
                     
                     # all MATminanom interatction
                     
                     MATmaxanom*MATminanom + MAPanom*MATminanom + BAL*MATminanom + damage*MATminanom + si*MATminanom +
                     PHYSIO*MATminanom + BA*MATminanom + RD*MATminanom + elev*MATminanom + Ndep*MATminanom +
                     SPCD.BA*MATminanom + non_SPCD.BA*MATminanom + prop.focal.ba*MATminanom+ 
                     # all MATmaxanom interatction
                     
                     MAPanom*MATmaxanom + BAL*MATmaxanom + damage*MATmaxanom + si*MATmaxanom +
                     PHYSIO*MATmaxanom + BA*MATmaxanom + RD*MATmaxanom + elev*MATmaxanom + Ndep*MATmaxanom +
                     SPCD.BA*MATmaxanom + non_SPCD.BA*MATmaxanom + prop.focal.ba*MATmaxanom, data = covariate.data, family = "binomial")
    
    
    
    glm.30 <-  glm(M ~ dia.diff + DIA + exp(DIA)+ aspect + slope + MATmax 
                   + MATmin+ MAP + MATmaxanom 
                   + MATminanom + MAPanom + BAL + damage + si +  PHYSIO + BA + RD + elev + Ndep + 
                     SPCD.BA + non_SPCD.BA + prop.focal.ba +
                     # annual growth interactcovarions
                     DIA*dia.diff + aspect*dia.diff  + slope*dia.diff  + MATmax*dia.diff  
                   + MATmin*dia.diff  + MAP*dia.diff  + MATmaxanom*dia.diff  
                   + MATminanom*dia.diff  + MAPanom*dia.diff  + BAL*dia.diff  + damage*dia.diff  + si*dia.diff +
                     +  PHYSIO*dia.diff + BA*dia.diff + RD*dia.diff + elev*dia.diff + Ndep*dia.diff +
                     SPCD.BA*dia.diff + non_SPCD.BA*dia.diff + prop.focal.ba*dia.diff +
                     # all diameter interactions
                     aspect*DIA  + slope*DIA  + MATmax*DIA  
                   + MATmin*DIA  + MAP*DIA  + MATmaxanom*DIA  
                   + MATminanom*DIA  + MAPanom*DIA  + BAL*DIA  + damage*DIA  + si*DIA + 
                     PHYSIO*DIA + BA*DIA + RD*DIA + elev*DIA + Ndep*DIA +
                     SPCD.BA*DIA + non_SPCD.BA*DIA + prop.focal.ba*DIA + 
                     # all aspect interactions
                     slope*aspect  + MATmax*aspect  
                   + MATmin*aspect  + MAP*aspect  + MATmaxanom*aspect  
                   + MATminanom*aspect  + MAPanom*aspect  + BAL*aspect  + damage*aspect  + si*aspect +
                     PHYSIO*aspect + BA*aspect + RD*aspect + elev*aspect + Ndep*aspect+
                     SPCD.BA*aspect + non_SPCD.BA*aspect + prop.focal.ba*aspect +
                     # all slope interactions
                     MATmax*slope  
                   + MATmin*slope  + MAP*slope  + MATmaxanom*slope  
                   + MATminanom*slope  + MAPanom*slope  + BAL*slope  + damage*slope  + si*slope + 
                     PHYSIO*slope + BA*slope + RD*slope + elev*slope + Ndep*slope +
                     SPCD.BA*slope + non_SPCD.BA*slope + prop.focal.ba*slope+ 
                     # all MATmax interactions
                     
                     MATmin*MATmax  + MAP*MATmax  + MATmaxanom*MATmax  
                   + MATminanom*MATmax  + MAPanom*MATmax  + BAL*MATmax  + damage*MATmax  + si*MATmax +
                     PHYSIO*MATmax + BA*MATmax + RD*MATmax + elev*MATmax +  Ndep*MATmax +
                     SPCD.BA*MATmax + non_SPCD.BA*MATmax + prop.focal.ba*MATmax+ 
                     # all MATmin interactions
                     
                     MAP*MATmin + MATminanom*MATmin
                   + MATminanom*MATmin + MAPanom*MATmin + BAL*MATmin + damage*MATmin + si*MATmin +
                     PHYSIO*MATmin + BA*MATmin + RD*MATmin + elev*MATmin + Ndep*MATmin +
                     SPCD.BA*MATmin + non_SPCD.BA*MATmin + prop.focal.ba*MATmin+ 
                     # all MAP interatction
                     MATminanom*MAP +
                     MATmaxanom*MAP + MAPanom*MAP + BAL*MAP + damage*MAP + si*MAP +
                     PHYSIO*MAP + BA*MAP + RD*MAP + elev*MAP + Ndep*MAP +
                     SPCD.BA*MAP + non_SPCD.BA*MAP + prop.focal.ba*MAP + 
                     
                     
                     # all MATminanom interatction
                     
                     MATmaxanom*MATminanom + MAPanom*MATminanom + BAL*MATminanom + damage*MATminanom + si*MATminanom +
                     PHYSIO*MATminanom + BA*MATminanom + RD*MATminanom + elev*MATminanom + Ndep*MATminanom +
                     SPCD.BA*MATminanom + non_SPCD.BA*MATminanom + prop.focal.ba*MATminanom+ 
                     # all MATmaxanom interatction
                     
                     MAPanom*MATmaxanom + BAL*MATmaxanom + damage*MATmaxanom + si*MATmaxanom +
                     PHYSIO*MATmaxanom + BA*MATmaxanom + RD*MATmaxanom + elev*MATmaxanom + Ndep*MATmaxanom +
                     SPCD.BA*MATmaxanom + non_SPCD.BA*MATmaxanom + prop.focal.ba*MATmaxanom +
                     # all MAPanom interatction
                     
                     BAL*MAPanom + damage*MAPanom + si*MAPanom +
                     PHYSIO*MAPanom + BA*MAPanom + RD*MAPanom +elev*MAPanom + Ndep*MAPanom+
                     SPCD.BA*MAPanom + non_SPCD.BA*MAPanom + prop.focal.ba*MAPanom , data = covariate.data, family = "binomial")
    
    glm.31 <-  glm(M ~ dia.diff + DIA + exp(DIA)+ aspect + slope + MATmax 
                   + MATmin+ MAP + MATmaxanom 
                   + MATminanom + MAPanom + BAL + damage + si +  PHYSIO + BA + RD + elev + Ndep + 
                     SPCD.BA + non_SPCD.BA + prop.focal.ba +
                     # annual growth interactcovarions
                     DIA*dia.diff + aspect*dia.diff  + slope*dia.diff  + MATmax*dia.diff  
                   + MATmin*dia.diff  + MAP*dia.diff  + MATmaxanom*dia.diff  
                   + MATminanom*dia.diff  + MAPanom*dia.diff  + BAL*dia.diff  + damage*dia.diff  + si*dia.diff +
                     +  PHYSIO*dia.diff + BA*dia.diff + RD*dia.diff + elev*dia.diff + Ndep*dia.diff +
                     SPCD.BA*dia.diff + non_SPCD.BA*dia.diff + prop.focal.ba*dia.diff +
                     # all diameter interactions
                     aspect*DIA  + slope*DIA  + MATmax*DIA  
                   + MATmin*DIA  + MAP*DIA  + MATmaxanom*DIA  
                   + MATminanom*DIA  + MAPanom*DIA  + BAL*DIA  + damage*DIA  + si*DIA + 
                     PHYSIO*DIA + BA*DIA + RD*DIA + elev*DIA + Ndep*DIA +
                     SPCD.BA*DIA + non_SPCD.BA*DIA + prop.focal.ba*DIA + 
                     # all aspect interactions
                     slope*aspect  + MATmax*aspect  
                   + MATmin*aspect  + MAP*aspect  + MATmaxanom*aspect  
                   + MATminanom*aspect  + MAPanom*aspect  + BAL*aspect  + damage*aspect  + si*aspect +
                     PHYSIO*aspect + BA*aspect + RD*aspect + elev*aspect + Ndep*aspect+
                     SPCD.BA*aspect + non_SPCD.BA*aspect + prop.focal.ba*aspect +
                     # all slope interactions
                     MATmax*slope  
                   + MATmin*slope  + MAP*slope  + MATmaxanom*slope  
                   + MATminanom*slope  + MAPanom*slope  + BAL*slope  + damage*slope  + si*slope + 
                     PHYSIO*slope + BA*slope + RD*slope + elev*slope + Ndep*slope +
                     SPCD.BA*slope + non_SPCD.BA*slope + prop.focal.ba*slope+ 
                     # all MATmax interactions
                     
                     MATmin*MATmax  + MAP*MATmax  + MATmaxanom*MATmax  
                   + MATminanom*MATmax  + MAPanom*MATmax  + BAL*MATmax  + damage*MATmax  + si*MATmax +
                     PHYSIO*MATmax + BA*MATmax + RD*MATmax + elev*MATmax +  Ndep*MATmax +
                     SPCD.BA*MATmax + non_SPCD.BA*MATmax + prop.focal.ba*MATmax+ 
                     # all MATmin interactions
                     
                     MAP*MATmin + MATminanom*MATmin
                   + MATminanom*MATmin + MAPanom*MATmin + BAL*MATmin + damage*MATmin + si*MATmin +
                     PHYSIO*MATmin + BA*MATmin + RD*MATmin + elev*MATmin + Ndep*MATmin +
                     SPCD.BA*MATmin + non_SPCD.BA*MATmin + prop.focal.ba*MATmin+ 
                     # all MAP interatction
                     MATminanom*MAP +
                     MATmaxanom*MAP + MAPanom*MAP + BAL*MAP + damage*MAP + si*MAP +
                     PHYSIO*MAP + BA*MAP + RD*MAP + elev*MAP + Ndep*MAP +
                     SPCD.BA*MAP + non_SPCD.BA*MAP + prop.focal.ba*MAP + 
                     
                     
                     # all MATminanom interatction
                     
                     MATmaxanom*MATminanom + MAPanom*MATminanom + BAL*MATminanom + damage*MATminanom + si*MATminanom +
                     PHYSIO*MATminanom + BA*MATminanom + RD*MATminanom + elev*MATminanom + Ndep*MATminanom +
                     SPCD.BA*MATminanom + non_SPCD.BA*MATminanom + prop.focal.ba*MATminanom+ 
                     # all MATmaxanom interatction
                     
                     MAPanom*MATmaxanom + BAL*MATmaxanom + damage*MATmaxanom + si*MATmaxanom +
                     PHYSIO*MATmaxanom + BA*MATmaxanom + RD*MATmaxanom + elev*MATmaxanom + Ndep*MATmaxanom +
                     SPCD.BA*MATmaxanom + non_SPCD.BA*MATmaxanom + prop.focal.ba*MATmaxanom +
                     # all MAPanom interatction
                     
                     BAL*MAPanom + damage*MAPanom + si*MAPanom +
                     PHYSIO*MAPanom + BA*MAPanom + RD*MAPanom +elev*MAPanom + Ndep*MAPanom+
                     SPCD.BA*MAPanom + non_SPCD.BA*MAPanom + prop.focal.ba*MAPanom  + 
                     # all BAL interatction
                     
                     damage*BAL + si*BAL + 
                     PHYSIO*BAL + BA*BAL + RD*BAL +elev*BAL + Ndep*BAL +
                     SPCD.BA*BAL + non_SPCD.BA*BAL + prop.focal.ba*BAL, data = covariate.data, family = "binomial")
    
    glm.32 <-  glm(M ~ dia.diff + DIA + exp(DIA)+ aspect + slope + MATmax 
                   + MATmin+ MAP + MATmaxanom 
                   + MATminanom + MAPanom + BAL + damage + si +  PHYSIO + BA + RD + elev + Ndep + 
                     SPCD.BA + non_SPCD.BA + prop.focal.ba +
                     # annual growth interactcovarions
                     DIA*dia.diff + aspect*dia.diff  + slope*dia.diff  + MATmax*dia.diff  
                   + MATmin*dia.diff  + MAP*dia.diff  + MATmaxanom*dia.diff  
                   + MATminanom*dia.diff  + MAPanom*dia.diff  + BAL*dia.diff  + damage*dia.diff  + si*dia.diff +
                     +  PHYSIO*dia.diff + BA*dia.diff + RD*dia.diff + elev*dia.diff + Ndep*dia.diff +
                     SPCD.BA*dia.diff + non_SPCD.BA*dia.diff + prop.focal.ba*dia.diff +
                     # all diameter interactions
                     aspect*DIA  + slope*DIA  + MATmax*DIA  
                   + MATmin*DIA  + MAP*DIA  + MATmaxanom*DIA  
                   + MATminanom*DIA  + MAPanom*DIA  + BAL*DIA  + damage*DIA  + si*DIA + 
                     PHYSIO*DIA + BA*DIA + RD*DIA + elev*DIA + Ndep*DIA +
                     SPCD.BA*DIA + non_SPCD.BA*DIA + prop.focal.ba*DIA + 
                     # all aspect interactions
                     slope*aspect  + MATmax*aspect  
                   + MATmin*aspect  + MAP*aspect  + MATmaxanom*aspect  
                   + MATminanom*aspect  + MAPanom*aspect  + BAL*aspect  + damage*aspect  + si*aspect +
                     PHYSIO*aspect + BA*aspect + RD*aspect + elev*aspect + Ndep*aspect+
                     SPCD.BA*aspect + non_SPCD.BA*aspect + prop.focal.ba*aspect +
                     # all slope interactions
                     MATmax*slope  
                   + MATmin*slope  + MAP*slope  + MATmaxanom*slope  
                   + MATminanom*slope  + MAPanom*slope  + BAL*slope  + damage*slope  + si*slope + 
                     PHYSIO*slope + BA*slope + RD*slope + elev*slope + Ndep*slope +
                     SPCD.BA*slope + non_SPCD.BA*slope + prop.focal.ba*slope+ 
                     # all MATmax interactions
                     
                     MATmin*MATmax  + MAP*MATmax  + MATmaxanom*MATmax  
                   + MATminanom*MATmax  + MAPanom*MATmax  + BAL*MATmax  + damage*MATmax  + si*MATmax +
                     PHYSIO*MATmax + BA*MATmax + RD*MATmax + elev*MATmax +  Ndep*MATmax +
                     SPCD.BA*MATmax + non_SPCD.BA*MATmax + prop.focal.ba*MATmax+ 
                     # all MATmin interactions
                     
                     MAP*MATmin + MATminanom*MATmin
                   + MATminanom*MATmin + MAPanom*MATmin + BAL*MATmin + damage*MATmin + si*MATmin +
                     PHYSIO*MATmin + BA*MATmin + RD*MATmin + elev*MATmin + Ndep*MATmin +
                     SPCD.BA*MATmin + non_SPCD.BA*MATmin + prop.focal.ba*MATmin+ 
                     # all MAP interatction
                     MATminanom*MAP +
                     MATmaxanom*MAP + MAPanom*MAP + BAL*MAP + damage*MAP + si*MAP +
                     PHYSIO*MAP + BA*MAP + RD*MAP + elev*MAP + Ndep*MAP +
                     SPCD.BA*MAP + non_SPCD.BA*MAP + prop.focal.ba*MAP + 
                     
                     
                     # all MATminanom interatction
                     
                     MATmaxanom*MATminanom + MAPanom*MATminanom + BAL*MATminanom + damage*MATminanom + si*MATminanom +
                     PHYSIO*MATminanom + BA*MATminanom + RD*MATminanom + elev*MATminanom + Ndep*MATminanom +
                     SPCD.BA*MATminanom + non_SPCD.BA*MATminanom + prop.focal.ba*MATminanom+ 
                     # all MATmaxanom interatction
                     
                     MAPanom*MATmaxanom + BAL*MATmaxanom + damage*MATmaxanom + si*MATmaxanom +
                     PHYSIO*MATmaxanom + BA*MATmaxanom + RD*MATmaxanom + elev*MATmaxanom + Ndep*MATmaxanom +
                     SPCD.BA*MATmaxanom + non_SPCD.BA*MATmaxanom + prop.focal.ba*MATmaxanom +
                     # all MAPanom interatction
                     
                     BAL*MAPanom + damage*MAPanom + si*MAPanom +
                     PHYSIO*MAPanom + BA*MAPanom + RD*MAPanom +elev*MAPanom + Ndep*MAPanom+
                     SPCD.BA*MAPanom + non_SPCD.BA*MAPanom + prop.focal.ba*MAPanom  + 
                     # all BAL interatction
                     
                     damage*BAL + si*BAL + 
                     PHYSIO*BAL + BA*BAL + RD*BAL +elev*BAL + Ndep*BAL +
                     SPCD.BA*BAL + non_SPCD.BA*BAL + prop.focal.ba*BAL+
                     # all remainingdamage interatction
                     
                     si*damage + PHYSIO*damage + BA*damage + RD*damage + elev*damage + Ndep*damage +
                     SPCD.BA*damage + non_SPCD.BA*damage + prop.focal.ba*damage, data = covariate.data, family = "binomial")
    
    glm.33 <-  glm(M ~ dia.diff + DIA + exp(DIA)+ aspect + slope + MATmax 
                   + MATmin+ MAP + MATmaxanom 
                   + MATminanom + MAPanom + BAL + damage + si +  PHYSIO + BA + RD + elev + Ndep + 
                     SPCD.BA + non_SPCD.BA + prop.focal.ba +
                     # annual growth interactcovarions
                     DIA*dia.diff + aspect*dia.diff  + slope*dia.diff  + MATmax*dia.diff  
                   + MATmin*dia.diff  + MAP*dia.diff  + MATmaxanom*dia.diff  
                   + MATminanom*dia.diff  + MAPanom*dia.diff  + BAL*dia.diff  + damage*dia.diff  + si*dia.diff +
                     +  PHYSIO*dia.diff + BA*dia.diff + RD*dia.diff + elev*dia.diff + Ndep*dia.diff +
                     SPCD.BA*dia.diff + non_SPCD.BA*dia.diff + prop.focal.ba*dia.diff +
                     # all diameter interactions
                     aspect*DIA  + slope*DIA  + MATmax*DIA  
                   + MATmin*DIA  + MAP*DIA  + MATmaxanom*DIA  
                   + MATminanom*DIA  + MAPanom*DIA  + BAL*DIA  + damage*DIA  + si*DIA + 
                     PHYSIO*DIA + BA*DIA + RD*DIA + elev*DIA + Ndep*DIA +
                     SPCD.BA*DIA + non_SPCD.BA*DIA + prop.focal.ba*DIA + 
                     # all aspect interactions
                     slope*aspect  + MATmax*aspect  
                   + MATmin*aspect  + MAP*aspect  + MATmaxanom*aspect  
                   + MATminanom*aspect  + MAPanom*aspect  + BAL*aspect  + damage*aspect  + si*aspect +
                     PHYSIO*aspect + BA*aspect + RD*aspect + elev*aspect + Ndep*aspect+
                     SPCD.BA*aspect + non_SPCD.BA*aspect + prop.focal.ba*aspect +
                     # all slope interactions
                     MATmax*slope  
                   + MATmin*slope  + MAP*slope  + MATmaxanom*slope  
                   + MATminanom*slope  + MAPanom*slope  + BAL*slope  + damage*slope  + si*slope + 
                     PHYSIO*slope + BA*slope + RD*slope + elev*slope + Ndep*slope +
                     SPCD.BA*slope + non_SPCD.BA*slope + prop.focal.ba*slope+ 
                     # all MATmax interactions
                     
                     MATmin*MATmax  + MAP*MATmax  + MATmaxanom*MATmax  
                   + MATminanom*MATmax  + MAPanom*MATmax  + BAL*MATmax  + damage*MATmax  + si*MATmax +
                     PHYSIO*MATmax + BA*MATmax + RD*MATmax + elev*MATmax +  Ndep*MATmax +
                     SPCD.BA*MATmax + non_SPCD.BA*MATmax + prop.focal.ba*MATmax+ 
                     # all MATmin interactions
                     
                     MAP*MATmin + MATminanom*MATmin
                   + MATminanom*MATmin + MAPanom*MATmin + BAL*MATmin + damage*MATmin + si*MATmin +
                     PHYSIO*MATmin + BA*MATmin + RD*MATmin + elev*MATmin + Ndep*MATmin +
                     SPCD.BA*MATmin + non_SPCD.BA*MATmin + prop.focal.ba*MATmin+ 
                     # all MAP interatction
                     MATminanom*MAP +
                     MATmaxanom*MAP + MAPanom*MAP + BAL*MAP + damage*MAP + si*MAP +
                     PHYSIO*MAP + BA*MAP + RD*MAP + elev*MAP + Ndep*MAP +
                     SPCD.BA*MAP + non_SPCD.BA*MAP + prop.focal.ba*MAP + 
                     
                     
                     # all MATminanom interatction
                     
                     MATmaxanom*MATminanom + MAPanom*MATminanom + BAL*MATminanom + damage*MATminanom + si*MATminanom +
                     PHYSIO*MATminanom + BA*MATminanom + RD*MATminanom + elev*MATminanom + Ndep*MATminanom +
                     SPCD.BA*MATminanom + non_SPCD.BA*MATminanom + prop.focal.ba*MATminanom+ 
                     # all MATmaxanom interatction
                     
                     MAPanom*MATmaxanom + BAL*MATmaxanom + damage*MATmaxanom + si*MATmaxanom +
                     PHYSIO*MATmaxanom + BA*MATmaxanom + RD*MATmaxanom + elev*MATmaxanom + Ndep*MATmaxanom +
                     SPCD.BA*MATmaxanom + non_SPCD.BA*MATmaxanom + prop.focal.ba*MATmaxanom +
                     # all MAPanom interatction
                     
                     BAL*MAPanom + damage*MAPanom + si*MAPanom +
                     PHYSIO*MAPanom + BA*MAPanom + RD*MAPanom +elev*MAPanom + Ndep*MAPanom+
                     SPCD.BA*MAPanom + non_SPCD.BA*MAPanom + prop.focal.ba*MAPanom  + 
                     # all BAL interatction
                     
                     damage*BAL + si*BAL + 
                     PHYSIO*BAL + BA*BAL + RD*BAL +elev*BAL + Ndep*BAL +
                     SPCD.BA*BAL + non_SPCD.BA*BAL + prop.focal.ba*BAL+
                     # all remainingdamage interatction
                     
                     si*damage + PHYSIO*damage + BA*damage + RD*damage + elev*damage + Ndep*damage +
                     SPCD.BA*damage + non_SPCD.BA*damage + prop.focal.ba*damage+
                     
                     PHYSIO*si + BA*si + RD*si + elev*si + Ndep*si +
                     SPCD.BA*si + non_SPCD.BA*si + prop.focal.ba*si, data = covariate.data, family = "binomial")
    
    glm.34 <-  glm(M ~ dia.diff + DIA + exp(DIA)+ aspect + slope + MATmax 
                   + MATmin+ MAP + MATmaxanom 
                   + MATminanom + MAPanom + BAL + damage + si +  PHYSIO + BA + RD + elev + Ndep + 
                     SPCD.BA + non_SPCD.BA + prop.focal.ba +
                     # annual growth interactcovarions
                     DIA*dia.diff + aspect*dia.diff  + slope*dia.diff  + MATmax*dia.diff  
                   + MATmin*dia.diff  + MAP*dia.diff  + MATmaxanom*dia.diff  
                   + MATminanom*dia.diff  + MAPanom*dia.diff  + BAL*dia.diff  + damage*dia.diff  + si*dia.diff +
                     +  PHYSIO*dia.diff + BA*dia.diff + RD*dia.diff + elev*dia.diff + Ndep*dia.diff +
                     SPCD.BA*dia.diff + non_SPCD.BA*dia.diff + prop.focal.ba*dia.diff +
                     # all diameter interactions
                     aspect*DIA  + slope*DIA  + MATmax*DIA  
                   + MATmin*DIA  + MAP*DIA  + MATmaxanom*DIA  
                   + MATminanom*DIA  + MAPanom*DIA  + BAL*DIA  + damage*DIA  + si*DIA + 
                     PHYSIO*DIA + BA*DIA + RD*DIA + elev*DIA + Ndep*DIA +
                     SPCD.BA*DIA + non_SPCD.BA*DIA + prop.focal.ba*DIA + 
                     # all aspect interactions
                     slope*aspect  + MATmax*aspect  
                   + MATmin*aspect  + MAP*aspect  + MATmaxanom*aspect  
                   + MATminanom*aspect  + MAPanom*aspect  + BAL*aspect  + damage*aspect  + si*aspect +
                     PHYSIO*aspect + BA*aspect + RD*aspect + elev*aspect + Ndep*aspect+
                     SPCD.BA*aspect + non_SPCD.BA*aspect + prop.focal.ba*aspect +
                     # all slope interactions
                     MATmax*slope  
                   + MATmin*slope  + MAP*slope  + MATmaxanom*slope  
                   + MATminanom*slope  + MAPanom*slope  + BAL*slope  + damage*slope  + si*slope + 
                     PHYSIO*slope + BA*slope + RD*slope + elev*slope + Ndep*slope +
                     SPCD.BA*slope + non_SPCD.BA*slope + prop.focal.ba*slope+ 
                     # all MATmax interactions
                     
                     MATmin*MATmax  + MAP*MATmax  + MATmaxanom*MATmax  
                   + MATminanom*MATmax  + MAPanom*MATmax  + BAL*MATmax  + damage*MATmax  + si*MATmax +
                     PHYSIO*MATmax + BA*MATmax + RD*MATmax + elev*MATmax +  Ndep*MATmax +
                     SPCD.BA*MATmax + non_SPCD.BA*MATmax + prop.focal.ba*MATmax+ 
                     # all MATmin interactions
                     
                     MAP*MATmin + MATminanom*MATmin
                   + MATminanom*MATmin + MAPanom*MATmin + BAL*MATmin + damage*MATmin + si*MATmin +
                     PHYSIO*MATmin + BA*MATmin + RD*MATmin + elev*MATmin + Ndep*MATmin +
                     SPCD.BA*MATmin + non_SPCD.BA*MATmin + prop.focal.ba*MATmin+ 
                     # all MAP interatction
                     MATminanom*MAP +
                     MATmaxanom*MAP + MAPanom*MAP + BAL*MAP + damage*MAP + si*MAP +
                     PHYSIO*MAP + BA*MAP + RD*MAP + elev*MAP + Ndep*MAP +
                     SPCD.BA*MAP + non_SPCD.BA*MAP + prop.focal.ba*MAP + 
                     
                     
                     # all MATminanom interatction
                     
                     MATmaxanom*MATminanom + MAPanom*MATminanom + BAL*MATminanom + damage*MATminanom + si*MATminanom +
                     PHYSIO*MATminanom + BA*MATminanom + RD*MATminanom + elev*MATminanom + Ndep*MATminanom +
                     SPCD.BA*MATminanom + non_SPCD.BA*MATminanom + prop.focal.ba*MATminanom+ 
                     # all MATmaxanom interatction
                     
                     MAPanom*MATmaxanom + BAL*MATmaxanom + damage*MATmaxanom + si*MATmaxanom +
                     PHYSIO*MATmaxanom + BA*MATmaxanom + RD*MATmaxanom + elev*MATmaxanom + Ndep*MATmaxanom +
                     SPCD.BA*MATmaxanom + non_SPCD.BA*MATmaxanom + prop.focal.ba*MATmaxanom +
                     # all MAPanom interatction
                     
                     BAL*MAPanom + damage*MAPanom + si*MAPanom +
                     PHYSIO*MAPanom + BA*MAPanom + RD*MAPanom +elev*MAPanom + Ndep*MAPanom+
                     SPCD.BA*MAPanom + non_SPCD.BA*MAPanom + prop.focal.ba*MAPanom  + 
                     # all BAL interatction
                     
                     damage*BAL + si*BAL + 
                     PHYSIO*BAL + BA*BAL + RD*BAL +elev*BAL + Ndep*BAL +
                     SPCD.BA*BAL + non_SPCD.BA*BAL + prop.focal.ba*BAL+
                     # all remainingdamage interatction
                     
                     si*damage + PHYSIO*damage + BA*damage + RD*damage + elev*damage + Ndep*damage +
                     SPCD.BA*damage + non_SPCD.BA*damage + prop.focal.ba*damage+
                     
                     PHYSIO*si + BA*si + RD*si + elev*si + Ndep*si +
                     SPCD.BA*si + non_SPCD.BA*si + prop.focal.ba*si+
                     BA*PHYSIO + RD*PHYSIO + elev*PHYSIO + Ndep*PHYSIO+
                     SPCD.BA*PHYSIO + non_SPCD.BA*PHYSIO + prop.focal.ba*PHYSIO, data = covariate.data, family = "binomial")
    
    
    glm.35 <-  glm(M ~ dia.diff + DIA + exp(DIA)+ aspect + slope + MATmax 
                   + MATmin+ MAP + MATmaxanom 
                   + MATminanom + MAPanom + BAL + damage + si +  PHYSIO + BA + RD + elev + Ndep + 
                     SPCD.BA + non_SPCD.BA + prop.focal.ba +
                     # annual growth interactcovarions
                     DIA*dia.diff + aspect*dia.diff  + slope*dia.diff  + MATmax*dia.diff  
                   + MATmin*dia.diff  + MAP*dia.diff  + MATmaxanom*dia.diff  
                   + MATminanom*dia.diff  + MAPanom*dia.diff  + BAL*dia.diff  + damage*dia.diff  + si*dia.diff +
                     +  PHYSIO*dia.diff + BA*dia.diff + RD*dia.diff + elev*dia.diff + Ndep*dia.diff +
                     SPCD.BA*dia.diff + non_SPCD.BA*dia.diff + prop.focal.ba*dia.diff +
                     # all diameter interactions
                     aspect*DIA  + slope*DIA  + MATmax*DIA  
                   + MATmin*DIA  + MAP*DIA  + MATmaxanom*DIA  
                   + MATminanom*DIA  + MAPanom*DIA  + BAL*DIA  + damage*DIA  + si*DIA + 
                     PHYSIO*DIA + BA*DIA + RD*DIA + elev*DIA + Ndep*DIA +
                     SPCD.BA*DIA + non_SPCD.BA*DIA + prop.focal.ba*DIA + 
                     # all aspect interactions
                     slope*aspect  + MATmax*aspect  
                   + MATmin*aspect  + MAP*aspect  + MATmaxanom*aspect  
                   + MATminanom*aspect  + MAPanom*aspect  + BAL*aspect  + damage*aspect  + si*aspect +
                     PHYSIO*aspect + BA*aspect + RD*aspect + elev*aspect + Ndep*aspect+
                     SPCD.BA*aspect + non_SPCD.BA*aspect + prop.focal.ba*aspect +
                     # all slope interactions
                     MATmax*slope  
                   + MATmin*slope  + MAP*slope  + MATmaxanom*slope  
                   + MATminanom*slope  + MAPanom*slope  + BAL*slope  + damage*slope  + si*slope + 
                     PHYSIO*slope + BA*slope + RD*slope + elev*slope + Ndep*slope +
                     SPCD.BA*slope + non_SPCD.BA*slope + prop.focal.ba*slope+ 
                     # all MATmax interactions
                     
                     MATmin*MATmax  + MAP*MATmax  + MATmaxanom*MATmax  
                   + MATminanom*MATmax  + MAPanom*MATmax  + BAL*MATmax  + damage*MATmax  + si*MATmax +
                     PHYSIO*MATmax + BA*MATmax + RD*MATmax + elev*MATmax +  Ndep*MATmax +
                     SPCD.BA*MATmax + non_SPCD.BA*MATmax + prop.focal.ba*MATmax+ 
                     # all MATmin interactions
                     
                     MAP*MATmin + MATminanom*MATmin
                   + MATminanom*MATmin + MAPanom*MATmin + BAL*MATmin + damage*MATmin + si*MATmin +
                     PHYSIO*MATmin + BA*MATmin + RD*MATmin + elev*MATmin + Ndep*MATmin +
                     SPCD.BA*MATmin + non_SPCD.BA*MATmin + prop.focal.ba*MATmin+ 
                     # all MAP interatction
                     MATminanom*MAP +
                     MATmaxanom*MAP + MAPanom*MAP + BAL*MAP + damage*MAP + si*MAP +
                     PHYSIO*MAP + BA*MAP + RD*MAP + elev*MAP + Ndep*MAP +
                     SPCD.BA*MAP + non_SPCD.BA*MAP + prop.focal.ba*MAP + 
                     
                     
                     # all MATminanom interatction
                     
                     MATmaxanom*MATminanom + MAPanom*MATminanom + BAL*MATminanom + damage*MATminanom + si*MATminanom +
                     PHYSIO*MATminanom + BA*MATminanom + RD*MATminanom + elev*MATminanom + Ndep*MATminanom +
                     SPCD.BA*MATminanom + non_SPCD.BA*MATminanom + prop.focal.ba*MATminanom+ 
                     # all MATmaxanom interatction
                     
                     MAPanom*MATmaxanom + BAL*MATmaxanom + damage*MATmaxanom + si*MATmaxanom +
                     PHYSIO*MATmaxanom + BA*MATmaxanom + RD*MATmaxanom + elev*MATmaxanom + Ndep*MATmaxanom +
                     SPCD.BA*MATmaxanom + non_SPCD.BA*MATmaxanom + prop.focal.ba*MATmaxanom +
                     # all MAPanom interatction
                     
                     BAL*MAPanom + damage*MAPanom + si*MAPanom +
                     PHYSIO*MAPanom + BA*MAPanom + RD*MAPanom +elev*MAPanom + Ndep*MAPanom+
                     SPCD.BA*MAPanom + non_SPCD.BA*MAPanom + prop.focal.ba*MAPanom  + 
                     # all BAL interatction
                     
                     damage*BAL + si*BAL + 
                     PHYSIO*BAL + BA*BAL + RD*BAL +elev*BAL + Ndep*BAL +
                     SPCD.BA*BAL + non_SPCD.BA*BAL + prop.focal.ba*BAL+
                     # all remainingdamage interatction
                     
                     si*damage + PHYSIO*damage + BA*damage + RD*damage + elev*damage + Ndep*damage +
                     SPCD.BA*damage + non_SPCD.BA*damage + prop.focal.ba*damage+
                     
                     PHYSIO*si + BA*si + RD*si + elev*si + Ndep*si +
                     SPCD.BA*si + non_SPCD.BA*si + prop.focal.ba*si+
                     BA*PHYSIO + RD*PHYSIO + elev*PHYSIO + Ndep*PHYSIO+
                     SPCD.BA*PHYSIO + non_SPCD.BA*PHYSIO + prop.focal.ba*PHYSIO+ 
                     BA*RD + elev*RD + Ndep*RD +
                     SPCD.BA*RD + non_SPCD.BA*RD + prop.focal.ba*RD, data = covariate.data, family = "binomial")
    
    
    glm.36 <-  glm(M ~ dia.diff + DIA + exp(DIA)+ aspect + slope + MATmax 
                   + MATmin+ MAP + MATmaxanom 
                   + MATminanom + MAPanom + BAL + damage + si +  PHYSIO + BA + RD + elev + Ndep + 
                     SPCD.BA + non_SPCD.BA + prop.focal.ba +
                     # annual growth interactcovarions
                     DIA*dia.diff + aspect*dia.diff  + slope*dia.diff  + MATmax*dia.diff  
                   + MATmin*dia.diff  + MAP*dia.diff  + MATmaxanom*dia.diff  
                   + MATminanom*dia.diff  + MAPanom*dia.diff  + BAL*dia.diff  + damage*dia.diff  + si*dia.diff +
                     +  PHYSIO*dia.diff + BA*dia.diff + RD*dia.diff + elev*dia.diff + Ndep*dia.diff +
                     SPCD.BA*dia.diff + non_SPCD.BA*dia.diff + prop.focal.ba*dia.diff +
                     # all diameter interactions
                     aspect*DIA  + slope*DIA  + MATmax*DIA  
                   + MATmin*DIA  + MAP*DIA  + MATmaxanom*DIA  
                   + MATminanom*DIA  + MAPanom*DIA  + BAL*DIA  + damage*DIA  + si*DIA + 
                     PHYSIO*DIA + BA*DIA + RD*DIA + elev*DIA + Ndep*DIA +
                     SPCD.BA*DIA + non_SPCD.BA*DIA + prop.focal.ba*DIA + 
                     # all aspect interactions
                     slope*aspect  + MATmax*aspect  
                   + MATmin*aspect  + MAP*aspect  + MATmaxanom*aspect  
                   + MATminanom*aspect  + MAPanom*aspect  + BAL*aspect  + damage*aspect  + si*aspect +
                     PHYSIO*aspect + BA*aspect + RD*aspect + elev*aspect + Ndep*aspect+
                     SPCD.BA*aspect + non_SPCD.BA*aspect + prop.focal.ba*aspect +
                     # all slope interactions
                     MATmax*slope  
                   + MATmin*slope  + MAP*slope  + MATmaxanom*slope  
                   + MATminanom*slope  + MAPanom*slope  + BAL*slope  + damage*slope  + si*slope + 
                     PHYSIO*slope + BA*slope + RD*slope + elev*slope + Ndep*slope +
                     SPCD.BA*slope + non_SPCD.BA*slope + prop.focal.ba*slope+ 
                     # all MATmax interactions
                     
                     MATmin*MATmax  + MAP*MATmax  + MATmaxanom*MATmax  
                   + MATminanom*MATmax  + MAPanom*MATmax  + BAL*MATmax  + damage*MATmax  + si*MATmax +
                     PHYSIO*MATmax + BA*MATmax + RD*MATmax + elev*MATmax +  Ndep*MATmax +
                     SPCD.BA*MATmax + non_SPCD.BA*MATmax + prop.focal.ba*MATmax+ 
                     # all MATmin interactions
                     
                     MAP*MATmin + MATminanom*MATmin
                   + MATminanom*MATmin + MAPanom*MATmin + BAL*MATmin + damage*MATmin + si*MATmin +
                     PHYSIO*MATmin + BA*MATmin + RD*MATmin + elev*MATmin + Ndep*MATmin +
                     SPCD.BA*MATmin + non_SPCD.BA*MATmin + prop.focal.ba*MATmin+ 
                     # all MAP interatction
                     MATminanom*MAP +
                     MATmaxanom*MAP + MAPanom*MAP + BAL*MAP + damage*MAP + si*MAP +
                     PHYSIO*MAP + BA*MAP + RD*MAP + elev*MAP + Ndep*MAP +
                     SPCD.BA*MAP + non_SPCD.BA*MAP + prop.focal.ba*MAP + 
                     
                     
                     # all MATminanom interatction
                     
                     MATmaxanom*MATminanom + MAPanom*MATminanom + BAL*MATminanom + damage*MATminanom + si*MATminanom +
                     PHYSIO*MATminanom + BA*MATminanom + RD*MATminanom + elev*MATminanom + Ndep*MATminanom +
                     SPCD.BA*MATminanom + non_SPCD.BA*MATminanom + prop.focal.ba*MATminanom+ 
                     # all MATmaxanom interatction
                     
                     MAPanom*MATmaxanom + BAL*MATmaxanom + damage*MATmaxanom + si*MATmaxanom +
                     PHYSIO*MATmaxanom + BA*MATmaxanom + RD*MATmaxanom + elev*MATmaxanom + Ndep*MATmaxanom +
                     SPCD.BA*MATmaxanom + non_SPCD.BA*MATmaxanom + prop.focal.ba*MATmaxanom +
                     # all MAPanom interatction
                     
                     BAL*MAPanom + damage*MAPanom + si*MAPanom +
                     PHYSIO*MAPanom + BA*MAPanom + RD*MAPanom +elev*MAPanom + Ndep*MAPanom+
                     SPCD.BA*MAPanom + non_SPCD.BA*MAPanom + prop.focal.ba*MAPanom  + 
                     # all BAL interatction
                     
                     damage*BAL + si*BAL + 
                     PHYSIO*BAL + BA*BAL + RD*BAL +elev*BAL + Ndep*BAL +
                     SPCD.BA*BAL + non_SPCD.BA*BAL + prop.focal.ba*BAL+
                     # all remainingdamage interatction
                     
                     si*damage + PHYSIO*damage + BA*damage + RD*damage + elev*damage + Ndep*damage +
                     SPCD.BA*damage + non_SPCD.BA*damage + prop.focal.ba*damage+
                     
                     PHYSIO*si + BA*si + RD*si + elev*si + Ndep*si +
                     SPCD.BA*si + non_SPCD.BA*si + prop.focal.ba*si+
                     BA*PHYSIO + RD*PHYSIO + elev*PHYSIO + Ndep*PHYSIO+
                     SPCD.BA*PHYSIO + non_SPCD.BA*PHYSIO + prop.focal.ba*PHYSIO+ 
                     BA*RD + elev*RD + Ndep*RD +
                     SPCD.BA*RD + non_SPCD.BA*RD + prop.focal.ba*RD+
                     elev*BA + Ndep*BA + 
                     SPCD.BA*BA + non_SPCD.BA*BA + prop.focal.ba*BA, data = covariate.data, family = "binomial")
    
    
    glm.37 <-  glm(M ~ dia.diff + DIA + exp(DIA)+ aspect + slope + MATmax 
                   + MATmin+ MAP + MATmaxanom 
                   + MATminanom + MAPanom + BAL + damage + si +  PHYSIO + BA + RD + elev + Ndep + 
                     SPCD.BA + non_SPCD.BA + prop.focal.ba +
                     # annual growth interactcovarions
                     DIA*dia.diff + aspect*dia.diff  + slope*dia.diff  + MATmax*dia.diff  
                   + MATmin*dia.diff  + MAP*dia.diff  + MATmaxanom*dia.diff  
                   + MATminanom*dia.diff  + MAPanom*dia.diff  + BAL*dia.diff  + damage*dia.diff  + si*dia.diff +
                     +  PHYSIO*dia.diff + BA*dia.diff + RD*dia.diff + elev*dia.diff + Ndep*dia.diff +
                     SPCD.BA*dia.diff + non_SPCD.BA*dia.diff + prop.focal.ba*dia.diff +
                     # all diameter interactions
                     aspect*DIA  + slope*DIA  + MATmax*DIA  
                   + MATmin*DIA  + MAP*DIA  + MATmaxanom*DIA  
                   + MATminanom*DIA  + MAPanom*DIA  + BAL*DIA  + damage*DIA  + si*DIA + 
                     PHYSIO*DIA + BA*DIA + RD*DIA + elev*DIA + Ndep*DIA +
                     SPCD.BA*DIA + non_SPCD.BA*DIA + prop.focal.ba*DIA + 
                     # all aspect interactions
                     slope*aspect  + MATmax*aspect  
                   + MATmin*aspect  + MAP*aspect  + MATmaxanom*aspect  
                   + MATminanom*aspect  + MAPanom*aspect  + BAL*aspect  + damage*aspect  + si*aspect +
                     PHYSIO*aspect + BA*aspect + RD*aspect + elev*aspect + Ndep*aspect+
                     SPCD.BA*aspect + non_SPCD.BA*aspect + prop.focal.ba*aspect +
                     # all slope interactions
                     MATmax*slope  
                   + MATmin*slope  + MAP*slope  + MATmaxanom*slope  
                   + MATminanom*slope  + MAPanom*slope  + BAL*slope  + damage*slope  + si*slope + 
                     PHYSIO*slope + BA*slope + RD*slope + elev*slope + Ndep*slope +
                     SPCD.BA*slope + non_SPCD.BA*slope + prop.focal.ba*slope+ 
                     # all MATmax interactions
                     
                     MATmin*MATmax  + MAP*MATmax  + MATmaxanom*MATmax  
                   + MATminanom*MATmax  + MAPanom*MATmax  + BAL*MATmax  + damage*MATmax  + si*MATmax +
                     PHYSIO*MATmax + BA*MATmax + RD*MATmax + elev*MATmax +  Ndep*MATmax +
                     SPCD.BA*MATmax + non_SPCD.BA*MATmax + prop.focal.ba*MATmax+ 
                     # all MATmin interactions
                     
                     MAP*MATmin + MATminanom*MATmin
                   + MATminanom*MATmin + MAPanom*MATmin + BAL*MATmin + damage*MATmin + si*MATmin +
                     PHYSIO*MATmin + BA*MATmin + RD*MATmin + elev*MATmin + Ndep*MATmin +
                     SPCD.BA*MATmin + non_SPCD.BA*MATmin + prop.focal.ba*MATmin+ 
                     # all MAP interatction
                     MATminanom*MAP +
                     MATmaxanom*MAP + MAPanom*MAP + BAL*MAP + damage*MAP + si*MAP +
                     PHYSIO*MAP + BA*MAP + RD*MAP + elev*MAP + Ndep*MAP +
                     SPCD.BA*MAP + non_SPCD.BA*MAP + prop.focal.ba*MAP + 
                     
                     
                     # all MATminanom interatction
                     
                     MATmaxanom*MATminanom + MAPanom*MATminanom + BAL*MATminanom + damage*MATminanom + si*MATminanom +
                     PHYSIO*MATminanom + BA*MATminanom + RD*MATminanom + elev*MATminanom + Ndep*MATminanom +
                     SPCD.BA*MATminanom + non_SPCD.BA*MATminanom + prop.focal.ba*MATminanom+ 
                     # all MATmaxanom interatction
                     
                     MAPanom*MATmaxanom + BAL*MATmaxanom + damage*MATmaxanom + si*MATmaxanom +
                     PHYSIO*MATmaxanom + BA*MATmaxanom + RD*MATmaxanom + elev*MATmaxanom + Ndep*MATmaxanom +
                     SPCD.BA*MATmaxanom + non_SPCD.BA*MATmaxanom + prop.focal.ba*MATmaxanom +
                     # all MAPanom interatction
                     
                     BAL*MAPanom + damage*MAPanom + si*MAPanom +
                     PHYSIO*MAPanom + BA*MAPanom + RD*MAPanom +elev*MAPanom + Ndep*MAPanom+
                     SPCD.BA*MAPanom + non_SPCD.BA*MAPanom + prop.focal.ba*MAPanom  + 
                     # all BAL interatction
                     
                     damage*BAL + si*BAL + 
                     PHYSIO*BAL + BA*BAL + RD*BAL +elev*BAL + Ndep*BAL +
                     SPCD.BA*BAL + non_SPCD.BA*BAL + prop.focal.ba*BAL+
                     # all remainingdamage interatction
                     
                     si*damage + PHYSIO*damage + BA*damage + RD*damage + elev*damage + Ndep*damage +
                     SPCD.BA*damage + non_SPCD.BA*damage + prop.focal.ba*damage+
                     
                     PHYSIO*si + BA*si + RD*si + elev*si + Ndep*si +
                     SPCD.BA*si + non_SPCD.BA*si + prop.focal.ba*si+
                     BA*PHYSIO + RD*PHYSIO + elev*PHYSIO + Ndep*PHYSIO+
                     SPCD.BA*PHYSIO + non_SPCD.BA*PHYSIO + prop.focal.ba*PHYSIO+ 
                     BA*RD + elev*RD + Ndep*RD +
                     SPCD.BA*RD + non_SPCD.BA*RD + prop.focal.ba*RD+
                     elev*BA + Ndep*BA + 
                     SPCD.BA*BA + non_SPCD.BA*BA + prop.focal.ba*BA+ 
                     elev*Ndep + SPCD.BA*elev + non_SPCD.BA*elev + prop.focal.ba*elev, data = covariate.data, family = "binomial")
    
    glm.38 <-  glm(M ~ dia.diff + DIA + exp(DIA)+ aspect + slope + MATmax 
                   + MATmin+ MAP + MATmaxanom 
                   + MATminanom + MAPanom + BAL + damage + si +  PHYSIO + BA + RD + elev + Ndep + 
                     SPCD.BA + non_SPCD.BA + prop.focal.ba +
                     # annual growth interactcovarions
                     DIA*dia.diff + aspect*dia.diff  + slope*dia.diff  + MATmax*dia.diff  
                   + MATmin*dia.diff  + MAP*dia.diff  + MATmaxanom*dia.diff  
                   + MATminanom*dia.diff  + MAPanom*dia.diff  + BAL*dia.diff  + damage*dia.diff  + si*dia.diff +
                     +  PHYSIO*dia.diff + BA*dia.diff + RD*dia.diff + elev*dia.diff + Ndep*dia.diff +
                     SPCD.BA*dia.diff + non_SPCD.BA*dia.diff + prop.focal.ba*dia.diff +
                     # all diameter interactions
                     aspect*DIA  + slope*DIA  + MATmax*DIA  
                   + MATmin*DIA  + MAP*DIA  + MATmaxanom*DIA  
                   + MATminanom*DIA  + MAPanom*DIA  + BAL*DIA  + damage*DIA  + si*DIA + 
                     PHYSIO*DIA + BA*DIA + RD*DIA + elev*DIA + Ndep*DIA +
                     SPCD.BA*DIA + non_SPCD.BA*DIA + prop.focal.ba*DIA + 
                     # all aspect interactions
                     slope*aspect  + MATmax*aspect  
                   + MATmin*aspect  + MAP*aspect  + MATmaxanom*aspect  
                   + MATminanom*aspect  + MAPanom*aspect  + BAL*aspect  + damage*aspect  + si*aspect +
                     PHYSIO*aspect + BA*aspect + RD*aspect + elev*aspect + Ndep*aspect+
                     SPCD.BA*aspect + non_SPCD.BA*aspect + prop.focal.ba*aspect +
                     # all slope interactions
                     MATmax*slope  
                   + MATmin*slope  + MAP*slope  + MATmaxanom*slope  
                   + MATminanom*slope  + MAPanom*slope  + BAL*slope  + damage*slope  + si*slope + 
                     PHYSIO*slope + BA*slope + RD*slope + elev*slope + Ndep*slope +
                     SPCD.BA*slope + non_SPCD.BA*slope + prop.focal.ba*slope+ 
                     # all MATmax interactions
                     
                     MATmin*MATmax  + MAP*MATmax  + MATmaxanom*MATmax  
                   + MATminanom*MATmax  + MAPanom*MATmax  + BAL*MATmax  + damage*MATmax  + si*MATmax +
                     PHYSIO*MATmax + BA*MATmax + RD*MATmax + elev*MATmax +  Ndep*MATmax +
                     SPCD.BA*MATmax + non_SPCD.BA*MATmax + prop.focal.ba*MATmax+ 
                     # all MATmin interactions
                     
                     MAP*MATmin + MATminanom*MATmin
                   + MATminanom*MATmin + MAPanom*MATmin + BAL*MATmin + damage*MATmin + si*MATmin +
                     PHYSIO*MATmin + BA*MATmin + RD*MATmin + elev*MATmin + Ndep*MATmin +
                     SPCD.BA*MATmin + non_SPCD.BA*MATmin + prop.focal.ba*MATmin+ 
                     # all MAP interatction
                     MATminanom*MAP +
                     MATmaxanom*MAP + MAPanom*MAP + BAL*MAP + damage*MAP + si*MAP +
                     PHYSIO*MAP + BA*MAP + RD*MAP + elev*MAP + Ndep*MAP +
                     SPCD.BA*MAP + non_SPCD.BA*MAP + prop.focal.ba*MAP + 
                     
                     
                     # all MATminanom interatction
                     
                     MATmaxanom*MATminanom + MAPanom*MATminanom + BAL*MATminanom + damage*MATminanom + si*MATminanom +
                     PHYSIO*MATminanom + BA*MATminanom + RD*MATminanom + elev*MATminanom + Ndep*MATminanom +
                     SPCD.BA*MATminanom + non_SPCD.BA*MATminanom + prop.focal.ba*MATminanom+ 
                     # all MATmaxanom interatction
                     
                     MAPanom*MATmaxanom + BAL*MATmaxanom + damage*MATmaxanom + si*MATmaxanom +
                     PHYSIO*MATmaxanom + BA*MATmaxanom + RD*MATmaxanom + elev*MATmaxanom + Ndep*MATmaxanom +
                     SPCD.BA*MATmaxanom + non_SPCD.BA*MATmaxanom + prop.focal.ba*MATmaxanom +
                     # all MAPanom interatction
                     
                     BAL*MAPanom + damage*MAPanom + si*MAPanom +
                     PHYSIO*MAPanom + BA*MAPanom + RD*MAPanom +elev*MAPanom + Ndep*MAPanom+
                     SPCD.BA*MAPanom + non_SPCD.BA*MAPanom + prop.focal.ba*MAPanom  + 
                     # all BAL interatction
                     
                     damage*BAL + si*BAL + 
                     PHYSIO*BAL + BA*BAL + RD*BAL +elev*BAL + Ndep*BAL +
                     SPCD.BA*BAL + non_SPCD.BA*BAL + prop.focal.ba*BAL+
                     # all remainingdamage interatction
                     
                     si*damage + PHYSIO*damage + BA*damage + RD*damage + elev*damage + Ndep*damage +
                     SPCD.BA*damage + non_SPCD.BA*damage + prop.focal.ba*damage+
                     
                     PHYSIO*si + BA*si + RD*si + elev*si + Ndep*si +
                     SPCD.BA*si + non_SPCD.BA*si + prop.focal.ba*si+
                     BA*PHYSIO + RD*PHYSIO + elev*PHYSIO + Ndep*PHYSIO+
                     SPCD.BA*PHYSIO + non_SPCD.BA*PHYSIO + prop.focal.ba*PHYSIO+ 
                     BA*RD + elev*RD + Ndep*RD +
                     SPCD.BA*RD + non_SPCD.BA*RD + prop.focal.ba*RD+
                     elev*BA + Ndep*BA + 
                     SPCD.BA*BA + non_SPCD.BA*BA + prop.focal.ba*BA+ 
                     elev*Ndep + SPCD.BA*elev + non_SPCD.BA*elev + prop.focal.ba*elev +
                     # SPCD.BA interactions
                     non_SPCD.BA*SPCD.BA + prop.focal.ba*SPCD.BA 
                   , data = covariate.data, family = "binomial")
    
    
    glm.39 <-  glm(M ~ dia.diff + DIA + exp(DIA)+ aspect + slope + MATmax 
                   + MATmin+ MAP + MATmaxanom 
                   + MATminanom + MAPanom + BAL + damage + si +  PHYSIO + BA + RD + elev + Ndep + 
                     SPCD.BA + non_SPCD.BA + prop.focal.ba +
                     # annual growth interactcovarions
                     DIA*dia.diff + aspect*dia.diff  + slope*dia.diff  + MATmax*dia.diff  
                   + MATmin*dia.diff  + MAP*dia.diff  + MATmaxanom*dia.diff  
                   + MATminanom*dia.diff  + MAPanom*dia.diff  + BAL*dia.diff  + damage*dia.diff  + si*dia.diff +
                     +  PHYSIO*dia.diff + BA*dia.diff + RD*dia.diff + elev*dia.diff + Ndep*dia.diff +
                     SPCD.BA*dia.diff + non_SPCD.BA*dia.diff + prop.focal.ba*dia.diff +
                     # all diameter interactions
                     aspect*DIA  + slope*DIA  + MATmax*DIA  
                   + MATmin*DIA  + MAP*DIA  + MATmaxanom*DIA  
                   + MATminanom*DIA  + MAPanom*DIA  + BAL*DIA  + damage*DIA  + si*DIA + 
                     PHYSIO*DIA + BA*DIA + RD*DIA + elev*DIA + Ndep*DIA +
                     SPCD.BA*DIA + non_SPCD.BA*DIA + prop.focal.ba*DIA + 
                     # all aspect interactions
                     slope*aspect  + MATmax*aspect  
                   + MATmin*aspect  + MAP*aspect  + MATmaxanom*aspect  
                   + MATminanom*aspect  + MAPanom*aspect  + BAL*aspect  + damage*aspect  + si*aspect +
                     PHYSIO*aspect + BA*aspect + RD*aspect + elev*aspect + Ndep*aspect+
                     SPCD.BA*aspect + non_SPCD.BA*aspect + prop.focal.ba*aspect +
                     # all slope interactions
                     MATmax*slope  
                   + MATmin*slope  + MAP*slope  + MATmaxanom*slope  
                   + MATminanom*slope  + MAPanom*slope  + BAL*slope  + damage*slope  + si*slope + 
                     PHYSIO*slope + BA*slope + RD*slope + elev*slope + Ndep*slope +
                     SPCD.BA*slope + non_SPCD.BA*slope + prop.focal.ba*slope+ 
                     # all MATmax interactions
                     
                     MATmin*MATmax  + MAP*MATmax  + MATmaxanom*MATmax  
                   + MATminanom*MATmax  + MAPanom*MATmax  + BAL*MATmax  + damage*MATmax  + si*MATmax +
                     PHYSIO*MATmax + BA*MATmax + RD*MATmax + elev*MATmax +  Ndep*MATmax +
                     SPCD.BA*MATmax + non_SPCD.BA*MATmax + prop.focal.ba*MATmax+ 
                     # all MATmin interactions
                     
                     MAP*MATmin + MATminanom*MATmin
                   + MATminanom*MATmin + MAPanom*MATmin + BAL*MATmin + damage*MATmin + si*MATmin +
                     PHYSIO*MATmin + BA*MATmin + RD*MATmin + elev*MATmin + Ndep*MATmin +
                     SPCD.BA*MATmin + non_SPCD.BA*MATmin + prop.focal.ba*MATmin+ 
                     # all MAP interatction
                     MATminanom*MAP +
                     MATmaxanom*MAP + MAPanom*MAP + BAL*MAP + damage*MAP + si*MAP +
                     PHYSIO*MAP + BA*MAP + RD*MAP + elev*MAP + Ndep*MAP +
                     SPCD.BA*MAP + non_SPCD.BA*MAP + prop.focal.ba*MAP + 
                     
                     
                     # all MATminanom interatction
                     
                     MATmaxanom*MATminanom + MAPanom*MATminanom + BAL*MATminanom + damage*MATminanom + si*MATminanom +
                     PHYSIO*MATminanom + BA*MATminanom + RD*MATminanom + elev*MATminanom + Ndep*MATminanom +
                     SPCD.BA*MATminanom + non_SPCD.BA*MATminanom + prop.focal.ba*MATminanom+ 
                     # all MATmaxanom interatction
                     
                     MAPanom*MATmaxanom + BAL*MATmaxanom + damage*MATmaxanom + si*MATmaxanom +
                     PHYSIO*MATmaxanom + BA*MATmaxanom + RD*MATmaxanom + elev*MATmaxanom + Ndep*MATmaxanom +
                     SPCD.BA*MATmaxanom + non_SPCD.BA*MATmaxanom + prop.focal.ba*MATmaxanom +
                     # all MAPanom interatction
                     
                     BAL*MAPanom + damage*MAPanom + si*MAPanom +
                     PHYSIO*MAPanom + BA*MAPanom + RD*MAPanom +elev*MAPanom + Ndep*MAPanom+
                     SPCD.BA*MAPanom + non_SPCD.BA*MAPanom + prop.focal.ba*MAPanom  + 
                     # all BAL interatction
                     
                     damage*BAL + si*BAL + 
                     PHYSIO*BAL + BA*BAL + RD*BAL +elev*BAL + Ndep*BAL +
                     SPCD.BA*BAL + non_SPCD.BA*BAL + prop.focal.ba*BAL+
                     # all remainingdamage interatction
                     
                     si*damage + PHYSIO*damage + BA*damage + RD*damage + elev*damage + Ndep*damage +
                     SPCD.BA*damage + non_SPCD.BA*damage + prop.focal.ba*damage+
                     
                     PHYSIO*si + BA*si + RD*si + elev*si + Ndep*si +
                     SPCD.BA*si + non_SPCD.BA*si + prop.focal.ba*si+
                     BA*PHYSIO + RD*PHYSIO + elev*PHYSIO + Ndep*PHYSIO+
                     SPCD.BA*PHYSIO + non_SPCD.BA*PHYSIO + prop.focal.ba*PHYSIO+ 
                     BA*RD + elev*RD + Ndep*RD +
                     SPCD.BA*RD + non_SPCD.BA*RD + prop.focal.ba*RD+
                     elev*BA + Ndep*BA + 
                     SPCD.BA*BA + non_SPCD.BA*BA + prop.focal.ba*BA+ 
                     elev*Ndep + SPCD.BA*elev + non_SPCD.BA*elev + prop.focal.ba*elev +
                     # SPCD.BA interactions
                     non_SPCD.BA*SPCD.BA + prop.focal.ba*SPCD.BA #+
                   # # non_SPCD.BA interactions
                   # non_SPCD.BA*prop.focal.ba
                   , data = covariate.data, family = "binomial")
    
    # mcfaddens rsquared
    list.mods <- list(glm.A, glm.B, glm.C, glm.D, glm.E, 
                      glm.F, glm.G, glm.H, glm.I, glm.J, 
                      glm.K, glm.L,glm.M, glm.N, glm.O, 
                      glm.1, glm.2, glm.3, glm.4, glm.5, 
                      glm.6, glm.7, glm.8, glm.9, glm.10, 
                      glm.11, glm.12, glm.13, glm.14, glm.15, 
                      glm.16, glm.17, glm.18, glm.19, glm.20, 
                      glm.20.a, glm.20.b, glm.20.c,
                      glm.21, glm.22, glm.23, glm.24, glm.25, 
                      glm.26, glm.27, glm.28, glm.29, glm.30, 
                      glm.31, glm.32, glm.33, glm.34, glm.35, 
                      glm.36, glm.37, glm.38)
    
    # get convergence list
    Convergence.list <- lapply(list.mods, FUN = function(x){x$converged}) 
    convergence.df <- do.call(rbind, Convergence.list)
    
    # get AICs
    AICS.list <- lapply(list.mods, FUN = function(x){x$aic})
    aics.df <- do.call(rbind, AICS.list)
    
    # get mcfaddens rsquared to look at model fit
    McFadden.rsq <- lapply(list.mods, FUN =  function(x){pscl::pR2(x)["McFadden"]})
    McFadden.rsq.df <- do.call(rbind, McFadden.rsq)
    
    
    # variable importance
    Var.importance.list <- lapply(list.mods, FUN = function(x){caret::varImp(x)})
    
    # AUC of predicted test data:
    library(ROCR)
    newdata <- test.covariate.data
    
    # m <- list.mods[[1]]
    #m <- glm.26
    get_AUC <- function(m, newdata = covariate.data){
      p <- predict(m, newdata=newdata, type="response")
      pr <- prediction(p, newdata$M)
      prf <- performance(pr, measure = "tpr", x.measure = "fpr")
      plot(prf)
      auc <- performance(pr, measure = "auc")
      auc <- auc@y.values[[1]]
      auc
    }
    AUC.lists <- list()
    for (p in 1:length(list.mods)){
      AUC.lists[[p]]<- get_AUC(m = list.mods[[p]], test.covariate.data)
    }
    AUC.lists <- lapply (list.mods, FUN = function(x){get_AUC(x, test.covariate.data)})
    AUC.df <- do.call(rbind, AUC.lists)
    
    model.diag <- data.frame(SPCD =  SPCD.df[i,]$SPCD,
                             Species = nspp[1:17, ] %>% filter(SPCD %in% SPCD.df[i,]$SPCD) %>% dplyr::select(COMMON),
                             model = paste0("model ", 1:56), 
                             model.no = as.numeric(1:56),
                             remper.correction = remper.cor.vector[j],
                             AUC = AUC.df[,1],
                             McFadden.Rsq = McFadden.rsq.df[,1], 
                             AIC = aics.df[,1], 
                             converged = convergence.df[,1])
    plot(model.diag$model.no, model.diag$AUC)
    plot(model.diag$model.no, model.diag$McFadden.Rsq)
    
    saveRDS(model.diag, paste0("SPCD_glm_output/GLM_model_diag_SPCD_", SPCD.df[i,]$SPCD, "_remp_", remper.cor.vector[j], ".RDS") )
    saveRDS(Var.importance.list, paste0("SPCD_glm_output/GLM_variable_importance_list_SPCD_", SPCD.df[i,]$SPCD, "_remp_", remper.cor.vector[j], ".RDS") )
    
    
    #--------------------------------------------------------------------------------------
    # PLOT UP VARIABLE IMPORTANCE FOR THE BEST FIT MODEL 
    #--------------------------------------------------------------------------------------
    AIC.best <- model.diag %>% filter(converged == TRUE)%>% mutate(minAIC = min(AIC))%>%  filter(AIC == minAIC)
    AIC.best$model.no
    
    AIC.best.varimp <- Var.importance.list[[AIC.best$model.no]]
    AIC.best.varimp$VARS <- rownames(AIC.best.varimp)
    AIC.best.ordered <- AIC.best.varimp %>% arrange(Overall)
    AIC.best.ordered$VARS <- factor(AIC.best.ordered$VARS, levels = unique(AIC.best.ordered$VARS))
    
    ggplot(AIC.best.ordered, aes(x = VARS, y = Overall))+geom_bar(stat = "identity")+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))+
      ylab("Variable Importance")+xlab("")+
      ggtitle(paste0("Variable Importance, ", nspp[1:17, ] %>% filter(SPCD %in% SPCD.df[i,]$SPCD) %>% dplyr::select(COMMON) , " model ", AIC.best$model.no))
    ggsave(filename = paste0("SPCD_glm_output/GLM_AIC_best_VARIMP_SPCD_", SPCD.df[i,]$SPCD, "_remp_", remper.cor.vector[j], ".png"), height = 5, width = 12)
    
    
    AIC.best$model.no
    
    AUC.best <- model.diag %>% filter(converged == TRUE)%>% mutate(maxAUC = max(AUC))%>%  filter(AUC == maxAUC)
    AUC.best$model.no
    
    AUC.best.varimp <- Var.importance.list[[AUC.best$model.no]]
    AUC.best.varimp$VARS <- rownames(AUC.best.varimp)
    AUC.best.ordered <- AUC.best.varimp %>% arrange(Overall)
    AUC.best.ordered$VARS <- factor(AUC.best.ordered$VARS, levels = unique(AUC.best.ordered$VARS))
    
    ggplot(AUC.best.ordered, aes(x = VARS, y = Overall))+geom_bar(stat = "identity")+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))+
      ylab("Variable Importance")+xlab("")+
      ggtitle(paste0("Variable Importance, ", nspp[1:17, ] %>% filter(SPCD %in% SPCD.df[i,]$SPCD) %>% dplyr::select(COMMON) , " model ", AUC.best$model.no))
    ggsave(filename = paste0("SPCD_glm_output/GLM_AUC_best_VARIMP_SPCD_", SPCD.df[i,]$SPCD, "_remp_", remper.cor.vector[j], ".png"), height = 5, width = 12)
    
    
    
    Rsq.best <- model.diag %>% filter(converged == TRUE)%>% mutate(maxRsq = max(McFadden.Rsq))%>%  filter(McFadden.Rsq == maxRsq)
    Rsq.best$model.no
    
    Rsq.best.varimp <- Var.importance.list[[Rsq.best$model.no]]
    Rsq.best.varimp$VARS <- rownames(Rsq.best.varimp)
    Rsq.best.ordered <- Rsq.best.varimp %>% arrange(Overall)
    Rsq.best.ordered$VARS <- factor(Rsq.best.ordered$VARS, levels = unique(Rsq.best.ordered$VARS))
    
    ggplot(Rsq.best.ordered, aes(x = VARS, y = Overall))+geom_bar(stat = "identity")+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))+
      ylab("Variable Importance")+xlab("")+
      ggtitle(paste0("Variable Importance, ", nspp[1:17, ] %>% filter(SPCD %in% SPCD.df[i,]$SPCD) %>% dplyr::select(COMMON) , " model ", Rsq.best$model.no))
    ggsave(filename = paste0("SPCD_glm_output/GLM_Rsq_best_VARIMP_SPCD_", SPCD.df[i,]$SPCD, "_remp_", remper.cor.vector[j], ".png"), height = 5, width = 12)
    
    ########################################################################
    # VIF for the species covariates
    #
    vif_values <-  vif(glm.20)
    VIF.cov <- data.frame(VIF = vif_values, 
                          covariates = names(vif_values))
    
    ggplot(VIF.cov, aes(x = covariates, y = VIF))+geom_bar(stat = "identity")+
      theme(axis.text.x = element_text(angle = 90))
    
    
    # Creating a correlation matrix
    corr <- round(cor(covariate.data[,3:ncol(covariate.data)]), 1)
    library(ggcorrplot)
    # generate correlation plots here:
    ggcorrplot(corr, hc.order = TRUE, type = "lower",
               lab = TRUE)
    ggsave(filename = paste0("SPCD_glm_output/GLM_Correlation_matrix", SPCD.df[i,]$SPCD, "_remp_", remper.cor.vector[j], ".png"), height = 12, width = 12)
    saveRDS(corr, paste0("SPCD_glm_output/GLM_Correlation_matrix_SPCD_",SPCD.id,"_predictors.rds") )
  }
}



# make a big table with the model results:
model.diags <- list()
for(i in 1:length(SPCD.df[,]$SPCD)){
  SPCD.id <- SPCD.df[i,]$SPCD
  # if(SPCD.id == 621){}else{
  model.diags[[i]] <- readRDS(paste0("SPCD_glm_output/GLM_model_diag_SPCD_", SPCD.df[i,]$SPCD, "_remp_", remper.cor.vector[j], ".RDS") )
  #}
}
model.diag <- do.call(rbind, model.diags)
model.diag %>% group_by(SPCD, COMMON)|> gt()

ggplot(model.diag %>% filter(converged == TRUE), aes(x = model.no, y = AUC, fill = COMMON, group = model.no))+geom_bar(stat= "identity",position = position_dodge2())#+position_dodge()
ggsave("SPCD_glm_output/GLM_all_species_AUC.png", height = 5, width = 8)

ggplot(model.diag %>% filter(converged == TRUE), aes(x = model.no, y = McFadden.Rsq))+geom_bar(stat= "identity") + facet_wrap(~COMMON, scales = "free_y")
ggsave("SPCD_glm_output/GLM_all_species_McFaddenRsq.png", height = 5, width = 8)

ggplot(model.diag %>% filter(converged == TRUE), aes(x = model.no, y = AUC))+geom_bar(stat= "identity") + facet_wrap(~COMMON, scales = "free_y")
ggsave("SPCD_glm_output/GLM_all_species_AUC.png", height = 5, width = 8)

ggplot(model.diag %>% filter(converged == TRUE), aes(x = model.no, y = AIC))+geom_bar(stat= "identity") + facet_wrap(~COMMON, scales = "free_y")
ggsave("SPCD_glm_output/GLM_all_species_AIC.png", height = 5, width = 8)

ggplot(model.diag %>% filter(converged ==TRUE), aes(AUC, AIC,  label = as.character(model.no)))+geom_text()+facet_wrap(~SPCD, scales = "free_y")
ggsave("SPCD_glm_output/GLM_all_species_AIC_AUC.png", height = 5, width = 8)

ggplot(model.diag %>% filter(converged ==TRUE), aes(AUC, McFadden.Rsq,  label = as.character(model.no)))+geom_text(size = 2)+facet_wrap(~SPCD, scales = "free_y")
ggsave("SPCD_glm_output/GLM_all_species_AUC_Rsq.png", height = 5, width = 8)

# make a table explaining the models:


glm.model.table <- data.frame(model.no = 1:49, 
                              description = c("MAP", 
                                              "MATmaxanom", 
                                              "MATminanom", 
                                              "MAPanom", 
                                              "BAL", 
                                              "percent damage", 
                                              "site index", 
                                              "Physiographic class", 
                                              "BA", 
                                              "Relative Density",
                                              "elevation", 
                                              "N depostion (wet + dry)", 
                                              "diameter difference", 
                                              
                                              ## sequentially adding in each variable
                                              "diameter difference + Diameter", 
                                              "exp(Diameter)", 
                                              "aspect", 
                                              "slope", 
                                              "MATmax", 
                                              "MATmin", 
                                              "MAP", 
                                              "MATmax anomaly", 
                                              "MATmin anomaly", 
                                              "MAP anomaly", 
                                              "BAL", 
                                              "percent damage", 
                                              "site index", 
                                              "physiographic class", 
                                              "BA", 
                                              "Relative Density", 
                                              "Elevation", 
                                              "N deposition", 
                                              ## adding in interactions
                                              "growth interactions", 
                                              "diameter interactions", 
                                              "aspect interactions", 
                                              " slope interactions", 
                                              "MATmax interactions", 
                                              "MATmin interactions", 
                                              "MAP interactions", 
                                              "MATmax anomaly interactions", 
                                              "MATmin anomaly interactions", 
                                              "MAP anomaly interactions", 
                                              "BAL interactions", 
                                              "percent damage interactions", 
                                              "site index interactions", 
                                              "Physiographic interactions", 
                                              "Relative Density interactions", 
                                              "elevation interactions",
                                              "Basal Area interactions", 
                                              "N dep interactions"
                              ), 
                              model.type = c(rep("single variable", 13), 
                                             rep("adding on to growth effect", 18), 
                                             rep("adding interaction terms", 18
                                             )))
model.diag <- left_join(model.diag, glm.model.table)

ggplot(model.diag %>% filter(converged == TRUE), aes(x = model.no, y = AUC, fill = COMMON, group = model.no))+geom_bar(stat= "identity",position = position_dodge2())#+position_dodge()
ggsave("SPCD_glm_output/GLM_all_species_AUC.png", height = 5, width = 8)

ggplot(model.diag %>% filter(converged == TRUE), aes(x = model.no, y = McFadden.Rsq, fill = model.type))+geom_bar(stat= "identity") + facet_wrap(~COMMON, scales = "free_y")
ggsave("SPCD_glm_output/GLM_all_species_McFaddenRsq.png", height = 5, width = 8)

ggplot(model.diag %>% filter(converged == TRUE), aes(x = model.no, y = AUC, fill = model.type))+geom_bar(stat= "identity") + facet_wrap(~COMMON, scales = "free_y")
ggsave("SPCD_glm_output/GLM_all_species_AUC.png", height = 5, width = 8)

ggplot(model.diag %>% filter(converged == TRUE), aes(x = model.no, y = AIC, fill = model.type))+geom_bar(stat= "identity") + facet_wrap(~COMMON, scales = "free_y")
ggsave("SPCD_glm_output/GLM_all_species_AIC.png", height = 5, width = 8)

ggplot(model.diag %>% filter(converged ==TRUE), aes(AUC, AIC,  label = as.character(model.no)))+geom_text()+facet_wrap(~SPCD, scales = "free_y")
ggsave("SPCD_glm_output/GLM_all_species_AIC_AUC.png", height = 5, width = 8)

ggplot(model.diag %>% filter(converged ==TRUE), aes(AUC, McFadden.Rsq,  label = as.character(model.no)))+geom_text(size = 2)+facet_wrap(~SPCD, scales = "free_y")
ggsave("SPCD_glm_output/GLM_all_species_AUC_Rsq.png", height = 5, width = 8)

write.csv(glm.model.table, "GLM_table.csv", quote = TRUE)

# get all of the correlated variables for all species--are there common ones that we should eliminate
inter.cov.cors <- list()
for (i in 1:17){
  #read in correlation matrix
  spcd.correlations <- readRDS(paste0("SPCD_glm_output/GLM_Correlation_matrix_SPCD_",nspp[i,]$SPCD,"_predictors.rds"))
  spcd.correlations[lower.tri(spcd.correlations)] <- NA # set the lower triangle to NA to filter out
  
  cor_df <- reshape2::melt(spcd.correlations, varnames = c("Variable1", "Variable2"), value.name = "Correlation") %>% 
                      filter(!is.na(Correlation) & !Correlation == 1) # filter out all NA values and the diagnal
  
  # how many variable combinations have high correlations?
  cor_df %>% filter(abs(Correlation) >= 0.5)
  
  # save in a list with the rest of the species
  inter.cov.corrs[[i]] <- cor_df %>% mutate(SPCD = nspp[i,]$SPCD, 
                    Species = nspp[i,]$Species) %>% select(SPCD, Species, Variable1, Variable2, Correlation) %>% 
                                     mutate(high.corr.F = ifelse(abs(Correlation) >= 0.5, "High", "Low"))
  
}

inter.cov.corrs.df <- do.call(rbind, inter.cov.corrs)

inter.cov.corrs.df %>% filter(abs(Correlation) >= 0.5)

# candidate variables for removal in the models
inter.cov.corrs.df %>% filter(high.corr.F == "High")  %>% group_by(Variable1) %>% summarise(n()) %>% arrange(desc(`n()`))
inter.cov.corrs.df %>% filter(high.corr.F == "High")  %>% group_by(Variable2) %>% summarise(n()) %>% arrange(desc(`n()`))

# what does the correlation look like for all the data?
# to do this we need to load in all the species dataframes
covariate.all <- list()
test.covariate.all <- list()
for(i in 1:length(SPCD.df$SPCD)){
  SPCD.id <- SPCD.df[i,]$SPCD
  common.name <- nspp[1:17,] %>% filter(SPCD %in% SPCD.id) %>% dplyr::select(COMMON)
  
  #if(SPCD.id == 621){
  # cat("Not running for yellow poplar, not enough data")
  #}else{
  

    
    remper.correction <- 0.5
    load(paste0("SPCD_GLM_standata/SPCD_", SPCD.id, "remper_correction_", remper.correction, ".Rdata")) # load the species code data
    covariate.data <- data.frame(M = mod.data$y, 
                                 SPCD = SPCD.id,
                                 annual.growth = mod.data$annual_growth, 
                                 dia.diff = mod.data$DIA.diff,
                                 DIA = mod.data$DIA, 
                                 si = mod.data$si, 
                                 slope = mod.data$slope, 
                                 aspect = mod.data$aspect, 
                                 MATmax = mod.data$MATmax, 
                                 MATmin= mod.data$MATmin, 
                                 MAP = mod.data$MAP, 
                                 BAL = mod.data$BAL, 
                                 damage = mod.data$damage, 
                                 MAPanom = mod.data$MAPanom, 
                                 MATmaxanom = mod.data$MATmaxanom, 
                                 MATminanom = mod.data$MATminanom, 
                                 PHYSIO = mod.data$PHYSIO, 
                                 BA = mod.data$BA, 
                                 RD = mod.data$RD, 
                                 elev = mod.data$elev, 
                                 Ndep = mod.data$Ndep, 
                                 SPCD.BA = mod.data$SPCD.BA,
                                 non_SPCD.BA = mod.data$non_SPCD.BA.scaled,
                                 prop.focal.ba = mod.data$prop.focal.ba.scaled 
    )
    
    test.covariate.data <- data.frame(M = test.data$M,
                                      annual.growth = mod.data$annual_growth_test, 
                                      SPCD = SPCD.id, 
                                      dia.diff = mod.data$DIA.diff_test,
                                      DIA = mod.data$DIA_test, 
                                      si = mod.data$si_test, 
                                      slope = mod.data$slope_test, 
                                      aspect = mod.data$aspect_test, 
                                      MATmax = mod.data$MATmax_test, 
                                      MATmin= mod.data$MATmin_test, 
                                      MAP = mod.data$MAP_test, 
                                      BAL = mod.data$BAL_test, 
                                      damage = mod.data$damage_test, 
                                      MAPanom = mod.data$MAPanom_test, 
                                      MATmaxanom = mod.data$MATmaxanom_test, 
                                      MATminanom = mod.data$MATminanom_test, 
                                      PHYSIO = mod.data$PHYSIO_test, 
                                      BA = mod.data$PHYSIO_test, 
                                      RD = mod.data$RD_test, 
                                      elev = mod.data$elev_test, 
                                      Ndep = mod.data$Ndep_test, 
                                      SPCD.BA = mod.data$SPCD.BA_test,
                                      non_SPCD.BA = mod.data$non_SPCD.BA.scaled_test,
                                      prop.focal.ba = mod.data$prop.focal.ba.scaled_test )
    
    covariate.all[[i]] <- covariate.data
    test.covariate.all[[i]] <- test.covariate.data
    
}
covariate.all.df <- do.call(rbind, covariate.all)
test.covariate.all.df <- do.call(rbind, test.covariate.all)


# Creating a correlation matrix
all.correlations <- round(cor(covariate.all.df[,4:ncol(covariate.all.df)]), 1)

# generate correlation plots here:
ggcorrplot(all.correlations, hc.order = TRUE, type = "lower",
           lab = TRUE)
ggsave(filename = paste0("SPCD_glm_output/GLM_Correlation_matrix_allspecies_remp_", remper.cor.vector[j], ".png"), height = 12, width = 12)
saveRDS(all.correlations, paste0("SPCD_glm_output/GLM_Correlation_matrix_SPCD_all_species_predictors.rds") )


# find covariate combinations for all the species that could be problematic:

all.correlations[lower.tri(all.correlations)] <- NA # set the lower triangle to NA to filter out

cor_df <- reshape2::melt(all.correlations, varnames = c("Variable1", "Variable2"), value.name = "Correlation") %>% 
  filter(!is.na(Correlation) & !Correlation == 1) # filter out all NA values and the diagnal

# how many variable combinations have high correlations?
cor_df %>% filter(abs(Correlation) >= 0.5)
# Candidates for removal and # of problem combinations they would elminate:
# 1. MATmin (2)
# 2. MATminanom (1)
# 3. RD (4)
# 3. SPCD.BA (2-3)


# save in a list with the rest of the species
inter.cov.corrs.comb <- cor_df %>% select(Variable1, Variable2, Correlation) %>% 
  mutate(high.corr.F = ifelse(abs(Correlation) >= 0.5, "High", "Low"))

# candidate variables for removal in the models
inter.cov.corrs.comb %>% filter(high.corr.F == "High")  %>% group_by(Variable1) %>% summarise(n()) %>% arrange(desc(`n()`))
inter.cov.corrs.comb %>% filter(high.corr.F == "High")  %>% group_by(Variable2) %>% summarise(n()) %>% arrange(desc(`n()`))

#####################################################################################
# redo corrations with out the problem variables
#####################################################################################
covariate.reduced.df <- covariate.all.df %>% select(-MATmin, -MATminanom, -RD, -SPCD.BA)
# Creating a correlation matrix
reduced.correlations <- round(cor(covariate.reduced.df[,4:ncol(covariate.reduced.df)]), 1)

cor_matrix <- Hmisc::rcorr(as.matrix(covariate.reduced.df[,4:ncol(covariate.reduced.df)]), type="pearson")
# most are significant, so I will focus on coefficient values


# generate correlation plots here:
ggcorrplot(reduced.correlations, hc.order = TRUE, type = "lower",lab = TRUE)
ggsave(filename = paste0("SPCD_glm_output/GLM_reduced_Correlation_matrix_allspecies_remp_", remper.cor.vector[j], ".png"), height = 12, width = 12)
saveRDS(all.correlations, paste0("SPCD_glm_output/GLM_reduced_Correlation_matrix_SPCD_all_species_predictors.rds") )


# find covariate combinations for all the species that could be problematic:
reduced.correlations[lower.tri(reduced.correlations)] <- NA # set the lower triangle to NA to filter out

red.cor_df <- reshape2::melt(reduced.correlations, varnames = c("Variable1", "Variable2"), value.name = "Correlation") %>% 
  filter(!is.na(Correlation) & !Correlation == 1) # filter out all NA values and the diagnal

# Now only one variable remains with higher correlations MATmax and Ndep for the full dataset
red.cor_df %>% filter(abs(Correlation) >= 0.5)

# lets see if there are any issues with some species that remain:
problem.variables.corr <- c("MATmin", "MATminanom", "RD", "SPCD.BA")
# candidate variables for removal in the models
inter.cov.corrs.df %>% 
  filter(!Variable1 %in% problem.variables) %>% 
  filter(!Variable2 %in% problem.variables) %>%
  filter(abs(Correlation) >= 0.6)

problem.variables.corr.species <- c("MATmin", "MATminanom", "RD", "SPCD.BA", # the same as above
                                    "non_SPCD.BA", "prop.focal.ba", "si", "elev") # but also more competition variables

# candidate variables for removal in the models
inter.cov.corrs.df %>% 
  filter(!Variable1 %in% problem.variables.corr.species) %>% 
  filter(!Variable2 %in% problem.variables.corr.species) %>%
  filter(abs(Correlation) >= 0.6)


######################################################################################
# GLMs for the non-correlated variables
######################################################################################

for(i in 1:length(SPCD.df$SPCD)){
  SPCD.id <- SPCD.df[i,]$SPCD
  common.name <- nspp[1:17,] %>% filter(SPCD %in% SPCD.id) %>% dplyr::select(COMMON)
  
  #if(SPCD.id == 621){
  # cat("Not running for yellow poplar, not enough data")
  #}else{
  
  for (j in 1:length(remper.cor.vector)){
    cat(paste("running glm mortality model for SPCD", SPCD.df[i,]$SPCD, common.name$COMMON, " remper correction", remper.cor.vector[j]))
    
    remper.correction <- remper.cor.vector[j]
    load(paste0("SPCD_GLM_standata/SPCD_", SPCD.id, "remper_correction_", remper.correction, ".Rdata")) # load the species code data
    covariate.data <- data.frame(M = mod.data$y, 
                                 annual.growth = mod.data$annual_growth, 
                                 dia.diff = mod.data$DIA.diff,
                                 DIA = mod.data$DIA, 
                                 si = mod.data$si, 
                                 slope = mod.data$slope, 
                                 aspect = mod.data$aspect, 
                                 MATmax = mod.data$MATmax, 
                                 MATmin= mod.data$MATmin, 
                                 MAP = mod.data$MAP, 
                                 BAL = mod.data$BAL, 
                                 damage = mod.data$damage, 
                                 MAPanom = mod.data$MAPanom, 
                                 MATmaxanom = mod.data$MATmaxanom, 
                                 MATminanom = mod.data$MATminanom, 
                                 PHYSIO = mod.data$PHYSIO, 
                                 BA = mod.data$BA, 
                                 RD = mod.data$RD, 
                                 elev = mod.data$elev, 
                                 Ndep = mod.data$Ndep, 
                                 SPCD.BA = mod.data$SPCD.BA,
                                 non_SPCD.BA = mod.data$non_SPCD.BA.scaled,
                                 prop.focal.ba = mod.data$prop.focal.ba.scaled 
    ) %>% select(!problem.variables.corr.species)
    
    test.covariate.data <- data.frame(M = test.data$M,
                                      annual.growth = mod.data$annual_growth_test, 
                                      dia.diff = mod.data$DIA.diff_test,
                                      DIA = mod.data$DIA_test, 
                                      si = mod.data$si_test, 
                                      slope = mod.data$slope_test, 
                                      aspect = mod.data$aspect_test, 
                                      MATmax = mod.data$MATmax_test, 
                                      MATmin= mod.data$MATmin_test, 
                                      MAP = mod.data$MAP_test, 
                                      BAL = mod.data$BAL_test, 
                                      damage = mod.data$damage_test, 
                                      MAPanom = mod.data$MAPanom_test, 
                                      MATmaxanom = mod.data$MATmaxanom_test, 
                                      MATminanom = mod.data$MATminanom_test, 
                                      PHYSIO = mod.data$PHYSIO_test, 
                                      BA = mod.data$PHYSIO_test, 
                                      RD = mod.data$RD_test, 
                                      elev = mod.data$elev_test, 
                                      Ndep = mod.data$Ndep_test, 
                                      SPCD.BA = mod.data$SPCD.BA_test,
                                      non_SPCD.BA = mod.data$non_SPCD.BA.scaled_test,
                                      prop.focal.ba = mod.data$prop.focal.ba.scaled_test )%>% 
      select(!problem.variables.corr.species)
    
    
    
    glm.A <-  glm(M ~ MAP , data = covariate.data, family = "binomial")
    glm.B <-  glm(M ~ MATmaxanom, data = covariate.data, family = "binomial")
    
    glm.C <-  glm(M ~ MAPanom, data = covariate.data, family = "binomial")
    glm.D <-  glm(M ~ BAL, data = covariate.data, family = "binomial")
    glm.E <-  glm(M ~ damage , data = covariate.data, family = "binomial")
    
    glm.F <-  glm(M ~ PHYSIO, data = covariate.data, family = "binomial")
    glm.G <-  glm(M ~ BA, data = covariate.data, family = "binomial")
    
    glm.H <-  glm(M ~ Ndep, data = covariate.data, family = "binomial")
    
    
    
    
    glm.1 <-  glm(M ~ dia.diff, data = covariate.data, family = "binomial")
    glm.2 <-  glm(M ~ dia.diff + DIA, data = covariate.data, family = "binomial")
    glm.3 <-  glm(M ~ dia.diff + DIA + exp(DIA), data = covariate.data, family = "binomial")
    
    glm.4 <-  glm(M ~ dia.diff + DIA + exp(DIA)+ aspect, data = covariate.data, family = "binomial")
    glm.5 <-  glm(M ~ dia.diff + DIA + exp(DIA)+ aspect + slope , data = covariate.data, family = "binomial")
    glm.6 <-  glm(M ~ dia.diff + DIA + exp(DIA) + aspect + slope + MATmax , data = covariate.data, family = "binomial")
    
    
    glm.7 <-  glm(M ~ dia.diff + DIA + exp(DIA)+ aspect + slope + MATmax 
                  + MATmax , data = covariate.data, family = "binomial")
    glm.8 <-  glm(M ~ dia.diff + DIA + exp(DIA) + aspect + slope + MATmax 
                  + MAP , data = covariate.data, family = "binomial")
    
    glm.9 <-  glm(M ~ dia.diff + DIA + exp(DIA)+ aspect + slope + MATmax 
                  + MAP + MATmaxanom , data = covariate.data, family = "binomial")
    
    glm.10 <-  glm(M ~ dia.diff + DIA + exp(DIA)+ aspect + slope + MATmax 
                   + MAP + MATmaxanom 
                   , data = covariate.data, family = "binomial")
    
    
    glm.11 <-  glm(M ~ dia.diff + DIA + exp(DIA)+ aspect + slope + MATmax 
                   + MAP + MATmaxanom 
                   + MAPanom , data = covariate.data, family = "binomial")
    
    glm.12 <-  glm(M ~ dia.diff + DIA + exp(DIA)+ aspect + slope + MATmax 
                   + MAP + MATmaxanom 
                   + MAPanom + BAL , data = covariate.data, family = "binomial")
    
    glm.13<-  glm(M ~ dia.diff + DIA + exp(DIA)+ aspect + slope + MATmax 
                  + MAP + MATmaxanom 
                  + MAPanom + BAL + damage, data = covariate.data, family = "binomial")
    
    
    glm.14 <-  glm(M ~ dia.diff + DIA + exp(DIA)+ aspect + slope + MATmax 
                   + MAP + MATmaxanom 
                   + MAPanom + BAL + damage  + PHYSIO, data = covariate.data, family = "binomial")
    
    glm.15 <-  glm(M ~ dia.diff + DIA + exp(DIA)+ aspect + slope + MATmax 
                   + MAP + MATmaxanom 
                   + MAPanom + BAL + damage  +  PHYSIO + BA , data = covariate.data, family = "binomial")
    
    
    glm.16 <-  glm(M ~ dia.diff + DIA + exp(DIA)+ aspect + slope + MATmax 
                   + MAP + MATmaxanom 
                   + MAPanom + BAL + damage +  PHYSIO + BA , data = covariate.data, family = "binomial")
    
    glm.17 <-  glm(M ~ dia.diff + DIA + exp(DIA)+ aspect + slope + MATmax 
                   + MAP + MATmaxanom 
                   + MAPanom + BAL + damage  +  PHYSIO + BA + Ndep, data = covariate.data, family = "binomial")
    
    
    
    
    glm.18 <-  glm(M ~ dia.diff + DIA + exp(DIA)+ aspect + slope + MATmax 
                   + MAP + MATmaxanom 
                   + MAPanom + BAL + damage  +  PHYSIO + BA  + Ndep + 
                     #SPCD.BA + non_SPCD.BA + prop.focal.ba +
                     # annual growth interactcovarions
                     DIA*dia.diff + aspect*dia.diff  + slope*dia.diff  + MATmax*dia.diff  
                   + MAP*dia.diff  + MATmaxanom*dia.diff  
                   + MAPanom*dia.diff  + BAL*dia.diff  + damage*dia.diff   +
                     +  PHYSIO*dia.diff + BA*dia.diff +  Ndep*dia.diff 
                   , data = covariate.data, family = "binomial")
    
    
    
    # all growth + diameter interaction terms
    glm.19 <-  glm(M ~ dia.diff + DIA + exp(DIA)+ aspect + slope + MATmax 
                   + MAP + MATmaxanom 
                   + MAPanom + BAL + damage  +  PHYSIO + BA  + Ndep + 
                     #SPCD.BA + non_SPCD.BA + prop.focal.ba +
                     # annual growth interactcovarions
                     DIA*dia.diff + aspect*dia.diff  + slope*dia.diff  + MATmax*dia.diff  
                   + MAP*dia.diff  + MATmaxanom*dia.diff  
                   + MAPanom*dia.diff  + BAL*dia.diff  + damage*dia.diff   +
                     +  PHYSIO*dia.diff + BA*dia.diff +  Ndep*dia.diff+
                     # all diameter interactions
                     aspect*DIA  + slope*DIA  + MATmax*DIA  
                   + MAP*DIA  + MATmaxanom*DIA  
                   + MAPanom*DIA  + BAL*DIA  + damage*DIA  +  
                     PHYSIO*DIA + BA*DIA  + Ndep*DIA , data = covariate.data, family = "binomial")
    
    glm.20 <-  glm(M ~ dia.diff + DIA + exp(DIA)+ aspect + slope + MATmax 
                   + MAP + MATmaxanom 
                   + MAPanom + BAL + damage  +  PHYSIO + BA  + Ndep + 
                     #SPCD.BA + non_SPCD.BA + prop.focal.ba +
                     # annual growth interactcovarions
                     DIA*dia.diff + aspect*dia.diff  + slope*dia.diff  + MATmax*dia.diff  
                   + MAP*dia.diff  + MATmaxanom*dia.diff  
                   + MAPanom*dia.diff  + BAL*dia.diff  + damage*dia.diff   
                   +  PHYSIO*dia.diff + BA*dia.diff +  Ndep*dia.diff+
                     # all diameter interactions
                     aspect*DIA  + slope*DIA  + MATmax*DIA  
                   + MAP*DIA  + MATmaxanom*DIA  
                   + MAPanom*DIA  + BAL*DIA  + damage*DIA  +  
                     PHYSIO*DIA + BA*DIA  + Ndep*DIA + 
                     # all aspect interactions
                     slope*aspect  + MATmax*aspect  
                   + MAP*aspect  + MATmaxanom*aspect  
                   + MAPanom*aspect  + BAL*aspect  + damage*aspect  + 
                     PHYSIO*aspect + BA*aspect   + Ndep*aspect, data = covariate.data, family = "binomial")
    
    
    
    glm.21 <-  glm(M ~ dia.diff + DIA + exp(DIA)+ aspect + slope + MATmax 
                   + MAP + MATmaxanom 
                   + MAPanom + BAL + damage  +  PHYSIO + BA  + Ndep + 
                     #SPCD.BA + non_SPCD.BA + prop.focal.ba +
                     # annual growth interactcovarions
                     DIA*dia.diff + aspect*dia.diff  + slope*dia.diff  + MATmax*dia.diff  
                   + MAP*dia.diff  + MATmaxanom*dia.diff  
                   + MAPanom*dia.diff  + BAL*dia.diff  + damage*dia.diff   
                   +  PHYSIO*dia.diff + BA*dia.diff +  Ndep*dia.diff+
                     # all diameter interactions
                     aspect*DIA  + slope*DIA  + MATmax*DIA  
                   + MAP*DIA  + MATmaxanom*DIA  
                   + MAPanom*DIA  + BAL*DIA  + damage*DIA  +  
                     PHYSIO*DIA + BA*DIA  + Ndep*DIA + 
                     # all aspect interactions
                     slope*aspect  + MATmax*aspect  
                   + MAP*aspect  + MATmaxanom*aspect  
                   + MAPanom*aspect  + BAL*aspect  + damage*aspect  + 
                     PHYSIO*aspect + BA*aspect   + Ndep*aspect+
                     # all slope interactions
                     MATmax*slope  
                   + MAP*slope  + MATmaxanom*slope  
                   + MAPanom*slope  + BAL*slope  + damage*slope  +
                     PHYSIO*slope + BA*slope  + Ndep*slope , data = covariate.data, family = "binomial")
    
    glm.22 <-  glm(M ~ dia.diff + DIA + exp(DIA)+ aspect + slope + MATmax 
                   + MAP + MATmaxanom 
                   + MAPanom + BAL + damage  +  PHYSIO + BA  + Ndep + 
                     #SPCD.BA + non_SPCD.BA + prop.focal.ba +
                     # annual growth interactcovarions
                     DIA*dia.diff + aspect*dia.diff  + slope*dia.diff  + MATmax*dia.diff  
                   + MAP*dia.diff  + MATmaxanom*dia.diff  
                   + MAPanom*dia.diff  + BAL*dia.diff  + damage*dia.diff   
                   +  PHYSIO*dia.diff + BA*dia.diff +  Ndep*dia.diff+
                     # all diameter interactions
                     aspect*DIA  + slope*DIA  + MATmax*DIA  
                   + MAP*DIA  + MATmaxanom*DIA  
                   + MAPanom*DIA  + BAL*DIA  + damage*DIA  +  
                     PHYSIO*DIA + BA*DIA  + Ndep*DIA + 
                     # all aspect interactions
                     slope*aspect  + MATmax*aspect  
                   + MAP*aspect  + MATmaxanom*aspect  
                   + MAPanom*aspect  + BAL*aspect  + damage*aspect  + 
                     PHYSIO*aspect + BA*aspect   + Ndep*aspect+
                     # all slope interactions
                     MATmax*slope  
                   + MAP*slope  + MATmaxanom*slope  
                   + MAPanom*slope  + BAL*slope  + damage*slope  +
                     PHYSIO*slope + BA*slope  + Ndep*slope +
                     # all MATmax interactions
                     
                     MAP*MATmax  + MATmaxanom*MATmax  
                   + MAPanom*MATmax  + BAL*MATmax  + damage*MATmax  + 
                     PHYSIO*MATmax + BA*MATmax +   Ndep*MATmax , data = covariate.data, family = "binomial")
    
    
    
    
    
    glm.23 <-  glm(M ~ dia.diff + DIA + exp(DIA)+ aspect + slope + MATmax 
                   + MAP + MATmaxanom 
                   + MAPanom + BAL + damage  +  PHYSIO + BA  + Ndep + 
                     #SPCD.BA + non_SPCD.BA + prop.focal.ba +
                     # annual growth interactcovarions
                     DIA*dia.diff + aspect*dia.diff  + slope*dia.diff  + MATmax*dia.diff  
                   + MAP*dia.diff  + MATmaxanom*dia.diff  
                   + MAPanom*dia.diff  + BAL*dia.diff  + damage*dia.diff   
                   +  PHYSIO*dia.diff + BA*dia.diff +  Ndep*dia.diff+
                     # all diameter interactions
                     aspect*DIA  + slope*DIA  + MATmax*DIA  
                   + MAP*DIA  + MATmaxanom*DIA  
                   + MAPanom*DIA  + BAL*DIA  + damage*DIA  +  
                     PHYSIO*DIA + BA*DIA  + Ndep*DIA + 
                     # all aspect interactions
                     slope*aspect  + MATmax*aspect  
                   + MAP*aspect  + MATmaxanom*aspect  
                   + MAPanom*aspect  + BAL*aspect  + damage*aspect  + 
                     PHYSIO*aspect + BA*aspect   + Ndep*aspect+
                     # all slope interactions
                     MATmax*slope  
                   + MAP*slope  + MATmaxanom*slope  
                   + MAPanom*slope  + BAL*slope  + damage*slope  +
                     PHYSIO*slope + BA*slope  + Ndep*slope +
                     # all MATmax interactions
                     
                     MAP*MATmax  + MATmaxanom*MATmax  
                   + MAPanom*MATmax  + BAL*MATmax  + damage*MATmax  + 
                     PHYSIO*MATmax + BA*MATmax +   Ndep*MATmax+ 
                     
                     # all MAP interatction
                     
                     MATmaxanom*MAP + MAPanom*MAP + BAL*MAP + damage*MAP + 
                     PHYSIO*MAP + BA*MAP +  Ndep*MAP , data = covariate.data, family = "binomial")
    
    
    
    glm.24 <-  glm(M ~ dia.diff + DIA + exp(DIA)+ aspect + slope + MATmax 
                   + MAP + MATmaxanom 
                   + MAPanom + BAL + damage  +  PHYSIO + BA  + Ndep + 
                     #SPCD.BA + non_SPCD.BA + prop.focal.ba +
                     # annual growth interactcovarions
                     DIA*dia.diff + aspect*dia.diff  + slope*dia.diff  + MATmax*dia.diff  
                   + MAP*dia.diff  + MATmaxanom*dia.diff  
                   + MAPanom*dia.diff  + BAL*dia.diff  + damage*dia.diff   
                   +  PHYSIO*dia.diff + BA*dia.diff +  Ndep*dia.diff+
                     # all diameter interactions
                     aspect*DIA  + slope*DIA  + MATmax*DIA  
                   + MAP*DIA  + MATmaxanom*DIA  
                   + MAPanom*DIA  + BAL*DIA  + damage*DIA  +  
                     PHYSIO*DIA + BA*DIA  + Ndep*DIA + 
                     # all aspect interactions
                     slope*aspect  + MATmax*aspect  
                   + MAP*aspect  + MATmaxanom*aspect  
                   + MAPanom*aspect  + BAL*aspect  + damage*aspect  + 
                     PHYSIO*aspect + BA*aspect   + Ndep*aspect+
                     # all slope interactions
                     MATmax*slope  
                   + MAP*slope  + MATmaxanom*slope  
                   + MAPanom*slope  + BAL*slope  + damage*slope  +
                     PHYSIO*slope + BA*slope  + Ndep*slope +
                     # all MATmax interactions
                     
                     MAP*MATmax  + MATmaxanom*MATmax  
                   + MAPanom*MATmax  + BAL*MATmax  + damage*MATmax  + 
                     PHYSIO*MATmax + BA*MATmax +   Ndep*MATmax+ 
                     
                     # all MAP interatction
                     
                     MATmaxanom*MAP + MAPanom*MAP + BAL*MAP + damage*MAP + 
                     PHYSIO*MAP + BA*MAP +  Ndep*MAP + 
                     # all MATmaxanom interatction
                     
                     MAPanom*MATmaxanom + BAL*MATmaxanom + damage*MATmaxanom + 
                     PHYSIO*MATmaxanom + BA*MATmaxanom + Ndep*MATmaxanom 
                   , data = covariate.data, family = "binomial")
    
    
    
    glm.25 <-  glm(M ~ dia.diff + DIA + exp(DIA)+ aspect + slope + MATmax 
                   + MAP + MATmaxanom 
                   + MAPanom + BAL + damage  +  PHYSIO + BA  + Ndep + 
                     #SPCD.BA + non_SPCD.BA + prop.focal.ba +
                     # annual growth interactcovarions
                     DIA*dia.diff + aspect*dia.diff  + slope*dia.diff  + MATmax*dia.diff  
                   + MAP*dia.diff  + MATmaxanom*dia.diff  
                   + MAPanom*dia.diff  + BAL*dia.diff  + damage*dia.diff   
                   +  PHYSIO*dia.diff + BA*dia.diff +  Ndep*dia.diff+
                     # all diameter interactions
                     aspect*DIA  + slope*DIA  + MATmax*DIA  
                   + MAP*DIA  + MATmaxanom*DIA  
                   + MAPanom*DIA  + BAL*DIA  + damage*DIA  +  
                     PHYSIO*DIA + BA*DIA  + Ndep*DIA + 
                     # all aspect interactions
                     slope*aspect  + MATmax*aspect  
                   + MAP*aspect  + MATmaxanom*aspect  
                   + MAPanom*aspect  + BAL*aspect  + damage*aspect  + 
                     PHYSIO*aspect + BA*aspect   + Ndep*aspect+
                     # all slope interactions
                     MATmax*slope  
                   + MAP*slope  + MATmaxanom*slope  
                   + MAPanom*slope  + BAL*slope  + damage*slope  +
                     PHYSIO*slope + BA*slope  + Ndep*slope +
                     # all MATmax interactions
                     
                     MAP*MATmax  + MATmaxanom*MATmax  
                   + MAPanom*MATmax  + BAL*MATmax  + damage*MATmax  + 
                     PHYSIO*MATmax + BA*MATmax +   Ndep*MATmax+ 
                     
                     # all MAP interatction
                     
                     MATmaxanom*MAP + MAPanom*MAP + BAL*MAP + damage*MAP + 
                     PHYSIO*MAP + BA*MAP +  Ndep*MAP + 
                     # all MATmaxanom interatction
                     
                     MAPanom*MATmaxanom + BAL*MATmaxanom + damage*MATmaxanom + 
                     PHYSIO*MATmaxanom + BA*MATmaxanom + Ndep*MATmaxanom +
                     # all MAPanom interatction
                     
                     BAL*MAPanom + damage*MAPanom + 
                     PHYSIO*MAPanom + BA*MAPanom  + Ndep*MAPanom , data = covariate.data, family = "binomial")
    
    glm.26 <-  glm(M ~ dia.diff + DIA + exp(DIA)+ aspect + slope + MATmax 
                   + MAP + MATmaxanom 
                   + MAPanom + BAL + damage  +  PHYSIO + BA  + Ndep + 
                     #SPCD.BA + non_SPCD.BA + prop.focal.ba +
                     # annual growth interactcovarions
                     DIA*dia.diff + aspect*dia.diff  + slope*dia.diff  + MATmax*dia.diff  
                   + MAP*dia.diff  + MATmaxanom*dia.diff  
                   + MAPanom*dia.diff  + BAL*dia.diff  + damage*dia.diff   
                   +  PHYSIO*dia.diff + BA*dia.diff +  Ndep*dia.diff+
                     # all diameter interactions
                     aspect*DIA  + slope*DIA  + MATmax*DIA  
                   + MAP*DIA  + MATmaxanom*DIA  
                   + MAPanom*DIA  + BAL*DIA  + damage*DIA  +  
                     PHYSIO*DIA + BA*DIA  + Ndep*DIA + 
                     # all aspect interactions
                     slope*aspect  + MATmax*aspect  
                   + MAP*aspect  + MATmaxanom*aspect  
                   + MAPanom*aspect  + BAL*aspect  + damage*aspect  + 
                     PHYSIO*aspect + BA*aspect   + Ndep*aspect+
                     # all slope interactions
                     MATmax*slope  
                   + MAP*slope  + MATmaxanom*slope  
                   + MAPanom*slope  + BAL*slope  + damage*slope  +
                     PHYSIO*slope + BA*slope  + Ndep*slope +
                     # all MATmax interactions
                     
                     MAP*MATmax  + MATmaxanom*MATmax  
                   + MAPanom*MATmax  + BAL*MATmax  + damage*MATmax  + 
                     PHYSIO*MATmax + BA*MATmax +   Ndep*MATmax+ 
                     
                     # all MAP interatction
                     
                     MATmaxanom*MAP + MAPanom*MAP + BAL*MAP + damage*MAP + 
                     PHYSIO*MAP + BA*MAP +  Ndep*MAP + 
                     # all MATmaxanom interatction
                     
                     MAPanom*MATmaxanom + BAL*MATmaxanom + damage*MATmaxanom + 
                     PHYSIO*MATmaxanom + BA*MATmaxanom + Ndep*MATmaxanom +
                     # all MAPanom interatction
                     
                     BAL*MAPanom + damage*MAPanom + 
                     PHYSIO*MAPanom + BA*MAPanom  + Ndep*MAPanom   + 
                     # all BAL interatction
                     
                     damage*BAL + 
                     PHYSIO*BAL + BA*BAL +  Ndep*BAL , data = covariate.data, family = "binomial")
    
    glm.27 <-  glm(M ~ dia.diff + DIA + exp(DIA)+ aspect + slope + MATmax 
                   + MAP + MATmaxanom 
                   + MAPanom + BAL + damage  +  PHYSIO + BA  + Ndep + 
                     #SPCD.BA + non_SPCD.BA + prop.focal.ba +
                     # annual growth interactcovarions
                     DIA*dia.diff + aspect*dia.diff  + slope*dia.diff  + MATmax*dia.diff  
                   + MAP*dia.diff  + MATmaxanom*dia.diff  
                   + MAPanom*dia.diff  + BAL*dia.diff  + damage*dia.diff   
                   +  PHYSIO*dia.diff + BA*dia.diff +  Ndep*dia.diff+
                     # all diameter interactions
                     aspect*DIA  + slope*DIA  + MATmax*DIA  
                   + MAP*DIA  + MATmaxanom*DIA  
                   + MAPanom*DIA  + BAL*DIA  + damage*DIA  +  
                     PHYSIO*DIA + BA*DIA  + Ndep*DIA + 
                     # all aspect interactions
                     slope*aspect  + MATmax*aspect  
                   + MAP*aspect  + MATmaxanom*aspect  
                   + MAPanom*aspect  + BAL*aspect  + damage*aspect  + 
                     PHYSIO*aspect + BA*aspect   + Ndep*aspect+
                     # all slope interactions
                     MATmax*slope  
                   + MAP*slope  + MATmaxanom*slope  
                   + MAPanom*slope  + BAL*slope  + damage*slope  +
                     PHYSIO*slope + BA*slope  + Ndep*slope +
                     # all MATmax interactions
                     
                     MAP*MATmax  + MATmaxanom*MATmax  
                   + MAPanom*MATmax  + BAL*MATmax  + damage*MATmax  + 
                     PHYSIO*MATmax + BA*MATmax +   Ndep*MATmax+ 
                     
                     # all MAP interatction
                     
                     MATmaxanom*MAP + MAPanom*MAP + BAL*MAP + damage*MAP + 
                     PHYSIO*MAP + BA*MAP +  Ndep*MAP + 
                     # all MATmaxanom interatction
                     
                     MAPanom*MATmaxanom + BAL*MATmaxanom + damage*MATmaxanom + 
                     PHYSIO*MATmaxanom + BA*MATmaxanom + Ndep*MATmaxanom +
                     # all MAPanom interatction
                     
                     BAL*MAPanom + damage*MAPanom + 
                     PHYSIO*MAPanom + BA*MAPanom  + Ndep*MAPanom   + 
                     # all BAL interatction
                     
                     damage*BAL + 
                     PHYSIO*BAL + BA*BAL +  Ndep*BAL+
                     # all remainingdamage interatction
                     
                     PHYSIO*damage + BA*damage +  Ndep*damage, data = covariate.data, family = "binomial")
    
    
    glm.28 <-  glm(M ~ dia.diff + DIA + exp(DIA)+ aspect + slope + MATmax 
                   + MAP + MATmaxanom 
                   + MAPanom + BAL + damage  +  PHYSIO + BA  + Ndep + 
                     #SPCD.BA + non_SPCD.BA + prop.focal.ba +
                     # annual growth interactcovarions
                     DIA*dia.diff + aspect*dia.diff  + slope*dia.diff  + MATmax*dia.diff  
                   + MAP*dia.diff  + MATmaxanom*dia.diff  
                   + MAPanom*dia.diff  + BAL*dia.diff  + damage*dia.diff   
                   +  PHYSIO*dia.diff + BA*dia.diff +  Ndep*dia.diff+
                     # all diameter interactions
                     aspect*DIA  + slope*DIA  + MATmax*DIA  
                   + MAP*DIA  + MATmaxanom*DIA  
                   + MAPanom*DIA  + BAL*DIA  + damage*DIA  +  
                     PHYSIO*DIA + BA*DIA  + Ndep*DIA + 
                     # all aspect interactions
                     slope*aspect  + MATmax*aspect  
                   + MAP*aspect  + MATmaxanom*aspect  
                   + MAPanom*aspect  + BAL*aspect  + damage*aspect  + 
                     PHYSIO*aspect + BA*aspect   + Ndep*aspect+
                     # all slope interactions
                     MATmax*slope  
                   + MAP*slope  + MATmaxanom*slope  
                   + MAPanom*slope  + BAL*slope  + damage*slope  +
                     PHYSIO*slope + BA*slope  + Ndep*slope +
                     # all MATmax interactions
                     
                     MAP*MATmax  + MATmaxanom*MATmax  
                   + MAPanom*MATmax  + BAL*MATmax  + damage*MATmax  + 
                     PHYSIO*MATmax + BA*MATmax +   Ndep*MATmax+ 
                     
                     # all MAP interatction
                     
                     MATmaxanom*MAP + MAPanom*MAP + BAL*MAP + damage*MAP + 
                     PHYSIO*MAP + BA*MAP +  Ndep*MAP + 
                     # all MATmaxanom interatction
                     
                     MAPanom*MATmaxanom + BAL*MATmaxanom + damage*MATmaxanom + 
                     PHYSIO*MATmaxanom + BA*MATmaxanom + Ndep*MATmaxanom +
                     # all MAPanom interatction
                     
                     BAL*MAPanom + damage*MAPanom + 
                     PHYSIO*MAPanom + BA*MAPanom  + Ndep*MAPanom   + 
                     # all BAL interatction
                     
                     damage*BAL + 
                     PHYSIO*BAL + BA*BAL +  Ndep*BAL+
                     # all remainingdamage interatction
                     
                     PHYSIO*damage + BA*damage +  Ndep*damage +
                     BA*PHYSIO  + Ndep*PHYSIO, data = covariate.data, family = "binomial")
    
    
    
    
    glm.29 <-  glm(M ~ dia.diff + DIA + exp(DIA)+ aspect + slope + MATmax 
                   + MAP + MATmaxanom 
                   + MAPanom + BAL + damage  +  PHYSIO + BA  + Ndep + 
                     #SPCD.BA + non_SPCD.BA + prop.focal.ba +
                     # annual growth interactcovarions
                     DIA*dia.diff + aspect*dia.diff  + slope*dia.diff  + MATmax*dia.diff  
                   + MAP*dia.diff  + MATmaxanom*dia.diff  
                   + MAPanom*dia.diff  + BAL*dia.diff  + damage*dia.diff   
                   +  PHYSIO*dia.diff + BA*dia.diff +  Ndep*dia.diff+
                     # all diameter interactions
                     aspect*DIA  + slope*DIA  + MATmax*DIA  
                   + MAP*DIA  + MATmaxanom*DIA  
                   + MAPanom*DIA  + BAL*DIA  + damage*DIA  +  
                     PHYSIO*DIA + BA*DIA  + Ndep*DIA + 
                     # all aspect interactions
                     slope*aspect  + MATmax*aspect  
                   + MAP*aspect  + MATmaxanom*aspect  
                   + MAPanom*aspect  + BAL*aspect  + damage*aspect  + 
                     PHYSIO*aspect + BA*aspect   + Ndep*aspect+
                     # all slope interactions
                     MATmax*slope  
                   + MAP*slope  + MATmaxanom*slope  
                   + MAPanom*slope  + BAL*slope  + damage*slope  +
                     PHYSIO*slope + BA*slope  + Ndep*slope +
                     # all MATmax interactions
                     
                     MAP*MATmax  + MATmaxanom*MATmax  
                   + MAPanom*MATmax  + BAL*MATmax  + damage*MATmax  + 
                     PHYSIO*MATmax + BA*MATmax +   Ndep*MATmax+ 
                     
                     # all MAP interatction
                     
                     MATmaxanom*MAP + MAPanom*MAP + BAL*MAP + damage*MAP + 
                     PHYSIO*MAP + BA*MAP +  Ndep*MAP + 
                     # all MATmaxanom interatction
                     
                     MAPanom*MATmaxanom + BAL*MATmaxanom + damage*MATmaxanom + 
                     PHYSIO*MATmaxanom + BA*MATmaxanom + Ndep*MATmaxanom +
                     # all MAPanom interatction
                     
                     BAL*MAPanom + damage*MAPanom + 
                     PHYSIO*MAPanom + BA*MAPanom  + Ndep*MAPanom   + 
                     # all BAL interatction
                     
                     damage*BAL + 
                     PHYSIO*BAL + BA*BAL +  Ndep*BAL+
                     # all remainingdamage interatction
                     
                     PHYSIO*damage + BA*damage +  Ndep*damage +
                     BA*PHYSIO  + Ndep*PHYSIO+
                     Ndep*BA , data = covariate.data, family = "binomial")
    
    
    
    
    # mcfaddens rsquared
    list.mods <- list(glm.A, glm.B, glm.C, glm.D, glm.E, 
                      glm.F, glm.G, glm.H,  
                      glm.1, glm.2, glm.3, glm.4, glm.5, 
                      glm.6, glm.7, glm.8, glm.9, glm.10, 
                      glm.11, glm.12, glm.13, glm.14, glm.15, 
                      glm.16, glm.17, glm.18, glm.19, glm.20, 
                      #glm.20.a, glm.20.b, glm.20.c,
                      glm.21, glm.22, glm.23, glm.24, glm.25, 
                      glm.26, glm.27, glm.28, glm.29)
    
    # get convergence list
    Convergence.list <- lapply(list.mods, FUN = function(x){x$converged}) 
    convergence.df <- do.call(rbind, Convergence.list)
    
    # get AICs
    AICS.list <- lapply(list.mods, FUN = function(x){x$aic})
    aics.df <- do.call(rbind, AICS.list)
    
    # get mcfaddens rsquared to look at model fit
    McFadden.rsq <- lapply(list.mods, FUN =  function(x){pscl::pR2(x)["McFadden"]})
    McFadden.rsq.df <- do.call(rbind, McFadden.rsq)
    
    
    # variable importance
    Var.importance.list <- lapply(list.mods, FUN = function(x){caret::varImp(x)})
    
    # AUC of predicted test data:
    library(ROCR)
    newdata <- test.covariate.data
    
    # m <- list.mods[[1]]
    #m <- glm.26
    get_AUC <- function(m, newdata = covariate.data){
      p <- predict(m, newdata=newdata, type="response")
      pr <- prediction(p, newdata$M)
      prf <- performance(pr, measure = "tpr", x.measure = "fpr")
      plot(prf)
      auc <- performance(pr, measure = "auc")
      auc <- auc@y.values[[1]]
      auc
    }
    AUC.lists <- list()
    for (p in 1:length(list.mods)){
      AUC.lists[[p]]<- get_AUC(m = list.mods[[p]], test.covariate.data)
    }
    AUC.lists <- lapply (list.mods, FUN = function(x){get_AUC(x, test.covariate.data)})
    AUC.df <- do.call(rbind, AUC.lists)
    
    model.diag <- data.frame(SPCD =  SPCD.df[i,]$SPCD,
                             Species = nspp[1:17, ] %>% filter(SPCD %in% SPCD.df[i,]$SPCD) %>% dplyr::select(COMMON),
                             model = paste0("model ", 1:37), 
                             model.no = as.numeric(1:37),
                             remper.correction = remper.cor.vector[j],
                             AUC = AUC.df[,1],
                             McFadden.Rsq = McFadden.rsq.df[,1], 
                             AIC = aics.df[,1], 
                             converged = convergence.df[,1])
    plot(model.diag$model.no, model.diag$AUC)
    plot(model.diag$model.no, model.diag$McFadden.Rsq)
    
    saveRDS(model.diag, paste0("SPCD_glm_output/GLM_reduced_model_diag_SPCD_", SPCD.df[i,]$SPCD, "_remp_", remper.cor.vector[j], ".RDS") )
    saveRDS(Var.importance.list, paste0("SPCD_glm_output/GLM_reduced_variable_importance_list_SPCD_", SPCD.df[i,]$SPCD, "_remp_", remper.cor.vector[j], ".RDS") )
    
    
    #--------------------------------------------------------------------------------------
    # PLOT UP VARIABLE IMPORTANCE FOR THE BEST FIT MODEL 
    #--------------------------------------------------------------------------------------
    AIC.best <- model.diag %>% filter(converged == TRUE)%>% mutate(minAIC = min(AIC))%>%  filter(AIC == minAIC)
    AIC.best$model.no
    
    AIC.best.varimp <- Var.importance.list[[AIC.best$model.no]]
    AIC.best.varimp$VARS <- rownames(AIC.best.varimp)
    AIC.best.ordered <- AIC.best.varimp %>% arrange(Overall)
    AIC.best.ordered$VARS <- factor(AIC.best.ordered$VARS, levels = unique(AIC.best.ordered$VARS))
    
    ggplot(AIC.best.ordered, aes(x = VARS, y = Overall))+geom_bar(stat = "identity")+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))+
      ylab("Variable Importance")+xlab("")+
      ggtitle(paste0("Variable Importance, ", nspp[1:17, ] %>% filter(SPCD %in% SPCD.df[i,]$SPCD) %>% dplyr::select(COMMON) , " model ", AIC.best$model.no))
    ggsave(filename = paste0("SPCD_glm_output/GLM_reduced_AIC_best_VARIMP_SPCD_", SPCD.df[i,]$SPCD, "_remp_", remper.cor.vector[j], ".png"), height = 5, width = 12)
    
    
    AIC.best$model.no
    
    AUC.best <- model.diag %>% filter(converged == TRUE)%>% mutate(maxAUC = max(AUC))%>%  filter(AUC == maxAUC) 
    if(length(AUC.best$SPCD) > 1){
      AUC.best <- AUC.best %>% filter(model.no == min(model.no))
    }
    AUC.best$model.no
    
    AUC.best.varimp <- Var.importance.list[[AUC.best$model.no]]
    AUC.best.varimp$VARS <- rownames(AUC.best.varimp)
    AUC.best.ordered <- AUC.best.varimp %>% arrange(Overall)
    AUC.best.ordered$VARS <- factor(AUC.best.ordered$VARS, levels = unique(AUC.best.ordered$VARS))
    
    ggplot(AUC.best.ordered, aes(x = VARS, y = Overall))+geom_bar(stat = "identity")+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))+
      ylab("Variable Importance")+xlab("")+
      ggtitle(paste0("Variable Importance, ", nspp[1:17, ] %>% filter(SPCD %in% SPCD.df[i,]$SPCD) %>% dplyr::select(COMMON) , " model ", AUC.best$model.no))
    ggsave(filename = paste0("SPCD_glm_output/GLM_reduced_AUC_best_VARIMP_SPCD_", SPCD.df[i,]$SPCD, "_remp_", remper.cor.vector[j], ".png"), height = 5, width = 12)
    
    
    
    Rsq.best <- model.diag %>% filter(converged == TRUE)%>% mutate(maxRsq = max(McFadden.Rsq))%>%  filter(McFadden.Rsq == maxRsq)
    Rsq.best$model.no
    
    Rsq.best.varimp <- Var.importance.list[[Rsq.best$model.no]]
    Rsq.best.varimp$VARS <- rownames(Rsq.best.varimp)
    Rsq.best.ordered <- Rsq.best.varimp %>% arrange(Overall)
    Rsq.best.ordered$VARS <- factor(Rsq.best.ordered$VARS, levels = unique(Rsq.best.ordered$VARS))
    
    ggplot(Rsq.best.ordered, aes(x = VARS, y = Overall))+geom_bar(stat = "identity")+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))+
      ylab("Variable Importance")+xlab("")+
      ggtitle(paste0("Variable Importance, ", nspp[1:17, ] %>% filter(SPCD %in% SPCD.df[i,]$SPCD) %>% dplyr::select(COMMON) , " model ", Rsq.best$model.no))
    ggsave(filename = paste0("SPCD_glm_output/GLM_reduced_Rsq_best_VARIMP_SPCD_", SPCD.df[i,]$SPCD, "_remp_", remper.cor.vector[j], ".png"), height = 5, width = 12)
    
    ########################################################################
    # VIF for the species covariates
    #
    vif_values <-  vif(glm.17)
    VIF.cov <- data.frame(VIF = vif_values, 
                          covariates = names(vif_values))
    
    ggplot(VIF.cov, aes(x = covariates, y = VIF))+geom_bar(stat = "identity")+
      theme(axis.text.x = element_text(angle = 90))
    
    
    # Creating a correlation matrix
    corr <- round(cor(covariate.data[,3:ncol(covariate.data)]), 1)
    library(ggcorrplot)
    # generate correlation plots here:
    ggcorrplot(corr, hc.order = TRUE, type = "lower",
               lab = TRUE)
    ggsave(filename = paste0("SPCD_glm_output/GLM_reduced_Correlation_matrix", SPCD.df[i,]$SPCD, "_remp_", remper.cor.vector[j], ".png"), height = 12, width = 12)
    saveRDS(corr, paste0("SPCD_glm_output/GLM_reduced_Correlation_matrix_SPCD_",SPCD.id,"_predictors.rds") )
  }
}

# make a big table with the model results:
model.diags <- list()
for(i in 1:length(SPCD.df[,]$SPCD)){
  SPCD.id <- SPCD.df[i,]$SPCD
  # if(SPCD.id == 621){}else{
  model.diags[[i]] <- readRDS(paste0("SPCD_glm_output/GLM_reduced_model_diag_SPCD_", SPCD.df[i,]$SPCD, "_remp_", remper.cor.vector[j], ".RDS") )
  #}
}
model.diag <- do.call(rbind, model.diags)
model.diag %>% group_by(SPCD, COMMON)|> gt()

ggplot(model.diag %>% filter(converged == TRUE), aes(x = model.no, y = AUC, fill = COMMON, group = model.no))+geom_bar(stat= "identity",position = position_dodge2())#+position_dodge()
ggsave("SPCD_glm_output/GLM_reduced_all_species_AUC.png", height = 5, width = 8)

ggplot(model.diag %>% filter(converged == TRUE), aes(x = model.no, y = McFadden.Rsq))+geom_bar(stat= "identity") + facet_wrap(~COMMON, scales = "free_y")
ggsave("SPCD_glm_output/GLM_reduced_all_species_McFaddenRsq.png", height = 5, width = 8)

ggplot(model.diag %>% filter(converged == TRUE), aes(x = model.no, y = AUC))+geom_bar(stat= "identity") + facet_wrap(~COMMON, scales = "free_y")
ggsave("SPCD_glm_output/GLM_reduced_all_species_AUC.png", height = 5, width = 8)

ggplot(model.diag %>% filter(converged == TRUE), aes(x = model.no, y = AIC))+geom_bar(stat= "identity") + facet_wrap(~COMMON, scales = "free_y")
ggsave("SPCD_glm_output/GLM_reduced_all_species_AIC.png", height = 5, width = 8)

ggplot(model.diag %>% filter(converged ==TRUE), aes(AUC, AIC,  label = as.character(model.no)))+geom_text()+facet_wrap(~SPCD, scales = "free_y")
ggsave("SPCD_glm_output/GLM_reduced_all_species_AIC_AUC.png", height = 5, width = 8)

ggplot(model.diag %>% filter(converged ==TRUE), aes(AUC, McFadden.Rsq,  label = as.character(model.no)))+geom_text(size = 2)+facet_wrap(~SPCD, scales = "free_y")
ggsave("SPCD_glm_output/GLM_reduced_all_species_AUC_Rsq.png", height = 5, width = 8)

# make a table explaining the models:


glm.model.table <- data.frame(model.no = 1:35, 
                              description = c("MAP", 
                                              "MATmaxanom", 
                                              #"MATminanom", 
                                              "MAPanom", 
                                              "BAL", 
                                              "percent damage", 
                                              #"site index", 
                                              "Physiographic class", 
                                              "BA", 
                                              #"Relative Density",
                                              #"elevation", 
                                              "N depostion (wet + dry)", 
                                              "diameter difference", 
                                              
                                              ## sequentially adding in each variable
                                              "diameter difference + Diameter", 
                                              "exp(Diameter)", 
                                              "aspect", 
                                              "slope", 
                                              "MATmax", 
                                              #"MATmin", 
                                              "MAP", 
                                              "MATmax anomaly", 
                                              #"MATmin anomaly", 
                                              "MAP anomaly", 
                                              "BAL", 
                                              "percent damage", 
                                              #"site index", 
                                              "physiographic class", 
                                              "BA", 
                                              #"Relative Density", 
                                              #"Elevation", 
                                              "N deposition", 
                                              ## adding in interactions
                                              "growth interactions", 
                                              "diameter interactions", 
                                              "aspect interactions", 
                                              "slope interactions", 
                                              "MATmax interactions", 
                                              #"MATmin interactions", 
                                              "MAP interactions", 
                                              "MATmax anomaly interactions", 
                                              #"MATmin anomaly interactions", 
                                              "MAP anomaly interactions", 
                                              "BAL interactions", 
                                              "percent damage interactions", 
                                              #"site index interactions", 
                                              "Physiographic interactions", 
                                              #"Relative Density interactions", 
                                              #"elevation interactions",
                                              "Basal Area interactions", 
                                              "N dep interactions"
                              ), 
                              model.type = c(rep("single variable", 9), 
                                             rep("adding on to growth effect",13), 
                                             rep("adding interaction terms", 13
                                             )))
model.diag <- left_join(model.diag, glm.model.table)

ggplot(model.diag %>% filter(converged == TRUE), aes(x = model.no, y = AUC, fill = COMMON, group = model.no))+geom_bar(stat= "identity",position = position_dodge2())#+position_dodge()
ggsave("SPCD_glm_output/GLM_reduced_all_species_AUC.png", height = 5, width = 8)

ggplot(model.diag %>% filter(converged == TRUE), aes(x = model.no, y = McFadden.Rsq, fill = model.type))+geom_bar(stat= "identity") + facet_wrap(~COMMON, scales = "free_y")
ggsave("SPCD_glm_output/GLM_reduced_all_species_McFaddenRsq.png", height = 5, width = 8)

ggplot(model.diag %>% filter(converged == TRUE), aes(x = model.no, y = AUC, fill = model.type))+geom_bar(stat= "identity") + facet_wrap(~COMMON, scales = "free_y")
ggsave("SPCD_glm_output/GLM_reduced_all_species_AUC.png", height = 5, width = 8)

ggplot(model.diag %>% filter(converged == TRUE), aes(x = model.no, y = AIC, fill = model.type))+geom_bar(stat= "identity") + facet_wrap(~COMMON, scales = "free_y")
ggsave("SPCD_glm_output/GLM_reduced_all_species_AIC.png", height = 5, width = 8)

ggplot(model.diag %>% filter(converged ==TRUE), aes(AUC, AIC,  label = as.character(model.no)))+geom_text()+facet_wrap(~SPCD, scales = "free_y")
ggsave("SPCD_glm_output/GLM_reduced_all_species_AIC_AUC.png", height = 5, width = 8)

ggplot(model.diag %>% filter(converged ==TRUE), aes(AUC, McFadden.Rsq,  label = as.character(model.no)))+geom_text(size = 2)+facet_wrap(~SPCD, scales = "free_y")
ggsave("SPCD_glm_output/GLM_reduced_all_species_AUC_Rsq.png", height = 5, width = 8)

write.csv(glm.model.table, "GLM_reduced_table.csv", quote = TRUE)


#####################################################################################################################################################
# Look at Random forest approach and see if we get similar answers for the covariates that are important
#####################################################################################################################################################
library(randomForest)
covariate.data.RF <- covariate.all.df %>% select(-annual.growth, -SPCD)
covariate.data.RF$PHYSIO <- as.factor(covariate.data.RF$PHYSIO)
covariate.data.RF$M <- as.factor(covariate.data.RF$M)
all.cov.model <- randomForest(
  formula = M ~ .,
  data = covariate.data.RF
)

plot(all.cov.model)
varImpPlot(all.cov.model)




covariate.data.RF.red <- covariate.all.df %>% 
  select(-annual.growth, -SPCD, -MATmin, -MATminanom, -RD, -SPCD.BA, -non_SPCD.BA)
covariate.data.RF.red$PHYSIO <- as.factor(covariate.data.RF.red$PHYSIO)
covariate.data.RF.red$M <- as.factor(covariate.data.RF.red$M)
red.cov.model <- randomForest(
  formula = M ~ .,
  data = covariate.data.RF.red
)

plot(red.cov.model)
varImpPlot(red.cov.model)





# predict held out observations --full model
test.covariate.data.RF <- test.covariate.all.df %>% select(-annual.growth, -SPCD)
test.covariate.data.RF$PHYSIO <- as.factor(test.covariate.data.RF$PHYSIO)
test.covariate.data.RF$M <- as.factor(test.covariate.data.RF$M)

oos.preds <- predict(all.cov.model, newdata=test.covariate.data.RF)
roc.test <- pROC::roc(test.covariate.data.RF$M, as.numeric(as.vector(oos.preds)))
oos.full.AUC <- pROC::auc(roc.test)


# predict held out observations --reduced model
test.covariate.data.RF.red <- test.covariate.all.df %>% 
  select(-annual.growth, -SPCD, -MATmin, -MATminanom, -RD, -SPCD.BA, -non_SPCD.BA, 
         -prop.focal.ba, -elev, -si, -PHYSIO)
test.covariate.data.RF.red$PHYSIO <- as.factor(test.covariate.data.RF.red$PHYSIO)
test.covariate.data.RF.red$M <- as.factor(test.covariate.data.RF.red$M)

oos.preds.red <- predict(red.cov.model, newdata=test.covariate.data.RF.red)
roc.test.red <- pROC::roc(test.covariate.data.RF.red$M, as.numeric(as.vector(oos.preds.red)))
oos.red.AUC <- pROC::auc(roc.test.red)

# reducing the covariates improves prediction

# now do RF models for each species
covariate.import <- covariate.import.red <- covariate.import.red.all <- AUC.list<- list()
for(i in 1:17){
  covariate.data.RF <- covariate.all.df %>% 
    filter (SPCD %in% nspp[i,]$SPCD)%>% 
    select(-annual.growth, -SPCD)
  
  
  covariate.data.RF$PHYSIO <- as.factor(covariate.data.RF$PHYSIO)
  levels(covariate.data.RF$PHYSIO) <- as.character(c(1:7))#levels(covariate.data.RF$PHYSIO) # some instances the training set has less levels than the testing
  
  covariate.data.RF$M <- as.factor(covariate.data.RF$M)
  all.cov.model <- randomForest(
    formula = M ~ .,
    data = covariate.data.RF
   
  )
  
  plot(all.cov.model)
  varImpPlot(all.cov.model)
  all.cov.model$importance
  
  
  # predict held out observations --full model
  test.covariate.data.RF <- test.covariate.all.df %>%  
    filter (SPCD %in% nspp[i,]$SPCD)%>% 
    select(-annual.growth, -SPCD)#%>% filter(PHYSIO < 7 & PHYSIO > 3)
  test.covariate.data.RF$PHYSIO <- as.factor(test.covariate.data.RF$PHYSIO)
  test.covariate.data.RF$M <- as.factor(test.covariate.data.RF$M)
  #PHYS.keeps <- levels(test.covariate.data.RF$PHYSIO)[levels(test.covariate.data.RF$PHYSIO) %in% levels(covariate.data.RF$PHYSIO)]
  
  levels(test.covariate.data.RF$PHYSIO) <- levels(covariate.data.RF$PHYSIO) # some instances the training set has less levels than the testing
  levels(test.covariate.data.RF$PHYSIO) <- as.character(c(1:7))#levels(covariate.data.RF$PHYSIO) # some instances the training set has less levels than the testing
  
  
  
  oos.preds <- predict(all.cov.model, newdata=test.covariate.data.RF)
  oos.full.AUC <- mltools::auc_roc( preds = as.numeric(as.vector(oos.preds)), 
                                actuals = as.numeric(as.character(test.covariate.data.RF$M)))
  #oos.full.AUC <- pROC::auc(roc.test)
  
  
  
  covariate.import[[i]] <- data.frame(SPCD = nspp[i,]$SPCD, 
             covariate = rownames(all.cov.model$importance), 
             MeanDecreaseGini = all.cov.model$importance[,1])
  
  
  # drop the correlated covariates
   covariate.data.RF.drop <- covariate.data.RF %>% 
     select( -MATmin, -MATminanom, -RD, -SPCD.BA, -non_SPCD.BA)
  
  
  dropped.cov.model <- randomForest(
    formula = M ~ .,
    data = covariate.data.RF.drop
    
  )
  
  plot(dropped.cov.model)
  varImpPlot(dropped.cov.model)
  
  
  # save the reduced covariate importances into a list:
  covariate.import.red[[i]] <- data.frame(SPCD = nspp[i,]$SPCD, 
                                      covariate = rownames(dropped.cov.model$importance), 
                                      MeanDecreaseGini = dropped.cov.model$importance[,1])
  
  # predict held out observations --dropped covariates model
  test.covariate.data.RF.red <- test.covariate.all.df %>%  
    filter (SPCD %in% nspp[i,]$SPCD)%>% 
    select(-annual.growth, -SPCD, -MATmin, -MATminanom, -RD, -SPCD.BA, -non_SPCD.BA)
  test.covariate.data.RF.red$PHYSIO <- as.factor(test.covariate.data.RF$PHYSIO)
  test.covariate.data.RF.red$M <- as.factor(test.covariate.data.RF$M)
  levels(test.covariate.data.RF.red$PHYSIO) <- levels(covariate.data.RF.drop$PHYSIO) # some instances the training set has less levels than the testing
  
  oos.preds.red <- predict(dropped.cov.model, newdata=test.covariate.data.RF.red)
  oos.red.AUC <- mltools::auc_roc( preds = as.numeric(as.vector(oos.preds.red)), 
                                    actuals = as.numeric(as.character(test.covariate.data.RF.red$M)))
  

  # drop the additional correlated covariates:
  problem.variables.corr.species <- c("MATmin", "MATminanom", "RD", "SPCD.BA", # the same as above
                                      "non_SPCD.BA", "prop.focal.ba", "si", "elev") # but also more competition variables
  
  covariate.data.RF.drop.all <- covariate.data.RF %>% 
    select( -MATmin, -MATminanom, -RD, -SPCD.BA, -non_SPCD.BA, 
            -prop.focal.ba, -si, -elev)
  
  
  all.dropped.cov.model <- randomForest(
    formula = M ~ .,
    data = covariate.data.RF.drop.all
    
  )
  
  plot(all.dropped.cov.model)
  varImpPlot(all.dropped.cov.model)
  
  
  # save the reduced covariate importances into a list:
  covariate.import.red.all[[i]] <- data.frame(SPCD = nspp[i,]$SPCD, 
                                          covariate = rownames(all.dropped.cov.model$importance), 
                                          MeanDecreaseGini = all.dropped.cov.model$importance[,1])
  
  # predict held out observations --dropped covariates model
  test.covariate.data.RF.red.all <- test.covariate.all.df %>%  
    filter (SPCD %in% nspp[i,]$SPCD)%>% 
    select(-annual.growth, -SPCD, -MATmin, -MATminanom, -RD, -SPCD.BA, -non_SPCD.BA,
           -prop.focal.ba, -si, -elev)
  
  test.covariate.data.RF.red.all$PHYSIO <- as.factor(test.covariate.data.RF.red.all$PHYSIO)
  test.covariate.data.RF.red.all$M <- as.factor(test.covariate.data.RF.red.all$M)
  levels(test.covariate.data.RF.red.all$PHYSIO) <- levels(covariate.data.RF.drop$PHYSIO) # some instances the training set has less levels than the testing
  
  oos.preds.red.all <- predict(all.dropped.cov.model, newdata=test.covariate.data.RF.red.all)
  oos.red.AUC.all <- mltools::auc_roc( preds = as.numeric(as.vector(oos.preds.red.all)), 
                                   actuals = as.numeric(as.character(test.covariate.data.RF.red.all$M)))
  
  
  
  
  
 # further reduce by removing the extra species correlated variables--
  AUC.list[[i]] <- data.frame(SPCD = nspp[i,]$SPCD, 
                              AUC.full = oos.full.AUC, 
                              AUC.reduced = oos.red.AUC, 
                              AUC.reduced.all = oos.red.AUC.all)
  
}
# gather full importances and plot
importances.spcd <- do.call(rbind, covariate.import)
ggplot(importances.spcd, aes(x = MeanDecreaseGini, y = covariate, color = SPCD))+geom_point()

# gather reduced importances and plot
reduced.mods.importances.spcd <- do.call(rbind, covariate.import.red)
ggplot(reduced.mods.importances.spcd, aes(x = MeanDecreaseGini, y = covariate, color = SPCD))+geom_point()

# gather reduced importances and plot
AUC.df <- do.call(rbind, AUC.list)
AUC.df %>% mutate(improved = ifelse(AUC.reduced > AUC.full, "Yes", "No"), 
                  improved.full = ifelse(AUC.reduced.all > AUC.full, "Yes", "No"))

AUC.df %>% mutate(gain = AUC.reduced - AUC.full, 
                  gain.full = AUC.reduced.all - AUC.full) 
# for most species, removing variables doesnt necessarily improve AUC scores, but it also doesnt drastically hurt them
