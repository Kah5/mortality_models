library(rstan)
library(MASS)
#library(here)
library(tidyverse)
#library(gt)
library(FIESTA)
library(dplyr)
library(mltools)
library(scales)
library(posterior)

options(mc.cores = parallel::detectCores())
# cleaned.data <- readRDS( "data-store/data/iplant/home/kellyheilman/mort_data/cleaned.data.mortality.TRplots.RDS")
# cleaned.data.no.ll <- cleaned.data %>% select(-LAT_FIADB, -LONG_FIADB)
# 
# cleaned.data <- cleaned.data %>% filter(!is.na(ba) & !is.na(slope) & ! is.na(physio) & !is.na(aspect))%>% 
#   dplyr::select(state, county, pltnum, cndtn, point, tree, PLOT.ID, cycle, spp, dbhcur, dbhold, damage, Species, SPCD,
#                                                remper, LAT_FIADB, LONG_FIADB, elev, DIA_DIFF, annual.growth, M, relative.growth, si, physio:RD) %>% distinct()
# 
# 
# nspp <- cleaned.data %>% group_by(SPCD) %>% summarise(n = n(), 
#                                                       pct = n/nrow(cleaned.data)) %>% arrange (desc(`pct`))
# 
# nspp$cumulative.pct <- cumsum(nspp$pct)
# 
# 
# 
# # link up to the species table:
# nspp$COMMON <- FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD),]$COMMON
# nspp$Species <- paste(FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD),]$GENUS, FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD),]$SPECIES)
# 
# #View(nspp)
# 
# nspp[1:17,]$COMMON
# 
# 
# nspp[1:17,] %>% mutate(pct = round(pct, 3), 
#                        cumulative.pct = round(cumulative.pct, 3)) %>% rename(`# of trees` = "n", 
#                        `% of trees` = "pct",
#                        `cumulative %` = "cumulative.pct", 
#                        `Common name` = "COMMON") %>%
#   dplyr::select(Species, `Common name`, SPCD, `# of trees`, `% of trees`, `cumulative %`)
# 
# 
# 
# cleaned.data$SPGRPCD <- FIESTA::ref_species[match(cleaned.data$SPCD, FIESTA::ref_species$SPCD),]$E_SPGRPCD
# 
# SPGRP.df <- FIESTA::ref_codes %>% filter(VARIABLE %in% "SPGRPCD") %>% filter(VALUE %in% unique(cleaned.data$SPGRPCD))
# cleaned.data$SPGRPNAME <- SPGRP.df[match(cleaned.data$SPGRPCD, SPGRP.df$VALUE),]$MEANING


#----------------------------------------------------------------------------------
# running stan models with most important variables
#----------------------------------------------------------------------------------


# for each species group, fit a model, plot the outputs, and save the results
# we source a function from another script
source("R/speciesModels/SPCD_run_stan.R")

# this runs a stan model and saves the outputs
#cleaned.data %>% group_by(SPCD) %>% summarise(n())
SPCD.df <- data.frame(SPCD = nspp[1:17, ]$SPCD, 
                      spcd.id = 1:17)
remper.cor.vector <- c(0.5)
#model.number <- 6
model.list <- 1:9
m <- 1


for(m in 3:9){ 
  
  model.number <- model.list[m]
  for(i in 9:1){# run for each of the 17 species
    common.name <- nspp[1:17, ] %>% filter(SPCD %in% SPCD.df[i,]$SPCD) %>% dplyr::select(COMMON)
    
    #for(m in 1:length(model.list)){  # run each of the 9 models
    
    
    for (j in 1:length(remper.cor.vector)){ # for the growth only model explore the consequences of other assumptions about remeasurement period
      cat(paste("running stan mortality model ", model.number, " for SPCD", SPCD.df[i,]$SPCD, common.name$COMMON, " remper correction", remper.cor.vector[j]))
      
      fit.1 <- SPCD_run_stan(SPCD.id = SPCD.df[i,]$SPCD,
                             model.no = model.number,
                             niter = 2000,
                             nchains = 3,
                             remper.correction = remper.cor.vector[j],
                             model.file = 'modelcode/mort_model_general.stan', 
                             output.folder = "SPCD_stanoutput_full_standardized_v3")
      SPCD.id <-  SPCD.df[i,]$SPCD
      
      # set up to make plots of the stan outputs 
      model.name <- paste0("mort_model_",model.number,"_SPCD_", SPCD.id, "_remper_correction_", remper.cor.vector[j])
      remp.cor <- remper.cor.vector[j]
      remper.correction <- remper.cor.vector[j]
      
      output.folder = "SPCD_stanoutput_full_standardized_v3/"
      
      source("R/speciesModels/SPCD_plot_stan.R")
      rm(fit.1)
    }
  }
}


# get the predicted AUC for each model 1-6:
for(i in 9:1){# run for each of the 17 species
  common.name <- nspp[1:17, ] %>% filter(SPCD %in% SPCD.df[i,]$SPCD) %>% dplyr::select(COMMON)
  
  for(m in 3:9){  # run each of the 9 models
    #for(m in 8:9){ 
    
    model.number <- model.list[m]
    for (j in 1:length(remper.cor.vector)){ # for the growth only model explore the consequences of other assumptions about remeasurement period
      cat(paste("running AUC stats for model ",model.number, " for SPCD", SPCD.df[i,]$SPCD, common.name$COMMON, " remper correction", remper.cor.vector[j]))
      
      SPCD.id <-  SPCD.df[i,]$SPCD
      
      # set up to make plots of the stan outputs 
      model.name <- paste0("mort_model_",model.number,"_SPCD_", SPCD.id, "_remper_correction_", remper.cor.vector[j])
      remp.cor <- remper.cor.vector[j]
      remper.correction <- remper.cor.vector[j]
      
      output.folder = "/home/rstudio/"
      
      source("R/speciesModels/SPCD_AUC_stan.R")
      rm(fit.1)
    }
  }
}

# now transfer all the files to the outputs directory:

# copy the data-store files
system(paste("cp -r",  "SPCD_stanoutput_full_standardized_v3/",
             "data-store/data/output/"))



# ##########################################################################################
# # for model 6 vary the remper length
# ##########################################################################################

# 
# # this runs a stan model and saves the outputs
# cleaned.data.full %>% group_by(SPCD) %>% summarise(n())
# SPCD.df <- data.frame(SPCD = nspp[1:17, ]$SPCD, 
#                       spcd.id = 1:17)
# remper.cor.vector <- c(0.1, 0.3, 0.7, 0.9)
# #model.number <- 6
# model.list <- 6
# 
# for(i in 1:17){# run for each of the 17 species
#   common.name <- nspp[1:17, ] %>% filter(SPCD %in% SPCD.df[i,]$SPCD) %>% dplyr::select(COMMON)
#   
#   for(m in 1:length(model.list)){  # run each of the 9 models
#     model.number <- model.list[m]
#     for (j in 1:length(remper.cor.vector)){ # for the gowoth only model explore the consequences of other assumptions about remeasurement period
#       cat(paste("running stan mortality model ",model.number, " for SPCD", SPCD.df[i,]$SPCD, common.name$COMMON, " remper correction", remper.cor.vector[j]))
#       
#       fit.1 <- SPCD_run_stan(SPCD.id = SPCD.df[i,]$SPCD,
#                              model.no = model.number,
#                              niter = 3000,
#                              nchains = 3,
#                              remper.correction = remper.cor.vector[j],
#                              model.file = 'modelcode/mort_model_general.stan', 
#                              output.folder = "C:/Users/KellyHeilman/Box/01. kelly.heilman Workspace/mortality/Eastern-Mortality/mortality_models/")
#       SPCD.id <-  SPCD.df[i,]$SPCD
#       #saveRDS(fit.1, paste0("SPCD_stanoutput_full/samples/model_",model.number,"_SPCD_",SPCD.id, "_remper_correction_", remper.cor.vector[j], ".RDS"))
#       # save_diagnostics (stanfitobj = fit.1, nchains = 2, model.no = model.number, remper.correction = remper.cor.vector[j])
#       
#       model.name <- paste0("mort_model_",model.number,"_SPCD_", SPCD.id, "_remper_correction_", remper.cor.vector[j])
#       remp.cor <- remper.cor.vector[j]
#       remper.correction <- remper.cor.vector[j]
#       
#       output.folder = "C:/Users/KellyHeilman/Box/01. kelly.heilman Workspace/mortality/Eastern-Mortality/mortality_models/"
#       
#       source("R/SPCD_plot_stan.R")
#       rm(fit.1)
#     }
#   }
# }
# 
# 
# # plotting effects of mortality model by species to compare remper growth effects
# # for SPCD 316, plot all estimated values:


