library(rstan)
library(MASS)
library(here)
library(tidyverse)
library(gt)
library(FIESTA)
library(dplyr)
library(mltools)
library(splines)
library(gt)

options(mc.cores = parallel::detectCores())
cleaned.data <- readRDS( "data/cleaned.data.mortality.TRplots.RDS")
cleaned.data <- cleaned.data %>% filter(!is.na(ba) & !is.na(slope) & ! is.na(physio) & !is.na(aspect))%>% dplyr::select(state, county, pltnum, cndtn, point, tree, PLOT.ID, cycle, spp, dbhcur, dbhold, damage, Species, SPCD,
                                                                                                                        remper, LAT_FIADB, LONG_FIADB, elev, DIA_DIFF, annual.growth, M, relative.growth, si, physio:RD) %>% distinct()
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
nspp$Species <- paste(FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD),]$GENUS, FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD),]$SPECIES)

#View(nspp)

nspp[1:17,]$COMMON


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


View(cleaned.data %>% filter(SPCD %in% unique(nspp[1:17,]$SPCD))%>% group_by( SPGRPNAME, SPCD) %>% summarise(n()))
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
cleaned.data.full$status

mortality.summary <- cleaned.data.full %>% filter(SPCD %in% unique(nspp[1:17,]$SPCD)) %>%
  group_by(SPCD, M) %>% summarise(ntrees = n()) %>%
  spread(M, ntrees) %>% rename("live" = `0`, 
                               "dead" = `1`)
mortality.summary$`Common name` <- FIESTA::ref_species[match(mortality.summary$SPCD, FIESTA::ref_species$SPCD),]$COMMON
mortality.summary %>% dplyr::select(`Common name`, SPCD, live, dead) %>% mutate(total = live + dead) %>% ungroup() |> gt()
#-----------------------------------------------------------------------------------------
# Make the species level datasets for the top 15 species to run the model
#-----------------------------------------------------------------------------------------
length(unique(cleaned.data$SPCD))

nspp[1:17,]$COMMON
# save these as .RDA files so we can just load, run the model, and 
SPCD.id <- 316#unique(cleaned.data$SPCD)[25]
set.seed(22)
summary(cleaned.data.full)
cleaned.data.full$status

# Adapted this function to create several different matrices of data for the general STAN model
# each dataset has a different combination of variables and interaction effects
# options for model effects
# Option A:
# 1. Diameter
# 2. Damage
# 3. Annual growth
# 4. Diameter + Annual growth + Damage
# 5. All Fixed effects
# 6. All Fixed effects and all growth interactions
# 7. All Fixed effects and all growth + Diameter interactions
# 8. All Fixed effects and all growth + Diameter + damage interactions
# 9. All Fixed effects and all interactions

# Option B:
# 1. Annual growth 
# 2. Diameter + Annual growth
# 3. Diameter + Annual growth + competition variables (RD, BAL, damage)
# 4. Diameter + Annual growth + competition variables (RD, BAL, damage) + Climate variables
# 5. Diameter + Annual growth + competition variables (RD, BAL, damage) + Climate variables + site/soil effects + ndep
# 6. All Fixed effects and all growth + Diameter interactions
# 7. model 5 + competition interactions
# 8. model 6 + climate interactions
# 9. All Fixed effects and all interactions

stan.model.table <- data.frame(model = 1:9, 
                               `Covariates` = c("Annual growth",  
                                                "diameter + annual growth",
                                                "diameter + annual growth + competition variables", 
                                                "diameter + annual growth + competition variables  + climate variables", 
                                                "diameter + annual growth + competition variables + climate variables + site/soil effects + ndep",
                                                "All Fixed effects and all growth + diameter interactions",
                                                "model 5 + competition interactions",
                                                "model 6 + climate interactions",
                                                "All Fixed effects and all interactions"))
stan.model.table |> gt()

# make another table with covariates
Covariate.table <- read.csv("C:/Users/KellyHeilman/Box/01. kelly.heilman Workspace/mortality/Eastern-Mortality/mortality_manuscript/Covariate_descriptions_table.csv")
Covariate.table |> gt()
# script that generates all the testing and training datasets
# source("R/SPCD_stan_data.R")
# # write the data for all 26 different species groups:
# for(i in 1:length(unique(nspp[1:17,]$SPCD))){
#   cat(i)
#   SPCD.stan.data(SPCD.id = nspp[i,]$SPCD, remper.correction = 0.5, cleaned.data.full = cleaned.data.full)
#   # SPCD.stan.data(SPCD.id = nspp[i,]$SPCD, remper.correction = 0.9, cleaned.data.full = cleaned.data.full)
#   # SPCD.stan.data(SPCD.id = nspp[i,]$SPCD, remper.correction = 0.1, cleaned.data.full = cleaned.data.full)
#   # SPCD.stan.data(SPCD.id = nspp[i,]$SPCD, remper.correction = 0.3, cleaned.data.full = cleaned.data.full)
#   # SPCD.stan.data(SPCD.id = nspp[i,]$SPCD, remper.correction = 0.7, cleaned.data.full = cleaned.data.full)
# }
# 

#----------------------------------------------------------------------------------
# running stan models with most important variables
#----------------------------------------------------------------------------------


# for each species group, fit a model, plot the outputs, and save the results
# we source a function from another script
source("R/speciesBasisModels/SPCD_run_stan_basis.R")

# this runs a stan model and saves the outputs
# SPGRPCD 2 throws uncerialize socklist error--too big of a diataset to parallelisze?
cleaned.data.full %>% group_by(SPCD) %>% summarise(n())
SPCD.df <- data.frame(SPCD = nspp[1:17, ]$SPCD, 
                      spcd.id = 1:17)
remper.cor.vector <- c(0.5)
#model.number <- 6
model.list <- 3:5

# load("dbh_smooth_gams/BDBH_GAM_models_SPCD_318.Rdata")
# plot.gam(m.dbh.1)
# m.dbh.2$coefficients
# m.dbh.3$coefficients
# m.dbh.4$coefficients
# m.dbh.5$coefficients
# m.dbh.5$smooth[[1]]$knots
# 
# b$smooth[[1]]
for(i in 1:17){# run for each of the 17 species
  common.name <- nspp[1:17, ] %>% filter(SPCD %in% SPCD.df[i,]$SPCD) %>% dplyr::select(COMMON)
  
  for(m in 1:length(model.list)){  # run each of the 9 models
    model.number <- model.list[m]
    for (j in 1:length(remper.cor.vector)){ # for the growth only model explore the consequences of other assumptions about remeasurement period
      cat(paste("running stan mortality model ",model.number, " for SPCD", SPCD.df[i,]$SPCD, common.name$COMMON, " remper correction", remper.cor.vector[j]))
      
      fit.1 <- SPCD_run_stan_basis(SPCD.id = SPCD.df[i,]$SPCD,
                             model.no = model.number,
                             niter = 100,
                             nchains = 1,
                             remper.correction = remper.cor.vector[j],
                             model.file = 'modelcode/mort_model_basis.stan',
                             output.folder = "C:/Users/KellyHeilman/Box/01. kelly.heilman Workspace/mortality/Eastern-Mortality/mortality_models/")
      
      SPCD.id <-  SPCD.df[i,]$SPCD
      #saveRDS(fit.1, paste0("SPCD_stanoutput_full/samples/model_",model.number,"_SPCD_",SPCD.id, "_remper_correction_", remper.cor.vector[j], ".RDS"))
      # save_diagnostics (stanfitobj = fit.1, nchains = 2, model.no = model.number, remper.correction = remper.cor.vector[j])
      
      model.name <- paste0("basis_model_",model.number,"_SPCD_", SPCD.id, "_remper_correction_", remper.cor.vector[j])
      remp.cor <- remper.cor.vector[j]
      remper.correction <- remper.cor.vector[j]
      output.folder = "C:/Users/KellyHeilman/Box/01. kelly.heilman Workspace/mortality/Eastern-Mortality/mortality_models/"
      source("R/speciesBasisModels/SPCD_plot_stan_basis.R")
      rm(fit.1)
    }
  }
}

# read in the outputs from the diameter basis function models for model 6
model.number <- 6
SPCD.id
remper.cor.vector <- 0.5
etas.list <- list()
for(i in 1:17){
  SPCD.id <- nspp[i,]$SPCD
  common.name <- nspp[1:17, ] %>% filter(SPCD %in% SPCD.id) %>% dplyr::select(COMMON)
  
   
  fit.1 <- readRDS( paste0("SPCD_stanoutput_full/samples/basis_model_",model.number,"_SPCD_",SPCD.id, "_remper_correction_0.5.RDS"))
  fit.1 <- as_draws_matrix(fit.1)
  
  etas <- subset_draws(fit.1, variable = "eta", chain = 1:2)
  rm(fit.1)
  etas.quant <- summarise_draws(etas, median, ~quantile(.x, probs = c(0.025, 0.975)))%>% rename(`ci.lo` = `2.5%`, `ci.hi` = `97.5%`)
  
  etas.quant$SPCD <- SPCD.id
  etas.quant$Species <- common.name
  etas.quant$model <- 6
  etas.quant$`Size Spline` <- paste0("Size Beta ", 1:9)
  etas.list[[i]] <- etas.quant
  
}

etas.df <- do.call(rbind, etas.list)
etas.df$`significance` <- ifelse(etas.df$ci.lo < 0 & etas.df$ci.hi < 0, "significant", 
                                            ifelse(etas.df$ci.lo > 0 & etas.df$ci.hi > 0, "significant", "not overlapping zero"))

ggplot(data = etas.df, aes(x = `Size Spline`, y = median, color = significance))+geom_point()+
  geom_errorbar(data = etas.df, aes(x = `Size Spline`, ymin = ci.lo, ymax = ci.hi, color = significance), width = 0.5)+
  theme_bw()+facet_wrap(~Species$COMMON)+ylab("Size Spline Estimates")+xlab("Size Spline Betas")+
  theme(axis.text.x = element_text(vjust = 0.5, angle = 90), 
        legend.position = "none")+scale_color_manual(values = c("not overlapping zero"="grey46", "significant"="black"))

ggsave(filename = "model_summary_full/DBH_basis_eta_values_model6_all_species.png", 
       height = 6, width = 6.5)


# plot the effect * the 
get.basis.preds <- function(SPP.num){
    SPCD.id <- nspp[SPP.num,]$SPCD
    common.name <- nspp[1:17, ] %>% filter(SPCD %in% SPCD.id) %>% dplyr::select(COMMON)
    
    load(paste0("SPCD_standata_basis/SPCD_",SPCD.id,"remper_correction_0.5model_6.Rdata"))
    
    # get diameter means and sd values
    all.data.diam.scaler <- rbind(train.data, test.data)%>% summarise(mean.dbh = mean(dbhcur), 
                                                          sd.dbh = sd(dbhcur))
    
    # calculate the diameter range and basis functions
    dia.range <- seq(quantile(dia, c(.05,.95))[1], quantile(dia, c(.05,.95))[2], by = 0.5)
    knots = quantile(dia, c(.05, .2, .5, .8, .95))
    basis_range_dia <- bs(dia.range, knots=knots, degree = 3, intercept = TRUE)
    dia.range.real <- dia.range*all.data.diam.scaler$sd.dbh + all.data.diam.scaler$mean.dbh
    
    fit.1 <- readRDS( paste0("SPCD_stanoutput_full/samples/basis_model_6_SPCD_",SPCD.id, "_remper_correction_0.5.RDS"))
    fit.1 <- as_draws_matrix(fit.1)
    
    etas <- subset_draws(fit.1, variable = "eta", chain = 1:2)
    rm(fit.1, train.data, test.data)
    pred <- list()
    for(i in 1:nrow(basis_range_dia)){
      pred[[i]] <- etas%*%as.matrix(basis_range_dia)[i,]
    }
    pred.df <- do.call(cbind, pred)
    pred.m <- reshape2::melt(pred.df) %>% group_by(Var2) %>% summarise(median = quantile(value, 0.5), 
                                                                       ci.lo = quantile(value, 0.025), 
                                                                       ci.hi = quantile(value, 0.975), 
                                                                       ci.lo.75 = quantile(value, 0.125), 
                                                                       ci.hi.75 = quantile(value, 0.875))
    pred.m$`DIA scaled` <- dia.range
    pred.m$`Diameter` <- dia.range.real
    pred.m$Species <- common.name$COMMON
    pred.m$SPCD <- SPCD.id
    pred.m
}

basis.preds <- lapply(1:17, get.basis.preds)
basis.preds.all <- do.call(rbind, basis.preds)

ggplot(data = basis.preds.all, aes(x = `DIA scaled`, y = median))+geom_line()+
  geom_ribbon(data = basis.preds.all, aes(x = `DIA scaled`, ymin = ci.lo, ymax = ci.hi), alpha = 0.5, fill = "#a6bddb")+
  geom_ribbon(data = basis.preds.all, aes(x = `DIA scaled`, ymin = ci.lo.75, ymax = ci.hi.75), alpha = 0.5, fill = "#1c9099")+
  facet_wrap(~Species, scales = "free_y")+theme_bw()+ylab("Log-Predicted Effect of Diameter on Survival")+xlab("Scaled Diameter")
ggsave("SPCD_stanoutput_full/Basis_predicted_effects_relative_diameter_all_SPP_model_6.png",
  height = 6, width = 8)

ggplot(data = basis.preds.all, aes(x = `Diameter`, y = median))+geom_line()+
  geom_ribbon(data = basis.preds.all, aes(x = `Diameter`, ymin = ci.lo, ymax = ci.hi), alpha = 0.5, fill = "#a6bddb")+
  geom_ribbon(data = basis.preds.all, aes(x = `Diameter`, ymin = ci.lo.75, ymax = ci.hi.75), alpha = 0.5, fill = "#1c9099")+
  facet_wrap(~Species, scales = "free_y")+theme_bw()+ylab("Log-Predicted Effect of Diameter on Survival")+xlab("Diameter (in)")
ggsave("SPCD_stanoutput_full/Basis_predicted_effects_abs_diameter_in_all_SPP_model_6.png",
       height = 6, width = 8)

ggplot(data = basis.preds.all, aes(x = `Diameter`*2.54, y = median))+geom_line()+
  geom_ribbon(data = basis.preds.all, aes(x = `Diameter`*2.54, ymin = ci.lo, ymax = ci.hi), alpha = 0.5, fill = "#a6bddb")+
  geom_ribbon(data = basis.preds.all, aes(x = `Diameter`*2.54, ymin = ci.lo.75, ymax = ci.hi.75), alpha = 0.5, fill = "#1c9099")+
  facet_wrap(~Species, scales = "free_y")+theme_bw()+ylab("Log-Predicted Effect of Diameter on Survival")+xlab("Diameter (cm)")
ggsave("SPCD_stanoutput_full/Basis_predicted_effects_abs_diameter_cm_all_SPP_model_6.png",
       height = 6, width = 8)
