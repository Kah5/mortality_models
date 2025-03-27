library(rstan)
library(MASS)
library(here)
library(tidyverse)
library(gt)
library(FIESTA)
library(dplyr)
library(mltools)
library(scales)

options(mc.cores = parallel::detectCores())
cleaned.data <- readRDS( "data/cleaned.data.mortality.TRplots.RDS")
cleaned.data.no.ll <- cleaned.data %>% select(-LAT_FIADB, -LONG_FIADB)
saveRDS(cleaned.data.no.ll, "data/cleaned.data.mortality.TRplots_no_LL.RDS")

cleaned.data <- cleaned.data %>% filter(!is.na(ba) & !is.na(slope) & ! is.na(physio) & !is.na(aspect))%>% 
  dplyr::select(state, county, pltnum, cndtn, point, tree, PLOT.ID, cycle, spp, dbhcur, dbhold, damage, Species, SPCD,
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
plot.medians <- unique(
  cleaned.data %>% ungroup() %>%
    dplyr::select(
      PLOT.ID,
      si,
      ba,
      BA_total,
      SPCD_BA, 
      non_SPCD_BA, 
      
      slope,
      aspect,
      MAP,
      MATmin,
      MATmax,
      damage.total,
      elev,
      Ndep.remper.avg,
      physio,
      RD
    )
) %>%
  ungroup() %>% summarise(
    si.median = median(si, na.rm = TRUE),
    RD.median = median(RD, na.rm = TRUE),
    ba.median = median(ba, na.rm = TRUE),
    BA_tot.median = median(BA_total, na.rm =
                             TRUE),
    nonSPCD_BA_tot.median = median(non_SPCD_BA, na.rm = TRUE),
    SPCD_BA.median = median(SPCD_BA, na.rm =
                              TRUE),
    slope.median = median(slope, na.rm = TRUE),
    aspect.median = median(aspect, na.rm = TRUE),
    damage.median = median(damage.total, na.rm =
                             TRUE),
    elev.median = median(elev, na.rm = TRUE),
    Ndep.median = median(Ndep.remper.avg, na.rm =
                           TRUE),
    physio.median = median(physio, na.rm = TRUE),
    MAP.median = median(MAP, na.rm =
                          TRUE),
    MATmin.median = median(MATmin, na.rm = TRUE),
    MATmax.median = median(MATmax, na.rm =
                             TRUE),
    RD.sd = sd(RD, na.rm = TRUE),
    ba.sd = sd(ba, na.rm = TRUE),
    
    BA_tot.sd = sd(BA_total, na.rm =
                     TRUE),
    nonSPCD_BA_tot.sd = sd(non_SPCD_BA, na.rm = TRUE),
    SPCD_BA.sd = sd(SPCD_BA, na.rm =
                      TRUE),
    si.sd = sd(si, na.rm =
                 TRUE),
    slope.sd = sd(slope, na.rm = TRUE),
    aspect.sd = sd(aspect, na.rm = TRUE),
    damage.sd = sd(damage.total, na.rm =
                     TRUE),
    elev.sd = sd(elev, na.rm = TRUE),
    Ndep.sd = sd(Ndep.remper.avg, na.rm =
                   TRUE),
    physio.sd = sd(physio, na.rm = TRUE),
    MAP.sd = sd(MAP, na.rm = TRUE),
    MATmin.sd = sd(MATmin, na.rm =
                     TRUE),
    MATmax.sd = sd(MATmax, na.rm = TRUE)
  )

#View(cleaned.data %>% group_by(SPGRPCD, SPCD) %>% summarise(n()))

cleaned.data.full <- cleaned.data

mortality.summary <- cleaned.data.full %>% filter(SPCD %in% unique(nspp[1:17,]$SPCD)) %>%
  group_by(SPCD, M) %>% summarise(ntrees = n()) %>%
  spread(M, ntrees) %>% rename("live" = `0`, 
                               "dead" = `1`)
mortality.summary$`Common name` <- FIESTA::ref_species[match(mortality.summary$SPCD, FIESTA::ref_species$SPCD),]$COMMON
mortality.summary %>% dplyr::select(`Common name`, SPCD, live, dead) %>% mutate(total = live + dead) %>% ungroup() |> gt()
colnames(cleaned.data.full)

#-----------------------------------------------------------------------------------------
# Make the species level datasets for the top 15 species to run the model
#-----------------------------------------------------------------------------------------
length(unique(cleaned.data$SPCD))

nspp[1:17,]$COMMON
# save these as .RDA files so we can just load, run the model, and 
SPCD.id <- 316#unique(cleaned.data$SPCD)[25]
set.seed(22)
summary(cleaned.data.full)

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

stan.model.table <- data.frame(Model = 1:9, 
                               `Covariates` = c("annual growth",  
                                                "annual growth & diameter ",
                                                "annual growth, diameter  & competition variables", 
                                                "annual growth, diameter, competition variables, climate variables", 
                                                "annual growth, diameter, competition variables, climate variables, site conditions ",
                                                "all fixed effects and all growth + diameter interactions",
                                                "model 5 + competition interactions",
                                                "model 6 + climate interactions",
                                                "all fixed effects and all interactions"), 
                               `Fixed.parameters` = c(1, 2, 7, 13, 18, 51, 103, 148, 158))
stan.model.table |> gt() |> cols_label(Fixed.parameters = "# of fixed parameters") |> gtsave("C:/Users/KellyHeilman/Box/01. kelly.heilman Workspace/mortality/Eastern-Mortality/mortality_manuscript/figures/model_table_summary.png")

# make another table with covariates
Covariate.table <- read.csv("C:/Users/KellyHeilman/Box/01. kelly.heilman Workspace/mortality/Eastern-Mortality/mortality_manuscript/Covariate_descriptions_table.csv")
Covariate.table %>% rename(`Covariate Group` = "Covariate.group", 
                           `Spatial Scale` = "Spatial.Scale")|> gt()|> 
  gtsave("C:/Users/KellyHeilman/Box/01. kelly.heilman Workspace/mortality/Eastern-Mortality/mortality_manuscript/figures/covariate_table_summary.png")
# script that generates all the testing and training datasets
source("R/speciesModels/SPCD_stan_data.R")
# write the data for all 26 different species groups:
for(i in 1:length(unique(nspp[1:17,]$SPCD))){
  cat(i)
  SPCD.stan.data(SPCD.id = nspp[i,]$SPCD, remper.correction = 0.5, cleaned.data.full = cleaned.data.full)
  # SPCD.stan.data(SPCD.id = nspp[i,]$SPCD, remper.correction = 0.9, cleaned.data.full = cleaned.data.full)
  # SPCD.stan.data(SPCD.id = nspp[i,]$SPCD, remper.correction = 0.1, cleaned.data.full = cleaned.data.full)
  # SPCD.stan.data(SPCD.id = nspp[i,]$SPCD, remper.correction = 0.3, cleaned.data.full = cleaned.data.full)
  # SPCD.stan.data(SPCD.id = nspp[i,]$SPCD, remper.correction = 0.7, cleaned.data.full = cleaned.data.full)
}
