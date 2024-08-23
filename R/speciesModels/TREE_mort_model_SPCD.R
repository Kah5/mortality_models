library(rstan)
library(MASS)
library(here)
library(tidyverse)
library(gt)
library(FIESTA)
library(dplyr)
library(mltools)

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
Covariate.table %>% rename(`Covariate Group` = "Covariate.group", 
                           `Spatial Scale` = "Spatial.Scale")|> gt()
# script that generates all the testing and training datasets
source("R/speciesModels/SPCD_stan_data.R")
# write the data for all 26 different species groups:
for(i in 1:length(unique(nspp[1:17,]$SPCD))){
  cat(i)
  SPCD.stan.data(SPCD.id = nspp[i,]$SPCD, remper.correction = 0.5, cleaned.data.full = cleaned.data.full)
   SPCD.stan.data(SPCD.id = nspp[i,]$SPCD, remper.correction = 0.9, cleaned.data.full = cleaned.data.full)
   SPCD.stan.data(SPCD.id = nspp[i,]$SPCD, remper.correction = 0.1, cleaned.data.full = cleaned.data.full)
   SPCD.stan.data(SPCD.id = nspp[i,]$SPCD, remper.correction = 0.3, cleaned.data.full = cleaned.data.full)
   SPCD.stan.data(SPCD.id = nspp[i,]$SPCD, remper.correction = 0.7, cleaned.data.full = cleaned.data.full)
}


#----------------------------------------------------------------------------------
# running stan models with most important variables
#----------------------------------------------------------------------------------


# for each species group, fit a model, plot the outputs, and save the results
# we source a function from another script
source("R/speciesModels/SPCD_run_stan.R")

# this runs a stan model and saves the outputs
# SPGRPCD 2 throws uncerialize socklist error--too big of a diataset to parallelisze?
cleaned.data.full %>% group_by(SPCD) %>% summarise(n())
SPCD.df <- data.frame(SPCD = nspp[1:17, ]$SPCD, 
                      spcd.id = 1:17)
remper.cor.vector <- c(0.1, 0.3, 0.7, 0.9)
#model.number <- 6
model.list <- 6

for(i in 1:17){# run for each of the 17 species
  common.name <- nspp[1:17, ] %>% filter(SPCD %in% SPCD.df[i,]$SPCD) %>% dplyr::select(COMMON)
  
  for(m in 1:length(model.list)){  # run each of the 9 models
    model.number <- model.list[m]
    for (j in 1:length(remper.cor.vector)){ # for the gowoth only model explore the consequences of other assumptions about remeasurement period
      cat(paste("running stan mortality model ",model.number, " for SPCD", SPCD.df[i,]$SPCD, common.name$COMMON, " remper correction", remper.cor.vector[j]))
      
      fit.1 <- SPCD_run_stan(SPCD.id = SPCD.df[i,]$SPCD,
                             model.no = model.number,
                             niter = 3000,
                             nchains = 3,
                             remper.correction = remper.cor.vector[j],
                            model.file = 'modelcode/mort_model_general.stan' )
      SPCD.id <-  SPCD.df[i,]$SPCD
      #saveRDS(fit.1, paste0("SPCD_stanoutput_full/samples/model_",model.number,"_SPCD_",SPCD.id, "_remper_correction_", remper.cor.vector[j], ".RDS"))
     # save_diagnostics (stanfitobj = fit.1, nchains = 2, model.no = model.number, remper.correction = remper.cor.vector[j])
      
      model.name <- paste0("mort_model_",model.number,"_SPCD_", SPCD.id, "_remper_correction_", remper.cor.vector[j])
      remp.cor <- remper.cor.vector[j]
      remper.correction <- remper.cor.vector[j]
      source("R/SPCD_plot_stan.R")
      rm(fit.1)
    }
  }
}





# plotting effects of mortality model by species to compare remper growth effects
# for SPCD 316, plot all estimated values:

# run for model 6
betas.list <- list()
SPCD.id
remper.cor.vector <- c(0.1, 0.3, 0.5, 0.7, 0.9)
model.no <- 6

get.beta.covariates <- function(SPCD.id, remper.cor.vector, model.no){
  print(SPCD.id)
  for(j in 1:length(remper.cor.vector)){
    print (paste0( "remper correction vector ", remper.cor.vector[j]))
    fit <- readRDS(paste0("SPCD_stanoutput_full/samples/model_",model.no,"_SPCD_",SPCD.id, "_remper_correction_", remper.cor.vector[j], ".RDS"))
    fit_ssm_df <- as_draws_df(fit) # takes awhile to convert to df
    
    # get all the covariates using posterior package
    betas.quant <- subset_draws(fit_ssm_df, variable = "u_beta") %>% summarise_draws(median, ~quantile(., probs = c(0.025, 0.975))) %>%
      rename(`ci.lo` = "2.5%", `ci.hi` = "97.5%") %>% 
      mutate(remper.cor = remper.cor.vector[j], 
             SPCD = SPCD.id)
    
    saveRDS(betas.quant, paste0("SPCD_stanoutput_full/betas/model_",model.no,"_SPCD_",SPCD.id, "_remper_correction_", remper.cor.vector[j], ".RDS"))
    
    # betas.list[[j]] <- betas.quant
    
    rm(fit)
    #betas.list
  }
  # betas.df <- do.call(rbind, betas.list)
}

big.betas <- list()
for(i in 1:17){
  get.beta.covariates(SPCD.id = SPCD.df[i,]$SPCD, remper.cor.vector = remper.cor.vector, model.no = 6)
}


# read in all the betas:
all.files <- paste0("SPCD_stanoutput_full/betas/", list.files("SPCD_stanoutput_full/betas/"))
big.betas <- lapply(all.files, readRDS)

betas.all.df <- do.call(rbind, big.betas)
betas.all.df %>% filter(variable %in% "u_beta[1]") %>% group_by(SPCD) %>% summarise(n())
betas.all.df$remper.cor <- as.character(betas.all.df$remper.cor)
betas.all.df$Species <- FIESTA::ref_species[match(betas.all.df$SPCD, ref_species$SPCD),]$COMMON

ggplot(data = na.omit(betas.all.df) %>% filter(variable %in% "u_beta[1]"), aes(x = remper.cor, y = median, fill = remper.cor))+geom_bar(stat = "identity", position = position_dodge())+
  geom_errorbar(data = na.omit(betas.all.df)  %>% filter(variable %in% "u_beta[1]"), aes(x = remper.cor, ymin = ci.lo, ymax = ci.hi, group = remper.cor), width = 0.1, position=position_dodge())+
  geom_abline(aes(slope = 0, intercept = 0), color = "grey", linetype = "dashed")+facet_wrap(~Species)+
  theme_bw(base_size = 12)+
  theme( axis.text.x = element_text(angle = 45, hjust = 1))+ylab("Effect of average annual growth on probability of survival")+xlab("Mortality year assumption (proportion of remper)")
ggsave("SPCD_stanoutput_full/images/sensitivity_avg_growth_to_remper.png")
  
ggplot(data = na.omit(betas.all.df) %>% filter(!SPCD %in% 621), aes(x = variable, y = median, color = remper.cor, group = variable))+geom_point(position=position_dodge(width = 2))+
  geom_errorbar(data = na.omit(betas.all.df) %>% filter(!SPCD %in% 621), aes(x = variable, ymin = ci.lo, ymax = ci.hi, color = remper.cor, group = variable), width = 0.1, position=position_dodge(width = 2))+
  geom_abline(aes(slope = 0, intercept = 0), color = "grey", linetype = "dashed")+facet_wrap(~SPCD)+
  theme_bw(base_size = 12)+
  theme( axis.text.x = element_text(angle = 45, hjust = 1))+ylab("Effect on probability of mortality")+xlab("Covariate")




# summarise all effects by spcies
plot_spp_effects <- function(SPP){
  ggplot(data = na.omit(betas.all.df) %>% filter(SPCD %in% SPP), aes(x = remper.cor, y = median, color = remper.cor, group = variable))+geom_point(position=position_dodge(1))+
    geom_errorbar(data = na.omit(betas.all.df) %>% filter(SPCD %in% SPP), aes(x = remper.cor, ymin = ci.lo, ymax = ci.hi, color = remper.cor, group = variable), width = 0.1, position=position_dodge())+
    geom_abline(aes(slope = 0, intercept = 0), color = "red", linetype = "dashed")+facet_wrap(~variable, scales = "free_y")+
    theme_bw(base_size = 12)+ggtitle(paste0( unique(betas.all.df %>% filter(SPCD %in% SPP) %>% select(Species))))+
    theme( axis.text.x = element_text(angle = 45, hjust = 1))+ylab("Effect on probability of mortality")+xlab("porportion of remper")
  ggsave(paste0("SPCD_stanoutput_full/images/all_model_effects_remper_cor_SPCD_", SPP, ".png"), height = 5, width = 8)
}
plot_spp_effects(SPP = 318)
lapply(unique(betas.all.df$SPCD), plot_spp_effects)

# generate parameter-specific tables of parameter effects

growth.table <- betas.all.df %>% filter(variable %in% "beta_growth") %>% 
  mutate(sig = ifelse(ci.lo < 0 & ci.hi < 0, "-", 
                      ifelse(ci.lo > 0 & ci.hi > 0, "+", "o")))%>%
  group_by(SPCD) %>% select(Species, SPCD, remper.cor, sig) %>% 
  spread(remper.cor, sig)%>% ungroup()|> gt() |>tab_header(
    title = "Effects of growth on p(mort)") |> tab_spanner(label = "% of remper", columns = starts_with("0."))

table.beta.param <- function(beta.param){
  clim.norms.table <- betas.all.df %>% filter(variable %in% beta.param) %>% 
    mutate(sig = factor(ifelse(ci.lo < 0 & ci.hi < 0, "-", 
                               ifelse(ci.lo > 0 & ci.hi > 0, "+", "o"))))%>%
    #mutate()
    group_by(variable, SPCD) %>% select(Species, SPCD, variable, remper.cor, sig) %>% 
    spread(remper.cor, sig)%>% ungroup()|> gt(groupname_col = "variable") |>tab_header(
      title = "Effects on p(mort)") |> tab_spanner(label = "% of remper", columns = starts_with("0."))|>
    cols_align("left", variable)|> tab_options(
      summary_row.background.color = "gray95",
      row_group.background.color = "#FFEFDB",
      row_group.as_column = TRUE
    ) |> 
    opt_align_table_header(align = "left")|>
    data_color(
      columns = starts_with("0."),
      method = "factor",
      levels = c("o", "+", "-"),
      palette = "viridis",
      domain = c(150E6, 170E6),
      reverse = TRUE
    )
  
  clim.norms.table
}

b.names <- unique(betas.all.df$variable)
# need to list these all out because you need chrome or chromium to save these tables?
table.beta.param(b.names[1])
table.beta.param(b.names[2])
table.beta.param(b.names[3])
table.beta.param(b.names[4])
table.beta.param(b.names[5])
table.beta.param(b.names[6])
table.beta.param(b.names[7])
table.beta.param(b.names[8])
table.beta.param(b.names[9])
table.beta.param(b.names[10])
table.beta.param(b.names[11])
table.beta.param(b.names[12])
table.beta.param(b.names[13])
table.beta.param(b.names[14])



# run the marginal posteriors for all species and generate plots
#source("R/plot_marginal_posteriors_SPCD.R")

# generate predictions for the test datasets
# how to evaluate models?
# predicted pmort vs observed mortality frequency?
# classification errors?

