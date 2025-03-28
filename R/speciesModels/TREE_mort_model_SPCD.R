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


for(m in 1:9){ 
  
  model.number <- model.list[m]
  for(i in 17:1){# run for each of the 17 species
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
for(i in 1:17){# run for each of the 17 species
  common.name <- nspp[1:17, ] %>% filter(SPCD %in% SPCD.df[i,]$SPCD) %>% dplyr::select(COMMON)
  
  for(m in c(1)){  # run each of the 9 models
    #for(m in 8:9){ 
    
    model.number <- model.list[m]
    for (j in 1:length(remper.cor.vector)){ # for the growth only model explore the consequences of other assumptions about remeasurement period
      cat(paste("running AUC stats for model ",model.number, " for SPCD", SPCD.df[i,]$SPCD, common.name$COMMON, " remper correction", remper.cor.vector[j]))
      
      SPCD.id <-  SPCD.df[i,]$SPCD
      
      # set up to make plots of the stan outputs 
      model.name <- paste0("mort_model_",model.number,"_SPCD_", SPCD.id, "_remper_correction_", remper.cor.vector[j])
      remp.cor <- remper.cor.vector[j]
      remper.correction <- remper.cor.vector[j]
      
      output.folder = "C:/Users/KellyHeilman/Box/01. kelly.heilman Workspace/mortality/Eastern-Mortality/mortality_models/"
      
      source("R/speciesModels/SPCD_AUC_stan.R")
      rm(fit.1)
    }
  }
}



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


##############################################################################3
# get the species-level covariates and make a big figure
# run for model 6
betas.list <- list()
SPCD.id
remper.cor.vector <- c(0.5)
model.no <- 6
output.folder = "C:/Users/KellyHeilman/Box/01. kelly.heilman Workspace/mortality/Eastern-Mortality/mortality_models/"

get.beta.covariates <- function(SPCD.id, remper.cor.vector, model.no){
  print(SPCD.id)
  for(j in 1:length(remper.cor.vector)){
    print (paste0( "remper correction vector ", remper.cor.vector[j]))
    fit <- readRDS(paste0(output.folder,"/SPCD_stanoutput_full_standardized/samples/model_",model.no,"_SPCD_",SPCD.id, "_remper_correction_", remper.cor.vector[j], ".RDS"))
    fit_ssm_df <- as_draws_df(fit) # takes awhile to convert to df
    
    # get all the covariates using posterior package
    betas.quant <- subset_draws(fit_ssm_df, variable = "u_beta") %>% summarise_draws(median, ~quantile(., probs = c(0.025, 0.975))) %>%
      rename(`ci.lo` = "2.5%", `ci.hi` = "97.5%") %>% 
      mutate(remper.cor = remper.cor.vector[j], 
             SPCD = SPCD.id)
    
    saveRDS(betas.quant, paste0(output.folder, "/SPCD_stanoutput_full_standardized/betas/model_",model.no,"_SPCD_",SPCD.id, "_remper_correction_", remper.cor.vector[j], ".RDS"))
    
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

# read in all the betas from just model 6, 0.5:
# read in all the betas:
all.files.05 <- paste0(output.folder, "SPCD_stanoutput_full_standardized/betas/", list.files(paste0(output.folder,"SPCD_stanoutput_full_standardized/betas/"), "model_6"))
big.betas <- lapply(all.files.05, readRDS)
betas.all.df <- do.call(rbind, big.betas)
betas.all.df$remper.cor <- as.character(betas.all.df$remper.cor)
betas.all.df$Species <- FIESTA::ref_species[match(betas.all.df$SPCD, ref_species$SPCD),]$COMMON

# get model parameter names
model.number <- 6
SPCD.id <- 318
remper.correction <- 0.5
load(paste0("SPCD_standata_general_full_standardized/SPCD_",SPCD.id, "remper_correction_", remper.correction,"model_",model.number, ".Rdata")) # load the species code data5
param.names <- data.frame(Parameter = colnames(mod.data$xM), 
                          variable = paste0("u_beta[", 1:51,"]"))

betas.all.df <- left_join(betas.all.df, param.names)
betas.all.df
spp.traits <- read.csv(paste0(output.folder, "NinemetsSpeciesTraits.csv"))
betas.all.df <- left_join(betas.all.df, spp.traits)

ggplot(data = na.omit(betas.all.df) %>% filter(variable %in% "u_beta[1]"), aes(x = Species, y = median, color = SFTWD_HRDWD))+geom_point()+
  geom_errorbar(data = na.omit(betas.all.df)  %>% filter(variable %in% "u_beta[1]"), aes(x = Species, ymin = ci.lo, ymax = ci.hi, color = SFTWD_HRDWD), width = 0.1)+
  geom_abline(aes(slope = 0, intercept = 0), color = "grey", linetype = "dashed")+
  theme_bw(base_size = 12)+
  theme( axis.text.x = element_text(angle = 45, hjust = 1))+ylab("Effect of average annual growth on probability of survival")

ggplot(data = na.omit(betas.all.df) %>% filter(variable %in% "u_beta[2]"), aes(x = Species, y = median))+geom_point()+
  geom_errorbar(data = na.omit(betas.all.df)  %>% filter(variable %in% "u_beta[2]"), aes(x = Species, ymin = ci.lo, ymax = ci.hi), width = 0.1)+
  geom_abline(aes(slope = 0, intercept = 0), color = "grey", linetype = "dashed")+
  theme_bw(base_size = 12)+
  theme( axis.text.x = element_text(angle = 45, hjust = 1))+ylab("Effect of Diameter on probability of survival")

ggplot(betas.all.df, aes(x = ShadeTol, y = median))+geom_point()+stat_smooth(method = "lm")+ facet_wrap(~Parameter, scales = "free_y")
ggsave(height = 6, width = 10, filename = "summaryfigures/betas_Shade_tolerance_model_6.png")
ggplot(betas.all.df, aes(x = DroughtTol, y = median))+geom_point()+stat_smooth(method = "lm")+ facet_wrap(~Parameter, scales = "free_y")
ggsave(height = 6, width = 10, filename = "summaryfigures/betas_Drought_tolerance_model_6.png")

ggplot(betas.all.df, aes(x = FloodTol, y = median))+geom_point()+stat_smooth(method = "lm")+ facet_wrap(~Parameter, scales = "free_y")
ggsave(height = 6, width = 10, filename = "summaryfigures/betas_Flood_tolerance_model_6.png")

ggplot(betas.all.df, aes(x = WOOD_SPGR_GREENVOL_DRYWT, y = median))+geom_point()+stat_smooth(method = "lm")+ facet_wrap(~Parameter, scales = "free_y")
ggsave(height = 6, width = 10, filename = "summaryfigures/betas_vol_drywt_model_6.png")

# lets just do for the fixed effects
# ordered by shade tolerance
betas.quant <- betas.all.df %>% arrange(by = ShadeTol) 
betas.quant$Species <- factor(betas.quant$Species, levels = unique(betas.quant$Species))
ggplot(data = na.omit(betas.quant) , aes(x = Species, y = median))+geom_point()+
  geom_errorbar(data = na.omit(betas.quant) , aes(x = Species, ymin = ci.lo, ymax = ci.hi), width = 0.1)+
  geom_abline(aes(slope = 0, intercept = 0), color = "grey", linetype = "dashed")+
  theme_bw(base_size = 12)+
  theme( axis.text.x = element_text(angle = 45, hjust = 1))+facet_wrap(~Parameter, scales = "free_y")

# ordered by drought tolerance
betas.quant <- betas.all.df %>% arrange(by = DroughtTol) 
betas.quant$Species <- factor(betas.quant$Species, levels = unique(betas.quant$Species))
ggplot(data = na.omit(betas.quant) , aes(x = Species, y = median))+geom_point()+
  geom_errorbar(data = na.omit(betas.quant) , aes(x = Species, ymin = ci.lo, ymax = ci.hi), width = 0.1)+
  geom_abline(aes(slope = 0, intercept = 0), color = "grey", linetype = "dashed")+
  theme_bw(base_size = 12)+
  theme( axis.text.x = element_text(angle = 45, hjust = 1))+facet_wrap(~Parameter, scales = "free_y")
ggsave(height = 5, width = 10, "summaryfigures/species_fixed_response_model_6_by_drought_tolerance.png")

# ordered by flood tolerance
betas.quant <- betas.all.df %>% arrange(by = FloodTol) 
betas.quant$Species <- factor(betas.quant$Species, levels = unique(betas.quant$Species))
ggplot(data = na.omit(betas.quant) , aes(x = Species, y = median))+geom_point()+
  geom_errorbar(data = na.omit(betas.quant) , aes(x = Species, ymin = ci.lo, ymax = ci.hi), width = 0.1)+
  geom_abline(aes(slope = 0, intercept = 0), color = "grey", linetype = "dashed")+
  theme_bw(base_size = 12)+
  theme( axis.text.x = element_text(angle = 45, hjust = 1))+facet_wrap(~Parameter, scales = "free_y")

# ordered by WOOD_SPGR_GREENVOL_DRYWT
betas.quant <- betas.all.df %>% arrange(by = WOOD_SPGR_GREENVOL_DRYWT) 
betas.quant$Species <- factor(betas.quant$Species, levels = unique(betas.quant$Species))
ggplot(data = na.omit(betas.quant) , aes(x = Species, y = median))+geom_point()+
  geom_errorbar(data = na.omit(betas.quant) , aes(x = Species, ymin = ci.lo, ymax = ci.hi), width = 0.1)+
  geom_abline(aes(slope = 0, intercept = 0), color = "grey", linetype = "dashed")+
  theme_bw(base_size = 12)+
  theme( axis.text.x = element_text(angle = 45, hjust = 1))+facet_wrap(~Parameter, scales = "free_y")



# lets just do for the fixed effects
# ordered by shade tolerance
betas.quant <- betas.all.df %>% arrange(by = ShadeTol)%>% filter(Parameter %in% param.names$Parameter[1:17])
betas.quant$Species <- factor(betas.quant$Species, levels = unique(betas.quant$Species))
ggplot(data = na.omit(betas.quant) , aes(x = Species, y = median))+geom_point()+
  geom_errorbar(data = na.omit(betas.quant) , aes(x = Species, ymin = ci.lo, ymax = ci.hi), width = 0.1)+
  geom_abline(aes(slope = 0, intercept = 0), color = "grey", linetype = "dashed")+
  theme_bw(base_size = 12)+
  
  theme( axis.text.x = element_text(angle = 45, hjust = 1))+facet_wrap(~Parameter, scales = "free_y")
ggsave(height = 5, width = 10, "summaryfigures/species_fixed_response_model_6_by_shade_tolerance.png")

# variables linked to Shade tolerence 
# BAL, Damage?, MAP, PHSYIO, RD

betas.quant$`significance` <- ifelse(betas.quant$ci.lo < 0 & betas.quant$ci.hi < 0, "significant", 
                                     ifelse(betas.quant$ci.lo > 0 & betas.quant$ci.hi > 0, "significant", "not overlapping zero"))

betas.quant %>% #filter(Parameter %in% "BAL.scaled") %>% 
  group_by(Species)  %>% 
  mutate(test = ifelse( median > 0 & significance %in% "significant", fontawesome::fa('plus-circle', fill = "#018571") |> html(), 
                        ifelse(median < 0 & significance %in% "significant", fontawesome::fa('minus-circle', fill = "#a6611a") |> html(),
                               ifelse( median > 0 & !significance %in% "significant", fontawesome::fa('plus-circle', fill = "#80cdc1") |> html(),
                                       
                                       fontawesome::fa('minus-circle', fill = "#dfc27d") |> html()))) ) %>% 
  select(Species, Parameter, test) %>% 
  spread( Species, test) %>% ungroup() |> gt() |> fmt_markdown()

parameter.group <- data.frame(Parameter = unique(betas.quant$Parameter), 
                              parameter.group = c("Growth", "Diameter", 
                                                  "Competition", "Competition", 
                                                  "Competition", "Competition", 
                                                  "Average Climate", "Average Climate", "Average Climate",
                                                  "Climate Anomoly","Climate Anomoly","Climate Anomoly",
                                                  "Site Characteristics","Site Characteristics","Site Characteristics","Site Characteristics","Site Characteristics" ))
betas.quant %>% #filter(Parameter %in% "BAL.scaled") %>% 
  group_by(Species)  %>% 
  mutate(test = ifelse( median > 0 & significance %in% "significant", fontawesome::fa('plus-circle', fill = "#018571") |> html(), 
                        ifelse(median < 0 & significance %in% "significant", fontawesome::fa('minus-circle', fill = "#a6611a") |> html(),
                               ifelse( median > 0 & !significance %in% "significant", fontawesome::fa('circle', fill = "#80cdc1") |> html(),
                                       
                                       fontawesome::fa('circle', fill = "#dfc27d") |> html()))) ) %>% 
  select(Species, Parameter, test) %>% 
  #filter(Parameter %in% c("BAL.scaled", "damage.scaled", "MAP.scaled")) %>% 
  spread( Species, test) %>% left_join(.,parameter.group) %>% 
  ungroup() %>% group_by(parameter.group) |> gt() |> fmt_markdown() |>
  tab_spanner(
    label = "shade intolerant",
    columns = 2:9) |>
  tab_spanner(
    label = "shade tolerant",
    columns = 10:13)|>
  tab_spanner(
    label = "very shade tolerant",
    columns = 14:18) |> gtsave(filename = "summaryfigures/species_model_6_summary_by_shade_tolerance_circles.html")


# only the values sorded by shade tolerance

betas.quant %>% #filter(Parameter %in% "BAL.scaled") %>% 
  group_by(Species)  %>% 
  mutate(test = ifelse( median > 0 & significance %in% "significant", fontawesome::fa('plus-circle', fill = "#018571") |> html(), 
                        ifelse(median < 0 & significance %in% "significant", fontawesome::fa('minus-circle', fill = "#a6611a") |> html(),
                               ifelse( median > 0 & !significance %in% "significant", fontawesome::fa('circle', fill = "#80cdc1") |> html(),
                                       
                                       fontawesome::fa('circle', fill = "#dfc27d") |> html()))) ) %>% 
  select(Species, Parameter, test) %>% 
  filter(Parameter %in% c("BAL.scaled", "damage.scaled", "MAP.scaled")) %>% 
  spread( Species, test) %>% left_join(.,parameter.group) %>% 
  ungroup() %>% group_by(parameter.group) |> gt() |> fmt_markdown() |>
  tab_spanner(
    label = "shade intolerant",
    columns = 2:9) |>
  tab_spanner(
    label = "shade tolerant",
    columns = 10:13)|>
  tab_spanner(
    label = "very shade tolerant",
    columns = 14:18) |> gtsave(filename = "summaryfigures/species_model_6_summary_by_shade_tolerance_circles_limited.html")


rep.scale = seq(-1, 1, by = 0.1)
sparkline <- "brown"
final_value <- "brown"
range_low <- "brown"
range_high <- "brown"
type <- "brown"

spkl_palette_decidious <- c(sparkline, final_value, range_low, range_high, type)

sparkline <- "forestgreen"
final_value <- "forestgreen"
range_low <- "forestgreen"
range_high <- "forestgreen"
type <- "forestgreen"

spkl_palette_conifer <- c(sparkline, final_value, range_low, range_high, type)


betas.quant %>% group_by(Parameter, Species) %>% summarise(predicted.values = list(exp(median*rep.scale)/(1+exp(median*rep.scale))), .groups = "drop") %>%
  group_by(Parameter) %>% spread(Species, predicted.values) %>% 
  #select( S, annual.growth.scaled, BAL.scaled, MAP.scaled) %>%
  ungroup() %>% left_join(.,parameter.group) %>% 
  ungroup()%>% 
  
  group_by(parameter.group) |> 
  gt() |>  
  gtExtras::gt_plt_sparkline(`yellow-poplar`,  fig_dim = c(5, 10),  label = FALSE, palette = spkl_palette_decidious) |>
  gtExtras::gt_plt_sparkline(`white ash`, label = FALSE, fig_dim = c(5, 10), palette = spkl_palette_decidious) |>
  gtExtras::gt_plt_sparkline(`black cherry`,  label = FALSE,fig_dim = c(5, 10), palette = spkl_palette_decidious)|>
  gtExtras::gt_plt_sparkline(`black oak`,  label = FALSE, fig_dim = c(5, 10), palette = spkl_palette_decidious) |>
  gtExtras::gt_plt_sparkline(`hickory spp.`,  label = FALSE, fig_dim = c(5, 10), palette = spkl_palette_decidious) |>
  gtExtras::gt_plt_sparkline(`northern red oak`,  label = FALSE, fig_dim = c(5, 10), palette = spkl_palette_decidious)|>
  gtExtras::gt_plt_sparkline(`white oak`,  label = FALSE, fig_dim = c(5, 10), palette = spkl_palette_decidious) |>
  gtExtras::gt_plt_sparkline(`chestnut oak`,  label = FALSE, fig_dim = c(5, 10), palette = spkl_palette_decidious) |>
  gtExtras::gt_plt_sparkline(`yellow birch`,  label = FALSE, fig_dim = c(5, 10), palette = spkl_palette_decidious) |>
  gtExtras::gt_plt_sparkline(`eastern white pine`,  label = FALSE, fig_dim = c(5, 10),palette = spkl_palette_conifer) |>
  gtExtras::gt_plt_sparkline(`red maple`,  label = FALSE, fig_dim = c(5, 10), palette = spkl_palette_decidious) |>
  gtExtras::gt_plt_sparkline(`northern white-cedar`,  label = FALSE, fig_dim = c(5, 10),palette = spkl_palette_conifer)|>
  gtExtras::gt_plt_sparkline(`red spruce`,  label = FALSE, fig_dim = c(5, 10),palette = spkl_palette_conifer) |>
  gtExtras::gt_plt_sparkline(`American beech`, label = FALSE, fig_dim = c(5, 10), palette = spkl_palette_decidious) |>
  gtExtras::gt_plt_sparkline(`sugar maple`,  label = FALSE, fig_dim = c(5, 10), palette = spkl_palette_decidious)|>
  gtExtras::gt_plt_sparkline(`eastern hemlock`,  label = FALSE, fig_dim = c(5, 10),palette = spkl_palette_conifer) |>
  gtExtras::gt_plt_sparkline(`balsam fir`,  label = FALSE, fig_dim = c(5, 10), palette = spkl_palette_conifer)|>
  tab_spanner(
    label = "shade intolerant",
    columns = 2:9) |>
  tab_spanner(
    label = "shade tolerant",
    columns = 10:13)|>
  tab_spanner(
    label = "very shade tolerant",
    columns = 14:18)|> gtsave(filename = "summaryfigures/species_model_6_summary_by_shade_tolerance.html")

# select only the variables that seem to be correlated to shade tolerance:

betas.quant %>% group_by(Parameter, Species) %>% summarise(predicted.values = list(exp(median*rep.scale)/(1+exp(median*rep.scale))), .groups = "drop") %>%
  group_by(Parameter) %>% spread(Species, predicted.values) %>% 
  #select( S, annual.growth.scaled, BAL.scaled, MAP.scaled) %>%
  ungroup() %>% left_join(.,parameter.group) %>% 
  ungroup()%>% filter(Parameter %in% c("BAL.scaled", "RD.scaled", "MAP.scaled", "Physio.scaled", "damage.scaled")) %>%
  group_by(parameter.group) |> 
  gt() |>  
  gtExtras::gt_plt_sparkline(`yellow-poplar`,  fig_dim = c(5, 10),  label = FALSE, palette = spkl_palette_decidious) |>
  gtExtras::gt_plt_sparkline(`white ash`, label = FALSE, fig_dim = c(5, 10), palette = spkl_palette_decidious) |>
  gtExtras::gt_plt_sparkline(`black cherry`,  label = FALSE,fig_dim = c(5, 10), palette = spkl_palette_decidious)|>
  gtExtras::gt_plt_sparkline(`black oak`,  label = FALSE, fig_dim = c(5, 10), palette = spkl_palette_decidious) |>
  gtExtras::gt_plt_sparkline(`hickory spp.`,  label = FALSE, fig_dim = c(5, 10), palette = spkl_palette_decidious) |>
  gtExtras::gt_plt_sparkline(`northern red oak`,  label = FALSE, fig_dim = c(5, 10), palette = spkl_palette_decidious)|>
  gtExtras::gt_plt_sparkline(`white oak`,  label = FALSE, fig_dim = c(5, 10), palette = spkl_palette_decidious) |>
  gtExtras::gt_plt_sparkline(`chestnut oak`,  label = FALSE, fig_dim = c(5, 10), palette = spkl_palette_decidious) |>
  gtExtras::gt_plt_sparkline(`yellow birch`,  label = FALSE, fig_dim = c(5, 10), palette = spkl_palette_decidious) |>
  gtExtras::gt_plt_sparkline(`eastern white pine`,  label = FALSE, fig_dim = c(5, 10),palette = spkl_palette_conifer) |>
  gtExtras::gt_plt_sparkline(`red maple`,  label = FALSE, fig_dim = c(5, 10), palette = spkl_palette_decidious) |>
  gtExtras::gt_plt_sparkline(`northern white-cedar`,  label = FALSE, fig_dim = c(5, 10),palette = spkl_palette_conifer)|>
  gtExtras::gt_plt_sparkline(`red spruce`,  label = FALSE, fig_dim = c(5, 10),palette = spkl_palette_conifer) |>
  gtExtras::gt_plt_sparkline(`American beech`, label = FALSE, fig_dim = c(5, 10), palette = spkl_palette_decidious) |>
  gtExtras::gt_plt_sparkline(`sugar maple`,  label = FALSE, fig_dim = c(5, 10), palette = spkl_palette_decidious)|>
  gtExtras::gt_plt_sparkline(`eastern hemlock`,  label = FALSE, fig_dim = c(5, 10),palette = spkl_palette_conifer) |>
  gtExtras::gt_plt_sparkline(`balsam fir`,  label = FALSE, fig_dim = c(5, 10), palette = spkl_palette_conifer)|>
  tab_spanner(
    label = "shade intolerant",
    columns = 2:9) |>
  tab_spanner(
    label = "shade tolerant",
    columns = 10:13)|>
  tab_spanner(
    label = "very shade tolerant",
    columns = 14:18)|> gtsave(filename = "summaryfigures/species_model_6_summary_by_shade_tolerance_four_vars.html")


betas.quant %>% #filter(Parameter %in% "BAL.scaled") %>% 
  group_by(Species)  %>% 
  mutate(test = ifelse( median > 0 & significance %in% "significant", fontawesome::fa('plus-circle', fill = "#018571") |> html(), 
                        ifelse(median < 0 & significance %in% "significant", fontawesome::fa('minus-circle', fill = "#a6611a") |> html(),
                               ifelse( median > 0 & !significance %in% "significant", fontawesome::fa('circle', fill = "#80cdc1") |> html(),
                                       
                                       fontawesome::fa('circle', fill = "#dfc27d") |> html()))) ) %>% 
  select(Species, Parameter, test) %>% ungroup()%>%
  spread( Parameter, test) %>% #left_join(.,parameter.group) %>% 
  ungroup() |> gt() |> fmt_markdown()|> gtsave(filename = "summaryfigures/species_model_6_summary_by_shade_tolerance_plus_circles.html")

#-----------------------------------------------------------------------------------------
# ordered by drought tolerance
#-----------------------------------------------------------------------------------------
betas.quant <- betas.all.df %>% arrange(by = DroughtTol) %>% filter(Parameter %in% param.names$Parameter[1:18])
betas.quant$Species <- factor(betas.quant$Species, levels = unique(betas.quant$Species))
ggplot(data = na.omit(betas.quant) , aes(x = Species, y = median))+geom_point()+
  geom_errorbar(data = na.omit(betas.quant) , aes(x = Species, ymin = ci.lo, ymax = ci.hi), width = 0.1)+
  geom_abline(aes(slope = 0, intercept = 0), color = "grey", linetype = "dashed")+
  theme_bw(base_size = 12)+
  theme( axis.text.x = element_text(angle = 45, hjust = 1))+facet_wrap(~Parameter, scales = "free_y")
ggsave(height = 5, width = 10, "summaryfigures/species_fixed_response_model_6_by_drought_tolerance.png")
# variables linked to Drought tolerence 
# Growth?, BAL, damage, DIA_growth interaction (more negative w/ drought tolerance), MAP, Slope?

betas.quant$`significance` <- ifelse(betas.quant$ci.lo < 0 & betas.quant$ci.hi < 0, "significant", 
                                     ifelse(betas.quant$ci.lo > 0 & betas.quant$ci.hi > 0, "significant", "not overlapping zero"))

betas.quant %>% #filter(Parameter %in% "BAL.scaled") %>% 
  group_by(Species)  %>% 
  mutate(test = ifelse( median > 0 & significance %in% "significant", fontawesome::fa('plus-circle', fill = "forestgreen") |> html(), 
                        ifelse(median < 0 & significance %in% "significant", fontawesome::fa('minus-circle', fill = "brown") |> html(),
                               ifelse( median > 0 & !significance %in% "significant", fontawesome::fa('plus-circle', fill = "lightgrey") |> html(),
                                       
                                       fontawesome::fa('minus-circle', fill = "lightgrey") |> html()))) ) %>% 
  select(Species, Parameter, test) %>% 
  spread( Species, test) %>% ungroup() |> gt() |> fmt_markdown()|>
  gtsave(filename = "summaryfigures/species_model_6_summary_by_drought_tolerance_circle.html")

# color by the magnitude of the effect

# create a table with the hexcolors
hex.table <- data.frame(bin = c("(-7,-5]", 
                                "(-5,-3]", 
                                "(-3,-1]", 
                                "(-1,0]", 
                                "(0,1]",
                                "(1,3]",
                                "(3,5]",
                                "(5,6]"), 
                        hexcolor = c("#8c510a",
                                     "#bf812d",
                                     "#dfc27d",
                                     "#f6e8c3",
                                     "#c7eae5",
                                     "#80cdc1",
                                     "#35978f",
                                     "#01665e"))


betas.spread <- data.frame(betas.quant %>% select(Parameter, median, Species) %>% spread(Parameter, median))
rownames(betas.spread) <- betas.spread$Species
res.hc <- hclust(dist( betas.spread[, 2:ncol(betas.spread)]),  method = "ward.D2")
fviz_dend(res.hc, cex = 0.5, k = 4, palette = "jco") 
fviz_dend(res.hc, cex = 0.5, k = 5, palette = "jco")

library(pheatmap)
pheatmap(t(betas.spread[, 2:ncol(betas.spread)]), cutree_cols = 4)

kmeans(mydata, 3, nstart = 25)
betas.quant  %>% mutate(bin = cut(median, breaks = c(-7,-5,-3,-1,0,1,3,5,6))) %>%
  left_join(., hex.table) %>%  group_by(Species, Parameter) %>%
  mutate(test = ifelse( median > 0 & significance %in% "significant", fontawesome::fa('plus-circle', fill = hexcolor) |> html(), 
                        ifelse(median < 0 & significance %in% "significant", fontawesome::fa('minus-circle', fill = hexcolor) |> html(),
                               ifelse( median > 0 & !significance %in% "significant", fontawesome::fa('plus-circle', fill = "lightgrey") |> html(),
                                       
                                       fontawesome::fa('minus-circle', fill = "lightgrey") |> html()))) ) %>% 
  select(Species, Parameter, test) %>% ungroup() %>%  
  # reorder factors
  mutate(Parameter = factor(Parameter, levels = c(# mostly positive
    "annual.growth.scaled", 
    
    # mostly negative
    "DIA_scaled", 
    "damage.scaled", 
    "ppt.anom", 
    "ba.scaled",
    
    # increasingly positive with DT
    "MAP.scaled", 
    "physio.scaled",
    "elev.scaled",
    "RD.scaled", 
    "non_SPCD.BA.scaled",
    "slope.scaled",
    "MATmax.scaled",
    
    # increasingly negative with DT
    "Ndep.scaled",
    "DIA_scaled_growth.int",
    "BAL.scaled",
    "tmin.anom", 
    "tmax.anom", 
    "MATmin.scaled",
    "aspect.scaled"))) %>% 
  group_by(Parameter) %>% 
  
  
  spread( Species, test) %>% ungroup() |> 
  gt() |> 
  opt_css(
    css = "
    #mygt .gt_col_heading {
      text-align: center;
      transform: rotate(-90deg);
      font-weight: bold;
    }
    "
  )|>
  fmt_markdown(columns = 2:18) |>
  gtsave(filename = "summaryfigures/species_model_6_summary_by_drought_tolerance_circle_scaled_color.html")

# MAP.scaled, Physio.scaled, tmax.anom
betas.quant %>% filter(Parameter %in% c("annual.growth.scaled", "damage.scaled", "elev.scaled", "MAP.scaled", 
                                        "MATmin.scaled", "Ndep.scaled", "Physio.scaled", "slope.scaled")) %>% 
  group_by(Species)  %>% 
  mutate(test = ifelse( median > 0 & significance %in% "significant", fontawesome::fa('plus-circle', fill = "forestgreen") |> html(), 
                        ifelse(median < 0 & significance %in% "significant", fontawesome::fa('minus-circle', fill = "brown") |> html(),
                               ifelse( median > 0 & !significance %in% "significant", fontawesome::fa('plus-circle', fill = "lightgrey") |> html(),
                                       
                                       fontawesome::fa('minus-circle', fill = "lightgrey") |> html()))) ) %>% 
  select(Species, Parameter, test) %>% 
  spread( Species, test) %>% ungroup() |> gt() |> fmt_markdown()|>
  gtsave(filename = "summaryfigures/species_model_6_summary_by_drought_tolerance_circle_limited.html")

unique(betas.quant[,c("Species","DroughtTol")])

betas.quant %>% group_by(Parameter, Species) %>% summarise(predicted.values = list(exp(median*rep.scale)/(1+exp(median*rep.scale))), .groups = "drop") %>%
  group_by(Parameter) %>% spread(Species, predicted.values) %>% 
  #select( S, annual.growth.scaled, BAL.scaled, MAP.scaled) %>%
  ungroup() %>% left_join(.,parameter.group) %>% 
  ungroup()%>% 
  
  group_by(parameter.group) |> 
  gt() |>  
  gtExtras::gt_plt_sparkline(`yellow-poplar`,  fig_dim = c(5, 10),  label = FALSE, palette = spkl_palette_decidious) |>
  gtExtras::gt_plt_sparkline(`white ash`, label = FALSE, fig_dim = c(5, 10), palette = spkl_palette_decidious) |>
  gtExtras::gt_plt_sparkline(`black cherry`,  label = FALSE,fig_dim = c(5, 10), palette = spkl_palette_decidious)|>
  gtExtras::gt_plt_sparkline(`black oak`,  label = FALSE, fig_dim = c(5, 10), palette = spkl_palette_decidious) |>
  gtExtras::gt_plt_sparkline(`hickory spp.`,  label = FALSE, fig_dim = c(5, 10), palette = spkl_palette_decidious) |>
  gtExtras::gt_plt_sparkline(`northern red oak`,  label = FALSE, fig_dim = c(5, 10), palette = spkl_palette_decidious)|>
  gtExtras::gt_plt_sparkline(`white oak`,  label = FALSE, fig_dim = c(5, 10), palette = spkl_palette_decidious) |>
  gtExtras::gt_plt_sparkline(`chestnut oak`,  label = FALSE, fig_dim = c(5, 10), palette = spkl_palette_decidious) |>
  gtExtras::gt_plt_sparkline(`yellow birch`,  label = FALSE, fig_dim = c(5, 10), palette = spkl_palette_decidious) |>
  gtExtras::gt_plt_sparkline(`eastern white pine`,  label = FALSE, fig_dim = c(5, 10),palette = spkl_palette_conifer) |>
  gtExtras::gt_plt_sparkline(`red maple`,  label = FALSE, fig_dim = c(5, 10), palette = spkl_palette_decidious) |>
  gtExtras::gt_plt_sparkline(`northern white-cedar`,  label = FALSE, fig_dim = c(5, 10),palette = spkl_palette_conifer)|>
  gtExtras::gt_plt_sparkline(`red spruce`,  label = FALSE, fig_dim = c(5, 10),palette = spkl_palette_conifer) |>
  gtExtras::gt_plt_sparkline(`American beech`, label = FALSE, fig_dim = c(5, 10), palette = spkl_palette_decidious) |>
  gtExtras::gt_plt_sparkline(`sugar maple`,  label = FALSE, fig_dim = c(5, 10), palette = spkl_palette_decidious)|>
  gtExtras::gt_plt_sparkline(`eastern hemlock`,  label = FALSE, fig_dim = c(5, 10),palette = spkl_palette_conifer) |>
  gtExtras::gt_plt_sparkline(`balsam fir`,  label = FALSE, fig_dim = c(5, 10), palette = spkl_palette_conifer)|>
  tab_spanner(
    label = "drought intolerant",
    columns = 2:5) |>
  tab_spanner(
    label = "moderate drought tolerance",
    columns = 6:12)|>
  tab_spanner(
    label = "drought tolerant",
    columns = 13:16)|>
  tab_spanner(
    label = "very drought tolerant",
    columns = 17:18)|>gtsave(filename = "summaryfigures/species_model_6_summary_by_drought_tolerance.html")


betas.quant %>% group_by(Parameter, Species) %>% summarise(predicted.values = list(exp(median*rep.scale)/(1+exp(median*rep.scale))), .groups = "drop") %>%
  group_by(Parameter) %>% spread(Species, predicted.values) %>% 
  filter( Parameter %in% c("annual.growth.scaled", "BAL.scaled", "MAP.scaled", 
                           "MATmax.scaled", "MATmin.scaled",
                           "Physio.scaled", "tmax.anom")) %>%
  ungroup() %>% left_join(.,parameter.group) %>% 
  ungroup()%>% 
  
  group_by(parameter.group) |> 
  gt() |>  
  gtExtras::gt_plt_sparkline(`yellow-poplar`,  fig_dim = c(5, 10),  label = FALSE, palette = spkl_palette_decidious) |>
  gtExtras::gt_plt_sparkline(`white ash`, label = FALSE, fig_dim = c(5, 10), palette = spkl_palette_decidious) |>
  gtExtras::gt_plt_sparkline(`black cherry`,  label = FALSE,fig_dim = c(5, 10), palette = spkl_palette_decidious)|>
  gtExtras::gt_plt_sparkline(`black oak`,  label = FALSE, fig_dim = c(5, 10), palette = spkl_palette_decidious) |>
  gtExtras::gt_plt_sparkline(`hickory spp.`,  label = FALSE, fig_dim = c(5, 10), palette = spkl_palette_decidious) |>
  gtExtras::gt_plt_sparkline(`northern red oak`,  label = FALSE, fig_dim = c(5, 10), palette = spkl_palette_decidious)|>
  gtExtras::gt_plt_sparkline(`white oak`,  label = FALSE, fig_dim = c(5, 10), palette = spkl_palette_decidious) |>
  gtExtras::gt_plt_sparkline(`chestnut oak`,  label = FALSE, fig_dim = c(5, 10), palette = spkl_palette_decidious) |>
  gtExtras::gt_plt_sparkline(`yellow birch`,  label = FALSE, fig_dim = c(5, 10), palette = spkl_palette_decidious) |>
  gtExtras::gt_plt_sparkline(`eastern white pine`,  label = FALSE, fig_dim = c(5, 10),palette = spkl_palette_conifer) |>
  gtExtras::gt_plt_sparkline(`red maple`,  label = FALSE, fig_dim = c(5, 10), palette = spkl_palette_decidious) |>
  gtExtras::gt_plt_sparkline(`northern white-cedar`,  label = FALSE, fig_dim = c(5, 10),palette = spkl_palette_conifer)|>
  gtExtras::gt_plt_sparkline(`red spruce`,  label = FALSE, fig_dim = c(5, 10),palette = spkl_palette_conifer) |>
  gtExtras::gt_plt_sparkline(`American beech`, label = FALSE, fig_dim = c(5, 10), palette = spkl_palette_decidious) |>
  gtExtras::gt_plt_sparkline(`sugar maple`,  label = FALSE, fig_dim = c(5, 10), palette = spkl_palette_decidious)|>
  gtExtras::gt_plt_sparkline(`eastern hemlock`,  label = FALSE, fig_dim = c(5, 10),palette = spkl_palette_conifer) |>
  gtExtras::gt_plt_sparkline(`balsam fir`,  label = FALSE, fig_dim = c(5, 10), palette = spkl_palette_conifer)|>
  tab_spanner(
    label = "drought intolerant",
    columns = 2:5) |>
  tab_spanner(
    label = "moderate drought tolerance",
    columns = 6:12)|>
  tab_spanner(
    label = "drought tolerant",
    columns = 13:16)|>
  tab_spanner(
    label = "very drought tolerant",
    columns = 17:18)|>gtsave(filename = "summaryfigures/species_model_6_summary_by_drought_tolerance_limited.html")


# ordered by flood tolerance
betas.quant <- betas.all.df %>% arrange(by = FloodTol) %>% filter(Parameter %in% param.names$Parameter[1:18])
betas.quant$Species <- factor(betas.quant$Species, levels = unique(betas.quant$Species))
ggplot(data = na.omit(betas.quant) , aes(x = Species, y = median))+geom_point()+
  geom_errorbar(data = na.omit(betas.quant) , aes(x = Species, ymin = ci.lo, ymax = ci.hi), width = 0.1)+
  geom_abline(aes(slope = 0, intercept = 0), color = "grey", linetype = "dashed")+
  theme_bw(base_size = 12)+
  theme( axis.text.x = element_text(angle = 45, hjust = 1))+facet_wrap(~Parameter, scales = "free_y")
# variables linked to flod tolerence 
# Elevation?


# ordered by WOOD_SPGR_GREENVOL_DRYWT
betas.quant <- betas.all.df %>% arrange(by = WOOD_SPGR_GREENVOL_DRYWT) %>% filter(Parameter %in% param.names$Parameter[1:18])
betas.quant$Species <- factor(betas.quant$Species, levels = unique(betas.quant$Species))
ggplot(data = na.omit(betas.quant) , aes(x = Species, y = median))+geom_point()+
  geom_errorbar(data = na.omit(betas.quant) , aes(x = Species, ymin = ci.lo, ymax = ci.hi), width = 0.1)+
  geom_abline(aes(slope = 0, intercept = 0), color = "grey", linetype = "dashed")+
  theme_bw(base_size = 12)+
  theme( axis.text.x = element_text(angle = 45, hjust = 1))+facet_wrap(~Parameter, scales = "free_y")


# read in all the betas:
all.files <- paste0("SPCD_stanoutput_full/betas/", list.files(paste0(output.folder,"SPCD_stanoutput_full/betas/")))
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



