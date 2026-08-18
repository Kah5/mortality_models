# `bayesian_model_averaging.R` ---
# This script relies outputs saved from prior files:
# marginal effect generated quantities
# `cmdstan_model_assessment.R` for the ELPD_LOO_df object that has the LOO estimated model weights
# helper functions

# ----------------------------------------------------------------------
# Setting up ---
# ----------------------------------------------------------------------
library(cmdstanr)
library(posterior)
library(dplyr)
library(FIESTA)
library(qs2)
library(ggplot2)
library(tidyverse)

# load some helper functions, including main_effects, int_labels, build_X_grid()
source("R/modelAssessment/bma_prediction_helpers.R")  

# set up output directories
output.dir = "C:/Users/KellyHeilman/Box/01. kelly.heilman Workspace/mortality/Eastern-Mortality/mortality_models/"
# directory with the predictive draws
pred_draws_dir <- paste0(output.dir, "SPCD_stanoutput_cmdstan/model_avg_draws")    

# read in species numbers
spp.table <- read.csv(file = paste0(output.dir, "/data/Hierarchical_obs_model_7.csv"))
spcd_ids <- spp.table$SPCD
n_spp_total <- length(spcd_ids) 

# names of the main effects
main_effects <- c("DIA_DIFF_scaled", "DIA_scaled", "ba.scaled", "BAL.scaled",
                  "damage.scaled", "MATmax.scaled", "MAP.scaled", "ppt.anom",
                  "tmax.anom", "slope.scaled", "aspect.scaled", "Ndep.scaled")



# ----------------------------------------------------------------------
# Load all species model weights from LOO_ELPD.df (generated in cmdstan_model_assessment) ---
# ----------------------------------------------------------------------
LOO_ELPD.df <- readRDS( paste0(output.dir, "SPCD_stanoutput_cmdstan/summary/All_LOO_comparisons.rds"))

# function to filter the weights for each species
load_species_weights <- function(spcd, 
                                 weighting = "stacking_wts")# can be c("pbma_wts","pbma_BB_wts", "stacking_wts"))
{
LOO_ELPD.df %>% filter(SPCD %in% spcd) %>% select(SPCD, model,  !!sym(weighting))  
}

# ----------------------------------------------------------------------
# Set up indices to draw predictions from all models--
# ----------------------------------------------------------------------
# Because weighting is done by sampling N*w draws from the posterior predictions of pSurvival for each predictor value 
# we need to make sure we are subsampling or weighting the posterior predictions using the same random draw number for each predictor
build_species_draw_plan <- function(spcd, weights_i, N = 4000) {
  
  weights_i[, weighting]
  
  # keep the models that have weight > 0
  weights <- weights_i[weights_i[,3] > 0,]
  model_names <- weights$model
  w <- weights[,3] / sum(weights[,3]) # resum
  
  # check how many draws each model has (they should all be the same, but this is to be save)

  n_draws_per_model <- vapply(model_names, function(wname) {
    
    structure_str <- if (grepl("^Hierarchical_model", wname)) "hierarchical" else "species_only"
    
    # get the model number
    k <- as.integer(gsub("\\D", "", wname))
    f <- file.path(pred_draws_dir, sprintf("SPCD_%s_model_%s_%s.qs", spcd, k, structure))
    if (!file.exists(f)) next
    res <- qs2::qs_read(f)
    
    nrow(res$prob_draws)
  }, integer(1))
  
  # get an index for which model to sample for each of the 4000 posteriors: 
  sampled_model <- sample(model_names, N, replace = TRUE, prob = w)
  # get the index within each original model draws to keep consistent across variables
  sampled_draw  <- vapply(sampled_model, function(m) sample(n_draws_per_model[m], 1), integer(1))
  
  list(model = sampled_model, draw = sampled_draw, N = N)
}

# mix_predictions function applies the plan (from build_species_draw_plan) to get a weighted sample for model stacked or average posterior
mix_predictions <- function(pred_list, plan) { # renamed to mix_predictions from apply_draw_plan
  
  # get the number of x values 
  grid_n <- ncol(pred_list[[1]])
  
  out <- matrix(NA_real_, plan$N, grid_n)
  for (n in seq_len(plan$N)) {
    out[n, ] <- pred_list[[plan$model[n]]][plan$draw[n], ]
  }
  
  # keep a record of which model this draw came from
  list(draws = out, model = plan$model)
}


# -----------------------------------------------------------------------
# MAIN EFFECTS: Get the marginal effect curves using full-Bayesian Model Averaging ---
# for each species x covariation, across all 18 models
# -----------------------------------------------------------------------
spcd = 12
focal_var = "DIA_scaled"
N = 4000 # number of total samples

# function to get the draws for a single focal variable
build_bma_curve <- function(spcd, focal_var, plan, weighting) {
  pred_list <- list()
  grid_val <- NULL
  
  for (structure in c("species_only", "hierarchical")) {
    for (k in 1:9) {
      f <- file.path(pred_draws_dir, sprintf("SPCD_%s_model_%s_%s.qs", spcd, k, structure))
      if (!file.exists(f)) next
      res <- qs2::qs_read(f)
      
      # get the index with the samples for the focal var, this is used to index res$grid_val (the variation in xM for that focal_var), and the draws
      rows <- which(res$focal_var == focal_var)
      
      if (length(rows) == 0) next  # this shouln't happen, but good guardrail
      
      wname <- paste0(if (structure == "hierarchical") "Hierarchical_model_" else "Species_model_", k)
      if (!(wname %in% plan$model)) next  # zero-weight models never get sampled; skip loading
      
      pred_list[[wname]] <- res$prob_draws[, rows, drop = FALSE]
     # w_list[wname] <- weights_i[weights_i$model %in% wname, weighting]
      grid_val <- res$grid_val[rows]  # should be identical across models
      
      
    }
  }
  
  if (length(pred_list) == 0) {
    stop("No saved draws found for SPCD ", spcd, " / ", focal_var,
         " -- did you run RUN_MAIN_LOOP <- TRUE in run_bma_generated_quantities.R?")
  }
  
  mixed <- mix_predictions(pred_list, plan)
 
  list(draws = mixed$draws, 
      model_source = mixed$model, 
      grid_val = grid_val, 
      spcd = spcd, 
      focal_var = focal_var, 
      weighting = weighting)
}

# this runs bma for all the predictor variables for a given species
build_bma_curves_for_species <- function(spcd, 
                                         weighting = "stacking_wts",
                                         focal_vars = main_effects, 
                                         N = 4000) {
  
  
  weights_i <- load_species_weights(spcd, weighting)
  plan <- build_species_draw_plan(spcd, weights_i, N = N)
  
  # apply the curve over all of the species
  curves <- lapply(focal_vars, function(fv) build_bma_curve(spcd, fv, plan, weighting))
  names(curves) <- focal_vars
  curves
}




# Summary function to compile 95% CI of model averages 
# NOTE: bma_curve needs to be the the name of the list element for a single predictor variable
summarize_marginal_draws <- function(bma_curve) {
  
  # add code to convert draws into annualized probability of mortality:
  annual_pMort <- 1 - bma_curve$draws
  annual_pSurv <- bma_curve$draws
  
  decadal_pMort <- 1 - bma_curve$draws^10
  decadal_pSurv <- bma_curve$draws^10
  
  
  # create summaries
  data.frame(predictor = bma_curve$focal_var,
             weighting = bma_curve$weighting,
             grid_val = bma_curve$grid_val,
             
             # get summaries of annual survival probability
             annual_pSurv.median = apply(annual_pSurv, 2, quantile, probs = 0.5),
             annual_pSurv.05.ci.lo = apply(annual_pSurv, 2, quantile, probs = 0.05), 
             annual_pSurv.95.ci.hi = apply(annual_pSurv, 2, quantile, probs = 0.95), 
             annual_pSurv.10.ci.lo = apply(annual_pSurv, 2, quantile, probs = 0.10), 
             annual_pSurv.90.ci.hi = apply(annual_pSurv, 2, quantile, probs = 0.90), 
             annual_pSurv.25.ci.lo = apply(annual_pSurv, 2, quantile, probs = 0.25), 
             annual_pSurv.75.ci.hi = apply(annual_pSurv, 2, quantile, probs = 0.75),
             
             
             # get summaries of annual survival probability
             annual_pMort.median = apply(annual_pMort, 2, quantile, probs = 0.5),
             annual_pMort.05.ci.lo = apply(annual_pMort, 2, quantile, probs = 0.05), 
             annual_pMort.95.ci.hi = apply(annual_pMort, 2, quantile, probs = 0.95),
             
             annual_pMort.10.ci.lo = apply(annual_pMort, 2, quantile, probs = 0.10), 
             annual_pMort.90.ci.hi = apply(annual_pMort, 2, quantile, probs = 0.90), 
             annual_pMort.25.ci.lo = apply(annual_pMort, 2, quantile, probs = 0.25), 
             annual_pMort.75.ci.hi = apply(annual_pMort, 2, quantile, probs = 0.75),             
             
             # get summaries of decadal survival probability
             decadal_pSurv.median = apply(decadal_pSurv, 2, quantile, probs = 0.5),
             decadal_pSurv.05.ci.lo = apply(decadal_pSurv, 2, quantile, probs = 0.05), 
             decadal_pSurv.95.ci.hi = apply(decadal_pSurv, 2, quantile, probs = 0.95),
             
             decadal_pSurv.10.ci.lo = apply(decadal_pSurv, 2, quantile, probs = 0.10), 
             decadal_pSurv.90.ci.hi = apply(decadal_pSurv, 2, quantile, probs = 0.90), 
             decadal_pSurv.25.ci.lo = apply(decadal_pSurv, 2, quantile, probs = 0.25), 
             decadal_pSurv.75.ci.hi = apply(decadal_pSurv, 2, quantile, probs = 0.75),
             
             # get summaries of decadal survival probability
             decadal_pMort.median = apply(decadal_pMort, 2, quantile, probs = 0.5),
             decadal_pMort.05.ci.lo = apply(decadal_pMort, 2, quantile, probs = 0.05), 
             decadal_pMort.95.ci.hi = apply(decadal_pMort, 2, quantile, probs = 0.95),
             decadal_pMort.10.ci.lo = apply(decadal_pMort, 2, quantile, probs = 0.10), 
             decadal_pMort.90.ci.hi = apply(decadal_pMort, 2, quantile, probs = 0.90), 
             decadal_pMort.25.ci.lo = apply(decadal_pMort, 2, quantile, probs = 0.25), 
             decadal_pMort.75.ci.hi = apply(decadal_pMort, 2, quantile, probs = 0.75)         
            )
}

# function to all draws in long format  
marginal_draws_long <- function(bma_curve) {
  
  # convert p_survival to long format
  psurv.long <- data.frame(draw_source = bma_curve$model_source,
                           bma_curve$draws)%>% 
    rownames_to_column(var = "draw")%>%
    pivot_longer(cols = starts_with("X"), 
                 names_to = "val_name", 
                 values_to = "annual_pSurv")
  
  # add species cod, the focal variable, the weighting, and join to pSurv_long
  data.frame(spcd = bma_curve$spcd,
             predictor = bma_curve$focal_var,
             weighting = bma_curve$weighting,
             grid_val = bma_curve$grid_val,
             val_name = colnames(data.frame(bma_curve$draws))) %>%
    left_join(., psurv.long, by = "val_name") %>% 
    
    # calculate additional probabilities
    mutate(annual_pMort = 1 - annual_pSurv, 
           decadal_pSurv = annual_pSurv^10, 
           decadal_pMort = 1 - annual_pSurv^10)
    
    
  
}

# Example output is a nested list of the species, with the predictor variables as the names
curve_SPCD12 <- build_bma_curves_for_species(spcd = 12,
                                      weighting = "stacking_wts", # can be c("pbma_wts","pbma_BB_wts", "stacking_wts")
                                      focal_vars = main_effects, # the name of the model predictor 
                                      N = 4000)

# get a summary for DIA_diff only
summary_df <- summarize_marginal_draws(bma_curve = curve_SPCD12$DIA_DIFF_scaled)


ggplot(summary_df)+geom_line(aes(x = grid_val, y = annual_pMort.median))+
  geom_ribbon(aes(x = grid_val, ymin =annual_pMort.05.ci.lo, ymax = annual_pMort.95.ci.hi), alpha = 0.25, fill = "red")+
  geom_ribbon(aes(x = grid_val, ymin =annual_pMort.25.ci.lo, ymax = annual_pMort.75.ci.hi), alpha = 0.65, fill = "red")



## get draws for DIA_diff only we can peak at which models are crontibuting the most
marg_draws_df <- marginal_draws_long(bma_curve = curve_SPCD12$DIA_DIFF_scaled)

ggplot(marg_draws_df)+geom_line(aes(x = grid_val, y = decadal_pMort, group = draw, color = draw_source), alpha = 0.05)+
  theme_bw(base_size = 12)+
  xlab(paste0(unique(marg_draws_df$predictor)))+facet_wrap(~draw_source)

  
# 'output_BMA_marginal_effects' function to get draws or draw summaries this to all species and main effects--------------
output_BMA_marginal_effects <- function(spcd_idx,  
         weighting = "stacking_wts", 
         main_effects = main_effects, 
         N = 4000, 
         return_format = "draws"){ # return_format can be "draws" or "summary"
  
  cat(paste("using", weighting, "for BMA on marginal effects for", spcd_idx, "\n"))
    # run all curves for a single species
    curve_SPCD <- build_bma_curves_for_species(spcd = spcd_idx,
                                                 weighting = "stacking_wts", # can be c("pbma_wts","pbma_BB_wts", "stacking_wts")
                                                 focal_vars = main_effects, # the name of the model predictor 
                                                 N = 4000)
    
    # use lapply to get a summary for all effects for this spcd_idx
    main_marg_summary <- do.call(rbind, lapply(main_effects, function(fv_idx){
      summarize_marginal_draws(bma_curve = curve_SPCD[[fv_idx]])}))%>%
      mutate(SPCD = spcd_idx) %>% select(SPCD, everything())
    main_marg_summary$COMMON_NAME <- ref_species[match(main_marg_summary$SPCD, ref_species$SPCD),]$COMMON_NAME
    
    # ggplot(data = main_marg_summary)+
    #   geom_ribbon(aes(x = grid_val, ymin = decadal_pMort.05.ci.lo, ymax = decadal_pMort.95.ci.hi, fill = COMMON_NAME), alpha = 0.5)+
    #   geom_line(aes(x = grid_val, y = decadal_pMort.median, color = COMMON_NAME))+
    #   facet_wrap(~predictor)
    
    
    ## get draws output for all of the main effects for this species
    main_marg_draws <- do.call(rbind, lapply(main_effects, function(fv_idx){
      marginal_draws_long(bma_curve = curve_SPCD[[fv_idx]])}))%>%
      mutate(SPCD = spcd_idx) %>% select(SPCD, everything())
    main_marg_draws$COMMON_NAME <- ref_species[match(main_marg_draws$SPCD, ref_species$SPCD),]$COMMON_NAME
    
    # # example plot of all draws colored by the model source--takes awhile due to object size!
    # ggplot(main_marg_draws)+
    #   geom_line(aes(x = grid_val, y = decadal_pMort, group = draw, color = draw_source), alpha = 0.05)+
    #   theme_bw(base_size = 12)+
    #   xlab(paste0(unique(marg_draws_df$predictor)))+
    #   facet_wrap(~predictor)+
    #   guides(color = guide_legend(override.aes = list(linewidth = 3, alpha = 1)))
    # save the objects:
    # if the directory doesnt exists, create it
    if (!dir.exists(paste0(output.dir, "SPCD_stanoutput_cmdstan/model_avg_draws/",weighting, "/"))) {
      dir.create(paste0(output.dir, "SPCD_stanoutput_cmdstan/model_avg_draws/",weighting, "/"), recursive = TRUE)
    }
    
    qs2::qs_save(main_marg_draws, paste0(output.dir, "SPCD_stanoutput_cmdstan/model_avg_draws/",weighting, "/Draws_marginal_main_", spcd_idx, "_", weighting, ".qs"))
    qs2::qs_save(main_marg_summary, paste0(output.dir, "SPCD_stanoutput_cmdstan/model_avg_draws/",weighting, "/Summary_marginal_main_", spcd_idx, "_", weighting, ".qs"))

if(return_format %in% "draws"){
  return(main_marg_draws)
} else{
  return(main_marg_summary)
}   
    
}


# apply function over all the species and get the summaries
# get summaries as the outputs...
marg_summary_all_list <- lapply(spcd_ids, FUN = function(species_code){
  output_BMA_marginal_effects(spcd_idx = species_code, 
                              weighting = "stacking_wts", 
                              main_effects = main_effects, 
                              N = 4000, 
                              return_format = "summary")
})

marg_summary_all_df <- do.call(rbind, marg_summary_all_list)
marg_summary_all_df$predictor <- factor(marg_summary_all_df$predictor, levels = main_effects)

# plot up model average marginal effects of 10-year mortality probabilities

# set the species order using the factors:
SP.TRAITS <- read.csv("data/NinemetsSpeciesTraits.csv") %>% filter(COMMON_NAME %in% unique(spp.table$COMMON))
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


# plot up the 50% CI only across the values:
ggplot(data = marg_summary_all_df)+
  #geom_ribbon(aes(x = grid_val, ymin = decadal_pMort.05.ci.lo, ymax = decadal_pMort.95.ci.hi, fill = COMMON_NAME), alpha = 0.25)+
  geom_ribbon(aes(x = grid_val, ymin = decadal_pMort.25.ci.lo, ymax = decadal_pMort.75.ci.hi, fill = COMMON_NAME), alpha = 0.5)+
  geom_line(aes(x = grid_val, y = decadal_pMort.median, color = COMMON_NAME))+
  facet_wrap(~predictor, scales = "free_y")+
  species_fill + species_color

# omit balsam fir
ggplot(data = marg_summary_all_df %>% filter(! COMMON_NAME %in% "balsam fir"))+
  geom_ribbon(aes(x = grid_val, ymin = decadal_pMort.05.ci.lo, ymax = decadal_pMort.95.ci.hi, fill = COMMON_NAME), alpha = 0.25)+
  geom_ribbon(aes(x = grid_val, ymin = decadal_pMort.25.ci.lo, ymax = decadal_pMort.75.ci.hi, fill = COMMON_NAME), alpha = 0.75)+
  geom_line(aes(x = grid_val, y = decadal_pMort.median, color = COMMON_NAME))+
  facet_wrap(~predictor, scales = "free_y")+
  species_fill + species_color

# just plot medians 
ggplot(data = marg_summary_all_df)+
 # geom_ribbon(aes(x = grid_val, ymin = decadal_pMort.05.ci.lo, ymax = decadal_pMort.95.ci.hi, fill = COMMON_NAME), alpha = 0.25)+
  #geom_ribbon(aes(x = grid_val, ymin = decadal_pMort.25.ci.lo, ymax = decadal_pMort.75.ci.hi, fill = COMMON_NAME), alpha = 0.75)+
  geom_line(aes(x = grid_val, y = decadal_pMort.median, color = COMMON_NAME), size = 1.5)+
  facet_wrap( vars(predictor), scales = "free_y")+
  theme_bw()+
  theme(panel.grid = element_blank(), strip.text.y = element_text(angle = 0))+
  #species_fill + 
  species_color+
  ylab("Predicted 10-year mortality probabilities (model stacking)")

# flip rows and columns
ggplot(data = marg_summary_all_df)+
  geom_ribbon(aes(x = grid_val, ymin = decadal_pMort.05.ci.lo, ymax = decadal_pMort.95.ci.hi, fill = COMMON_NAME), alpha = 0.25)+
  geom_ribbon(aes(x = grid_val, ymin = decadal_pMort.25.ci.lo, ymax = decadal_pMort.75.ci.hi, fill = COMMON_NAME), alpha = 0.75)+
  geom_line(aes(x = grid_val, y = decadal_pMort.median, color = COMMON_NAME))+
  facet_grid(vars(predictor), vars(COMMON_NAME),  scales = "free_y")+
  theme_bw()+
  theme(panel.grid = element_blank(), strip.text.y = element_text(angle = 0))+
  species_fill + species_color+
  ylab("Predicted 10-year mortality probabilities (model stacking)")


# the three northern conifers have large uncertainties--
ggplot(data = marg_summary_all_df %>% filter(! COMMON_NAME %in% c("balsam fir", "red spruce", "northern white-cedar")))+
geom_ribbon(aes(x = grid_val, ymin = decadal_pMort.05.ci.lo, ymax = decadal_pMort.95.ci.hi, fill = COMMON_NAME), alpha = 0.25)+
  geom_ribbon(aes(x = grid_val, ymin = decadal_pMort.25.ci.lo, ymax = decadal_pMort.75.ci.hi, fill = COMMON_NAME), alpha = 0.75)+
  geom_line(aes(x = grid_val, y = decadal_pMort.median, color = COMMON_NAME))+
  facet_grid(vars(predictor), vars(COMMON_NAME), scales = "free_y")+
  theme_bw()+
  theme(panel.grid = element_blank(), strip.text.y = element_text(angle = 0), 
        strip.text.x = element_text(angle = 90))+
  species_fill + species_color+
  ylab("Predicted 10-year mortality probabilities (model stacking)")+
  xlab("Scaled Covariate")

ggplot(data = marg_summary_all_df %>% filter(COMMON_NAME %in% c("balsam fir", "red spruce", "northern white-cedar")))+
  geom_ribbon(aes(x = grid_val, ymin = decadal_pMort.05.ci.lo, ymax = decadal_pMort.95.ci.hi, fill = COMMON_NAME), alpha = 0.25)+
  geom_ribbon(aes(x = grid_val, ymin = decadal_pMort.25.ci.lo, ymax = decadal_pMort.75.ci.hi, fill = COMMON_NAME), alpha = 0.75)+
  geom_line(aes(x = grid_val, y = decadal_pMort.median, color = COMMON_NAME))+
  facet_grid(vars(predictor), vars(COMMON_NAME))+
  theme_bw()+
  theme(panel.grid = element_blank(), strip.text.y = element_text(angle = 0),strip.text.x = element_text(angle = 90))+
  species_fill + species_color+
  ylab("Predicted 10-year mortality probabilities (model stacking)")

ggplot(data = marg_summary_all_df %>% filter(COMMON_NAME %in% c("northern red oak", "chestnut oak", "white oak")))+
  geom_ribbon(aes(x = grid_val, ymin = decadal_pMort.05.ci.lo, ymax = decadal_pMort.95.ci.hi, fill = COMMON_NAME), alpha = 0.25)+
  geom_ribbon(aes(x = grid_val, ymin = decadal_pMort.25.ci.lo, ymax = decadal_pMort.75.ci.hi, fill = COMMON_NAME), alpha = 0.75)+
  geom_line(aes(x = grid_val, y = decadal_pMort.median, color = COMMON_NAME))+
  facet_grid(vars(predictor), vars(COMMON_NAME), scales = "free_y")+
  theme_bw()+
  theme(panel.grid = element_blank(), strip.text.y = element_text(angle = 0),strip.text.x = element_text(angle = 90))+
  species_fill + species_color+
  ylab("Predicted 10-year mortality probabilities (model stacking)")

ggplot(data = marg_summary_all_df %>% filter(COMMON_NAME %in% c("red maple","sugar maple", "American beech")))+
  geom_ribbon(aes(x = grid_val, ymin = decadal_pMort.05.ci.lo, ymax = decadal_pMort.95.ci.hi, fill = COMMON_NAME), alpha = 0.25)+
  geom_ribbon(aes(x = grid_val, ymin = decadal_pMort.25.ci.lo, ymax = decadal_pMort.75.ci.hi, fill = COMMON_NAME), alpha = 0.75)+
  geom_line(aes(x = grid_val, y = decadal_pMort.median, color = COMMON_NAME))+
  facet_grid(vars(predictor), vars(COMMON_NAME), scales = "free_y")+
  theme_bw()+
  theme(panel.grid = element_blank(), strip.text.y = element_text(angle = 0),strip.text.x = element_text(angle = 90))+
  species_fill + species_color+
  ylab("Predicted 10-year mortality probabilities (model stacking)")

ggplot(data = marg_summary_all_df %>% filter(COMMON_NAME %in% c("hickory spp.","white ash", "black cherry")))+
  geom_ribbon(aes(x = grid_val, ymin = decadal_pMort.05.ci.lo, ymax = decadal_pMort.95.ci.hi, fill = COMMON_NAME), alpha = 0.25)+
  geom_ribbon(aes(x = grid_val, ymin = decadal_pMort.25.ci.lo, ymax = decadal_pMort.75.ci.hi, fill = COMMON_NAME), alpha = 0.75)+
  geom_line(aes(x = grid_val, y = decadal_pMort.median, color = COMMON_NAME))+
  facet_grid(vars(predictor), vars(COMMON_NAME))+
  theme_bw()+
  theme(panel.grid = element_blank(), strip.text.y = element_text(angle = 0),strip.text.x = element_text(angle = 90))+
  species_fill + species_color+
  ylab("Predicted 10-year mortality probabilities (model stacking)")

ggplot(data = marg_summary_all_df %>% filter(COMMON_NAME %in% c("eastern hemlock","American beech", "yellow-poplar")))+
  geom_ribbon(aes(x = grid_val, ymin = decadal_pMort.05.ci.lo, ymax = decadal_pMort.95.ci.hi, fill = COMMON_NAME), alpha = 0.25)+
  geom_ribbon(aes(x = grid_val, ymin = decadal_pMort.25.ci.lo, ymax = decadal_pMort.75.ci.hi, fill = COMMON_NAME), alpha = 0.75)+
  geom_line(aes(x = grid_val, y = decadal_pMort.median, color = COMMON_NAME))+
  facet_grid(vars(predictor), vars(COMMON_NAME))+
  theme_bw()+
  theme(panel.grid = element_blank(), strip.text.y = element_text(angle = 0),strip.text.x = element_text(angle = 90))+
  species_fill + species_color+
  ylab("Predicted 10-year mortality probabilities (model stacking)")


# do the same thing but get a list of draws (this will be large! ~ 2GB)
marg_draws_all_list <- lapply(spcd_ids, FUN = function(species_code){
  output_BMA_marginal_effects(spcd_idx = species_code, 
                              weighting = "stacking_wts", 
                              main_effects = main_effects, 
                              N = 4000, 
                              return_format = "draws")
})

# expanding into a very big dataframe
marg_draws_all_df <- do.call(rbind, marg_draws_all_list)

# arrange the species--
marg_draws_all_df$COMMON_NAME <- factor(marg_draws_all_df$COMMON_NAME, levels = c(SP.TRAITS$COMMON_NAME))

# make draws summary plots for different effects to identify which model s
marg_draws_all_df %>% filter(predictor %in% c("DIA_DIFF_scaled", "DIA_scaled"))|>
ggplot()+
  geom_line(aes(x = grid_val, y = decadal_pMort, group = draw, color = draw_source), alpha = 0.05)+
  theme_bw(base_size = 12)+
  xlab("Scaled Predictor")+
  ylab("10-year predicted mortality probability")+
  facet_grid(vars(predictor), vars(COMMON_NAME), scales = "free_y")+
  theme(panel.grid = element_blank(),
        strip.text.y = element_text(angle = 0), 
        strip.text.x = element_text(angle = 90), 
        legend.position = "bottom", 
        legend.direction = "horizontal")+
  guides( color = guide_legend(override.aes = list(linewidth = 2, alpha = 1)))+
    labs(color = "Weighted Draw Source")

ggsave(paste0(output.dir, "SPCD_stanoutput_cmdstan/summary/Size_effects_pMort10_draws.png"), 
       height = 6, width = 16)

# do this for the plot neighborhood effects:
c("DIA_DIFF_scaled", "DIA_scaled", "ba.scaled", "BAL.scaled",
  "damage.scaled", "MATmax.scaled", "MAP.scaled", "ppt.anom",
  "tmax.anom", "slope.scaled", "aspect.scaled", "Ndep.scaled")

marg_draws_all_df %>% filter(predictor %in% c("ba.scaled", "BAL.scaled", "damage.scaled"))|>
  ggplot()+
  geom_line(aes(x = grid_val, y = decadal_pMort, group = draw, color = draw_source), alpha = 0.05)+
  theme_bw(base_size = 12)+
  xlab("Scaled Predictor")+
  ylab("10-year predicted mortality probability")+
  facet_grid(vars(predictor), vars(COMMON_NAME), scales = "free_y")+
  theme(panel.grid = element_blank(),
        strip.text.y = element_text(angle = 0), 
        strip.text.x = element_text(angle = 90), 
        legend.position = "bottom", 
        legend.direction = "horizontal")+
  guides( color = guide_legend(override.aes = list(linewidth = 2, alpha = 1)))+
  labs(color = "Weighted Draw Source")

ggsave(paste0(output.dir, "SPCD_stanoutput_cmdstan/summary/Neighborhood_effects_pMort10_draws.png"), 
       height = 8, width = 16)

# for the climate effects
marg_draws_all_df %>% filter(predictor %in% c("MATmax.scaled", "MAP.scaled", 
                                              "ppt.anom",
                                              "tmax.anom"))|>
  ggplot()+
  geom_line(aes(x = grid_val, y = decadal_pMort, group = draw, color = draw_source), alpha = 0.05)+
  theme_bw(base_size = 12)+
  xlab("Scaled Predictor")+
  ylab("10-year predicted mortality probability")+
  facet_grid(vars(predictor), vars(COMMON_NAME), scales = "free_y")+
  theme(panel.grid = element_blank(),
        strip.text.y = element_text(angle = 0), 
        strip.text.x = element_text(angle = 90), 
        legend.position = "bottom", 
        legend.direction = "horizontal")+
  guides( color = guide_legend(override.aes = list(linewidth = 2, alpha = 1)))+
  labs(color = "Weighted Draw Source")

ggsave(paste0(output.dir, "SPCD_stanoutput_cmdstan/summary/Climate_effects_pMort10_draws.png"), 
       height = 9, width = 16)



# for the site condition effects
marg_draws_all_df %>% filter(predictor %in% c("Ndep.scaled","slope.scaled", "aspect.scaled"))|>
  ggplot()+
  geom_line(aes(x = grid_val, y = decadal_pMort, group = draw, color = draw_source), alpha = 0.05)+
  theme_bw(base_size = 12)+
  xlab("Scaled Predictor")+
  ylab("10-year predicted mortality probability")+
  facet_grid(vars(predictor), vars(COMMON_NAME), scales = "free_y")+
  theme(panel.grid = element_blank(),
        strip.text.y = element_text(angle = 0), 
        strip.text.x = element_text(angle = 90), 
        legend.position = "bottom", 
        legend.direction = "horizontal")+
  guides( color = guide_legend(override.aes = list(linewidth = 2, alpha = 1)))+
  labs(color = "Weighted Draw Source")

ggsave(paste0(output.dir, "SPCD_stanoutput_cmdstan/summary/Site_condition_effects_pMort10_draws.png"), 
       height = 8, width = 16)

rm(marg_draws_all_df, decadal_pMort, decadal_pSurv)
rm(draws, curves_316, curves, hier_results, bma_curve, curve_Ndep_SPCD316, curve_SPCD, curve_SPCD12)
rm(psurv.long)
rm(pred_list, res, res_sp)
rm(main_marg_summary, main_marg_draws, marg_draws_list_all, marg_draws_df, marge_summary_all_list)
gc()

#################################################################################
# TODO: Placeholder for the interaction effects-------
#################################################################################




#################################################################################
# Bayesian model averaging (stacking) for yhat, yrep, mhat, mrep, betas, and alphas -------
#################################################################################

# we can use similar appraoch to marginal effect averaging above
post_draws_dir <- paste0(output.dir, "SPCD_stanoutput_cmdstan/predicted_mort/")

# set up functions to get the species paths
# NOTE: species and hierarchical models had inconsistent naming of the generated quantities, 
# created a lookup table to translate these names, so when we request "y_rep" in the following functions, 
# we always get the out of sample survival assignment

# y_rep and y_hat naming are flipped for the species_only
# species: "y_rep" = in sample and "y_hat" = out of sample
# hierarchical: "y_hat" = in sample and "y_rep" = out of sample

y_type_lookup <- data.frame(
  actual_y_type = c("y_rep", "y_hat", "pSurv_rep", "pSurv_hat"),
  heir_y_type = c("y_rep", "y_hat", "pSurv_rep", "pSurv_hat"),
  
  SPP_y_type = c("y_hat", "y_rep", "pSurv_hat", "pSurv_rep"),# the mixed up labels for species model
  data_type = c("out-of-sample", "in-sample", "out-of-sample", "in-sample"),
  parameter = c("Survival Assignment", "Survival Assignment", "annual survival probability", "annual survival probability")
)


# function to read the species paths
species_path <- function(y_type, spcd, k) {
  # change the naming of the species files to read in the correct file:
  y_type_real <- y_type_lookup[match(y_type, y_type_lookup$actual_y_type),]$SPP_y_type
  
  # get the filename
  paste0(post_draws_dir, y_type_real, "_samps_mort_model_", k, "_SPCD_", spcd, "_remper_correction_0.5_niter_1000_nchain_4.qs")
}

# function to read the hierarchical paths
hier_path     <- function(y_type, spcd, k) {
  # just to be consistent with species paths
  y_type_real <- y_type_lookup[match(y_type, y_type_lookup$actual_y_type),]$heir_y_type
  paste0(post_draws_dir, y_type_real, "_samps_SPCD_", spcd, "_hierarchical_mort_model_",k,"_niter_1000_nchain_4.qs")
}


# get posterior matrix actually reads and outputs these paths
get_posterior_matrix <- function(structure,  spcd, y_type, k, spp_index = NULL) {
 
   if (structure == "species_only") {
    qs2::qs_read(species_path( y_type, spcd, k))          # draws x N_species, already per-species
  } else {
    full <- qs2::qs_read(hier_path(y_type, spcd, k))            # draws x N_total (all species stacked)
    full
  }
}


get_posterior_matrix(structure = "species_only",  
                     spcd = "261", 
                     y_type = "y_hat",
                     k = 1, 
                     spp_index = NULL) %>% ncol()

get_posterior_matrix(structure = "hierarchical",  
                     spcd = "261", 
                     y_type = "y_hat",
                     k = 1, 
                     spp_index = NULL) %>% ncol()

post_draw_plan <- function(spcd, weights_i, N = 4000) {
  
  weights_i[, weighting]
  
  # keep the models that have weight > 0
  weights <- weights_i[weights_i[,3] > 0,]
  model_names <- weights$model
  w <- weights[,3] / sum(weights[,3]) # resum
  
  

  # get an index for which model to sample for each of the 4000 posteriors: 
  sampled_model <- sample(model_names, N, replace = TRUE, prob = w)
  # get the index within each original model draws to keep consistent across variables
  # assuming that N = the number of draws across all the models
  sampled_draw  <- vapply(sampled_model, function(m) sample(N, 1), integer(1))
  
  list(model = sampled_model, draw = sampled_draw, N = N)
}




# this is giving an error at stopifnot(length(unique(n_obs)) == 1)  # same set of trees across all model structures
# I think this is because y_rep and y_hat naming are flipped for hierarchical vs species models
stack_y_species <- function(spcd, 
                        y_type = c("y_rep"),
                        N = 4000, 
                        weighting = "stacking_wts",
                        y_plan = NULL) {
 
  # get the weighting and indices to sample from for each species
  weights_i <- load_species_weights(spcd, weighting)
  plan <- if (is.null(y_plan)) post_draw_plan(spcd, weights_i, N) else y_plan
  
  # get the full predicted list from all the models that have non-zero weights for this species
  pred_list <- list()
  for (structure in c("species_only", "hierarchical")) {
    for (k in 1:9) {
      # wname = naming structure of the model in the weight dataframe
      wname <- paste0(if (structure == "hierarchical") "Hierarchical_model_" else "Species_model_", k)
      if (!(wname %in% plan$model)) next  # zero-weight models were never sampled; skip loading
      
      # get the matrix of y values
      mat <- tryCatch(get_posterior_matrix(structure, spcd, y_type, k), error = function(e) NULL)
      if (is.null(mat)) next
      
      # get the max number of iterations needed from this model and check that we have that
      needed <- max(plan$draw[plan$model == wname])
      if (nrow(mat) < needed) {
        stop(wname, "'s saved ", type, " has only ", nrow(mat), " draws, but the plan ",
             "needs draw index ", needed, ".")
      }
      pred_list[[wname]] <- mat # 
    }
  }
  
  # check the number of observations for this focal prediction parameter
  n_obs <- vapply(pred_list, ncol, integer(1))
  stopifnot(length(unique(n_obs)) == 1)  # throw error if we dont have the same number trees across all model structures
  
  
  # use indices in the plan to output the predicted y or psurv values
  mixed <- mix_predictions(pred_list, plan)  # draws x n_obs matrix of BMA-mixed simulated y (0/1)
  
  # get the number of tree observations for the species
  tree_N <- unique(n_obs)
  var_names <- paste0(y_type, "_",1:tree_N)
  draws_df <- as.data.frame(mixed$draws)
  colnames(draws_df) <- var_names  
  draws_matrix <- as_draws_matrix(draws_df)
  
  list(draws = draws_matrix, 
       model_source = mixed$model, 
       nObs = tree_N, 
       spcd = spcd, 
       y_type = y_type, 
       weighting = weighting)
}

# example usage for pSurv_hat for one species
psurv_261 <- stack_y_species(spcd = 261, 
                y_type = "pSurv_hat", 
                N = 4000, 
                y_plan = NULL)

#psurv_261$draws %>% as_draws_df() %>% summarise_draws()


# create a function to process/summarise the stacked posteriors
# 1. Save the weighted samples pSurv_hat, pSurv_rep, y_hat, y_rep for each species as draws_df for later use
# 2. Calculate is and oos AUC scores, confusion matrices for each species
# 3. Summarise the pSurv_hat, pSurv_rep, y_hat, y_rep for each species


spcd = 261
weighting = "stacking_wts"
N = 4000
process_bma_posterior_pSurv_y <- function(spcd, weighting, N){
      # start function
      
    cat(paste0("getting weighted posteriors for spcd ", spcd, "\n"))
  # get species in and out of sample y data to compare to:
      load(paste0("SPCD_standata_general_full_standardized_v3/", "SPCD_", spcd, "remper_correction_0.5model_9.Rdata"))
      #length(mod.data$y)
      
      
      # set up plan for each species once, then feed it inot the stack_y_species
      weights_i <- load_species_weights(spcd, weighting)
      plan <- post_draw_plan(spcd, weights_i, N)
      
      # get the weighted probabilities and predicted classifications
      psurv_hat <- stack_y_species(spcd = spcd, 
                                   y_type = "pSurv_hat", 
                                   N = 4000, 
                                   y_plan = plan)
      psurv_rep <- stack_y_species(spcd = spcd, 
                                   y_type = "pSurv_rep", 
                                   N = 4000, 
                                   y_plan = plan)
      
      y_hat <- stack_y_species(spcd = spcd, 
                               y_type = "y_hat", 
                               N = 4000, 
                               y_plan = plan)
      
      y_rep <- stack_y_species(spcd = spcd, 
                                   y_type = "y_rep", 
                                   N = 4000, 
                                   y_plan = plan)
      
      
      cat(paste0("saving weighted posterior draws for ", spcd, "\n"))
      # get just the draws
      y_rep_samps <- y_rep$draws
      y_hat_samps <- y_hat$draws
      
      pSurv_rep_samps <- psurv_rep$draws
      pSurv_hat_samps <- psurv_hat$draws
      
      qs2::qs_save(y_rep_samps,     paste0(output.dir, "SPCD_stanoutput_cmdstan/predicted_mort/y_rep_samps_SPCD_", spcd, "_", weighting, ".qs"))
      qs2::qs_save(y_hat_samps,     paste0(output.dir, "SPCD_stanoutput_cmdstan/predicted_mort/y_hat_samps_SPCD_", spcd, "_", weighting, ".qs"))
      qs2::qs_save(pSurv_rep_samps, paste0(output.dir, "SPCD_stanoutput_cmdstan/predicted_mort/pSurv_rep_samps_SPCD_", spcd, "_", weighting, ".qs"))
      qs2::qs_save(pSurv_hat_samps, paste0(output.dir, "SPCD_stanoutput_cmdstan/predicted_mort/pSurv_hat_samps_SPCD_", spcd, "_", weighting, ".qs"))
      
      # convert remper probabilities to annualized probabilities of survival:
      Remper_matrix <- matrix(mod.data$Remper, nrow = nrow(pSurv_hat_samps), ncol = ncol(pSurv_hat_samps), byrow = TRUE)
      #length(unique(Remper_matrix[,1])) ==1 # check that the byrow is right
      
      pSannual_hat <- pSurv_hat_samps^(1/Remper_matrix) 
      
      Remperoos_matrix <- matrix(mod.data$Remperoos, nrow = nrow(pSurv_rep_samps), ncol = ncol(pSurv_rep_samps), byrow = TRUE)
      #length(unique(Remper_matrix[,1])) ==1 # check that the byrow is right
      
      pSannual_rep <- pSurv_rep_samps^(1/Remperoos_matrix) 
      
      qs2::qs_save(pSannual_rep, paste0(output.dir, "SPCD_stanoutput_cmdstan/predicted_mort/pSurv_annual_rep_samps_SPCD_", spcd, "_", weighting, ".qs"))
      qs2::qs_save(pSannual_hat, paste0(output.dir, "SPCD_stanoutput_cmdstan/predicted_mort/pSurv_annual_hat_samps_SPCD_", spcd, "_", weighting, ".qs"))
      
  # does the predicted count of y_hat match that in the data?
      is.ppc.bar <- ppc_bars(y = mod.data$y, as_draws_matrix(y_hat_samps))+xlab("Survival Classifcation (In-sample)")
      oos.ppc.bar <- ppc_bars(y = mod.data$ytest, as_draws_matrix(y_rep_samps))+xlab("Survival Classifcation (out-of-sample)")
      combined.bar <- cowplot::plot_grid( is.ppc.bar, oos.ppc.bar, align = "hv")
      cowplot::save_plot(paste0(output.dir, "SPCD_stanoutput_cmdstan/images/post_pred_y_SPCD_", spcd, "_", weighting, ".png"), 
                         plot = combined.bar, 
                         bg = "white")
      rm(is.ppc.bar, oos.ppc.bar, combined.bar)
      
      
      # --- AUC + confusion matrix over posterior draws ------------------------
      actuals     <- mod.data$y
      actuals.oos <- mod.data$ytest
      cat(paste0("estimating AUC and confusion matrix on weighted posterior draws for ", spcd, "\n"))
      
     # library(pROC)
      

      # calculate AUC for each draw:------
      AUC.is.samples.df  <- apply(pSurv_hat_samps, 1, function(p) as.numeric(pROC::auc(actuals, p, quiet = TRUE)))
      AUC.oos.samples.df <- apply(pSurv_rep_samps, 1, function(p) as.numeric(pROC::auc(actuals.oos, p, quiet = TRUE)))
      #rm(pSurv_hat_samps, pSurv_rep_samps, y_hat.quant, y_rep.quant, pSurv_hat.quant, pSurv_rep.quant); gc()
      
      preds.is.class <- y_hat_samps == 1
      confusion.is_draws <- data.frame(
        TP_draws = rowSums(preds.is.class[, actuals == 1, drop = FALSE]),
        FP_draws = rowSums(preds.is.class[, actuals == 0, drop = FALSE]),
        TN_draws = rowSums(!preds.is.class[, actuals == 0, drop = FALSE]),
        FN_draws = rowSums(!preds.is.class[, actuals == 1, drop = FALSE])
      ) %>%
        mutate(`True survival rate` = TP_draws / (TP_draws + FN_draws),
               `True mortality rate` = TN_draws / (TN_draws + FP_draws),
               model.number = "BMA", 
               type = "in-sample",
               model.type = weighting, 
               SPCD = spcd)
      
      preds.oos.class <- y_rep_samps == 1
      confusion.oos_draws <- data.frame(
        TP_draws = rowSums(preds.oos.class[, actuals.oos == 1, drop = FALSE]),
        FP_draws = rowSums(preds.oos.class[, actuals.oos == 0, drop = FALSE]),
        TN_draws = rowSums(!preds.oos.class[, actuals.oos == 0, drop = FALSE]),
        FN_draws = rowSums(!preds.oos.class[, actuals.oos == 1, drop = FALSE])
      ) %>%
        mutate(`True survival rate` = TP_draws / (TP_draws + FN_draws),
               `True mortality rate` = TN_draws / (TN_draws + FP_draws),
               model.number = "BMA", 
               type = "out-of-sample",
               model.type = "weighting", 
               SPCD = spcd)
      
      confusion.is_draws$AUC  <- AUC.is.samples.df
      confusion.oos_draws$AUC <- AUC.oos.samples.df
      AUC.confusion_draws <- rbind(confusion.is_draws, confusion.oos_draws)
      
      qs2::qs_save(AUC.confusion_draws,
                   paste0(output.dir, "SPCD_stanoutput_cmdstan/AUC/AUC_draws_SPCD_", spcd, "_", weighting, ".qs"))
}

# run for all of the species
process_bma_posterior_pSurv_y(spcd = 261,
                              weighting = "stacking_wts",
                              N = 4000)

for(s in 17:1){
  
  process_bma_posterior_pSurv_y(spcd = spp.table$SPCD[s],
                                weighting = "stacking_wts",
                                N = 4000)
}

#---------------------------------------------------------------------------------------------------------
# State, County, Overall scaling of weighted posterior probability of mortality to get mortality rates---
SPCD.id <- 97
weighting =  "stacking_wts"
output.dir = "C:/Users/KellyHeilman/Box/01. kelly.heilman Workspace/mortality/Eastern-Mortality/mortality_models/"
model.type = "Model Average"

# get plot-level expansion factors
PLOT <- read_delim(paste0(output.dir,"data/formatted_older_matching_plts_PLOT.txt"))

plot_expansions <- PLOT %>% select(PLOT.ID, state, county, pltnum, cndtn, cycle, expacr, expvol) %>% distinct()
plot(plot_expansions$expacr, plot_expansions$expvol)


# Function: calculate_state_county_rates_BMA----
# very similar to `calculate_state_county_rates()` in R/modelAssessment/cmdstan_model_assessment.R, but renames weighted files
calculate_state_county_rates_BMA <- function(SPCD.id, weighting, model.type){
  
  cat("\n",paste0("Scaling county and state mortality rates from posteriors: ", SPCD.id,", ",model.type,": ",weighting))
  
        # read in pSurv weighted posteriors
        pSurv_rep_samps <- qs2::qs_read(paste0(output.dir, "SPCD_stanoutput_cmdstan/predicted_mort/pSurv_rep_samps_SPCD_", SPCD.id, "_", weighting, ".qs"))
        pSurv_hat_samps<- qs2::qs_read(paste0(output.dir, "SPCD_stanoutput_cmdstan/predicted_mort/pSurv_hat_samps_SPCD_", SPCD.id, "_", weighting, ".qs"))
        
        # read in teh model data for this species
        load(file = paste0("SPCD_standata_general_full_standardized_v3/SPCD_",SPCD.id,"remper_correction_0.5model_9.Rdata"))
        
        # set up lookup tables to link training and testing pSurv indices to tree.id, state, county, Survival observations, volfac expansion factors
        state_county_lookup_train <- train.data %>% mutate(tree.id = 1:length(train.data$state))%>% 
          left_join(.,plot_expansions, by = join_by(state, county, pltnum, cndtn, PLOT.ID, cycle))%>%
          select(tree.id, remper, volfac, S, pltnum, state, county, expacr) %>% 
          mutate(data.type = "in-sample", 
                 ST_CTY = paste0(state, "_", county))
        
        state_county_lookup_test <- test.data %>% mutate(tree.id = 1:length(test.data$state))%>% 
          left_join(.,plot_expansions, by = join_by(state, county, pltnum, cndtn, PLOT.ID, cycle))%>%
          select(tree.id, remper, volfac, S, pltnum, state, county, expacr) %>% 
          mutate(data.type = "out-of-sample", 
                 ST_CTY = paste0(state, "_", county))
        
        
        # set up matrices for matrix multiplication:
        # remper-number of years between remeasurements
        # volfac - trees per acre conversion factor for the tree
        # expacre - acres represented by the plot in state-level reporting
        
        Remper_matrix <- matrix(mod.data$Remper, nrow = nrow(pSurv_hat_samps), ncol = ncol(pSurv_hat_samps), byrow = TRUE)
        Remperoos_matrix <- matrix(mod.data$Remperoos, nrow = nrow(pSurv_rep_samps), ncol = ncol(pSurv_rep_samps), byrow = TRUE)
        
        volfac_hat_matrix <- matrix(train.data$volfac, nrow = nrow(pSurv_hat_samps), ncol = ncol(pSurv_hat_samps), byrow = TRUE)
        volfac_rep_matrix <- matrix(test.data$volfac, nrow = nrow(pSurv_rep_samps), ncol = ncol(pSurv_rep_samps), byrow = TRUE)
        
        # get expacr matrices:
        Expacr_rep_matrix <- matrix(state_county_lookup_test$expacr, 
                                    nrow = nrow(pSurv_rep_samps), ncol = ncol(pSurv_rep_samps), byrow = TRUE)
        
        Expacr_hat_matrix  <- matrix(state_county_lookup_train$expacr, 
                                     nrow = nrow(pSurv_hat_samps), ncol = ncol(pSurv_hat_samps), byrow = TRUE)
        
        # general mortality rate calculations:
        
        # Pred_Mort_rate = Predicted deaths/Exposure
        # Obs_Mort_rate =  Observed deaths/Exposure
        # Predicted deaths over the remper at the tree-scale: 
              # E_mort = (1-pSurv_remper)*volfac or (1-(pSurv_annual^remper))*volfac, if pSurv is annual
              # units = trees/acre
        
        # Observed deaths over the remper at the tree-scale: 
        # Deaths_i = (1-S)*volfac 
        # units = trees/acre
        
        # Exposure at the tree-scale (assume same for predicted and observed): 
              # Exposure = remper*volfac 
              # units = tree-years/acre
        
        # calculate the posterior expected # of trees dead based on tree-level volfac
        Post_Emort_remp_volfac_rep <- (1-pSurv_rep_samps)*volfac_rep_matrix
        Post_Emort_remp_volfac_hat <- (1-pSurv_hat_samps)*volfac_hat_matrix
        
        # calculate matrices of observed exposure for these predictions (#trees/acre * n years) tree-years
        Exposure_mat_rep <-  volfac_rep_matrix*Remperoos_matrix
        Exposure_mat_hat <-  volfac_hat_matrix*Remper_matrix
        
        
        # aggregate by states and counties:----------
        # for aggregation to higher domains (state, county, etc), we need to account for differences in plot designs
        # using expacre to put everything on a total basis
        
        # Pred_Mort_rate = Predicted deaths/Exposure
        # Obs_Mort_rate =  Observed deaths/Exposure
        
        # Predicted deaths over the remper at the tree-scale: 
        # E_mort = (1-pSurv_remper)*volfac*expacre or (1-(pSurv_annual^remper))*volfac*expacre, if pSurv is annual
        # units = trees
        
        # Observed deaths over the remper at the tree-scale: 
        # Deaths_i = (1-S)*volfac*expacre 
        # units = trees
        
        # Exposure at the tree-scale (assume same for predicted and observed): 
        # Exposure = remper*volfac 
        # units = tree-years/acre
        # Exposure_acre = remper*volfac*expacre 
        # units = tree-years/acre
        # 
        
        # weight posterior mortality rates by expacr for county and state scales
        # calculate the posterior expected # of trees dead based on tree-level volfac
        Post_Emort_remp_volfac_expacre_rep <- Expacr_rep_matrix*((1-pSurv_rep_samps)*volfac_rep_matrix)
        Post_Emort_remp_volfac_expacre_hat <- Expacr_hat_matrix*(1-pSurv_hat_samps)*volfac_hat_matrix
        
        # calculate matrices of observed exposure for these predictions (#trees/acre * n years) tree-years
        Exposure_expacre_mat_rep <-   Expacr_rep_matrix*(volfac_rep_matrix*Remperoos_matrix)
        Exposure_expacre_mat_hat <-   Expacr_hat_matrix*(volfac_hat_matrix*Remper_matrix)
        
        
        # overall estimates of observed and predicted mortality rates ---  
        # note this is a bit redundant here because we are not indexing by any group, but we do in functions below, so this is just consistent
        # calculate tree-level mortality rates
        # for out of sample data:
        st_Pmort_rep <- Post_Emort_remp_volfac_rep
        st_Expos_rep <- Exposure_mat_rep
        st_obs_train <- state_county_lookup_test%>% 
          mutate(Exposure_i = volfac*remper)%>% # tree obseved exposure
          mutate(Deaths_i = volfac*(1-S))%>% # tree observed deaths
          mutate(Exposure_EXPN_i = Exposure_i*expacr, # tree expanded exposure
                 Deaths_EXPN_i = Deaths_i*expacr) # tree expanded deaths
        
        # for in sample data:
        st_Pmort_hat <- Post_Emort_remp_volfac_hat
        st_Expos_hat <- Exposure_mat_hat
        st_obs_test <- state_county_lookup_train %>% 
          mutate(Exposure_i = volfac*remper)%>%
          mutate(Deaths_i = volfac*(1-S)) %>% 
          mutate(Exposure_EXPN_i = Exposure_i*expacr, 
                 Deaths_EXPN_i = Deaths_i*expacr)
        
        
        # do the same with expacre weighted values:
        st_Pmort_rep_EXPN <- Post_Emort_remp_volfac_expacre_rep # expanded predicted mortality
        st_Expos_rep_EXPN <- Exposure_expacre_mat_rep # expanded exposure matrix
        
        
        st_Pmort_hat_EXPN <- Post_Emort_remp_volfac_expacre_hat
        st_Expos_hat_EXPN <- Exposure_expacre_mat_hat
        
        # # tree-level mortality rates--
        # we could output this, but we omitting for now
        # tree_E_M_expn_i <- cbind(st_Pmort_hat_EXPN, st_Pmort_rep_EXPN)
        # tree_Expos_expn_i <- cbind(st_Expos_hat_EXPN, st_Expos_rep_EXPN)
        # tree_Pred_mort_expn <- tree_E_M_expn_i/tree_Expos_expn_i
        # 
        # 
        # tree_obs_Expos_expn_i = c(st_obs_train$Exposure_EXPN_i, st_obs_test$Exposure_EXPN_i)# calculating the "exposure for each tree" remper* volfac, in tree/acre-years
        # tree_obs_Deaths_expn_i= c(st_obs_train$Deaths_EXPN_i, st_obs_test$Deaths_EXPN_i) # calculate the number of deaths per tree observed over remper
        # 
        # tree_obs_mort_expn <- tree_obs_Deaths_expn_i/tree_obs_Expos_expn_i
        # 
        # tree_mort_rate_summary <- data.frame(Obs_mort_expn = tree_obs_mort_expn, 
        #                                      Pred_mort_expn.mean = colMeans(tree_Pred_mort_expn))
        # 
        # 
        
        # Regional mortality rates:
        OVERALL_mort_rate <- data.frame(E_mort = rowSums(cbind(st_Pmort_hat, st_Pmort_rep)), 
                                        Exposure = rowSums(cbind(st_Expos_hat, st_Expos_rep)),
                                        E_mort_expn = rowSums(cbind(st_Pmort_hat_EXPN, st_Pmort_rep_EXPN)), 
                                        Exposure_expn = rowSums(cbind(st_Expos_hat_EXPN, st_Expos_rep_EXPN)), 
                                        
                                        
                                        Exposure_obs = sum(st_obs_test$Exposure_i, st_obs_train$Exposure_i),# calculating the "exposure for each tree" remper* volfac, in tree/acre-years
                                        Deaths_obs = sum(st_obs_test$Deaths_i, st_obs_train$Deaths_i),
                                        
                                        Exposure_expn_obs = sum(st_obs_test$Exposure_EXPN_i, st_obs_train$Exposure_EXPN_i),# calculating the "exposure for each tree" remper* volfac, in tree/acre-years
                                        Deaths_expn_obs = sum(st_obs_test$Deaths_EXPN_i, st_obs_train$Deaths_EXPN_i),#, # calculate the number of deaths per tree observed over remper
                                        total_obs = length(c(st_obs_test$Deaths_i, st_obs_train$Deaths_i)),
                                        
                                        SPCD = SPCD.id, 
                                        model.number = ifelse(weighting %in% "stacking_wts", "Stacked", weighting), 
                                        model.type = model.type) %>%
          mutate(Obs_mort_rate =  Deaths_obs/Exposure_obs,
                 Pred_mort_rate = E_mort/Exposure, 
                 Obs_mort_rate_expn = Deaths_expn_obs/Exposure_expn_obs, 
                 Pred_mort_rate_expn = E_mort_expn/Exposure_expn)
        
        #ggplot(OVERALL_mort_rate)+geom_histogram(aes(Pred_mort_rate_expn))+geom_vline(aes(xintercept = Obs_mort_rate_expn))
        
        
        # calculate state-level estimates of observed and predicted mortality rates
        
        # unique states
        state_ids <- unique(c(state_county_lookup_train$state, state_county_lookup_test$state))
        
        state_mort_rate_list <- lapply(state_ids, FUN = function(state_cd){
          
          # get index for the focal state
          st_index_train <- state_county_lookup_train$state == state_cd
          st_index_test <- state_county_lookup_test$state == state_cd
          
          # use index to get only E_mort, exposure, and observation data from the state of interest
          st_Pmort_rep <- Post_Emort_remp_volfac_rep[,st_index_test]
          st_Expos_rep <- Exposure_mat_rep[,st_index_test]
          st_obs_train <- state_county_lookup_test[st_index_test,]%>% 
            mutate(Exposure_i = volfac*remper)%>%
            mutate(Deaths_i = volfac*(1-S))%>% 
            mutate(Exposure_EXPN_i = Exposure_i*expacr, 
                   Deaths_EXPN_i = Deaths_i*expacr)
          
          st_Pmort_hat <- Post_Emort_remp_volfac_hat[,st_index_train]
          st_Expos_hat <- Exposure_mat_hat[,st_index_train]
          st_obs_test <- state_county_lookup_train[st_index_train,] %>% 
            mutate(Exposure_i = volfac*remper)%>%
            mutate(Deaths_i = volfac*(1-S)) %>% 
            mutate(Exposure_EXPN_i = Exposure_i*expacr, 
                   Deaths_EXPN_i = Deaths_i*expacr)
          
          
          # do the same with expacre weighted values:
          st_Pmort_rep_EXPN <- Post_Emort_remp_volfac_expacre_rep[,st_index_test]
          st_Expos_rep_EXPN <- Exposure_expacre_mat_rep[,st_index_test]
          
          
          st_Pmort_hat_EXPN <- Post_Emort_remp_volfac_expacre_hat[,st_index_train]
          st_Expos_hat_EXPN <- Exposure_expacre_mat_hat[,st_index_train]
          
          
          # combine together:
          st_mort_rate <- data.frame(E_mort = rowSums(cbind(st_Pmort_hat, st_Pmort_rep)), 
                                     Exposure = rowSums(cbind(st_Expos_hat, st_Expos_rep)),
                                     E_mort_expn = rowSums(cbind(st_Pmort_hat_EXPN, st_Pmort_rep_EXPN)), 
                                     Exposure_expn = rowSums(cbind(st_Expos_hat_EXPN, st_Expos_rep_EXPN)), 
                                     
                                     state = state_cd,
                                     Exposure_obs = sum(st_obs_test$Exposure_i, st_obs_train$Exposure_i),# calculating the "exposure for each tree" remper* volfac, in tree/acre-years
                                     Deaths_obs = sum(st_obs_test$Deaths_i, st_obs_train$Deaths_i),
                                     
                                     Exposure_expn_obs = sum(st_obs_test$Exposure_EXPN_i, st_obs_train$Exposure_EXPN_i),# calculating the "exposure for each tree" remper* volfac, in tree/acre-years
                                     Deaths_expn_obs = sum(st_obs_test$Deaths_EXPN_i, st_obs_train$Deaths_EXPN_i),#, # calculate the number of deaths per tree observed over remper
                                     total_obs = length(c(st_obs_test$Deaths_i, st_obs_train$Deaths_i)),
                                     
                                     SPCD = SPCD.id, 
                                     model.number = ifelse(weighting %in% "stacking_wts", "Stacked", weighting),  
                                     model.type = model.type) %>%
            mutate(Obs_mort_rate =  Deaths_obs/Exposure_obs,
                   Pred_mort_rate = E_mort/Exposure, 
                   Obs_mort_rate_expn = Deaths_expn_obs/Exposure_expn_obs, 
                   Pred_mort_rate_expn = E_mort_expn/Exposure_expn)
          return(st_mort_rate)
        })
        
        state_mort_df <- do.call(rbind,state_mort_rate_list) 
        
        # get county level predicted mortality summaries:
        STCTY_ids <- unique(c(state_county_lookup_train$ST_CTY, state_county_lookup_test$ST_CTY))
        # ST_CT_id <- STCTY_ids[1]
        
        county_mort_rate_list <- lapply(STCTY_ids, FUN = function(ST_CT_id){
          
          # get index for the focal state
          st_index_train <- state_county_lookup_train$ST_CTY == ST_CT_id
          st_index_test <- state_county_lookup_test$ST_CTY == ST_CT_id
          
          # calculate tree level expected mortality and exposure 
          st_Pmort_rep <- Post_Emort_remp_volfac_rep[,st_index_test]
          st_Expos_rep <- Exposure_mat_rep[,st_index_test]
          st_obs_train <- state_county_lookup_test[st_index_test,]%>% 
            mutate(Exposure_i = volfac*remper)%>%
            mutate(Deaths_i = volfac*(1-S))%>% 
            mutate(Exposure_EXPN_i = Exposure_i*expacr, 
                   Deaths_EXPN_i = Deaths_i*expacr)
          
          st_Pmort_hat <- Post_Emort_remp_volfac_hat[,st_index_train]
          st_Expos_hat <- Exposure_mat_hat[,st_index_train]
          st_obs_test <- state_county_lookup_train[st_index_train,] %>% 
            mutate(Exposure_i = volfac*remper)%>%
            mutate(Deaths_i = volfac*(1-S)) %>% 
            mutate(Exposure_EXPN_i = Exposure_i*expacr, 
                   Deaths_EXPN_i = Deaths_i*expacr)
          
          
          # do the same with expacre weighted values:
          st_Pmort_rep_EXPN <- Post_Emort_remp_volfac_expacre_rep[,st_index_test]
          st_Expos_rep_EXPN <- Exposure_expacre_mat_rep[,st_index_test]
          
          
          st_Pmort_hat_EXPN <- Post_Emort_remp_volfac_expacre_hat[,st_index_train]
          st_Expos_hat_EXPN <- Exposure_expacre_mat_hat[,st_index_train]
          
          
          
          # combine together:
          st_mort_rate <- data.frame(E_mort = rowSums(cbind(st_Pmort_hat, st_Pmort_rep)), 
                                     Exposure = rowSums(cbind(st_Expos_hat, st_Expos_rep)), 
                                     E_mort_expn = rowSums(cbind(st_Pmort_hat_EXPN, st_Pmort_rep_EXPN)), 
                                     Exposure_expn = rowSums(cbind(st_Expos_hat_EXPN, st_Expos_rep_EXPN)), 
                                     
                                     
                                     ST_CTY = ST_CT_id,
                                     
                                     Exposure_obs = sum(st_obs_test$Exposure_i, st_obs_train$Exposure_i),# calculating the "exposure for each tree" remper* volfac, in tree/acre-years
                                     Deaths_obs = sum(st_obs_test$Deaths_i, st_obs_train$Deaths_i),#, # calculate the number of deaths per tree observed over remper
                                     Exposure_expn_obs = sum(st_obs_test$Exposure_EXPN_i, st_obs_train$Exposure_EXPN_i),# calculating the "exposure for each tree" remper* volfac, in tree/acre-years
                                     Deaths_expn_obs = sum(st_obs_test$Deaths_EXPN_i, st_obs_train$Deaths_EXPN_i),#, # calculate the number of deaths per tree observed over remper
                                     
                                     total_obs = length(c(st_obs_test$Deaths_i, st_obs_train$Deaths_i)),
                                     SPCD = SPCD.id, 
                                     model.number = ifelse(weighting %in% "stacking_wts", "Stacked", weighting), 
                                     model.type = model.type) %>%
            
            mutate(Obs_mort_rate =  Deaths_obs/Exposure_obs,
                   Pred_mort_rate = E_mort/Exposure,
                   Obs_mort_rate_expn = Deaths_expn_obs/Exposure_expn_obs, 
                   Pred_mort_rate_expn = E_mort_expn/Exposure_expn)
          return(st_mort_rate)
        })
        
        county_mort_df <- do.call(rbind, county_mort_rate_list) 
        
        county_mort_summary <- county_mort_df %>%  
          group_by(SPCD, model.number, model.type, ST_CTY)%>%
          
          summarise(obs_M_median = median(Obs_mort_rate, na.rm =TRUE), 
                    n_obs = median(total_obs, na.rm =TRUE),
                    pred_M_median = median(Pred_mort_rate, na.rm =TRUE), 
                    pred_M_5.ci.lo = quantile(Pred_mort_rate, 0.05, na.rm =TRUE), 
                    pred_M_95.ci.hi = quantile(Pred_mort_rate, 0.95, na.rm =TRUE),
                    pred_M_25.ci.lo = quantile(Pred_mort_rate, 0.25, na.rm =TRUE), 
                    pred_M_75.ci.hi = quantile(Pred_mort_rate, 0.75, na.rm =TRUE), 
                    
                    obs_M_expn_median = median(Obs_mort_rate_expn, na.rm =TRUE), 
                    
                    pred_M_expn_median = median(Pred_mort_rate_expn, na.rm =TRUE), 
                    pred_M_expn_5.ci.lo = quantile(Pred_mort_rate_expn, 0.05, na.rm =TRUE), 
                    pred_M_expn_95.ci.hi = quantile(Pred_mort_rate_expn, 0.95, na.rm =TRUE),
                    pred_M_expn_25.ci.lo = quantile(Pred_mort_rate_expn, 0.25, na.rm =TRUE), 
                    pred_M_expn_75.ci.hi = quantile(Pred_mort_rate_expn, 0.75, na.rm =TRUE), 
                    .groups = "drop")
        
        
        # save the county-level outputs for each species:
        # save the samples weighted by volface of each species in the county
        
        qs_save(county_mort_df, paste0(
          output.dir,
          "SPCD_stanoutput_cmdstan/predicted_mort/ST_CTY_mort_rate_samps_",SPCD.id,"_",weighting,".qs"
        ))
        qs_save(state_mort_df, paste0(
          output.dir,
          "SPCD_stanoutput_cmdstan/predicted_mort/State_mort_rate_samps_",SPCD.id,"_",model.type, "_",weighting,".qs"
        ))
        
        qs_save(OVERALL_mort_rate, paste0(
          output.dir,
          "SPCD_stanoutput_cmdstan/predicted_mort/Regional_mort_rate_samps_",SPCD.id,"_",model.type, "_",weighting,".qs"
        ))
        
    return(county_mort_summary)
}

stacked.ests <- do.call(rbind, lapply(length(spp.table$SPCD):1, FUN = function(x){
  calculate_state_county_rates_BMA(SPCD.id = spp.table$SPCD[x], 
                                   weighting =  "stacking_wts",
                                   model.type = "Model Average")
}))
stacked.ests$COMMON_NAME <- ref_species[match(stacked.ests$SPCD, ref_species$SPCD),]$COMMON_NAME
# plot predicted vs observed
stacked.ests %>% filter(n_obs > 25)|>
  ggplot()+geom_point(aes(x = obs_M_expn_median, y = pred_M_expn_median, color = COMMON_NAME))+
  geom_errorbar(aes(x = obs_M_expn_median, ymin = pred_M_expn_5.ci.lo, ymax = pred_M_expn_95.ci.hi, color = COMMON_NAME))+
  geom_errorbar(aes(x = obs_M_expn_median, ymin = pred_M_expn_25.ci.lo, ymax = pred_M_expn_75.ci.hi, color = COMMON_NAME))+
  theme_bw()

###################################################################################
# get weighted posterior draws for each beta

# read in betas for hiearachical models, save each species betas in a separate file
betas.heir <- list.files(paste0(output.dir, "SPCD_stanoutput_cmdstan/betas/"), pattern = "hierarchical_mort_model_", full.names = T)
hier_draws_path <- function(k)
  paste0(output.dir, "SPCD_stanoutput_cmdstan/betas/u_beta_alpha_samps_hierarchical_mort_model_",k,"_niter_1000_nchain_4.qs")  # ADAPT: your actual naming

hier_draws_spp_save_path <- function(spcd, k)
  paste0(output.dir, "SPCD_stanoutput_cmdstan/betas/u_beta_alpha_samps_hierarchical_mort_model_",k,"_SPCD_",spcd,"_niter_1000_nchain_4.qs")  # ADAPT: your actual naming


species_draws_path <- function(spcd, k){
  
  paste0(output.dir, "SPCD_stanoutput_cmdstan/betas/u_beta_alpha_samps_mort_model_",k,"_SPCD_",spcd,"_remper_correction_0.5_niter_1000_nchain_4.qs")
}
resave_hier_betas <- do.call(rbind, lapply(1:9, function(i){
  
  cat(paste0("model ", i, "\n"))
  ubetas <- qs2::qs_read(hier_draws_path(k = i) ) %>% subset_draws(., "u_beta")#%>% 
    #summarise_draws() %>% mutate(model.number = i) %>%
    # mutate(param = "u_beta")
    # 
  alphas <- qs2::qs_read(hier_draws_path(k = i) ) %>% subset_draws(., "alpha_SPP")
  
  for(s in 1:17){
    cat(s)
  species.ubeta_idx <- colnames(ubetas) %in% paste0("u_beta[", s, ",", 1:78, "]")
  species.alpha_idx <- colnames(alphas) %in% paste0("alpha_SPP[", s, "]")
  
  spp_alphas_betas <- bind_draws(alphas[, species.alpha_idx], ubetas[, species.ubeta_idx]) 
  
  colnames(spp_alphas_betas) <- c("alpha_SPP", paste0("u_beta[",1:(length(colnames(spp_alphas_betas)) -1),"]"))
  
  qs2::qs_save(object = spp_alphas_betas, 
               file = hier_draws_spp_save_path(spcd = spp.table[match(s, spp.table$SPP),]$SPCD, k = i))
  }
}))



# get posterior beta matrix function
get_posterior_betas_matrix <- function(structure,  spcd,  k, spp_index = NULL) {
  
  if (structure == "species_only") {
    qs2::qs_read(species_draws_path(spcd, k))          # draws x N_species, already per-species
  } else {
    qs2::qs_read(hier_draws_spp_save_path(spcd, k))            # draws x N_total (all species stacked)
    
  }
}

mix_beta_predictions <- function(pred_list, plan){
  
  # get the number of x values 
  grid_n <- 79
  
  out <- matrix(NA_real_, plan$N, grid_n)
  for (n in seq_len(plan$N)) {
    ncol_mod <-  ncol(pred_list[[plan$model[n]]])
    out[n, 1:ncol_mod] <- pred_list[[plan$model[n]]][plan$draw[n], 1:ncol_mod] # ensures we only add to the right columns
    
  }
  
  # keep a record of which model this draw came from
  list(draws = out, model = plan$model)
}

# now get a draw from the posterior betas based on the plan


stack_betas_alphas_species <- function(spcd, 
                            #y_type = c("y_rep"),
                            N = 4000, 
                            weighting = "stacking_wts",
                            y_plan = NULL) {
  
  # get the weighting and indices to sample from for each species
  weights_i <- load_species_weights(spcd, weighting)
  plan <- if (is.null(y_plan)) post_draw_plan(spcd, weights_i, N) else y_plan
  
  # get the full predicted list from all the models that have non-zero weights for this species
  pred_list <- list()
  for (structure in c("species_only", "hierarchical")) {
    for (k in 1:9) {
      # wname = naming structure of the model in the weight dataframe
      wname <- paste0(if (structure == "hierarchical") "Hierarchical_model_" else "Species_model_", k)
      if (!(wname %in% plan$model)) next  # zero-weight models were never sampled; skip loading
      
      # get the matrix of y values
      mat <- tryCatch(get_posterior_betas_matrix(structure, spcd, k), error = function(e) NULL)
      if (is.null(mat)) next
      
      # get the max number of iterations needed from this model and check that we have that
      needed <- max(plan$draw[plan$model == wname])
      if (nrow(mat) < needed) {
        stop(wname, "'s saved ", type, " has only ", nrow(mat), " draws, but the plan ",
             "needs draw index ", needed, ".")
      }
      pred_list[[wname]] <- mat # 
    }
  }
  #unique(plan$model) %in% names(pred_list) 
  list(pred_list = pred_list, 
       plan = plan)





  # use indices in the plan to output the predicted y or psurv values
  mixed <- mix_beta_predictions(pred_list, plan)  # takes the draws plan and gets the pred_list
  
  # get the number of tree observations for the species
  var_names <- c("alpha_SPP", paste0("u_beta[",1:78,"]"))
  draws_df <- as.data.frame(mixed$draws)
  colnames(draws_df) <- var_names  
  draws_matrix <- as_draws_matrix(draws_df)
  
  #draws_matrix$model_source <- mixed$model
  qs2::qs_save(draws_matrix,     paste0(output.dir, "SPCD_stanoutput_cmdstan/betas/u_betas_alpha_samps_SPCD_", spcd, "_", weighting, ".qs"))
  
 draws_df <-  draws_df %>% mutate(model_source = mixed$model)%>%
    mutate(spcd = spcd,
           weighting = weighting) %>%
    select(model_source, spcd, weighting, everything())
  
  qs2::qs_save(draws_matrix,     paste0(output.dir, "SPCD_stanoutput_cmdstan/betas/u_betas_alpha_samps_SPCD_", spcd, "_", weighting, ".qs"))
  qs2::qs_save(draws_df,     paste0(output.dir, "SPCD_stanoutput_cmdstan/betas/mixed_df_betas_samps_SPCD_", spcd, "_", weighting, ".qs"))
  
  #draws_df
  list(draws = draws_matrix,
       model_source = mixed$model,
       spcd = spcd,
       weighting = weighting)
}

# example usage for pSurv_hat for one species
betas_261 <- stack_betas_alphas_species(spcd = 261, 
                             N = 4000, 
                             y_plan = NULL)

betas_261$draws %>% as_draws_df() %>% summarise_draws()
summarize_beta_posteriors <- function(x) {
  c(median = median(x, na.rm =TRUE),
    ci.lo.2.5 = quantile(x, probs = c(0.025), na.rm =TRUE),
    ci.hi.97.5 = quantile(x, probs = c(0.975), na.rm =TRUE), 
    ci.lo.5 = quantile(x, probs = c(0.05), na.rm =TRUE),
    ci.hi.95 = quantile(x, probs = c(0.95), na.rm =TRUE), 
    ci.lo.25 = quantile(x, probs = c(0.25), na.rm =TRUE),
    ci.hi.75 = quantile(x, probs = c(0.75), na.rm =TRUE))
}

stacked_betas_list <- list()
for(s in 1:length(spcd_ids)){
  betas_spcd <- stack_betas_alphas_species(spcd = spcd_ids[s], 
                                          N = 4000, 
                                          y_plan = NULL)
  
  #betas_spcd$model_source[4]
  betas.NA.omited <- betas_spcd$draws %>% summarise_draws(., summarize_beta_posteriors) # just omits NA values
  
  beta_draws <- betas_spcd$draws
  beta_draws[is.na(beta_draws)] <- 0
  
  betas.NA.as.0 <- beta_draws %>% summarise_draws(., summarize_beta_posteriors) %>%
    mutate(spcd = betas_spcd$spcd, 
           weighting = betas_spcd$weighting)# just omits NA values
  colnames(betas.NA.as.0) <- c("variable", "median", "ci.lo.2.5", "ci.hi.97.5", "ci.lo.5", "ci.hi.95", "ci.lo.25", "ci.hi.75", 
                               "spcd", "weighting")
  stacked_betas_list[[s]] <- betas.NA.as.0
}
stacked_betas_df <- do.call(rbind, stacked_betas_list)

stacked_betas_df %>% filter(variable %in% "alpha_SPP")|>
  ggplot()+geom_pointrange(aes(x = as.character(spcd), y = median, ymin = ci.lo.5, ymax = ci.hi.95))+
  ylab("alpha_SPP stacked")

stacked_betas_df %>% filter(variable %in% "u_beta[1]")|>
  ggplot()+geom_pointrange(aes(x = as.character(spcd), y = median, ymin = `ci.lo.25%`, ymax = `ci.hi.75%`))+
  ylab("u_beta[1] effect")

# create a function to process/summarise the stacked posteriors
# 1. Save the weighted samples pSurv_hat, pSurv_rep, y_hat, y_rep for each species as draws_df for later use
# 2. Calculate is and oos AUC scores, confusion matrices for each species
# 3. Summarise the pSurv_hat, pSurv_rep, y_hat, y_rep for each species


spcd = 261
weighting = "stacking_wts"
N = 4000
process_bma_posterior_pSurv_y <- function(spcd, weighting, N){
  # start function
  
  cat(paste0("getting weighted posteriors for spcd ", spcd, "\n"))
  # get species in and out of sample y data to compare to:
  load(paste0("SPCD_standata_general_full_standardized_v3/", "SPCD_", spcd, "remper_correction_0.5model_9.Rdata"))
  #length(mod.data$y)
  
  
  # set up plan for each species once, then feed it inot the stack_y_species
  weights_i <- load_species_weights(spcd, weighting)
  plan <- post_draw_plan(spcd, weights_i, N)
  
  # get the weighted probabilities and predicted classifications
  psurv_hat <- stack_y_species(spcd = spcd, 
                               y_type = "pSurv_hat", 
                               N = 4000, 
                               y_plan = plan)
  psurv_rep <- stack_y_species(spcd = spcd, 
                               y_type = "pSurv_rep", 
                               N = 4000, 
                               y_plan = plan)
  
  y_hat <- stack_y_species(spcd = spcd, 
                           y_type = "y_hat", 
                           N = 4000, 
                           y_plan = plan)
  
  y_rep <- stack_y_species(spcd = spcd, 
                           y_type = "y_rep", 
                           N = 4000, 
                           y_plan = plan)
  
  
  cat(paste0("saving weighted posterior draws for ", spcd, "\n"))
  # get just the draws
  y_rep_samps <- y_rep$draws
  y_hat_samps <- y_hat$draws
  
  pSurv_rep_samps <- psurv_rep$draws
  pSurv_hat_samps <- psurv_hat$draws
  
  qs2::qs_save(y_rep_samps,     paste0(output.dir, "SPCD_stanoutput_cmdstan/predicted_mort/y_rep_samps_SPCD_", spcd, "_", weighting, ".qs"))
  qs2::qs_save(y_hat_samps,     paste0(output.dir, "SPCD_stanoutput_cmdstan/predicted_mort/y_hat_samps_SPCD_", spcd, "_", weighting, ".qs"))
  qs2::qs_save(pSurv_rep_samps, paste0(output.dir, "SPCD_stanoutput_cmdstan/predicted_mort/pSurv_rep_samps_SPCD_", spcd, "_", weighting, ".qs"))
  qs2::qs_save(pSurv_hat_samps, paste0(output.dir, "SPCD_stanoutput_cmdstan/predicted_mort/pSurv_hat_samps_SPCD_", spcd, "_", weighting, ".qs"))
  
  # convert remper probabilities to annualized probabilities of survival:
  Remper_matrix <- matrix(mod.data$Remper, nrow = nrow(pSurv_hat_samps), ncol = ncol(pSurv_hat_samps), byrow = TRUE)
  #length(unique(Remper_matrix[,1])) ==1 # check that the byrow is right
  
  pSannual_hat <- pSurv_hat_samps^(1/Remper_matrix) 
  
  Remperoos_matrix <- matrix(mod.data$Remperoos, nrow = nrow(pSurv_rep_samps), ncol = ncol(pSurv_rep_samps), byrow = TRUE)
  #length(unique(Remper_matrix[,1])) ==1 # check that the byrow is right
  
  pSannual_rep <- pSurv_rep_samps^(1/Remperoos_matrix) 
  
  qs2::qs_save(pSannual_rep, paste0(output.dir, "SPCD_stanoutput_cmdstan/predicted_mort/pSurv_annual_rep_samps_SPCD_", spcd, "_", weighting, ".qs"))
  qs2::qs_save(pSannual_hat, paste0(output.dir, "SPCD_stanoutput_cmdstan/predicted_mort/pSurv_annual_hat_samps_SPCD_", spcd, "_", weighting, ".qs"))
  
  # does the predicted count of y_hat match that in the data?
  is.ppc.bar <- ppc_bars(y = mod.data$y, as_draws_matrix(y_hat_samps))+xlab("Survival Classifcation (In-sample)")
  oos.ppc.bar <- ppc_bars(y = mod.data$ytest, as_draws_matrix(y_rep_samps))+xlab("Survival Classifcation (out-of-sample)")
  combined.bar <- cowplot::plot_grid( is.ppc.bar, oos.ppc.bar, align = "hv")
  cowplot::save_plot(paste0(output.dir, "SPCD_stanoutput_cmdstan/images/post_pred_y_SPCD_", spcd, "_", weighting, ".png"), 
                     plot = combined.bar, 
                     bg = "white")
  rm(is.ppc.bar, oos.ppc.bar, combined.bar)
  
  
  # --- AUC + confusion matrix over posterior draws ------------------------
  actuals     <- mod.data$y
  actuals.oos <- mod.data$ytest
  cat(paste0("estimating AUC and confusion matrix on weighted posterior draws for ", spcd, "\n"))
  
  # library(pROC)
  
  
  # calculate AUC for each draw:------
  AUC.is.samples.df  <- apply(pSurv_hat_samps, 1, function(p) as.numeric(pROC::auc(actuals, p, quiet = TRUE)))
  AUC.oos.samples.df <- apply(pSurv_rep_samps, 1, function(p) as.numeric(pROC::auc(actuals.oos, p, quiet = TRUE)))
  #rm(pSurv_hat_samps, pSurv_rep_samps, y_hat.quant, y_rep.quant, pSurv_hat.quant, pSurv_rep.quant); gc()
  
  preds.is.class <- y_hat_samps == 1
  confusion.is_draws <- data.frame(
    TP_draws = rowSums(preds.is.class[, actuals == 1, drop = FALSE]),
    FP_draws = rowSums(preds.is.class[, actuals == 0, drop = FALSE]),
    TN_draws = rowSums(!preds.is.class[, actuals == 0, drop = FALSE]),
    FN_draws = rowSums(!preds.is.class[, actuals == 1, drop = FALSE])
  ) %>%
    mutate(`True survival rate` = TP_draws / (TP_draws + FN_draws),
           `True mortality rate` = TN_draws / (TN_draws + FP_draws),
           model.number = "BMA", 
           type = "in-sample",
           model.type = weighting, 
           SPCD = spcd)
  
  preds.oos.class <- y_rep_samps == 1
  confusion.oos_draws <- data.frame(
    TP_draws = rowSums(preds.oos.class[, actuals.oos == 1, drop = FALSE]),
    FP_draws = rowSums(preds.oos.class[, actuals.oos == 0, drop = FALSE]),
    TN_draws = rowSums(!preds.oos.class[, actuals.oos == 0, drop = FALSE]),
    FN_draws = rowSums(!preds.oos.class[, actuals.oos == 1, drop = FALSE])
  ) %>%
    mutate(`True survival rate` = TP_draws / (TP_draws + FN_draws),
           `True mortality rate` = TN_draws / (TN_draws + FP_draws),
           model.number = "BMA", 
           type = "out-of-sample",
           model.type = "weighting", 
           SPCD = spcd)
  
  confusion.is_draws$AUC  <- AUC.is.samples.df
  confusion.oos_draws$AUC <- AUC.oos.samples.df
  AUC.confusion_draws <- rbind(confusion.is_draws, confusion.oos_draws)
  
  qs2::qs_save(AUC.confusion_draws,
               paste0(output.dir, "SPCD_stanoutput_cmdstan/AUC/AUC_draws_SPCD_", spcd, "_", weighting, ".qs"))
}

# run for all of the species
process_bma_posterior_pSurv_y(spcd = 261,
                              weighting = "stacking_wts",
                              N = 4000)


# outputs are saved for later use in manuscript figures document
