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
             
             # get summaries of annual survival probability
             annual_pMort.median = apply(annual_pMort, 2, quantile, probs = 0.5),
             annual_pMort.05.ci.lo = apply(annual_pMort, 2, quantile, probs = 0.05), 
             annual_pMort.95.ci.hi = apply(annual_pMort, 2, quantile, probs = 0.95), 
             
             
             # get summaries of decadal survival probability
             decadal_pSurv.median = apply(decadal_pSurv, 2, quantile, probs = 0.5),
             decadal_pSurv.05.ci.lo = apply(decadal_pSurv, 2, quantile, probs = 0.05), 
             decadal_pSurv.95.ci.hi = apply(decadal_pSurv, 2, quantile, probs = 0.95), 
             
             # get summaries of decadal survival probability
             decadal_pMort.median = apply(decadal_pMort, 2, quantile, probs = 0.5),
             decadal_pMort.05.ci.lo = apply(decadal_pMort, 2, quantile, probs = 0.05), 
             decadal_pMort.95.ci.hi = apply(decadal_pMort, 2, quantile, probs = 0.95)
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
  geom_ribbon(aes(x = grid_val, ymin =annual_pMort.05.ci.lo, ymax = annual_pMort.95.ci.hi), alpha = 0.5)



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
ggplot(data = marg_summary_all_df)+
  geom_ribbon(aes(x = grid_val, ymin = decadal_pMort.05.ci.lo, ymax = decadal_pMort.95.ci.hi, fill = COMMON_NAME), alpha = 0.5)+
  geom_line(aes(x = grid_val, y = decadal_pMort.median, color = COMMON_NAME))+
  facet_grid(vars(predictor), vars(COMMON_NAME))



# do the same thing but get a list of draws (this will be large!)
marg_draws_all_list <- lapply(spcd_ids, FUN = function(species_code){
  output_BMA_marginal_effects(spcd_idx = species_code, 
                              weighting = "stacking_wts", 
                              main_effects = main_effects, 
                              N = 4000, 
                              return_format = "draws")
})

marg_draws_all_df <- do.call(rbind, marg_draws_all_list)

ggplot(data = marg_summary_all_df)+
  geom_ribbon(aes(x = grid_val, ymin = decadal_pMort.05.ci.lo, ymax = decadal_pMort.95.ci.hi, fill = COMMON_NAME), alpha = 0.5)+
  geom_line(aes(x = grid_val, y = decadal_pMort.median, color = COMMON_NAME))+
  facet_grid(vars(predictor), vars(COMMON_NAME))



