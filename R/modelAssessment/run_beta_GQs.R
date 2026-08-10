# run_bma_generated_quantities.R
# Author: Kelly A Heilman
# This script uses cmdstan generate_quantities to estimate marginal responses to each main beta effect.
# Steps include:
# Building a covariate prediction grid 
# runs cmdstanr generate_quantities() using posterior saved draws for betas, alphas, etc
# saves posterior predictions of annual probability of survival in response to varying main effects
# output is used in bayesian model averaging

# Setting up ----
library(cmdstanr)
library(posterior)
library(dplyr)

source("R/modelAssessment/bma_prediction_helpers.R")  # main_effects, int_labels, build_X_grid()

output.dir = "C:/Users/KellyHeilman/Box/01. kelly.heilman Workspace/mortality/Eastern-Mortality/mortality_models/"

spp.table <- read.csv(file = paste0(output.dir, "/data/Hierarchical_obs_model_7.csv"))

spcd_ids <- spp.table$SPCD
n_spp_total       <- length(spcd_ids) 

remper_correction <- 0.5  
n_models          <- 9
grid_length        <- 25     # points per covariate grid
remper_value_pred  <- 1      # annual-scale survival prob; 


# Functions to return files for training data, and posterior beta and alpha draws-----------------
# species-only training data + draws (one file per species x model)
species_traindata_path <- function(spcd, k){
  sprintf("SPCD_standata_general_full_standardized_v3/SPCD_%sremper_correction_%smodel_%s.Rdata",
          spcd, remper_correction, k)
  }
species_traindata_path(spcd = spcd_ids[1], k = 1)
"C:/Users/KellyHeilman/Box/01. kelly.heilman Workspace/mortality/Eastern-Mortality/mortality_models/SPCD_stanoutput_cmdstan/betas/u_beta_alpha_samps_mort_model_1_SPCD_12_remper_correction_0.5_niter_1000_nchain_4.qs"

species_draws_path <- function(spcd, k){
  #qs2::qs_save(fit.1, paste0(output.dir,"SPCD_stanoutput_cmdstan/fittedmodels/", model.name, ".qs"))
  
  paste0(output.dir, "SPCD_stanoutput_cmdstan/betas/u_beta_alpha_samps_mort_model_",k,"_SPCD_",spcd,"_remper_correction_0.5_niter_1000_nchain_4.qs")
}

# hierarchical training data + draws (one shared file per model structure, all species combined)
hier_traindata_path <- function(k)
  sprintf("SPCD_standata_json/hierarchical_data_model_%s.json", k)  # ADAPT


hier_draws_path <- function(k)
 paste0(output.dir, "SPCD_stanoutput_cmdstan/betas/u_beta_alpha_samps_hierarchical_mort_model_",k,"_niter_1000_nchain_4.qs")  # ADAPT: your actual naming

output_dir <- paste0(output.dir, "SPCD_stanoutput_cmdstan/model_avg_draws")
dir.create(output_dir, showWarnings = FALSE)


# ----------------------------------------------------------------------
# 1. Compile the species and hierarchical generated-quantities stan programs---


gq_species_mod <- cmdstan_model("modelcode/mort_model_general_predict.stan")
gq_hier_mod    <- cmdstan_model("modelcode/hierarchical_QR_gen_quants.stan")

# ----------------------------------------------------------------------
# 2. Build a prediction grid over xM ranges to use in GQ predictions ---
# Vary main effects, and assume any models that are missing that effect have 0 as the xgrid
# note that this is missing the interaction terms

build_tagged_grid <- function(model_colnames, grid_length = 25,
                               focal_vars_requested = main_effects,
                               include_missing_as_flat = FALSE) {
  
  present_vars <- intersect(focal_vars_requested, model_colnames)
  missing_vars <- setdiff(focal_vars_requested, model_colnames)

  pieces <- lapply(present_vars, function(fv) {
    
   
    grid_vals <- seq(-2, 2, length.out = grid_length)
    
    Xg <- build_X_grid(fv, grid_vals, main_effects, int_labels)
    Xg <- Xg[, model_colnames, drop = FALSE]
    list(X = Xg, focal_var = fv, grid_val = grid_vals)
  })

  if (include_missing_as_flat && length(missing_vars) > 0) {
    flat_pieces <- lapply(missing_vars, function(fv) {
      grid_vals <- seq(-2, 2, length.out = grid_length)  # label only -- X doesn't vary
      baseline_row <- setNames(rep(0, length(model_colnames)), model_colnames)
      Xg <- matrix(rep(baseline_row, each = grid_length), nrow = grid_length,
                   dimnames = list(NULL, model_colnames))
      list(X = Xg, focal_var = fv, grid_val = grid_vals)
    })
    pieces <- c(pieces, flat_pieces)
  }

  X_all <- do.call(rbind, lapply(pieces, `[[`, "X"))
  tag_focal <- unlist(lapply(pieces, function(p) rep(p$focal_var, nrow(p$X))))
  tag_grid  <- unlist(lapply(pieces, `[[`, "grid_val"))

  list(X = X_all, focal_var = tag_focal, grid_val = tag_grid)
}



# ----------------------------------------------------------------------
# 3a. Function to run Species-only generated quantities ----

run_species_only_gq <- function(spcd, k, grid_obj, remper_value = remper_value_pred) {
  
  load(species_traindata_path(spcd, k))  # brings in `mod.data` (N, K, xM, y, Remper)
  
  
  draws <- qs2::qs_read(species_draws_path(spcd, k))  # cmdstanr draws object or CmdStanMCMC fit
  colnames(draws)
  
  Nrep <- nrow(grid_obj$X)
  data_list <- list(
    N = mod.data$N, 
    K = mod.data$K, 
    xM = mod.data$xM,
    y = mod.data$y, 
    Remper = mod.data$Remper,
    Nrep = Nrep, 
    xMrep = grid_obj$X,
    Remperoos = rep(remper_value, Nrep)
  )

  gq <- gq_species_mod$generate_quantities(fitted_params = as_draws_array(draws), 
                                           data = data_list,
                                           parallel_chains = 4)
  
  prob_draws <- gq$draws("mMhat", format = "matrix")  # draws x Nrep

  list(prob_draws = prob_draws, focal_var = grid_obj$focal_var,
       grid_val = grid_obj$grid_val, spcd = spcd, model = k, structure = "species_only")
}

# ----------------------------------------------------------------------
# 3b. Function to run Hierarchical generated quantities ---

# note this fuction isnt used anymore
run_hier_gq <- function(spp_index, spcd, k, grid_obj, remper_value = remper_value_pred) {
  
  hier_train <- fromJSON(hier_traindata_path(k))  # list: N, K, Nspp, SPP, xM, Remper
  draws <- qs2::qs_read(hier_draws_path(k)) %>% 
    subset_draws(., c("mu_alpha", "sigma_alpha", "alpha_SPP", "u_beta", "mu_beta"))
 

  Nrep <- nrow(grid_obj$X)
  data_list <- list(
    N = 17,#hier_train$N, 
    K = hier_train$K, 
    Nspp = hier_train$Nspp,
    SPP = rep(spp_index, 17), 
    xM = hier_train$xM[1:17,], 
    Remper = hier_train$Remper[1:17],
    Nrep = Nrep,
    SPPrep = rep(spp_index, Nrep), 
    xMrep = grid_obj$X,
    Remperoos = rep(remper_value, Nrep)
  )

  gq <- gq_hier_mod$generate_quantities(fitted_params = as_draws_array(draws), 
                                        data = data_list, 
                                        parallel_chains = 4)
  prob_draws <- gq$draws("mMrep", format = "matrix")  # draws x Nrep

  list(prob_draws = prob_draws, focal_var = grid_obj$focal_var,
       grid_val = grid_obj$grid_val, spcd = spcd, model = k, structure = "hierarchical")
}



run_hier_gq_raw <- function(k, K, xMrep, SPPrep, remper_value = remper_value_pred) {
  
  Nrep <- nrow(xMrep)
  data_list <- list(
    Nspp = n_spp_total, K = K,
    N = 0L, SPP = integer(0), xM = matrix(numeric(0), nrow = 0, ncol = K),
    Remper = numeric(0),
    Nrep = Nrep, SPPrep = SPPrep, xMrep = xMrep,
    Remperoos = rep(remper_value, Nrep)
  )
 
  
  draws <- qs2::qs_read(hier_draws_path(k)) %>% 
    subset_draws(., c("mu_alpha", "sigma_alpha", "alpha_SPP", "u_beta", "mu_beta"))
  
  
  
  gq <- gq_hier_mod$generate_quantities(fitted_params = draws, data = data_list)
  gq$draws("mMrep", format = "matrix")  # draws x Nrep
  
}

run_hier_gq_all_species <- function(k, 
                                    model_colnames, 
                                    spcd_ids,
                                    spp_indices = seq_along(spcd_ids),
                                    focal_vars_requested = main_effects,
                                    include_missing_as_flat = TRUE,
                                    remper_value = remper_value_pred) {
  
  K <- length(model_colnames)
  n_spp <- length(spcd_ids)
  stopifnot(length(spp_indices) == n_spp)
  
  grid_obj <- build_tagged_grid(model_colnames, grid_length,
                                focal_vars_requested = focal_vars_requested,
                                include_missing_as_flat = include_missing_as_flat)
  rows_per_spp <- nrow(grid_obj$X)
  
  X_rep_full   <- do.call(rbind, replicate(n_spp, grid_obj$X, simplify = FALSE))
  SPPrep_full  <- rep(spp_indices, each = rows_per_spp)  # true fitted-model index, not list position
  focal_full   <- rep(grid_obj$focal_var, times = n_spp)
  gridval_full <- rep(grid_obj$grid_val, times = n_spp)
  
  prob_draws_full <- run_hier_gq_raw(k, K, X_rep_full, SPPrep_full, remper_value)
  
  # ---- reorganize: split columns back out by species ----
  results <- vector("list", n_spp)
  names(results) <- as.character(spcd_ids)
  for (s in seq_len(n_spp)) {
    cols <- ((s - 1) * rows_per_spp + 1):(s * rows_per_spp)
    results[[s]] <- list(
      prob_draws = prob_draws_full[, cols, drop = FALSE],
      focal_var = focal_full[cols], grid_val = gridval_full[cols],
      spcd = spcd_ids[s], model = k, structure = "hierarchical"
    )
  }
  results
}

# ----------------------------------------------------------------------
# 4. Helper function to build the covariate grid and run the right gq function based on the model structure---


generate_prob_grid <- function(structure, spcd, k, spp_index = NULL, model_colnames,
                                focal_vars_requested = main_effects,
                                include_missing_as_flat = TRUE) {
  
  grid_obj <- build_tagged_grid(model_colnames, grid_length,
                                 focal_vars_requested = focal_vars_requested,
                                 include_missing_as_flat = include_missing_as_flat)

  if (structure == "species_only") {
    run_species_only_gq(spcd, k, grid_obj)
  } else if (structure == "hierarchical") {
    stopifnot(!is.null(spp_index))
    run_hier_gq(spp_index, spcd, k, grid_obj)
  } else {
    stop("structure must be 'species_only' or 'hierarchical'")
  }
}

# ----------------------------------------------------------------------
# 5. Run GQ to get marginal responses to betas from each species-level model ---

for (k in 1:n_models) {

  for (s in seq_along(spcd_ids)) {
  spcd <- spcd_ids[s]



    # --- species-only ---
    # column names come from the model's own saved training data,
    # so subsetting always matches exactly regardless of model size
    load(species_traindata_path(spcd, k))
    cols_k <- colnames(mod.data$xM)
    cols_k[cols_k %in% "aspect.int"] <- "Ndep.scaled_aspect.int" # fix column name match for model 9
    
    res_sp <- tryCatch(
      generate_prob_grid(structure = "species_only", spcd, k, model_colnames = cols_k),
      error = function(e) { message("species_only failed: SPCD ", spcd, " model ", k, ": ", e$message); NULL }
    )
    if (!is.null(res_sp)) {
      qs2::qs_save(res_sp, file.path(output_dir, sprintf("SPCD_%s_model_%s_species_only.qs", spcd, k)))
    }

}
}  

  
  
  
# ----------------------------------------------------------------------
# 5. Run GQ to get marginal responses to betas from each Hierarchcial model ---

  for (k in 1:n_models) {
    
    load(species_traindata_path(spcd_ids[1], k))  # any one species' file just to get this
    # model's column set (shared across species)
    cols_k <- colnames(mod.data$xM)
    cols_k[cols_k %in% "aspect.int"] <- "Ndep.scaled_aspect.int" # fix column name match for model 9
    
    
    hier_results <- tryCatch(
      run_hier_gq_all_species(k, model_colnames = cols_k, spcd_ids),
      error = function(e) { message("hierarchical batch failed: model ", k, ": ", e$message); NULL }
    )
    
    # save species_level results
    if (!is.null(hier_results)) {
      for (s in seq_along(spcd_ids)) {
        qs2::qs_save(hier_results[[s]],
                file.path(output_dir, sprintf("SPCD_%s_model_%s_hierarchical.qs", spcd_ids[s], k)))
      }
    }
  }
  

