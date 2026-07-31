# Author: Kelly A. Heilman
# joint_species_mortality_model_cmdstan.R
# About this script------------------------
# All cmdstan hierarchical models were fit in python, using cmdstanpy and saved as cmdstan .csv files
# Used cyverse DE
# stan code used in "modelcode/QR_reparam_hierarchical.stan"
# generated quantities from posterior parameter estimates stan "modelcode/hierarchical_QR_gen_quants.stan"
library(tidyverse)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(qs2)
library(jsonlite)
library(pROC)
library(loo)
library(data.table)
color_scheme_set("brightblue")

# set up input and output directories----

output.dir <- "C:/Users/KellyHeilman/Box/01. kelly.heilman Workspace/mortality/Eastern-Mortality/mortality_models/"
csv.subdir <- "SPCD_stanoutput_cmdstan/fittedmodels"   # <- parent folder holding one subfolder per model
json.dir   <- file.path("SPCD_standata_json") # <- persistent location of the model input jsons

nspp <- data.frame(SPCD = c(316, 318, 833, 832, 261, 531, 802, 129, 762, 12, 541, 97, 621, 400, 371, 241, 375))
nspp$Species <- paste(FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD), ]$GENUS,
                      FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD), ]$SPECIES)
nspp$COMMON <- FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD), ]$COMMON_NAME

spp.table <- data.frame(SPCD.id = nspp[1:17, ]$SPCD, spp = 1:17, COMMON = nspp[1:17, ]$COMMON)

# setting up remeasurement default length
remper.cor.vector <- c(0.5)
j <- 1
model.list <- 1:9 # number of models

# compile the Stan model used for prediction of generated quantities
gq.file <- file.path(getwd(), "modelcode", "hierarchical_QR_gen_quants.stan")
gq.mod  <- cmdstan_model(gq.file)


# 1. Set up functions to get each model csvs and only select what we need---


# each model's chain csvs live in their own folder 
# folder labels indicate the model number ("hierarchical_model_[model.number]")
# all models are labeled with save date time and a chain number, but
# function to list the model directory
get_model_dir <- function(model.number, output.dir, csv.subdir) {
  file.path(output.dir, csv.subdir, paste0("hierarchical_model_", model.number))
}
#get_model_dir(model.number = 1, output.dir = output.dir, csv.subdir = csv.subdir)

# function to list the csv files in the model directory
get_model_csv_files <- function(model.number, output.dir, csv.subdir) {
  model_dir <- get_model_dir(model.number, output.dir, csv.subdir)
  if (!dir.exists(model_dir)) {
    stop("No folder found for model ", model.number, " at: ", model_dir)
  }
  files <- Sys.glob(file.path(model_dir, "*.csv"))
  files <- files[!grepl("profile", files, ignore.case = TRUE)]  # profile csvs are separate, handled below
  if (length(files) == 0) {
    stop("No cmdstanpy output CSVs found for model ", model.number, " in: ", model_dir)
  }
  files
}

#get_model_csv_files(model.number = 1, output.dir = output.dir, csv.subdir = csv.subdir)



# Read the model input data from json locations (cmdstanpy did not save data) --

get_model_data <- function(model.number, json.dir) {
  json_path <- file.path(json.dir, paste0("hierarchical_data_model_", model.number, ".json"))
  if (!file.exists(json_path)) {
    stop("Data json not found for model ", model.number, " at: ", json_path)
  }
  fromJSON(json_path)
}

#get_model_data(model.number = 1, json.dir = json.dir)

# Function to get the model metadata only:
get_model_metadata <- function(csv_files) {
  read_cmdstan_csv(csv_files, variables = "", sampler_diagnostics = NULL)
}


# function to only read in the post_warmup_draws of named variables
read_draws <- function(csv_files, variables, format = "draws_matrix") {
  read_cmdstan_csv(csv_files, variables = variables, format = format)$post_warmup_draws
}


# 2.Functions to get sampler diagnostics and save a summary---
# note that we did not save the profile files with the csv files, so there wont actually be any profiles

get_model_profiles <- function(model.number, output.dir, csv.subdir) {
  
  model_dir <- get_model_dir(model.number, output.dir, csv.subdir)
  profile_files <- Sys.glob(file.path(model_dir, "*profile*.csv"))
  if (length(profile_files) == 0) {
    warning("No profiling CSVs found for model ", model.number, " -- skipping profile join.")
    return(NULL)
  }
  profile_list <- setNames(lapply(profile_files, data.table::fread), seq_along(profile_files))
  bind_rows(profile_list, .id = "chain")
}

# sampler diagnostics summary--used in save_model_diagnostics function
summarize_sampler_diagnostics <- function(csv_files, nchain) {
  info <- get_model_metadata(csv_files)
  diag_draws <- info$post_warmup_sampler_diagnostics
  diag_df <- posterior::as_draws_df(diag_draws)
  
  max_td <- suppressWarnings(as.numeric(info$metadata$max_treedepth))
  if (is.na(max_td)) max_td <- 10  # cmdstan default if max_treedepth was not listed in data
  
  diagnostics_summary <- diag_df %>%
    group_by(.chain) %>%
    summarise(
      tot.max.treedepth = sum(treedepth__ == max_td),
      n_divergent       = sum(divergent__),
      tot_leapfrog      = sum(n_leapfrog__),
      avg_accept_stat   = mean(accept_stat__),
      # manual E-BFMI 
      ebfmi             = mean(diff(energy__)^2) / stats::var(energy__)
    ) %>%
    rename(chain = .chain) %>%
    mutate(chain = as.character(chain))
  
  list(diagnostics_summary = diagnostics_summary, metadata = info$metadata)
}

save_model_diagnostics <- function(model.number, niter, nwarmup, nchain, output.dir, csv.subdir) {
  
  model.name <- paste0("hierarchical_mort_model_", model.number, "_niter_", niter, "_nchain_", nchain)
  csv_files  <- get_model_csv_files(model.number, output.dir, csv.subdir)
  
  cat(paste("\nReading diagnostics for hierarchical mortality model", model.number,
            "remper correction", remper.cor.vector[j], "\n"))
  
  diag <- summarize_sampler_diagnostics(csv_files, nchain)
  profile_df <- get_model_profiles(model.number, output.dir, csv.subdir)
  
  profile.stats <- if (!is.null(profile_df)) {
    profile_df %>%
      left_join(diag$diagnostics_summary, by = "chain") %>%
      mutate(model.no = model.number, model = model.name)
  } else {
    diag$diagnostics_summary %>% mutate(model.no = model.number, model = model.name)
  }
  
  # add a flag for max tree depth and divergent transitions
  cat( paste("Model",model.number, "hit maxtreedepth", sum(profile.stats$tot.max.treedepth), "times and had", sum(profile.stats$n_divergent), "divergent transitions" ))
  
  
  write.csv(profile.stats,
            paste0(output.dir, "SPCD_stanoutput_cmdstan/diagnostics/sample_diagnostics_", model.name, ".csv"),
            row.names = FALSE)
  
  rm(diag, profile_df, profile.stats); gc()
  invisible(csv_files)
}


# 3. special function for getting beta_alpha_draws ---

# specific to get the beta_alpha draws once
get_beta_alpha_draws <- function(csv_files, model.name, output.dir, save = TRUE, format = "draws_matrix") {
  vars <- c("alpha_SPP", "mu_alpha", "sigma_alpha", "u_beta", "mu_beta", "sigma_theta")
  beta_alpha_samps <- read_draws(csv_files, variables = vars, format = format)
  
  if (save) {
    qs2::qs_save(beta_alpha_samps,
                 paste0(output.dir, "SPCD_stanoutput_cmdstan/betas/u_beta_alpha_samps_", model.name, ".qs"))
  }
  beta_alpha_samps
}

# get only necessary parameters for generated quantities
get_gq_fitted_params_csv <- function(csv_files, model.number, model.name, output.dir, nchain) {
  
  all.vars <- get_model_metadata(csv_files)
  # all stan variables:
  # all.vars$metadata$stan_variables
  
  # get basenames of the variables
  base_names <- sub("\\[.*\\]$", "", all.vars$metadata$model_params)
  
  # keep only variables we will use
  exclude_vars <- c("mM", "log_lik", # both are of length >100,000 each, so we dont need
                    "z_theta", "sigma_theta", "mu_theta", "z_alpha_SPP","theta_beta") 
  
  
  GQ.vars <- all.vars$metadata$model_params[!base_names %in% exclude_vars]
  
  # get the draws_array for use:
  draws_for_GQ <- read_cmdstan_csv(csv_files, variables = GQ.vars, 
                                   format = "draws_array")
  return(draws_for_GQ)
  
}


# 4. function to run LOO stats and save `loo` objects----


get_species_loo <- function(model.number, niter, nwarmup, nchain, output.dir, csv.subdir, json.dir) {
  # set up species names and csv files
  model.name <- paste0("hierarchical_mort_model_", model.number, "_niter_", niter, "_nchain_", nchain)
  csv_files  <- get_model_csv_files(model.number, output.dir, csv.subdir)
  
  # read in model data from a local json directory
  mod.data <- get_model_data(model.number, json.dir)
  
  # only read in log_lik_draws, this will take awhile
  cat(paste0("Reading in log_liklihood draws: Model ", model.number))
  
  log_lik_draws <- read_draws(csv_files, variables = "log_lik", format = "draws_matrix")
  
  n_draws_per_chain <- nrow(log_lik_draws) / nchain
  chain_id <- rep(seq_len(nchain), each = n_draws_per_chain)
  
  for (s in seq_len(mod.data$Nspp)) {
    SPCD.id <- nspp[s, ]$SPCD
    cat(paste("loo results for SPCD", SPCD.id, "species number", s, "\n"))
    
    # get species index
    spec.idx <- which(mod.data$SPP == s)
    spp_log_lik <- log_lik_draws[, spec.idx, drop = FALSE]
    
    spp_r_eff_samps <- loo::relative_eff(exp(spp_log_lik), chain_id = chain_id, cores = 2)
    spp_loo_results <- loo::loo(spp_log_lik, r_eff = spp_r_eff_samps, cores = 2)
    
    qs2::qs_save(spp_loo_results,
                 paste0(output.dir, "SPCD_stanoutput_cmdstan/LOO/LOO_results_", model.name, "_SPCD_", SPCD.id, ".qs"))
    
    rm(spp_loo_results, spp_r_eff_samps, spp_log_lik); gc()
  }
  
  rm(log_lik_draws); gc()
  invisible(NULL)
}

# 5. Get generated quantities predictions species -----
summarize_posteriors <- function(x) {
  c(median = median(x),
    quantile(x, probs = c(0.025)),
    quantile(x, probs = c(0.975)),
    mean = mean(x), sd = sd(x), min = min(x), max = max(x),
    rhat = rhat_basic(x), ess_bulk = ess_bulk(x), ess_tail = ess_tail(x, na.rm = TRUE))
}

run_gen_quants <- function(model.number, niter, nwarmup, nchain, output.dir, csv.subdir,
                           json.dir, gq.mod, nparallel = 4) {
  
  model.name <- paste0("hierarchical_mort_model_", model.number, "_niter_", niter, "_nchain_", nchain)
  csv_files  <- get_model_csv_files(model.number, output.dir, csv.subdir)
  
  mod.data <- get_model_data(model.number, json.dir)
  
  
  # read the beta/alpha draws once per model & save:
  beta_alpha_samps <- get_beta_alpha_draws(csv_files, model.name, output.dir, save = TRUE)
  # plot up beta_alpha_samps traceplots:
  par.names <- colnames(beta_alpha_samps)[1:(ncol(beta_alpha_samps))]
  pdf( paste0(output.dir,"SPCD_stanoutput_cmdstan/images/traceplots_",model.name,".pdf"))
  #specify to save plots in 0x0 grid
  par(mfrow = c(8,3))
  for (p in 1:length(par.names)) {   
    print(mcmc_trace( beta_alpha_samps, pars = par.names[p]))
  }
  dev.off()
  
  
  
  cat(paste0("Extracting parameters needed for Generated Quantities: Model ", model.number))
  # get the draws array for generated quantities
  GQ_post_draws <- posterior::as_draws_array(
    get_gq_fitted_params_csv(csv_files = csv_files, 
                             model.number = model.number, 
                             model.name = model.name, 
                             output.dir = output.dir, 
                             nchain = nchain) 
  )
  
  
  
  # example for all the species alphas and all of the first betas
  # mcmc_trace(beta_alpha_samps, pars = c(paste0("alpha_SPP[",1:17,"]")))
  # mcmc_trace(beta_alpha_samps, pars = c(paste0("u_beta[",1:17,",",1,"]")))
  
  get_species_gen_quantities <- function(s) {
    SPCD.id <- nspp[s, ]$SPCD
    cat(paste("generating posterior predictions for SPCD", SPCD.id, "species number", s, "\n"))
    
    spec.idx     <- mod.data$SPP    %in% s
    spec_rep.idx <- mod.data$SPPrep %in% s
    
    spp.mod.data <- mod.data
    
    # need the full xMfit and N for recaluculating the QR decomposition
    # spp.mod.data$xM_fit     <- as.matrix(mod.data$xM) 
    # spp.mod.data$N_fit <- mod.data$N
    
    # set up species-specific data for gen quants
    spp.mod.data$SPP    <- mod.data$SPP[spec.idx]
    spp.mod.data$xM     <- as.matrix(mod.data$xM[spec.idx, ])
    spp.mod.data$Remper <- mod.data$Remper[spec.idx]
    spp.mod.data$N      <- length(spp.mod.data$SPP)
    spp.mod.data$y      <- mod.data$y[spec.idx]
    
    spp.mod.data$SPPrep    <- mod.data$SPPrep[spec_rep.idx]
    spp.mod.data$xMrep     <- as.matrix(mod.data$xMrep[spec_rep.idx, ])
    spp.mod.data$Remperoos <- mod.data$Remperoos[spec_rep.idx]
    spp.mod.data$Nrep      <- length(spp.mod.data$SPPrep)
    spp.mod.data$ytest     <- mod.data$ytest[spec_rep.idx]
    
    
    
    
    
    
    gen_quants <- gq.mod$generate_quantities(
      fitted_params  = GQ_post_draws,#csv_files,      # feed it cmdstancsv paths
      data           = spp.mod.data,
      seed           = 123,
      parallel_chains = nchain
    )
    
    
    
    
    y_rep_samps     <- gen_quants$draws(variables = "y_rep",  format = "draws_matrix")
    y_hat_samps     <- gen_quants$draws(variables = "y_hat",  format = "draws_matrix")
    pSurv_rep_samps <- gen_quants$draws(variables = "mMrep",  format = "draws_matrix")
    pSurv_hat_samps <- gen_quants$draws(variables = "mMhat",  format = "draws_matrix")
    
    qs2::qs_save(y_rep_samps,     paste0(output.dir, "SPCD_stanoutput_cmdstan/predicted_mort/y_rep_samps_SPCD_", SPCD.id, "_", model.name, ".qs"))
    qs2::qs_save(y_hat_samps,     paste0(output.dir, "SPCD_stanoutput_cmdstan/predicted_mort/y_hat_samps_SPCD_", SPCD.id, "_", model.name, ".qs"))
    qs2::qs_save(pSurv_rep_samps, paste0(output.dir, "SPCD_stanoutput_cmdstan/predicted_mort/pSurv_rep_samps_SPCD_", SPCD.id, "_", model.name, ".qs"))
    qs2::qs_save(pSurv_hat_samps, paste0(output.dir, "SPCD_stanoutput_cmdstan/predicted_mort/pSurv_hat_samps_SPCD_", SPCD.id, "_", model.name, ".qs"))
    
    rm(gen_quants); gc()
    
    # --- convergence / diagnostics table -----------------------------------
    u_betas_alpha.quant <- summarise_draws(beta_alpha_samps, summarize_posteriors)
    y_hat.quant     <- summarise_draws(y_hat_samps,     .cores = 4, summarize_posteriors)
    y_rep.quant     <- summarise_draws(y_rep_samps,     .cores = 4, summarize_posteriors)
    pSurv_hat.quant <- summarise_draws(pSurv_hat_samps, .cores = 4, summarize_posteriors)
    pSurv_rep.quant <- summarise_draws(pSurv_rep_samps, .cores = 4, summarize_posteriors)
    
    convergence.stats <- rbind(u_betas_alpha.quant, y_hat.quant, y_rep.quant,
                               pSurv_hat.quant, pSurv_rep.quant) %>%
      mutate(model.number = model.number, model.type = "Hierarchical",
             SPCD = SPCD.id, remper.correction = remper.cor.vector[j])
    
    write.csv(convergence.stats,
              paste0(output.dir, "SPCD_stanoutput_cmdstan/diagnostics/Rhats_ESS_quantiles_SPCD_", SPCD.id, "_", model.name, ".csv"))
    rm(convergence.stats, u_betas_alpha.quant); gc()
    
    # --- AUC + confusion matrix over posterior draws ------------------------
    actuals     <- spp.mod.data$y
    actuals.oos <- spp.mod.data$ytest
    
    AUC.is.samples.df  <- apply(pSurv_hat_samps, 1, function(p) as.numeric(pROC::auc(actuals, p, quiet = TRUE)))
    AUC.oos.samples.df <- apply(pSurv_rep_samps, 1, function(p) as.numeric(pROC::auc(actuals.oos, p, quiet = TRUE)))
    rm(pSurv_hat_samps, pSurv_rep_samps, y_hat.quant, y_rep.quant, pSurv_hat.quant, pSurv_rep.quant); gc()
    
    preds.is.class <- y_hat_samps == 1
    confusion.is_draws <- data.frame(
      TP_draws = rowSums(preds.is.class[, actuals == 1, drop = FALSE]),
      FP_draws = rowSums(preds.is.class[, actuals == 0, drop = FALSE]),
      TN_draws = rowSums(!preds.is.class[, actuals == 0, drop = FALSE]),
      FN_draws = rowSums(!preds.is.class[, actuals == 1, drop = FALSE])
    ) %>%
      mutate(`True survival rate` = TP_draws / (TP_draws + FN_draws),
             `True mortality rate` = TN_draws / (TN_draws + FP_draws),
             model.number = model.number, type = "in-sample",
             model.type = "Hierarchical", SPCD = SPCD.id)
    
    preds.oos.class <- y_rep_samps == 1
    confusion.oos_draws <- data.frame(
      TP_draws = rowSums(preds.oos.class[, actuals.oos == 1, drop = FALSE]),
      FP_draws = rowSums(preds.oos.class[, actuals.oos == 0, drop = FALSE]),
      TN_draws = rowSums(!preds.oos.class[, actuals.oos == 0, drop = FALSE]),
      FN_draws = rowSums(!preds.oos.class[, actuals.oos == 1, drop = FALSE])
    ) %>%
      mutate(`True survival rate` = TP_draws / (TP_draws + FN_draws),
             `True mortality rate` = TN_draws / (TN_draws + FP_draws),
             model.number = model.number, type = "out-of-sample",
             model.type = "Hierarchical", SPCD = SPCD.id)
    
    confusion.is_draws$AUC  <- AUC.is.samples.df
    confusion.oos_draws$AUC <- AUC.oos.samples.df
    AUC.confusion_draws <- rbind(confusion.is_draws, confusion.oos_draws)
    
    qs2::qs_save(AUC.confusion_draws,
                 paste0(output.dir, "SPCD_stanoutput_cmdstan/AUC/AUC_draws_SPCD_", SPCD.id, "_", model.name, ".qs"))
    
    
    # AUC.confusion_draws %>% group_by(SPCD, type) %>% summarise(median(AUC), 
    #                                                            quantile(AUC, 0.025), 
    #                                                            quantile(AUC, 0.975))
    rm(AUC.is.samples.df, AUC.oos.samples.df, actuals, actuals.oos,
       y_rep_samps, y_hat_samps, preds.is.class, preds.oos.class,
       confusion.is_draws, confusion.oos_draws, AUC.confusion_draws)
    gc()
    invisible(NULL)
  }
  
  cat(paste0("Generating Quantities by species: Model ", model.number))
  
  
  lapply(rev(seq_len(mod.data$Nspp)), get_species_gen_quantities)
  rm(beta_alpha_samps); gc()
  invisible(NULL)
}


# 6. Get model diagnostics for all the species -----


niter <- 1000; nwarmup <- 1000; nchain <- 4

for (m in model.list) {
  save_model_diagnostics(model.number = m, niter = niter, nwarmup = nwarmup,
                         nchain = nchain, output.dir = output.dir, csv.subdir = csv.subdir)
  
  get_species_loo(model.number = m, niter = niter, nwarmup = nwarmup,
                  nchain = nchain, output.dir = output.dir, csv.subdir = csv.subdir, json.dir = json.dir)
  
  run_gen_quants(model.number = m, niter = niter, nwarmup = nwarmup, nchain = nchain,
                 output.dir = output.dir, csv.subdir = csv.subdir, json.dir = json.dir,
                 gq.mod = gq.mod, nparallel = 4)
}
