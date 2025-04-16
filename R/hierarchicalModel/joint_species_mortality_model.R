# run the joint general species model for model 5

library(rstan)
library(MASS)
#library(here)
library(tidyverse)
#library(gt)
library(FIESTA)
library(dplyr)
library(posterior)
library(mltools)
output.folder <- "/home/rstudio/"



# # get the complete species list
nspp <- data.frame(SPCD = c(316, 318, 833, 832, 261, 531, 802, 129, 762,  12, 541,  97, 621, 400, 371, 241, 375))
nspp$Species <- paste(FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD),]$GENUS, FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD),]$SPECIES)
nspp$COMMON <- paste(FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD),]$GENUS, FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD),]$SPECIES)


# read in the test data for all the species
spp.table <- data.frame(SPCD.id = nspp[1:17,]$SPCD, 
                        spp = 1:17, 
                        COMMON = nspp[1:17,]$COMMON)
model.no <- 6

set.seed(22)
# run a giant for loop that runs models 1-9:
for(m in 8:9) {
  
  
  model.no <- m
  
  cat(paste0("running model number ", model.no))
  
  xM.list <-
    xMrep.list <-
    y.list <-
    y.test.list <-
    nSPP.list <-
    nSPP.rep.list <- remper.list <- remper.rep.list <- list()
  
  for (i in 1:17) {
    # SPCD.id
    SPCD.id <- spp.table[i, ]$SPCD.id
    load(
      paste0(
        "SPCD_standata_general_full_standardized_v3/SPCD_",
        SPCD.id,
        "remper_correction_0.5model_",
        model.no,
        ".Rdata"
      )
    ) # load the species code data
    
    
    
    xM.list[[i]] <- mod.data$xM
    xMrep.list[[i]] <- mod.data$xMrep
    y.list[[i]] <- mod.data$y
    y.test.list[[i]] <- mod.data$ytest
    nSPP.list[[i]] <- rep(i, length(mod.data$y))
    nSPP.rep.list[[i]] <- rep(i, length(mod.data$ytest))
    
    remper.list[[i]] <- train.data$remper
    remper.rep.list[[i]] <- test.data$remper
  }
  
  
  
  
  mod.data.full <-  list(
    xM = do.call(rbind, xM.list),
    xMrep = do.call(rbind, xMrep.list),
    y = unlist(y.list),
    ytest = unlist(y.test.list),
    SPP = unlist(nSPP.list),
    SPPrep = unlist(nSPP.rep.list),
    Remper = unlist(remper.list),
    Remperoos = unlist(remper.rep.list)
  )
  mod.data.full$K <- ncol(mod.data.full$xM)
  mod.data.full$N <- length(mod.data.full$y)
  
  mod.data.full$Nspp <- 17
  mod.data.full$Nrep <- length(mod.data.full$ytest)
  saveRDS (
    mod.data.full,
    paste0(
      output.folder,
      "SPCD_stanoutput_joint_v3/all_SPCD_model_",
      model.no,
      ".RDS"
    )
  )
  mod.data.full <-
    readRDS (
      paste0(
        output.folder,
        "SPCD_stanoutput_joint_v3/all_SPCD_model_",
        model.no,
        ".RDS"
      )
    )
  
  
  saveRDS(spp.table,
          paste0(output.folder, "SPCD_stanoutput_joint_v3/spp.table.rds"))
  spp.table <-
    readRDS(paste0(output.folder, "SPCD_stanoutput_joint_v3/spp.table.rds"))
  
  
  mort.spp <- data.frame(y = mod.data.full$y,
                         SPP = mod.data.full$SPP)
  mort.spp <- data.frame(cbind(mort.spp, mod.data.full$xM))
  
  #colnames(mort.spp) <- c("y", "SPP")
  mort.spp %>% group_by(SPP) %>% summarise(n())
  mort.spp %>% group_by(y, SPP) %>% slice_sample(n = 1000)
  
  # get the out of sample data
  mortrep.spp <- data.frame(ytest = mod.data.full$ytest,
                            SPPrep = mod.data.full$SPPrep)
  mortrep.spp <- data.frame(cbind(mortrep.spp, mod.data.full$xMrep))
  #colnames(mort.spp) <- c("y", "SPP")
  mortrep.spp %>% group_by(SPPrep) %>% summarise(n())
  mortrep.spp %>% group_by(ytest, SPPrep) %>% slice_sample(n = 1000)
  
  
  
  
  
  num_cores <-  parallel::detectCores()
  
  options(mc.cores = parallel::detectCores())
  
  # start time for the model:
  start.time <- Sys.time()
  fit.1 <- stan(
    file = "modelcode/reparam_hierarchical_model.stan" ,
    data = mod.data.full,
    seed = 22,
    init = 0,
    iter = 1000,
    chains = 2,
    verbose = FALSE,
    control =  list(max_treedepth = 15, adapt_delta = 0.99),
    #list(adapt_delta = 0.99, stepsize = 0.5, max_treedepth = 15),#, stepsize = 0.01, max_treedepth = 15),
    sample_file = paste0("hierarchical_model_", model.no),
    #adapt_delta = 0.99,
    # pars = c("alpha_species", "beta", "mu_alpha", "mu_beta"))
    pars = c(
      "alpha_SPP",
      "u_beta",
      # the species-specific params
      "mu_alpha",
      "mu_beta",
      "y_rep",
      "mMrep",
      ## out of sample predictions
      "y_hat",
      "mMhat"
    )
  ) ## out of sample predictions
  
  end.time <- Sys.time()
  # Calculate elapsed time in seconds
  elapsed_time <-
    as.numeric(difftime(end.time, start.time, units = "secs"))
  
  # Calculate core hours
  core_hours <-
    elapsed_time * num_cores / 3600  # convert seconds to hours
  remper.correction <- 0.5
  
  time.diag <- data.frame(SPP = mod.data.full$SPP,
                          y = mod.data.full$y) %>% group_by(SPP) %>% summarise(n()) %>% ungroup() %>%
    mutate(total = sum(`n()`)) %>% group_by(SPP) %>% mutate(weight = `n()` /
                                                              total) %>% ungroup() %>%
    mutate(
      total.core.hours = core_hours,
      total.elapsed.time = elapsed_time,
      cores = num_cores,
      model = model.no,
      remper = 0.5
    ) %>% group_by(SPP) %>%
    mutate(core.hours = weight * total.core.hours) %>% rename(`spp` = "SPP") %>%
    left_join(spp.table)
  
  
  
  write.csv(
    time.diag,
    paste0(
      output.folder,
      "SPCD_stanoutput_joint_v3/joint_model_time_diag_SPCD_joint_model_",
      model.no,
      "_remper_",
      remper.correction,
      ".csv"
    )
  )
  
  
  saveRDS(
    fit.1,
    paste0(
      output.folder,
      "SPCD_stanoutput_joint_v3/samples/model_",
      model.no,
      "all_SPCD_remper_correction_0.5.RDS"
    )
  )
  
  joint.samples <- as_draws_df(fit.1)
  
  alpha.p <- subset_draws(joint.samples, variable = "mu_alpha")
  alpha.spp <- subset_draws(joint.samples, variable = "alpha_SPP")
  
  saveRDS(
    alpha.p,
    paste0(
      output.folder,
      "SPCD_stanoutput_joint_v3/alpha.p_model_",
      model.no,
      "_1000samples.rds"
    )
  )
  saveRDS(
    alpha.spp,
    paste0(
      output.folder,
      "SPCD_stanoutput_joint_v3/alpha.spp_model_",
      model.no,
      "_1000samples.rds"
    )
  )
  
  beta.p <-
    subset_draws(joint.samples, variable = "mu_beta")#, chain = 1:2, iteration = 500:1500)
  bet0a.spp <-
    subset_draws(joint.samples, variable = "u_beta")#, chain = 1:2, iteration = 500:1500)
  saveRDS(
    beta.p,
    paste0(
      output.folder,
      "SPCD_stanoutput_joint_v3/beta_model_",
      model.no,
      "_1000samples.rds"
    )
  )
  saveRDS(
    bet0a.spp,
    paste0(
      output.folder,
      "SPCD_stanoutput_joint_v3/u_betas_model_",
      model.no,
      "_1000samples.rds"
    )
  )
  
  # sigmas <-
  #   subset_draws(joint.samples, variable = c("sigma_s", "sigma_aS"))#, chain = 1:2, iteration = 500:1500)
  # saveRDS(
  #   sigmas,
  #   paste0(
  #     output.folder,
  #     "SPCD_stanoutput_joint_v3/sigmas_model_",
  #     model.no,
  #     "_1000samples.rds"
  #   )
  #)
  
  yrep <-
    subset_draws(joint.samples, variable = "y_rep")#, chain = 1:2, iteration = 500:1500)
  yhat <-
    subset_draws(joint.samples, variable = "y_hat")#, chain = 1:2, iteration = 500:1500)
  saveRDS(
    yrep,
    paste0(
      output.folder,
      "SPCD_stanoutput_joint_v3/yrep_model_",
      model.no,
      "_1000samples.rds"
    )
  )
  saveRDS(
    yhat,
    paste0(
      output.folder,
      "SPCD_stanoutput_joint_v3/yhat_model_",
      model.no,
      "_1000samples.rds"
    )
  )
  
  
  mMrep <-
    subset_draws(joint.samples, variable = "mMrep")#, chain = 1:2, iteration = 500:1500)
  mMhat <-
    subset_draws(joint.samples, variable = "mMhat")#, chain = 1:2, iteration = 500:1500)
  saveRDS(
    mMrep,
    paste0(
      output.folder,
      "SPCD_stanoutput_joint_v3/mMrep_model_",
      model.no,
      "_1000samples.rds"
    )
  )
  saveRDS(
    mMhat,
    paste0(
      output.folder,
      "SPCD_stanoutput_joint_v3/mMhat_model_",
      model.no,
      "_1000samples.rds"
    )
  )
  
  # ## read in all the outputs to generate AUC scores
  # mMrep <-
  #   readRDS(
  #     paste0(
  #       output.folder,
  #       "SPCD_stanoutput_joint_v3/mMrep_model_",
  #       model.no,
  #       "_1000samples.rds"
  #     )
  #   )
  # mMhat <-
  #   readRDS(
  #     mMhat,
  #     paste0(
  #       output.folder,
  #       "SPCD_stanoutput_joint_v3/mMhat_model_",
  #       model.no,
  #       "_1000samples.rds"
  #     )
  #   )
  # 
  # yhat <-
  #   readRDS(
  #     paste0(
  #       output.folder,
  #       "SPCD_stanoutput_joint_v3/yhat_model_",
  #       model.no,
  #       "_1000samples.rds"
  #     )
  #   )
  # yrep <-
  #   readRDS(
  #     paste0(
  #       output.folder,
  #       "SPCD_stanoutput_joint_v3/yrep_model_",
  #       model.no,
  #       "_1000samples.rds"
  #     )
  #   )
  # 
  ## Traceplots
  # make and save the traceplots for the betas
  # 33 separate covariates, with 17 species and one population parameter
  
  for (j in 1:length(colnames(mod.data.full$xM))) {
    cov.params = c(paste0("u_beta[", 1:17, ",", j, "]"),
                   paste0("mu_beta[", j, "]"))
    plt.t <- traceplot(fit.1, pars = cov.params)
    ggsave(
      filename = paste0(
        "SPCD_stanoutput_joint_v3/images/traceplot_",
        colnames(mod.data.full$xM)[j],
        "_model_",
        model.no,
        ".png"
      ),
      plot = plt.t,
      width = 8,
      height = 6
    )
    
  }
  
  # plot up the alphas
  cov.params = c(paste0("alpha_SPP[", 1:17, "]"), "mu_alpha")
  plt.t <- traceplot(fit.1, pars = cov.params)
  ggsave(
    filename = paste0(
      "SPCD_stanoutput_joint_v3/images/traceplot_alphas_SPP_model_",
      model.no,
      ".png"
    ),
    plot = plt.t,
    width = 8,
    height = 6
  )
  
  
  
  #dev.off()
  
  
  beta.names <- data.frame(parameter = colnames(mod.data.full$xM),
                           param.no = 1:length(colnames(mod.data.full$xM)))
  
  nvariables <- length(names(fit.1))
  
  nvariables
  
  
  
  # compare the observed status to the predicted status
  yrep.estimates <-
    readRDS(
      paste0(
        output.folder,
        "SPCD_stanoutput_joint_v3/yrep_model_",
        model.no,
        "_1000samples.rds"
      )
    )#%>%
  #subset_draws(chain = 1:3, iteration = 1000:1500)
  
  yrep.quant <-
    summarise_draws(yrep.estimates, median, ~ quantile(.x, probs = c(0.025, 0.975))) %>% rename(`ci.lo` = `2.5%`, `ci.hi` = `97.5%`)
  
  yrep.quant$Mobs <- as.character(mod.data.full$ytest)
  yrep.quant$spp <- mod.data.full$SPPrep
  yrep.quant <- left_join(yrep.quant, spp.table)
  
  ggplot(yrep.quant, aes(as.character(Mobs), y = median, fill = Mobs)) + geom_violin() +
    ylab("Median predicted survival status") +
    xlab("Observed out-of-sample tree survival status") + scale_fill_manual(values =  c("0" = "#a6611a", "1" = "forestgreen")) +
    theme_bw() + theme(legend.position = "none") + facet_wrap( ~ COMMON)
  ggsave(
    height = 8,
    width = 8,
    units = "in",
    paste0(
      output.folder,
      "SPCD_stanoutput_joint_v3/images/Yrepsurv_oosample_median_vs_obs_violin_model_",
      model.no,
      "_species_joint_model_remper_corr_",
      remper.correction,
      ".png"
    )
  )
  
  ggplot(yrep.quant, aes(as.character(Mobs), y = ci.hi, fill = Mobs)) + geom_violin() +
    ylab("97.5% quantile of predicted survival status") +
    xlab("Observed out-of-sample tree survival status") + scale_fill_manual(values =  c("0" = "#a6611a", "1" = "forestgreen")) +
    theme_bw() + theme(legend.position = "none") + facet_wrap( ~ COMMON)
  ggsave(
    height = 8,
    width = 8,
    units = "in",
    paste0(
      output.folder,
      "SPCD_stanoutput_joint_v3/images/Yrepsurv_oossample_95pct_vs_obs_violin_model_",
      model.no,
      "_species_joint_model_remper_corr_",
      remper.correction,
      ".png"
    )
  )
  
  
  # compare the observed status to the predicted status
  yhat.estimates <-
    readRDS(
      paste0(
        output.folder,
        "SPCD_stanoutput_joint_v3/yhat_model_",
        model.no,
        "_1000samples.rds"
      )
    )#%>%
  # subset_draws(chain = 1:3, iteration = 1000:1500)
  yhat.quant <-
    summarise_draws(yhat.estimates, median, ~ quantile(.x, probs = c(0.025, 0.975))) %>% rename(`ci.lo` = `2.5%`, `ci.hi` = `97.5%`)
  
  yhat.quant$Mobs <- as.character(mod.data.full$y)
  yhat.quant$spp <- mod.data.full$SPP
  yhat.quant <- left_join(yhat.quant, spp.table)
  
  ggplot(yhat.quant, aes(as.character(Mobs), y = median, fill = Mobs)) + geom_violin() +
    ylab("Median predicted survival status") +
    xlab("Observed in-sample tree survival status") + scale_fill_manual(values =  c("0" = "#a6611a", "1" = "forestgreen")) +
    theme_bw() + theme(legend.position = "none") + facet_wrap( ~ COMMON)
  ggsave(
    height = 8,
    width = 8,
    units = "in",
    paste0(
      output.folder,
      "SPCD_stanoutput_joint_v3/images/Yhatsurv_insample_median_vs_obs_violin_model_",
      model.no,
      "_species_joint_model_remper_corr_",
      remper.correction,
      ".png"
    )
  )
  
  ggplot(yhat.quant, aes(as.character(Mobs), y = ci.hi, fill = Mobs)) + geom_violin() +
    ylab("97.5% quantile of predicted survival status") +
    xlab("Observed in-sample tree survival status") + scale_fill_manual(values =  c("0" = "#a6611a", "1" = "forestgreen")) +
    theme_bw() + theme(legend.position = "none") + facet_wrap( ~ COMMON)
  ggsave(
    height = 8,
    width = 8,
    units = "in",
    paste0(
      output.folder,
      "SPCD_stanoutput_joint_v3/images/Yhatsurv_insample_95pct_vs_obs_violin_model_",
      model.no,
      "_species_joint_model_remper_corr_",
      remper.correction,
      ".png"
    )
  )
  
  
  # survival probability for in-sample data
  psurv.estimates <-
    readRDS(
      paste0(
        output.folder,
        "SPCD_stanoutput_joint_v3/mMrep_model_",
        model.no,
        "_1000samples.rds"
      )
    )#%>%
  #subset_draws(chain = 1:3, iteration = 1000:1500)
  
  psurv.quant <-
    summarise_draws(psurv.estimates, median, ~ quantile(.x, probs = c(0.025, 0.975))) %>% rename(`ci.lo` = `2.5%`, `ci.hi` = `97.5%`)
  
  #hist(psurv.quant$median)
  
  psurv.quant$Mobs <- as.character(mod.data.full$ytest)
  psurv.quant$spp <- mod.data.full$SPPrep
  psurv.quant <- left_join(psurv.quant, spp.table)
  
  ll.test.pmort <- psurv.quant
  saveRDS(
    ll.test.pmort,
    paste0(
      output.folder,
      "SPCD_stanoutput_joint_v3/ll.train.pmort.RDS"
    )
  )
  
  rm(yhat.quant, rep.quant)
  # survival probability for held out-sample data
  psurv.hat.estimates <-
    readRDS(
      paste0(
        output.folder,
        "SPCD_stanoutput_joint_v3/mMhat_model_",
        model.no,
        "_1000samples.rds"
      )
    )#%>%
  #subset_draws(chain = 1:3, iteration = 1000:1500)
  
  psurv.hat.quant <-
    summarise_draws(psurv.hat.estimates, median, ~ quantile(.x, probs = c(0.025, 0.975))) %>% rename(`ci.lo` = `2.5%`, `ci.hi` = `97.5%`)
  
  
  psurv.hat.quant$Mobs <- as.character(mod.data.full$y)
  psurv.hat.quant$spp <- mod.data.full$SPP
  psurv.hat.quant <- left_join(psurv.hat.quant, spp.table)
  
  
  ll.train.pmort <- psurv.hat.quant
  
  saveRDS(ll.train.pmort,
          "SPCD_stanoutput_joint_v3/ll.train.pmort.RDS")
  
  ############################################################
  # get accuracy of prediction
  # Accuracy
  #ext_fit <- rstan::extract(fit.1)
  # accuracy.is <- mean(as.vector(yrep.quant$median) == mod.data.full$y)
  # accuracy.oos <- mean(as.vector(yhat.quant$median) == mod.data.full$ytest)
  
  # AUC using mltools auc_roc function
  # for in sample data
  actuals = mod.data.full$y
  preds = as.vector(psurv.hat.quant$median)
  auc.is <- auc_roc(preds, actuals)
  
  # get the range of responses for each sample:
  auc.is.list <- list()
  p.surv.hat <-
    psurv.hat.estimates %>% select(paste0("mMhat[", 1:mod.data.full$N, "]"))
  
  # need to split up by species
  
  auc.is.list <-
    lapply(
      1:nrow(psurv.hat.estimates),
      FUN = function(x) {
        spp.auc <- list()
        for (s in unique(mod.data.full$SPP)) {
          spp.idx <- mod.data.full$SPP %in% s
          spp.auc[[s]] <-
            auc_roc( as.vector(as.numeric(p.surv.hat[x, spp.idx])), mod.data.full$y[spp.idx])
        }
        spp.auc.df <- do.call(cbind, spp.auc)
        
        overall.auc <-
          auc_roc( as.vector(as.numeric(p.surv.hat[x, ])), mod.data.full$y)
        overall.auc
        
        all.auc.values <- cbind(overall.auc, spp.auc.df)
        colnames(all.auc.values) <- c("overall", 1:17)
        all.auc.values
      }
    )
  
  AUC.is.df <-
    do.call(rbind, auc.is.list) %>% reshape2::melt() %>% group_by(Var2) %>%
    summarise(
      median = median(value),
      auc.ci.lo = quantile(value, 0.025),
      auc.ci.hi = quantile(value, 0.975)
    ) %>%
    rename("spp" = "Var2")#
  
  
  AUC.is.df$spp <- c(18, 1:17)
  #full.auc.quants <- quantile(AUC.is.df[,1], c(0.025, 0.5, 0.975) )
  
  
  ## for out of sampled data
  actuals = mod.data.full$ytest
  preds = as.vector(psurv.quant$median)
  
  auc.oos <- auc_roc(preds, actuals)
  
  # get the range of responses for each sample:
  auc.oos.list <- list()
  p.surv.rep <-
    psurv.estimates %>% select(paste0("mMrep[", 1:mod.data.full$Nrep, "]"))
  
  
  auc.oos.list <- lapply(
    1:nrow(psurv.estimates),
    FUN = function(x) {
      spp.auc <- list()
      for (s in unique(mod.data.full$SPP)) {
        spp.idx <- mod.data.full$SPPrep %in% s
        spp.auc[[s]] <-
          auc_roc( as.vector(as.numeric(p.surv.rep[x, spp.idx])), mod.data.full$ytest[spp.idx])
      }
      spp.auc.df <- do.call(cbind, spp.auc)
      
      overall.auc <-
        auc_roc( as.vector(as.numeric(p.surv.rep[x, ])), mod.data.full$ytest)
      overall.auc
      
      all.auc.values <- cbind(overall.auc, spp.auc.df)
      colnames(all.auc.values) <- c("overall", 1:17)
      all.auc.values
      
      
    }
  )
  
  
  AUC.oos.df <-
    do.call(rbind, auc.oos.list) %>% reshape2::melt() %>% group_by(Var2) %>%
    summarise(
      median.oos = median(value),
      auc.oos.ci.lo = quantile(value, 0.025),
      auc.oos.ci.hi = quantile(value, 0.975)
    ) %>%
    rename("spp" = "Var2")
  AUC.oos.df$spp <- c(18, 1:17)
  
  #AUC.oos.df <- do.call(rbind, auc.oos.list)
  
  saveRDS(
    AUC.oos.df,
    paste0(
      output.folder,
      "SPCD_stanoutput_joint_v3/AUC_oos_with_uncertainty.rds"
    )
  )
  saveRDS(
    AUC.is.df,
    paste0(
      output.folder,
      "SPCD_stanoutput_joint_v3/AUC_is_with_uncertainty.rds"
    )
  )
  
  # save AUC samples:
  
  
  AUC.oos.samples <-
    do.call(rbind, auc.oos.list) %>% reshape2::melt() %>% group_by(Var2) 
  
  AUC.is.samples <-
    do.call(rbind, auc.is.list) %>% reshape2::melt() %>% group_by(Var2) 
  
  
  saveRDS(
    AUC.oos.samples,
    paste0(
      output.folder,
      "SPCD_stanoutput_joint_v3_model_",model.no,"/AUC_oos_samples.rds"
    )
  )
  saveRDS(
    AUC.is.samples,
    paste0(
      output.folder,
      "SPCD_stanoutput_joint_v3_model_",model.no,"/AUC_is_samples.rds"
    )
  )
  
  joint.table <- data.frame(SPCD.id = 1000,
                            spp = 18,
                            COMMON = "population")
  
  spp.joint.table <- rbind(spp.table, joint.table)
  
  AUC.is.df <- left_join(AUC.is.df, spp.joint.table)
  AUC.oos.df <- left_join(AUC.oos.df, spp.joint.table)
  
  model.assessment.df <- data.frame(
    SPCD = AUC.is.df$SPCD.id,
    COMMON = AUC.is.df$COMMON,
    model =  "hierarchical",
    remper.correction = remper.correction,
    
    auc.insample = auc.is,
    auc.insample.median = AUC.is.df$median,
    auc.insample.lo = AUC.is.df$auc.ci.lo,
    auc.insample.hi = AUC.is.df$auc.ci.hi,
    # out of sample
    auc.oosample = auc.oos,
    auc.oosample.median = AUC.oos.df$median.oos,
    auc.oosample.lo = AUC.oos.df$auc.oos.ci.lo,
    auc.oosample.hi = AUC.oos.df$auc.oos.ci.hi
  )
  
  
  
  write.csv(
    model.assessment.df ,
    paste0(
      output.folder,
      "SPCD_stanoutput_joint_v3/Accuracy_df_model_",
      model.no,
      "_remper_0.5_species_joint_model_remper_corr_",
      remper.correction,
      ".csv"
    ),
    row.names = FALSE
  )
  
  
  # now move all the model outputs into their own folder, and copy to cyverse analysis outputs
  # rename all of the outputs to a folder with the model number
  file.rename(
    "SPCD_stanoutput_joint_v3",
    paste0("SPCD_stanoutput_joint_v3_model_", model.no)
  )
  
  # copy to the data-store output
  system(paste(
    "cp -r",
    paste0("SPCD_stanoutput_joint_v3_model_", model.no),
    "data-store/data/output/"
  ))
  
  # remove the original files to free up space--permanately deletes
  unlink(paste0("SPCD_stanoutput_joint_v3_model_", model.no),
         recursive = T)
  
  
  # clean up the wd:
  rm(list=setdiff(ls(), c("nspp", "spp.table", "m", "output.folder")))
  
  # make sure that we recreate output folders and copy input folders:
  dir.create("SPCD_stanoutput_joint_v3")
  dir.create("SPCD_stanoutput_joint_v3/images")
  dir.create("SPCD_stanoutput_joint_v3/computational_resources")
  dir.create("SPCD_stanoutput_joint_v3/samples")
  dir.create("SPCD_stanoutput_joint_v3/predicted_mort")
  
}





# ####################################
# # redo the accuracy files here:
# ######################################
# # note that we dont need this long term--this was to fix outputs already saved with the incorrect AUC values
# output.folder <- "data-store/data/output/"
# for(m in 1:9) {
#   
#   
#   model.no <- m
#   
#   cat(paste0("running model number ", model.no))
#   
#   mod.data.full <-
#     readRDS (
#       paste0(
#         output.folder,
#         "data-store/data/output/SPCD_stanoutput_joint_v3_model_",model.no,"/all_SPCD_model_",
#         model.no,
#         ".RDS"
#       )
#     )
#   
#   # get the psurv.hat.quant and psurv.quant objects
#   psurv.hat.estimates <-
#     readRDS(
#       paste0(
#         output.folder,
#         "SPCD_stanoutput_joint_v3_model_",model.no,"/mMhat_model_",
#         model.no,
#         "_1000samples.rds"
#       )
#     )
#   
#   psurv.hat.quant <-
#     summarise_draws(psurv.hat.estimates, median, ~ quantile(.x, probs = c(0.025, 0.975))) %>% rename(`ci.lo` = `2.5%`, `ci.hi` = `97.5%`)
#   
#   
#   psurv.hat.quant$Mobs <- as.character(mod.data.full$y)
#   psurv.hat.quant$spp <- mod.data.full$SPP
#   psurv.hat.quant <- left_join(psurv.hat.quant, spp.table)
#   
#   # for in sample data
#   actuals = mod.data.full$y
#   preds = as.vector(psurv.hat.quant$median)
#   auc.is <- auc_roc(preds, actuals)
#   
#   # get the range of responses for each sample:
#   auc.is.list <- list()
#   p.surv.hat <-
#     psurv.hat.estimates %>% select(paste0("mMhat[", 1:mod.data.full$N, "]"))
#   
#   # need to split up by species
#   
#   auc.is.list <-
#     lapply(
#       1:nrow(psurv.hat.estimates),
#       FUN = function(x) {
#         spp.auc <- list()
#         for (s in unique(mod.data.full$SPP)) {
#           spp.idx <- mod.data.full$SPP %in% s
#           spp.auc[[s]] <-
#             auc_roc(as.vector(as.numeric(p.surv.hat[x, spp.idx])), mod.data.full$y[spp.idx])
#         }
#         spp.auc.df <- do.call(cbind, spp.auc)
#         
#         overall.auc <-
#           auc_roc( as.vector(as.numeric(p.surv.hat[x, ])), mod.data.full$y)
#         overall.auc
#         
#         all.auc.values <- cbind(overall.auc, spp.auc.df)
#         colnames(all.auc.values) <- c("overall", 1:17)
#         all.auc.values
#       }
#     )
#   
#   AUC.is.df <-
#     do.call(rbind, auc.is.list) %>% reshape2::melt() %>% group_by(Var2) %>%
#     summarise(
#       median = median(value),
#       auc.ci.lo = quantile(value, 0.025),
#       auc.ci.hi = quantile(value, 0.975)
#     ) %>%
#     rename("spp" = "Var2")#
#   
#   
#   AUC.is.df$spp <- c(18, 1:17)
#   
#   
#   ## for out of sampled data
#   psurv.estimates <-
#     readRDS(
#       paste0(
#         output.folder,
#         "SPCD_stanoutput_joint_v3_model_",model.no,"/mMrep_model_",
#         model.no,
#         "_1000samples.rds"
#       )
#     )
#   
#   psurv.quant <-
#     summarise_draws(psurv.estimates, median, ~ quantile(.x, probs = c(0.025, 0.975))) %>% rename(`ci.lo` = `2.5%`, `ci.hi` = `97.5%`)
#   
#   #hist(psurv.quant$median)
#   
#   psurv.quant$Mobs <- as.character(mod.data.full$ytest)
#   psurv.quant$spp <- mod.data.full$SPPrep
#   psurv.quant <- left_join(psurv.quant, spp.table)
#   
#   
#   actuals = mod.data.full$ytest
#   preds = as.vector(psurv.quant$median)
#   
#   auc.oos <- auc_roc(preds, actuals)
#   
#   # get the range of responses for each sample:
#   auc.oos.list <- list()
#   p.surv.rep <-
#     psurv.estimates %>% select(paste0("mMrep[", 1:mod.data.full$Nrep, "]"))
#   
#   
#   auc.oos.list <- lapply(
#     1:nrow(psurv.estimates),
#     FUN = function(x) {
#       spp.auc <- list()
#       for (s in unique(mod.data.full$SPP)) {
#         spp.idx <- mod.data.full$SPPrep %in% s
#         spp.auc[[s]] <-
#           auc_roc( as.vector(as.numeric(p.surv.rep[x, spp.idx])), mod.data.full$ytest[spp.idx])
#       }
#       spp.auc.df <- do.call(cbind, spp.auc)
#       
#       overall.auc <-
#         auc_roc( as.vector(as.numeric(p.surv.rep[x, ])), mod.data.full$ytest)
#       overall.auc
#       
#       all.auc.values <- cbind(overall.auc, spp.auc.df)
#       colnames(all.auc.values) <- c("overall", 1:17)
#       all.auc.values
#       
#       
#     }
#   )
#   
#   
#   AUC.oos.df <-
#     do.call(rbind, auc.oos.list) %>% reshape2::melt() %>% group_by(Var2) %>%
#     summarise(
#       median.oos = median(value),
#       auc.oos.ci.lo = quantile(value, 0.025),
#       auc.oos.ci.hi = quantile(value, 0.975)
#     ) %>%
#     rename("spp" = "Var2")
#   AUC.oos.df$spp <- c(18, 1:17)
#   
#   #AUC.oos.df <- do.call(rbind, auc.oos.list)
#   
#   saveRDS(
#     AUC.oos.df,
#     paste0(
#       output.folder,
#       "SPCD_stanoutput_joint_v3_model_",model.no,"/AUC_oos_with_uncertainty.rds"
#     )
#   )
#   saveRDS(
#     AUC.is.df,
#     paste0(
#       output.folder,
#       "SPCD_stanoutput_joint_v3_model_",model.no,"/AUC_is_with_uncertainty.rds"
#     )
#   )
#   
#   joint.table <- data.frame(SPCD.id = 1000,
#                             spp = 18,
#                             COMMON = "population")
#   
#   spp.joint.table <- rbind(spp.table, joint.table)
#   
#   AUC.is.df <- left_join(AUC.is.df, spp.joint.table)
#   AUC.oos.df <- left_join(AUC.oos.df, spp.joint.table)
#   
#   model.assessment.df <- data.frame(
#     SPCD = AUC.is.df$SPCD.id,
#     COMMON = AUC.is.df$COMMON,
#     model =  "hierarchical",
#     remper.correction = remper.correction,
#     
#     auc.insample = auc.is,
#     auc.insample.median = AUC.is.df$median,
#     auc.insample.lo = AUC.is.df$auc.ci.lo,
#     auc.insample.hi = AUC.is.df$auc.ci.hi,
#     # out of sample
#     auc.oosample = auc.oos,
#     auc.oosample.median = AUC.oos.df$median.oos,
#     auc.oosample.lo = AUC.oos.df$auc.oos.ci.lo,
#     auc.oosample.hi = AUC.oos.df$auc.oos.ci.hi
#   )
#   
#   
#   
#   write.csv(
#     model.assessment.df ,
#     paste0(
#       output.folder,
#       "SPCD_stanoutput_joint_v3_model_",model.no,"/Accuracy_df_model_",
#       model.no,
#       "_remper_0.5_species_joint_model_remper_corr_",
#       remper.correction,
#       ".csv"
#     ),
#     row.names = FALSE
#   )
#   
#   
#   
# }



# # make a big plot of the predictions
# 
# library(maps)
# library(mapdata)
# states <- map_data("state")
# #9=CT, 25=MA, 33=NH, 23=ME, 50=VT, 44=RI, 42=PA, 39=OH, 54=WV
# state_sub <- filter(states, region %in% c("connecticut","maine","new hampshire","vermont","new york", "new jersey",
#                                           "rhode island","pennsylvania","ohio","west virginia", "massachusetts", "virginia", "delaware", 
#                                           "north carolina", "kentucky", "tennessee", "michigan", "indiana", "district of columbia", "south carolina", 
#                                           "georgia", "maryland"))
# 
# canada <- map_data("worldHires", "Canada")
# 
# # plot distribution 
# ll.test.pmort  <- ll.test.pmort  %>% mutate(`p(mort)` = 1- median) %>%  mutate(Mort.quantiles = cut(`p(mort)`, 
#                                                                                                     breaks = c(0,0.01,0.05, 0.1, 0.2, 0.3, 0.40, 0.50, 1), 
#                                                                                                     include.lowest=TRUE))
# 
# ggplot() +
#   geom_polygon(data = canada, 
#                aes(x=long, y=lat, group = group), 
#                color = "black", fill = "white") +
#   geom_polygon(data = state_sub, 
#                aes(x=long, y=lat, group = group), 
#                color = "black", fill = "white") +
#   
#   geom_point(data = ll.test.pmort, aes( x = LONG_FIADB, y = LAT_FIADB, color = Mort.quantiles), size = 0.5)+theme_bw()+facet_wrap(~COMMON)+
#   scale_color_manual(values = c(
#     "[0,0.01]" = "#d0d1e6", 
#     "(0.01,0.05]" = "yellow",
#     "(0.05,0.1]" = "yellow2",
#     "(0.1,0.2]" = "#fee090", 
#     "(0.22,0.3]" = "#fec44f", 
#     "(0.3,0.4]" = "#fc8d59", 
#     "(0.4,0.5]" =  "#f46d43", 
#     "(0.5,1]" = "#d73027" 
#     
#     
#   ))+
#   coord_sf(xlim = c(-85, -68), ylim = c(35, 48))+theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'lightblue'))
# 
# ggsave(height = 10, width = 18, units = "in", "SPCD_stanoutput_joint/images/species_pmort_hat_distribution_maps.png")
# 
# # plot distribution 
# 
# # plot distribution 
# ll.train.pmort  <- ll.train.pmort  %>% mutate(Mort.quantiles = cut(`p(mort)`, 
#                                                                    breaks = c(0, 0.1, 0.2, 0.3, 0.40, 0.50, 1), 
#                                                                    include.lowest=TRUE))
# 
# ggplot() +
#   geom_polygon(data = canada, 
#                aes(x=long, y=lat, group = group), 
#                color = "black", fill = "white") +
#   geom_polygon(data = state_sub, 
#                aes(x=long, y=lat, group = group), 
#                color = "black", fill = "white") +
#   
#   geom_point(data = ll.train.pmort, aes( x = LONG_FIADB, y = LAT_FIADB, color = Mort.quantiles), size = 0.5)+theme_bw()+facet_wrap(~COMMON)+
#   scale_color_manual(values = c(
#     "[0,0.1]" = "#d0d1e6", 
#     "(0.1,0.2]" = "#fee090", 
#     "(0.22,0.3]" = "#fec44f", 
#     "(0.3,0.4]" = "#fc8d59", 
#     "(0.4,0.5]" =  "#f46d43", 
#     "(0.5,1]" = "#d73027" 
#     
#     
#   ))+
#   coord_sf(xlim = c(-85, -68), ylim = c(35, 48))+theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'lightblue'))
# 
# 
# 
# 
# ggsave(height = 10, width = 18, units = "in", "SPCD_stanoutput/images/species_pmort_insample_distribution_maps.png")
# 
# # plot distribution 
# 
# ggplot() +
#   geom_polygon(data = canada, 
#                aes(x=long, y=lat, group = group), 
#                color = "black", fill = "white") +
#   geom_polygon(data = state_sub, 
#                aes(x=long, y=lat, group = group), 
#                color = "black", fill = "white") +
#   
#   geom_point(data = ll.train.pmort %>% filter(spp %in% spp.table[i,]$SPCD.id), aes( x = LONG_FIADB, y = LAT_FIADB, color = `p(mort)`), size = 0.5)+theme_bw()+facet_wrap(~COMMON)+
#   scale_color_viridis()+
#   coord_sf(xlim = c(-85, -68), ylim = c(35, 48))+theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'lightblue'))
# 
# ggsave(height = 5, width = 9, units = "in", paste0("SPCD_stanoutput/images/species_",spp.table[i,]$SPCD.id,"_pmort_rep_distribution_maps.png"))
# 
# 
# 
# 
# ggplot(data = ll.test.pmort, aes( x = LONG_FIADB, y = LAT_FIADB, color = median))+
#   geom_point()+facet_wrap(~COMMON)
# 
# ggplot(data = ll.train.pmort, aes( x = LONG_FIADB, y = LAT_FIADB, color = median))+
#   geom_point()+facet_wrap(~COMMON)
# 
# 
# ##########################################################
# ###########
# # Use the Loo package to compute PSIS-LOO and check diagnositcs
# # this may take awhile..
# # Extract pointwise log-likelihood
# log_lik_1 <- extract_log_lik(fit.1, merge_chains = FALSE)
# 
# # provide relative effective sample sizes, to bettin estimate PSIS 
# r_eff <- relative_eff(exp(log_lik_1))
# 
# # preferably use more than 0 cores (as many cores as possible)
# # will use value of 'mc.cores' option if cores is not specified
# loo_1 <- loo(log_lik_1, r_eff = r_eff, save_psis = TRUE)
# print(loo_1)
# # save the loo_1 object
# 
# save(log_lik_1, r_eff, loo_1, file = paste0("SPCD_stanoutput_joint/LOO_model_",model.no,"remper.corr_", "_species_", SPCD.id ,"_remper_corr_",remper.correction, ".Rdata"))
# 
# # ###------------------------------------
# # psis <- loo_1$psis_object
# # keep_obs <- sample(1:length(mod.data$y), 100)
# # lw <- weights(psis)
# # 
# # 
# # ppc1 <- bayesplot::ppc_loo_intervals(mod.data$y, 
# #                                      yrep = as.matrix(yrep.estimates), 
# #                                      psis_object = psis, subset = keep_obs, order = "median") 
# # # ppc0 <- rstantools::ppc_loo_pit_overlay(mod.data$y, yrep =  as.matrix(yrep.estimates), lw = lw)
# # # ppc3 <- bayesplot::ppc_loo_pit_qq(mod.data$y, yrep =  as.matrix(yrep.estimates), lw = lw)
# 
