library(tidyverse)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(qs2)
library(jsonlite)
library(pROC)
color_scheme_set("brightblue")

#output.folder <- "/home/rstudio/"
output.dir = "C:/Users/KellyHeilman/Box/01. kelly.heilman Workspace/mortality/Eastern-Mortality/mortality_models/"


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

options(mc.cores = parallel::detectCores())
SPCD.df <- data.frame(SPCD = nspp[1:17, ]$SPCD, 
                      spcd.id = 1:17)
remper.cor.vector <- c(0.5)
#model.number <- 6
model.list <- 1:9


# compile the model once to save time:

hierarchical.updated.file <- file.path(getwd(), "modelcode", "test_reparam_hierarchical.stan")
hierarchical.updated.mod <- cmdstan_model(hierarchical.updated.file )

hierarchical.predict.file <- file.path(getwd(), "modelcode", "predict_hierarchical.stan")
hierarchical.predict <- cmdstan_model(hierarchical.predict.file)



m <- 1
j <- 1

niter <- 100
nwarmup <- 20
nchain <- 4
nparallel <- nchain
print.progress = 10

run.hierarchical.models <-  function( 
                                m, 
                                nparallel,
                                niter, 
                                nwarmup, 
                                nchain, 
                                output.dir, 
                                print.progress = 500){


model.number <- model.list[m]
cat(paste("running hierarchical mortality model ", model.number, " remper correction", remper.cor.vector[j]))
model.name <- paste0("hierarchical_mort_model_",model.number, "_niter_", niter, "_nchain_", nchain)


fit.1 <- hierarchical.updated.mod$sample(
  data = paste0("SPCD_standata_json/hierarchical_data_model_",model.number,".json"), # path to json data files
  seed = 123,
  chains = nchain,
  iter_warmup = nwarmup,
  iter_sampling = niter,
  parallel_chains = nchain,
  #adapt_delta = 0.95, # increased adapt_delta to reduce divergent transition warnings
  init = 0.5, # added to reduce the log(0) probability warnings during initial sampling
  refresh = print.progress # print update every 100 iters
  
)

#fit.1 <- qs2::qs_read(paste0(output.dir,"SPCD_stanoutput_cmdstan/fittedmodels/notpreloaded", model.name, ".qs"))



# this takes the longest to save
cat(paste("\n Saving stan mortality model fit "))
fit.1$sampler_diagnostics()
fit.1$time()
all_draws <- fit.1$draws()

#log_lik_samps <-  fit.1$draws(variables = c("log_lik"), format = "draws_matrix")
#qs2::qs_save(log_lik_samps, paste0(output.dir,"SPCD_stanoutput_cmdstan/LOO/log_lik_samps_", model.name, ".qs"))
qs2::qs_save(fit.1, paste0(output.dir,"SPCD_stanoutput_cmdstan/fittedmodels/",model.name, ".qs"))
#fit.1 <- qs2::qs_read( paste0(output.dir,"SPCD_stanoutput_cmdstan/fittedmodels/",model.name, ".qs"))

beta_alpha_samps <-  fit.1$draws(variables = 
                                   c("alpha_SPP",  "mu_alpha", "sigma_alpha",
                                     "u_beta", "mu_beta", "sigma_beta"), 
                                 format = "draws_matrix")
qs2::qs_save(beta_alpha_samps,paste0(output.dir,"SPCD_stanoutput_cmdstan/betas/u_beta_alpha_samps_", model.name, ".qs"))
#beta_alpha_samps <- qs2::qs_read(paste0(output.dir,"SPCD_stanoutput_cmdstan/betas/u_beta_alpha_samps_", model.name, ".qs"))



qs2::qs_save(fit.1, paste0(output.dir,"SPCD_stanoutput_cmdstan/fittedmodels/", model.name, ".qs"))
#fit.1 <- qs2::qs_read( paste0(output.dir,"SPCD_stanoutput_cmdstan/fittedmodels/", model.name, ".qs"))


# sampler diagnostics -----
# divergent transitions, time per chain, cores, etc
sampler_diag <- fit.1$time()$chains %>%
  mutate(total_allchains = fit.1$time()$total, 
         ncores = nchain,
         nchain = nchain, 
         niter = niter, 
         nwarmup = nwarmup,
         model.number = model.number, 
         model.type = "Hierarchical", 
         SPCD = "Hierarchical", 
         remper.correction = remper.cor.vector[j])%>%
  mutate(num_divergent = fit.1$diagnostic_summary()$num_divergent, 
         num_max_treedepth = fit.1$diagnostic_summary()$num_max_treedepth, 
         ebfmi = fit.1$diagnostic_summary()$ebfmi)

write.csv(sampler_diag, paste0(output.dir, 
                               "SPCD_stanoutput_cmdstan/diagnostics/sample_diagnostics_", 
                               model.name, ".csv"), row.names = FALSE)
rm(sampler_diag)
# Extract posterior draws and save separately ----
#cat(paste("\n Extracting posterior draws"))

# get loo results and save the log_lik, ypredictions, and mmhat and mmrep
# model fit statistics -----

#log_lik_samps <- qs2::qs_read("C:/Users/KellyHeilman/Box/01. kelly.heilman Workspace/mortality/Eastern-Mortality/mortality_models/SPCD_stanoutput_cmdstan/LOO/log_lik_samps_hierarchical_mort_model_1_remper_correction_0.5_niter_1000_nchain_4.qs")

# calculate the log_lik_samps by species for each model:---

# read in model data for comparison
mod.data <- fromJSON(fit.1$data_file())
#parameter_draws_1 <- beta_alpha_samps

# used for data argument to loo_i
mod.data_df <- data.frame(mod.data$y, 
                          mod.data$SPP, 
                          mod.data$Remper,
                          mod.data$xM)
colnames(mod.data_df) <- c("y", "SPP", "Remper", paste0("xM.", 1:mod.data$K))



llfun_logistic_remper <- function(data_i, 
                                  draws, 
                                  log = TRUE) {
  
  x_i <- as.matrix(data_i[,4:ncol(data_i)])
  species_i <- data_i$SPP
  remper_i <- data_i$Remper
  K <- 1:ncol(x_i)
  
  logit_pred <- draws[,paste0("alpha_SPP[",species_i, "]")] + draws[,paste0("u_beta[",species_i,",",K,"]")] %*% t(x_i)
  p_surv_remper <- (1/(1 + exp(-logit_pred)))^remper_i
  
  dbinom(x = data_i$y, 
         size = 1, 
         prob =p_surv_remper, 
         log = log)
}


for(s in 1:mod.data$Nspp){
  
  SPCD.id <- nspp[s,]$SPCD
  cat(paste("loo results for SPCD", SPCD.id, "species number", s, "\n"))
  
  mod.data_species <- mod.data_df %>% filter(SPP == 17)

  spp_r_eff_samps <- loo::relative_eff(llfun_logistic_remper, 
                             log = FALSE, # relative_eff wants likelihood not log-likelihood values
                             chain_id = rep(1:nchain, each = niter), 
                             data = mod.data_species, 
                             draws = beta_alpha_samps, 
                             cores = 2)
  
  
  spp_loo_results <-
    loo::loo(
      llfun_logistic_remper,
      
      cores = 2,
      # these next objects were computed above
      r_eff = spp_r_eff_samps, 
      draws = beta_alpha_samps,
      data = mod.data_species
    )
  qs2::qs_save(spp_loo_results, paste0(output.dir,"SPCD_stanoutput_cmdstan/LOO/LOO_results_", model.name,"_SPCD_", SPCD.id, ".qs"))
  # clean up
  rm(spp_loo_results, spp_r_eff_samps)
  
}


# 
# 
# 
# # compute relative efficiency (this is slow and optional but is recommended to allow 
# # for adjusting PSIS effective sample size based on MCMC effective sample size)
# 
# # the sameloo::loo_compare(loo_species_all, loo_ss_species_all)
# 
# 
# for(s in 1:mod.data$Nspp){
# 
#     SPCD.id <- nspp[s,]$SPCD
#     cat(paste("loo results for SPCD", SPCD.id, "species number", s, "\n"))
#     
#     spec.idx <- mod.data$SPP %in% s
#     
#     # subset the log_liklihood samples for each species
#     spp_log_lik_samps <- log_lik_samps[,spec.idx]
#     
#     # get r_eff for the species
#     spp_r_eff_samps <- loo::relative_eff(exp(spp_log_lik_samps), chain_id = rep(1:nchain, each = niter))
#     # get loo results and save
#     spp_loo_results <-  loo::loo(spp_log_lik_samps, r_eff = spp_r_eff_samps, cores = 2)
#     qs2::qs_save(spp_loo_results, paste0(output.dir,"SPCD_stanoutput_cmdstan/LOO/LOO_results_", model.name,"_SPCD_", SPCD.id, ".qs"))
#     # clean up
#     rm(spp_loo_results, spp_log_lik_samps, spp_r_eff_samps)
# }
# rm(log_lik_samps)
gc()





# PLOTS: parameters -----
# traceplots (betas and alphas)
par.names <- colnames(beta_alpha_samps)[1:(ncol(beta_alpha_samps))]
pdf( paste0(output.dir,"SPCD_stanoutput_cmdstan/images/traceplots_",model.name,".pdf"))
#specify to save plots in 0x0 grid
par(mfrow = c(8,3))
for (p in 1:length(par.names)) {   
  print(mcmc_trace( beta_alpha_samps, pars = par.names[p]))
}
dev.off()




# get summary of the beta parameters
summarize_posteriors <- function(x){
  c(
    median = median(x), 
    quantile(x, probs = c(0.025)), 
    quantile(x, probs = c(0.975)),
    mean = mean(x), 
    sd = sd(x), 
    min = min(x), 
    max = max(x), 
    rhat = rhat_basic(x), 
    ess_bulk = ess_bulk(x), 
    ess_tail = ess_tail(x, na.rm =TRUE)
  )
} 



u_betas_alpha.quant <- summarise_draws(beta_alpha_samps,  summarize_posteriors)
# just plot the species effects once overall
u_betas_alpha.quant <- u_betas_alpha.quant %>% rename("ci.lo" = `2.5%`, 
                                                      "ci.hi" = `97.5%`)
u_betas_alpha.quant$SPCD <- SPCD.id
u_betas_alpha.quant$Covariate <-  u_betas_alpha.quant$variable
# get overlapping zero to color the error bars
u_betas_alpha.quant$`significance` <- ifelse(u_betas_alpha.quant$ci.lo < 0 & u_betas_alpha.quant$ci.hi < 0, "significant", 
                                             ifelse(u_betas_alpha.quant$ci.lo > 0 & u_betas_alpha.quant$ci.hi > 0, "significant", "not overlapping zero"))

ggplot(data = na.omit(u_betas_alpha.quant), aes(x = Covariate, y = median, color = significance))+geom_point()+
  geom_errorbar(data = na.omit(u_betas_alpha.quant), aes(x = Covariate , ymin = ci.lo, ymax = ci.hi, color = significance), width = 0.1)+
  geom_abline(aes(slope = 0, intercept = 0), color = "grey", linetype = "dashed")+theme_bw(base_size = 10)+
  theme( axis.text.x = element_text(angle = 65, hjust = 1), panel.grid  = element_blank(), legend.position = "none")+ylab("Effect on survival")+xlab("Parameter")+
  scale_color_manual(values = c("not overlapping zero"="darkgrey", "significant"="black"))+
  ggtitle(paste0("Species Model Posterior Estimates\nSPCD ",SPCD.id, " model ", model.number ))

ggsave(height = 5, width = 10, units = "in",
       paste0(output.dir, "/images/Estimated_effects_on_survival_",model.name,".png"))


}



# test for model 1
# run.hierarchical.models(
#   m = 1, 
#   nparallel = nparallel,
#   niter = 1000, 
#   nwarmup = 500, 
#   nchain = 4, 
#   output.dir = output.dir, 
#   print.progress = 100)


  run.hierarchical.models(
    m = 7, 
    nparallel = 4,
    niter = 1000, 
    nwarmup = 500, 
    nchain = 4, 
    output.dir = output.dir, 
    print.progress = 100)



lapply(7:9, FUN = function(x){
  run.hierarchical.models(
    m = x, 
    nparallel = 4,
    niter = 1000, 
    nwarmup = 500, 
    nchain = 4, 
    output.dir = output.dir, 
    print.progress = 100)
})






# get the model output for model 1 and generate predictions:
m <- 1
nparallel = 4
niter = 1000 
nwarmup = 500 
nchain = 4

# options:
model.number <- model.list[m]
cat(paste("running hierarchical mortality model ", model.number, " remper correction", remper.cor.vector[j]))
model.name <- paste0("hierarchical_mort_model_",model.number, "_niter_", niter, "_nchain_", nchain)


fit.1 <- qs2::qs_read( paste0(output.dir,"SPCD_stanoutput_cmdstan/fittedmodels/",model.name, ".qs"))

# get the model data
mod.data <-  fromJSON(fit.1$data_file())
#parameter_draws_1 <- beta_alpha_samps

# used for data argument to loo_i
mod.data_df <- data.frame(mod.data$y, 
                          mod.data$SPP, 
                          mod.data$Remper,
                          mod.data$xM)
colnames(mod.data_df) <- c("y", "SPP", "Remper", paste0("xM.", 1:mod.data$K))





# use model$generate_quantities in cmdstan to generate prediction for each species---
# summary function used to summarise posteriors
summarize_posteriors <- function(x){
  c(
    median = median(x), 
    quantile(x, probs = c(0.025)), 
    quantile(x, probs = c(0.975)),
    mean = mean(x), 
    sd = sd(x), 
    min = min(x), 
    max = max(x), 
    rhat = rhat_basic(x), 
    ess_bulk = ess_bulk(x), 
    ess_tail = ess_tail(x, na.rm =TRUE)
  )
} 


#function to generate predicted y and p_surv from hierarchical posteriors:
# --Read in the posteriors from the hierarchical model number (m)
# --generates predictions & does AUC estimation for species (s)
# --lapply(species.num, predictions)
# --save predictions


Run_gen_quants <- function(m, # model number 
                           nparallel = 4, 
                           niter = 1000, 
                           nwarmup = 500, 
                           nchain = 4, 
                           output.dir = output.dir){
  
  # options:
  model.number <- model.list[m]
  cat(paste("running hierarchical mortality model ", model.number, " remper correction", remper.cor.vector[j]))
  model.name <- paste0("hierarchical_mort_model_",model.number, "_niter_", niter, "_nchain_", nchain)
  
  
  fit.1 <- qs2::qs_read( paste0(output.dir,"SPCD_stanoutput_cmdstan/fittedmodels/",model.name, ".qs"))
  
  # get the model data
  mod.data <-  fromJSON(fit.1$data_file())
  #parameter_draws_1 <- beta_alpha_samps
  
  # used for data argument to loo_i
  mod.data_df <- data.frame(mod.data$y, 
                            mod.data$SPP, 
                            mod.data$Remper,
                            mod.data$xM)
  colnames(mod.data_df) <- c("y", "SPP", "Remper", paste0("xM.", 1:mod.data$K))
  
  

  
get_species_gen_quantitites <- function(s){
    SPCD.id <- nspp[s,]$SPCD
    cat(paste("generating posterior predictions of SPCD", SPCD.id, "species number", s, "\n"))
    
    spec.idx <- mod.data$SPP %in% s
    spec_rep.idx <- mod.data$SPPrep %in% s
    
    spp.mod.data <- mod.data
    
    # get insample obs for species
    spp.mod.data$SPP <- mod.data$SPP[spec.idx]
    spp.mod.data$xM <- as.matrix(mod.data$xM[spec.idx,])
    spp.mod.data$Remper <- mod.data$Remper[spec.idx]
    spp.mod.data$N <- length(spp.mod.data$SPP)
    spp.mod.data$y <- mod.data$y[spec.idx]
    
    
    
    spp.mod.data$SPPrep <- mod.data$SPPrep[spec_rep.idx]
    spp.mod.data$xMrep <- as.matrix(mod.data$xMrep[spec_rep.idx,])
    spp.mod.data$Remperoos <- mod.data$Remperoos[spec_rep.idx]
    spp.mod.data$Nrep <- length(spp.mod.data$SPPrep)
    spp.mod.data$ytest <- mod.data$ytest[spec_rep.idx]
    
    #??generate_quantities()
    gen_quants <- hierarchical.predict$generate_quantities(
      fitted_params = fit.1, 
      
      data = spp.mod.data, # path to json data files
      seed = 123,
      parallel_chains = nparallel
    )
    
    
    cat(paste("\n Extracting posterior draws--this can take some time...\n"))
    y_rep_samps <- gen_quants$draws(variables = c("y_rep"), format = "draws_matrix")
    y_hat_samps <-  gen_quants$draws(variables = c("y_hat"), format = "draws_matrix")
   
    
    pSurv_rep_samps <-  gen_quants$draws(variables = c("mMrep"), format = "draws_matrix")
    pSurv_hat_samps <-  gen_quants$draws(variables = c("mMhat"), format = "draws_matrix")
    
    
    # save all to their own objects:
    
    qs2::qs_save(y_rep_samps,paste0(output.dir,"SPCD_stanoutput_cmdstan/predicted_mort/y_rep_samps_SPCD_",SPCD.id, "_", model.name, ".qs"))
    qs2::qs_save(y_hat_samps,paste0(output.dir,"SPCD_stanoutput_cmdstan/predicted_mort/y_hat_samps_SPCD_",SPCD.id, "_", model.name, ".qs"))
    
    qs2::qs_save(pSurv_rep_samps,paste0(output.dir,"SPCD_stanoutput_cmdstan/predicted_mort/pSurv_rep_samps_SPCD_",SPCD.id, "_",  model.name, ".qs"))
    qs2::qs_save(pSurv_hat_samps,paste0(output.dir,"SPCD_stanoutput_cmdstan/predicted_mort/pSurv_hat_samps_SPCD_",SPCD.id, "_",  model.name, ".qs"))
    
    
    #System.time(y_hat.quant <- gen_quants$summary(variables = c("y_rep"), summarize_posteriors))
    
    rm(gen_quants)
    gc()
    
    
    cat(paste("\n Getting model diagnostics and summaries"))
    # convergence statistics summaries ----
    beta_alpha_samps <-  fit.1$draws(variables = 
                                         c("alpha_SPP",  "mu_alpha", "sigma_alpha",
                                           "u_beta", "mu_beta", "sigma_beta"), 
                                       format = "draws_matrix")
    
    
    # get predicted draws for survival (0,1), annual and remper survival probabilities for in-sample ("hat") and held-out ("rep")
    u_betas_alpha.quant <- summarise_draws(beta_alpha_samps,  summarize_posteriors)
    
  
    
    system.time(y_hat.quant <- summarise_draws(y_hat_samps, .cores = 4,  summarize_posteriors))
    
    

  
   
    y_rep.quant <-  summarise_draws(y_rep_samps, .cores = 4,  summarize_posteriors)
    
    # pSurv_ predicted survival probabilities over the remper for in sample and out of sample
    pSurv_hat.quant <-  summarise_draws(pSurv_hat_samps, .cores = 4,  summarize_posteriors)
    pSurv_rep.quant <-  summarise_draws(pSurv_rep_samps, .cores = 4,  summarize_posteriors)
    
    
    
    
    convergence.stats <- rbind(u_betas_alpha.quant, 
                               y_hat.quant, 
                               y_rep.quant, 
                               pSurv_hat.quant, 
                               pSurv_rep.quant)%>%
      mutate(model.number = model.number, 
             model.type = "Hierarchical",
             SPCD = SPCD.id, 
             remper.correction = 0.5)
    
    
    write.csv(convergence.stats, 
              paste0(output.dir, "SPCD_stanoutput_cmdstan/diagnostics/Rhats_ESS_quantiles_SPCD_",SPCD.id, "_", model.name, ".csv"))
    
    
    rm(convergence.stats)
    gc()
    
    
    cat(paste("\n Calculating AUC scores"))
    # AUC scores
    # for in sample data
    actuals = spp.mod.data$y
    actuals.oos = spp.mod.data$ytest
    
    # use the posterior draws from pSurv to estimate draws of AUC responses for each sample:
    
    AUC.is.samples.df <- apply(pSurv_rep_samps , 
                               MARGIN = 1, function(prob){
                                 as.numeric(pROC::auc(actuals, prob, quiet = TRUE))
                               })
    
    AUC.oos.samples.df <- apply(pSurv_hat_samps, 
                                MARGIN = 1, function(prob){
                                  as.numeric(pROC::auc(actuals.oos, prob, quiet = TRUE))
                                })
    
    rm(pSurv_hat_samps, pSurv_rep_samps)
    gc()
    
    # get the confusion matrix over the draws--true postives/negatives, false positives/negatives
    # in-sample:
    #preds.is <-  #%>% select(-.chain, -.iteration, -.draw) %>% as.matrix()
    preds.is.class <- y_rep_samps == 1
    
    
    confusion.is_draws <-   data.frame(      
      TP_draws = rowSums(preds.is.class[,actuals == 1, drop = FALSE]),
      FP_draws = rowSums(preds.is.class[,actuals == 0, drop = FALSE]),
      
      TN_draws = rowSums(!preds.is.class[,actuals == 0, drop = FALSE]),
      FN_draws = rowSums(!preds.is.class[,actuals == 1, drop = FALSE])
    )%>%
      mutate(`True survival rate` = TP_draws/(TP_draws+FN_draws), 
             `True mortality rate` = TN_draws/(TN_draws+FP_draws), 
             model.number = model.number, 
             type = "in-sample", 
             model.type = "Hierarchical", 
             SPCD = SPCD.id
      )
    
    
    # out-of-sample:
    preds.oos.class <- y_hat_samps == 1
    
    
    
    confusion.oos_draws <-   data.frame(      
      TP_draws = rowSums(preds.oos.class[,actuals.oos == 1, drop = FALSE]),
      FP_draws = rowSums(preds.oos.class[,actuals.oos == 0, drop = FALSE]),
      
      TN_draws = rowSums(!preds.oos.class[,actuals.oos == 0, drop = FALSE]),
      FN_draws = rowSums(!preds.oos.class[,actuals.oos == 1, drop = FALSE])
    )%>%
      mutate(`True survival rate` = TP_draws/(TP_draws+FN_draws), 
             `True mortality rate` = TN_draws/(TN_draws+FP_draws), 
             model.number = model.number, 
             type = "out-of-sample", 
             model.type = "Hierarchical", 
             SPCD = SPCD.id
      )
    
    # combine the AUCs and confusion matrices together with draws;
    confusion.is_draws$AUC <- AUC.is.samples.df
    confusion.oos_draws$AUC <- AUC.oos.samples.df 
    
    AUC.confusion_draws <- rbind(confusion.is_draws, confusion.oos_draws)
    
    qs2::qs_save(AUC.confusion_draws, paste0(output.dir,"SPCD_stanoutput_cmdstan/AUC/AUC_draws_SPCD_",SPCD.id,"_", model.name, ".qs"))
    
    
    gc()

}

lapply(17:1, get_species_gen_quantitites)

}


# test for model 1
Run_gen_quants(
    m = 1, 
    nparallel = nparallel,
    niter = 1000, 
    nwarmup = 500, 
    nchain = 4, 
    output.dir = output.dir, 
    print.progress = 100)

# run for m = 2:6


# summarise the hierarchical models 1-6 so far by species----



