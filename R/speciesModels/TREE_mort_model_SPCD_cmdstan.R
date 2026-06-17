library(tidyverse)
library(cmdstanr)
library(posterior)
library(bayesplot)
color_scheme_set("brightblue")
#check_cmdstan_toolchain()
#cmdstan_path()

output.dir = "C:/Users/KellyHeilman/Box/01. kelly.heilman Workspace/mortality/Eastern-Mortality/mortality_models/"

nspp <- data.frame(SPCD = c(316, 318, 833, 832, 261, 531, 802, 129, 762,  12, 541,  97, 621, 400, 371, 241, 375))
nspp$Species <- paste(FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD),]$GENUS, FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD),]$SPECIES)
nspp$COMMON <- paste(FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD),]$GENUS, FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD),]$SPECIES)



options(mc.cores = parallel::detectCores())
SPCD.df <- data.frame(SPCD = nspp[1:17, ]$SPCD, 
                      spcd.id = 1:17)
remper.cor.vector <- c(0.5)
#model.number <- 6
model.list <- 1:9


# compile the model once to save time:
species.file <- file.path(getwd(), "modelcode", "mort_model_general.stan")
species.mod <- cmdstan_model(species.file)

i <- 17
m <- 6
j <- 1

niter <- 1000
nwarmup <- 500
nchain <- 4
nparallel <- 4

#TODO: multi-threading with reduce_sum

#for(i in 16:1){# run for each of the 17 species
  for(m in 1:9){ 

#for(m in 1:9){ 
  
  model.number <- model.list[m]

  common.name <- nspp[1:17, ] %>% filter(SPCD %in% SPCD.df[i,]$SPCD) %>% dplyr::select(COMMON)
  SPCD.id <- SPCD.df[i,]$SPCD
    #for(m in 1:length(model.list)){  # run each of the 9 models
    
    
   # for (j in 1:length(remper.cor.vector)){ # for the growth only model explore the consequences of other assumptions about remeasurement period
      cat(paste("running stan mortality model ", model.number, " for SPCD", SPCD.df[i,]$SPCD, common.name$COMMON, " remper correction", remper.cor.vector[j]))
     
      # 
       fit.1 <- species.mod$sample(
        data = paste0("SPCD_standata_json/SPCD_",SPCD.id,"remper_correction_",
                      remper.cor.vector[j],"model_",model.number,".json"), # path to json data files
        seed = 123,
        chains = nchain,
        iter_warmup = nwarmup,
        iter_sampling = niter,
        parallel_chains = nparallel,
        adapt_delta = 0.95,
        init = 0.5,
        refresh = 50 # print update every 500 iters

      )
     
       
      # save the fit object for later 
      model.name <- paste0("mort_model_",model.number,"_SPCD_", SPCD.id, 
                           "_remper_correction_", remper.cor.vector[j], "_niter_", niter, "_nchain_", nchain)
      remp.cor <- remper.cor.vector[j]
      remper.correction <- remper.cor.vector[j]
      
     fit.1$save_object(file = paste0(output.dir,"SPCD_stanoutput_cmdstan/", model.name, ".rds"))
      
      
      #fit.1 <- readRDS(paste0(output.dir,"SPCD_stanoutput_cmdstan/", model.name, ".rds"))
      # get diagnostics like R-hat and ESS
       
      
      # Extract posterior draws and save separately ----
      beta_alpha_samps <-  fit.1$draws(variables = c("alpha_SPP", "u_beta"), format = "df")
      
      # get loo results and save the log_lik, ypredictions, and mmhat and mmrep
    
      log_lik_samps <-  fit.1$draws(variables = c("log_lik"), format = "df")
      
      
      
      y_rep_samps <-  fit.1$draws(variables = c("y_rep"), format = "df")
      y_hat_samps <-  fit.1$draws(variables = c("y_hat"), format = "df")

      pSurv_rep_samps <-  fit.1$draws(variables = c("mMrep"), format = "df")
      pSurv_hat_samps <-  fit.1$draws(variables = c("mMhat"), format = "df")

      pSannual_rep_samps <-  fit.1$draws(variables = c("pSannualrep"), format = "df")
      pSannual_hat_samps <-  fit.1$draws(variables = c("pSannualhat"), format = "df")
      
      # save all to their own objects:
      saveRDS(log_lik_samps,paste0(output.dir,"SPCD_stanoutput_cmdstan/log_lik_samps_", model.name, ".rds"))
      saveRDS(beta_alpha_samps,paste0(output.dir,"SPCD_stanoutput_cmdstan/u_beta_alpha_samps_", model.name, ".rds"))
      
      saveRDS(y_rep_samps,paste0(output.dir,"SPCD_stanoutput_cmdstan/y_rep_samps_", model.name, ".rds"))
      saveRDS(y_hat_samps,paste0(output.dir,"SPCD_stanoutput_cmdstan/y_hat_samps_", model.name, ".rds"))
      
      saveRDS(pSurv_rep_samps,paste0(output.dir,"SPCD_stanoutput_cmdstan/pSurv_rep_samps_", model.name, ".rds"))
      saveRDS(pSurv_hat_samps,paste0(output.dir,"SPCD_stanoutput_cmdstan/pSurv_hat_samps_", model.name, ".rds"))
      
      saveRDS(pSannual_rep_samps,paste0(output.dir,"SPCD_stanoutput_cmdstan/pSannual_rep_samps_", model.name, ".rds"))
      saveRDS(pSannual_hat_samps,paste0(output.dir,"SPCD_stanoutput_cmdstan/pSannual_hat_samps_", model.name, ".rds"))
      
      
      # sampler diagnostics -----
          # divergent transitions, time per chain, cores, etc
      sampler_diag <- fit.1$time()$chains %>%
        mutate(total_allchains = fit.1$time()$total, 
               ncores = nparallel,
               nchain = nchain, 
               niter = niter, 
               nwarmup = nwarmup,
               model.number = model.number, 
               model.type = "Species", 
               SPCD = SPCD.id, 
               remper.correction = remper.cor.vector[j])%>%
        mutate(num_divergent = fit.1$diagnostic_summary()$num_divergent, 
               num_max_treedepth = fit.1$diagnostic_summary()$num_max_treedepth, 
               ebfmi = fit.1$diagnostic_summary()$ebfmi)
      
      write.csv(sampler_diag, paste0(output.dir, 
                                     "SPCD_stanoutput_cmdstan/sample_diagnostics_", 
                                     model.name, ".csv"), row.names = FALSE)
      
          #convergence statistics & summary for all parameters we do inference on:
      u_betas_alpha.quant <-  fit.1$summary(variables = c("alpha_SPP", "u_beta"), 
                                       median,           
                                       ~quantile(.x, probs = c(0.025, 0.975)), 
                                       mean, sd, min, max, 
                                       rhat, ess_bulk, ess_tail)     
      
      # get predicted draws for survival (0,1), annual and remper survival probabilities for in-sample ("hat") and held-out ("rep")
      y_hat.quant <-  fit.1$summary(variables = c("y_hat"), 
                                    median,           
                                    ~quantile(.x, probs = c(0.025, 0.975)), 
                                    mean, sd, min, max, 
                                    rhat, ess_bulk, ess_tail)
      
      y_rep.quant <-  fit.1$summary(variables = c( "y_rep"), 
                                    median,           
                                    ~quantile(.x, probs = c(0.025, 0.975)), 
                                    mean, sd, min, max, 
                                    rhat, ess_bulk, ess_tail)
      
      # pSurv_ predicted survival probabilities over the remper for in sample and out of sample
      pSurv_hat.quant <-  fit.1$summary(variables = c("mMhat"), 
                                        median,           
                                        ~quantile(.x, probs = c(0.025, 0.975)), 
                                        mean, sd, min, max, 
                                        rhat, ess_bulk, ess_tail) 
      pSurv_rep.quant <-  fit.1$summary(variables = c( "mMrep"), 
                                        median,           
                                        ~quantile(.x, probs = c(0.025, 0.975)), 
                                        mean, sd, min, max, 
                                        rhat, ess_bulk, ess_tail) 
      
      
      # pSannual_rep predicted annual survival probabilities for in sample and out of sample trees
      pSannual_hat.quant <-  fit.1$summary(variables = c("pSannualhat"), 
                                        median,           
                                        ~quantile(.x, probs = c(0.025, 0.975)), 
                                        mean, sd, min, max, 
                                        rhat, ess_bulk, ess_tail) 
      pSannual_rep.quant <-  fit.1$summary(variables = c( "pSannualrep"), 
                                           median,           
                                           ~quantile(.x, probs = c(0.025, 0.975)), 
                                           mean, sd, min, max, 
                                           rhat, ess_bulk, ess_tail) 
      
      
      convergence.stats <- rbind(u_betas_alphas, 
                                 y_hat.quant, 
                                 y_rep.quant, 
                                 pSurv_hat.quant, 
                                 pSurv_rep.quant, 
                                 pSannual_hat.quant, 
                                 pSannual_rep.quant)%>%
        mutate(model.number = model.number, 
               model.type = "Species",
               SPCD = SPCD.id, 
               remper.correction = remper.cor.vector[j])
     
      
      write.csv(convergence.stats, 
                paste0(output.dir, "SPCD_stanoutput_cmdstan/Rhats_ESS_quantiles_", model.name, ".csv"))
      
      
      # model fit statistics -----
          # loo results
      loo_results <-  fit.1$loo()
      saveRDS(loo_results, paste0(output.dir,"SPCD_stanoutput_cmdstan/LOO_results_", model.name, ".rds"))
      
      # read in model data for comparison
      mod.data <- fromJSON(fit.1$data_file())
      
          # AUC scores
      # for in sample data
      actuals = mod.data$y
      preds = as.vector(pSurv_rep.quant$median)
      auc.is <- pROC::auc( actuals, preds) %>% as.numeric()
      
      # assuming that a probability <= 0.8 results in mortality
      confusion.is <- data.frame(obs_outcome = actuals,
                                 prob = preds, 
                                 type = "in-sample") %>%
        mutate(pred_outcome = ifelse(prob > 0.8, 1, 0))%>%
          mutate(TP = ifelse(pred_outcome == 1 & obs_outcome == 1, 1, 0),
                 FP = ifelse(pred_outcome == 1 & obs_outcome == 0, 1, 0),
                 TN = ifelse(pred_outcome == 0 & obs_outcome == 0, 1, 0),
                 FN = ifelse(pred_outcome == 0 & obs_outcome == 1, 1, 0))%>%
        group_by(type)%>%
          summarise(`True Survival` = sum(TP),
                    `False Survival` = sum(FP),
                    `True Mortality` = sum(TN),
                    `False Mortality` = sum(FN))%>%
          mutate(`True survival rate` = `True Survival`/(`True Survival`+`False Mortality`),
                 `True mortality rate` = `True Mortality`/(`True Mortality`+`False Survival`))

      # for out of sample data
      actuals.oos = mod.data$ytest
      preds.oos = as.vector(pSurv_hat.quant$median)
      auc.oos <- pROC::auc( actuals.oos, preds.oos) %>% as.numeric()
      
      
      confusion.oos <- data.frame(obs_outcome = actuals.oos,
                                 prob = preds.oos, 
                                 type = "out-of-sample") %>%
        mutate(pred_outcome = ifelse(prob > 0.8, 1, 0))%>%
        mutate(TP = ifelse(pred_outcome == 1 & obs_outcome == 1, 1, 0),
               FP = ifelse(pred_outcome == 1 & obs_outcome == 0, 1, 0),
               TN = ifelse(pred_outcome == 0 & obs_outcome == 0, 1, 0),
               FN = ifelse(pred_outcome == 0 & obs_outcome == 1, 1, 0))%>%
        group_by(type)%>%
        summarise(`True Survival` = sum(TP),
                  `False Survival` = sum(FP),
                  `True Mortality` = sum(TN),
                  `False Mortality` = sum(FN))%>%
        mutate(`True survival rate` = `True Survival`/(`True Survival`+`False Mortality`),
               `True mortality rate` = `True Mortality`/(`True Mortality`+`False Survival`))
      
    
      # get the range of responses for each sample:
    
      AUC.is.samples.df <- apply(pSurv_rep_samps %>% select(-.chain, -.iteration, -.draw), 
                         MARGIN = 1, function(prob){
        as.numeric(pROC::auc(actuals, prob))
      })
      
      AUC.oos.samples.df <- apply(pSurv_hat_samps %>% select(-.chain, -.iteration, -.draw), 
                         MARGIN = 1, function(prob){
                           as.numeric(pROC::auc(actuals.oos, prob))
                         })
      
      # get the confusion matrix over the draws--true postives/negatives, false positives/negatives
      # in-sample:
      preds.is <- y_rep_samps %>% select(-.chain, -.iteration, -.draw) %>% as.matrix()
      preds.is.class <- preds.is == 1
      
      ySurv = actuals == 1
      yMort = actuals == 0
      
     
   confusion.is_draws <-   data.frame(      
      TP_draws = rowSums(preds.is.class[,ySurv, drop = FALSE]),
      FP_draws = rowSums(preds.is.class[,yMort, drop = FALSE]),
      
      TN_draws = rowSums(!preds.is.class[,yMort, drop = FALSE]),
      FN_draws = rowSums(!preds.is.class[,ySurv, drop = FALSE])
      )%>%
        mutate(`True survival rate` = TP_draws/(TP_draws+FN_draws), 
               `True mortality rate` = TN_draws/(TN_draws+FP_draws), 
               model.number = model.number, 
               type = "in-sample", 
               model.type = "Species", 
               SPCD = SPCD.id
               )
   
   
   # out-of-sample:
   preds.oos <- y_hat_samps %>% select(-.chain, -.iteration, -.draw) %>% as.matrix()
   preds.oos.class <- preds.oos == 1
   
   ySurv = actuals.oos == 1
   yMort = actuals.oos == 0
   
   
   confusion.oos_draws <-   data.frame(      
     TP_draws = rowSums(preds.oos.class[,ySurv, drop = FALSE]),
     FP_draws = rowSums(preds.oos.class[,yMort, drop = FALSE]),
     
     TN_draws = rowSums(!preds.oos.class[,yMort, drop = FALSE]),
     FN_draws = rowSums(!preds.oos.class[,ySurv, drop = FALSE])
   )%>%
     mutate(`True survival rate` = TP_draws/(TP_draws+FN_draws), 
            `True mortality rate` = TN_draws/(TN_draws+FP_draws), 
            model.number = model.number, 
            type = "out-of-sample", 
            model.type = "Species", 
            SPCD = SPCD.id
     )
      
     # hist(confusion.oos_draws$`True survival rate`)
     # hist(confusion.oos_draws$`True mortality rate`)
     
     
    # combine the AUCs and confusion matrices together with draws;
   confusion.is_draws$AUC <- AUC.is.samples.df
   confusion.oos_draws$AUC <- AUC.oos.samples.df 
   
   AUC.confusion_draws <- rbind(confusion.is_draws, confusion.oos_draws)
   
   saveRDS(AUC.confusion_draws, paste0(output.dir,"SPCD_stanoutput_cmdstan/AUC_draws_", model.name, ".rds"))
     
   # # save a summary with 95% CI  
   #  AUC.confusion_draws %>% group_by(model.number, type, SPCD)%>%
   #    summarise(AUC_median = median(AUC), 
   #              AUC_ci.lo = quantile(AUC, 0.025), 
   #              AUC_ci.hi = quantile(AUC, 0.975), 
   #              
   #              True_surv_rate = median(`True survival rate`), 
   #              True_surv_rate_ci.lo = quantile(`True survival rate`, 0.025),
   #              True_surv_rate_ci.hi = quantile(`True survival rate`, 0.975),
   #              
   #              True_mort_rate = median(`True mortality rate`), 
   #              True_mort_rate_ci.lo = quantile(`True mortality rate`, 0.025),
   #              True_mort_rate_ci.hi = quantile(`True mortality rate`, 0.975))
      
    
  
      # PLOTS: parameters -----
      # traceplots (betas and alphas)
   par.names <- colnames(beta_alpha_samps)[1:(ncol(beta_alpha_samps)-3)]
   pdf( paste0(output.dir,"SPCD_stanoutput_cmdstan/images/traceplots_",model.name,".pdf"))
   #specify to save plots in 0x0 grid
   par(mfrow = c(8,3))
   for (p in 1:length(par.names)) {   
     print(mcmc_trace( beta_alpha_samps, pars = par.names[p]))
   }
   dev.off()
   
       # correlation plots & summary to evaluate any non-identifiability in posteriors
   # post.correlations <- cor(beta_alpha_samps %>% select(-.chain, -.iteration, -.draw))
   # post.correlations[upper.tri(post.correlations)] <- NA
   # post.correlations %>% reshape2::melt(.)%>% filter(!is.na(value))%>% filter(!value == 1)%>%
   #   arrange(desc(value))%>% filter(value >0.5)
   # 
   # pairs(beta_alpha_samps[,par.names])#, 
   #            pars = c(par.names),
   #            off_diag_args = list(size = 0.75, alpha = 0.5))
      # dotplots of betas & alphas with real parameter names
   u_betas_alpha.quant <- u_betas_alpha.quant %>% rename("ci.lo" = `2.5%`, 
                                                         "ci.hi" = `97.5%`)
   u_betas_alpha.quant$SPCD <- SPCD.id
   u_betas_alpha.quant$Covariate <-  u_betas_alpha.quant$variable
   
   # reorder by the value of the covariate
  # u_betas_alpha.quant <- u_betas_alpha.quant %>% arrange(by = median)
   #u_betas_alpha.quant$Covariate <- factor(u_betas_alpha.quant$Covariate, levels = u_betas_alpha.quant$Covariate)
   
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
   
      
      # PLOTS: model fit -----
      # predicted vs observed
      # yrep
      
      
      # source("R/speciesModels/SPCD_plot_stan.R")
      rm(fit.1)
#     }
   }
#}

# do the summaries on diagnostics, time, etc
#fit.1 <- readRDS(paste0(output.dir,"SPCD_stanoutput_cmdstan/", model.name, ".rds"))



# get convergence stats---
convergence.files <- list.files(path = paste0(output.dir,"SPCD_stanoutput_cmdstan/"), 
                               pattern = paste0("Rhats_ESS"), full.names = TRUE)
convergence_all <- do.call(rbind, lapply(convergence.files, read.csv))

convergence.summary<- convergence_all %>%  #group_by(SPCD, model.type, model.number, chain_id)%>%
  mutate(chain_core_hours = total*(1/60)*(1/60)*(ncores/nchain)) %>%
  group_by(SPCD, model.type, model.number)%>%
  summarise(total_core_hours = sum(chain_core_hours), 
            total_iter = sum(niter), 
            total_warmup = sum(nwarmup))


diagnostic.files <- list.files(path = paste0(output.dir,"SPCD_stanoutput_cmdstan/"), pattern = paste0("diagnostics"), full.names = TRUE)
diags_all <- do.call(rbind, lapply(diagnostic.files, read.csv))
diag.summary<- diags_all %>%  #group_by(SPCD, model.type, model.number, chain_id)%>%
  mutate(chain_core_hours = total*(1/60)*(1/60)*(ncores/nchain)) %>%
  group_by(SPCD, model.type, model.number)%>%
  summarise(total_core_hours = sum(chain_core_hours), 
            total_iter = sum(niter), 
            total_warmup = sum(nwarmup), 
            
            total_divergent = sum(num_divergent),
            total_max_treedepth = sum(num_max_treedepth))
  

# get loo results and compare---
loo.files <- list.files(path = paste0(output.dir,"SPCD_stanoutput_cmdstan/"), pattern = paste0("LOO_results"), full.names = TRUE)
loo_results_all <- lapply(loo.files, readRDS)
loo::loo_compare(loo_results_all) # best fit based on loo elpd differences


# get auc results and compare---
AUC.files <- list.files(path = paste0(output.dir,"SPCD_stanoutput_cmdstan/"), pattern = paste0("AUC"), full.names = TRUE)
AUC_results_all <- do.call(rbind, lapply(AUC.files, function(file.name){
  AUC.confusion_draws <- readRDS(file.name)
   AUC.confusion_draws |> group_by(model.number, type, SPCD)%>%
     summarise(AUC_median = median(AUC),
               AUC_ci.lo = quantile(AUC, 0.025),
               AUC_ci.hi = quantile(AUC, 0.975),

               True_surv_rate = median(`True survival rate`),
               True_surv_rate_ci.lo = quantile(`True survival rate`, 0.025),
               True_surv_rate_ci.hi = quantile(`True survival rate`, 0.975),

               True_mort_rate = median(`True mortality rate`),
               True_mort_rate_ci.lo = quantile(`True mortality rate`, 0.025),
               True_mort_rate_ci.hi = quantile(`True mortality rate`, 0.975))
}))

ggplot(data = AUC_results_all)+
  geom_pointrange(aes(x = model.number, y = AUC_median, ymin = AUC_ci.lo, ymax = AUC_ci.hi))+
  facet_wrap(~type)

ggplot(data = AUC_results_all)+
  geom_pointrange(aes(x = model.number, y = True_surv_rate, ymin = True_surv_rate_ci.lo, ymax = True_surv_rate_ci.hi))+
  facet_wrap(~type)
 
ggplot(data = AUC_results_all)+
  geom_pointrange(aes(x = model.number, y = True_mort_rate, ymin = True_mort_rate_ci.lo, ymax = True_mort_rate_ci.hi))+
  facet_wrap(~type)
