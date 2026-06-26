library(tidyverse)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(qs2)
library(jsonlite)
library(pROC)
color_scheme_set("brightblue")
#check_cmdstan_toolchain()
#cmdstan_path()

output.dir = "C:/Users/KellyHeilman/Box/01. kelly.heilman Workspace/mortality/Eastern-Mortality/mortality_models/"

nspp <- data.frame(SPCD = c(316, 318, 833, 832, 261, 531, 802, 129, 762,  12, 541,  97, 621, 400, 371, 241, 375))
nspp$Species <- paste(FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD),]$GENUS, FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD),]$SPECIES)
nspp$COMMON_NAME <- FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD),]$COMMON_NAME



options(mc.cores = parallel::detectCores())
SPCD.df <- data.frame(SPCD = nspp[1:17, ]$SPCD, 
                      spcd.id = 1:17)
remper.cor.vector <- c(0.5)
#model.number <- 6
model.list <- 1:9


# compile the model once to save time:
species.file <- file.path(getwd(), "modelcode", "mort_model_general.stan")
species.mod <- cmdstan_model(species.file)

# compile the posterior prediction stan model (prediction only)
predict.species.file <- file.path(getwd(), "modelcode", "mort_model_general_predict.stan")
predict.species.mod <- cmdstan_model(predict.species.file)

i <- 13
m <- 1
j <- 1

niter <- 1000
nwarmup <- 500
nchain <- 4
nparallel <- 4

#TODO: multi-threading with reduce_sum

#for(i in 16:1){# run for each of the 17 species
  


 run.species.models <-  function(i, 
                                 m, 
                                 nparallel,
                                 niter, 
                                 nwarmup, 
                                 nchain, 
                                 output.dir, 
                                 print.progress = 500){
  
  model.number <- model.list[m]

  common.name <- nspp[1:17, ] %>% filter(SPCD %in% SPCD.df[i,]$SPCD) %>% dplyr::select(COMMON)
  SPCD.id <- SPCD.df[i,]$SPCD
    
    
    
   # for (j in 1:length(remper.cor.vector)){ # for the growth only model explore the consequences of other assumptions about remeasurement period
      cat(paste("\n Sampling stan mortality model ", model.number, " for SPCD", SPCD.df[i,]$SPCD, common.name$COMMON, " remper correction", remper.cor.vector[j]))
     
     
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
        refresh = print.progress # print update every 500 iters
      )
     
      
      # save the fit object for later 
      model.name <- paste0("mort_model_",model.number,"_SPCD_", SPCD.id, 
                           "_remper_correction_", remper.cor.vector[j], "_niter_", niter, "_nchain_", nchain)
      
      
      
      
      # this takes the longest to save
     #fit.1$save_object(file = paste0(output.dir,"SPCD_stanoutput_cmdstan/fittedmodels/", model.name, ".rds"))
      #fit.1$save_object(file = paste0(output.dir,"SPCD_stanoutput_cmdstan/fittedmodels/", model.name, ".qs"), format = "qs2")
    cat(paste("\n Saving stan mortality model fit "))
      
    qs2::qs_save(fit.1, paste0(output.dir,"SPCD_stanoutput_cmdstan/fittedmodels/", model.name, ".qs"))
    #fit.1 <- qs2::qs_read( paste0(output.dir,"SPCD_stanoutput_cmdstan/fittedmodels/", model.name, ".qs"))
    
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
                                    "SPCD_stanoutput_cmdstan/diagnostics/sample_diagnostics_", 
                                    model.name, ".csv"), row.names = FALSE)
      rm(sampler_diag)
      # Extract posterior draws and save separately ----
     #cat(paste("\n Extracting posterior draws"))
     
      # get loo results and save the log_lik, ypredictions, and mmhat and mmrep
      # model fit statistics -----
      # loo results
      loo_results <-  fit.1$loo(cores = 4)
      qs2::qs_save(loo_results, paste0(output.dir,"SPCD_stanoutput_cmdstan/LOO/LOO_results_", model.name, ".qs"))
      rm(loo_results)
      
      log_lik_samps <-  fit.1$draws(variables = c("log_lik"), format = "draws_matrix")
      qs2::qs_save(log_lik_samps,paste0(output.dir,"SPCD_stanoutput_cmdstan/LOO/log_lik_samps_", model.name, ".qs"))
      rm(log_lik_samps)
      
      beta_alpha_samps <-  fit.1$draws(variables = c("alpha_SPP", "u_beta"), format = "draws_matrix")
      qs2::qs_save(beta_alpha_samps,paste0(output.dir,"SPCD_stanoutput_cmdstan/betas/u_beta_alpha_samps_", model.name, ".qs"))
      
      
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
      
      
     
      
     
      
      
      # use model$generate_quantities in cmdstan to generate predictions
      gen_quants <- predict.species.mod$generate_quantities(
        fitted_params = fit.1, 
        data = paste0("SPCD_standata_json/SPCD_",SPCD.id,"remper_correction_",
                      remper.cor.vector[j],"model_",model.number,".json"), # path to json data files
        seed = 123,
        parallel_chains = nparallel
        
        
      )
      
      #cat(paste("\n Extracting posterior draws"))
      y_rep_samps <- gen_quants$draws(variables = c("y_rep"), format = "draws_matrix")
      y_hat_samps <-  gen_quants$draws(variables = c("y_hat"), format = "draws_matrix")

      pSurv_rep_samps <-  gen_quants$draws(variables = c("mMrep"), format = "draws_matrix")
      pSurv_hat_samps <-  gen_quants$draws(variables = c("mMhat"), format = "draws_matrix")

      # pSannual_rep_samps <-  gen_quants$draws(variables = c("pSannualrep"), format = "draws_matrix")
      # pSannual_hat_samps <-  gen_quants$draws(variables = c("pSannualhat"), format = "draws_matrix")
      
      # save all to their own objects:
      
      qs2::qs_save(y_rep_samps,paste0(output.dir,"SPCD_stanoutput_cmdstan/predicted_mort/y_rep_samps_", model.name, ".qs"))
      qs2::qs_save(y_hat_samps,paste0(output.dir,"SPCD_stanoutput_cmdstan/predicted_mort/y_hat_samps_", model.name, ".qs"))
      
      qs2::qs_save(pSurv_rep_samps,paste0(output.dir,"SPCD_stanoutput_cmdstan/predicted_mort/pSurv_rep_samps_", model.name, ".qs"))
      qs2::qs_save(pSurv_hat_samps,paste0(output.dir,"SPCD_stanoutput_cmdstan/predicted_mort/pSurv_hat_samps_", model.name, ".qs"))
      
      # qs2::qs_save(pSannual_rep_samps,paste0(output.dir,"SPCD_stanoutput_cmdstan/predicted_mort/pSannual_rep_samps_", model.name, ".qs"))
      # qs2::qs_save(pSannual_hat_samps,paste0(output.dir,"SPCD_stanoutput_cmdstan/predicted_mort/pSannual_hat_samps_", model.name, ".qs"))
      rm(gen_quants)
      gc()
     
     
      cat(paste("\n Getting model diagnostics and summaries"))
      # convergence statistics summaries ----
          
      
      # get predicted draws for survival (0,1), annual and remper survival probabilities for in-sample ("hat") and held-out ("rep")
      
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
      
      
      #beta_alpha_samps 
      u_betas_alpha.quant <- summarise_draws(beta_alpha_samps,  summarize_posteriors)
      
      #system.time(y_hat.quant2 <-  apply(y_hat_samps, MARGIN = 2, summarize_posteriors))
      system.time(y_hat.quant <- summarise_draws(y_hat_samps,  summarize_posteriors))
      
      
      y_rep.quant <-  summarise_draws(y_rep_samps,  summarize_posteriors)
      
      # pSurv_ predicted survival probabilities over the remper for in sample and out of sample
      pSurv_hat.quant <-  summarise_draws(pSurv_hat_samps,  summarize_posteriors)
      pSurv_rep.quant <-  summarise_draws(pSurv_rep_samps,  summarize_posteriors)
      
      
      # pSannual_rep predicted annual survival probabilities for in sample and out of sample trees
      # pSannual_hat.quant <-  gen_quants$summary(variables = c("pSannualhat"), 
      #                                   median,           
      #                                   ~quantile(.x, probs = c(0.025, 0.975)), 
      #                                   mean, sd, min, max, 
      #                                   rhat, ess_bulk, ess_tail) 
      # pSannual_rep.quant <-  gen_quants$summary(variables = c( "pSannualrep"), 
      #                                      median,           
      #                                      ~quantile(.x, probs = c(0.025, 0.975)), 
      #                                      mean, sd, min, max, 
      #                                      rhat, ess_bulk, ess_tail) 
      
      
      convergence.stats <- rbind(u_betas_alpha.quant, 
                                 y_hat.quant, 
                                 y_rep.quant, 
                                 pSurv_hat.quant, 
                                 pSurv_rep.quant)%>%#, 
                                 #pSannual_hat.quant, 
                                 #pSannual_rep.quant)%>%
        mutate(model.number = model.number, 
               model.type = "Species",
               SPCD = SPCD.id, 
               remper.correction = remper.cor.vector[j])
     
      
      write.csv(convergence.stats, 
                paste0(output.dir, "SPCD_stanoutput_cmdstan/diagnostics/Rhats_ESS_quantiles_", model.name, ".csv"))
      
      
      rm(convergence.stats)
      gc()
      
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
      
    
      # read in model data for comparison
      mod.data <- fromJSON(fit.1$data_file())
      
          # AUC scores
      # for in sample data
      actuals = mod.data$y
      #preds = as.vector(pSurv_rep.quant$median)
      #auc.is <- pROC::auc( actuals, preds, quiet = TRUE) %>% as.numeric()
      
      # assuming that a probability <= 0.8 results in mortality
      # confusion.is <- data.frame(obs_outcome = actuals,
      #                            prob = preds, 
      #                            type = "in-sample") %>%
      #   mutate(pred_outcome = ifelse(prob > 0.8, 1, 0))%>%
      #     mutate(TP = ifelse(pred_outcome == 1 & obs_outcome == 1, 1, 0),
      #            FP = ifelse(pred_outcome == 1 & obs_outcome == 0, 1, 0),
      #            TN = ifelse(pred_outcome == 0 & obs_outcome == 0, 1, 0),
      #            FN = ifelse(pred_outcome == 0 & obs_outcome == 1, 1, 0))%>%
      #   group_by(type)%>%
      #     summarise(`True Survival` = sum(TP),
      #               `False Survival` = sum(FP),
      #               `True Mortality` = sum(TN),
      #               `False Mortality` = sum(FN))%>%
      #     mutate(`True survival rate` = `True Survival`/(`True Survival`+`False Mortality`),
      #            `True mortality rate` = `True Mortality`/(`True Mortality`+`False Survival`))

      # for out of sample data
      actuals.oos = mod.data$ytest
      #preds.oos = as.vector(pSurv_hat.quant$median)
      #auc.oos <- pROC::auc( actuals.oos, preds.oos, quiet = TRUE) %>% as.numeric()
      
      
      # confusion.oos <- data.frame(obs_outcome = actuals.oos,
      #                            prob = preds.oos, 
      #                            type = "out-of-sample") %>%
      #   mutate(pred_outcome = ifelse(prob > 0.8, 1, 0))%>%
      #   mutate(TP = ifelse(pred_outcome == 1 & obs_outcome == 1, 1, 0),
      #          FP = ifelse(pred_outcome == 1 & obs_outcome == 0, 1, 0),
      #          TN = ifelse(pred_outcome == 0 & obs_outcome == 0, 1, 0),
      #          FN = ifelse(pred_outcome == 0 & obs_outcome == 1, 1, 0))%>%
      #   group_by(type)%>%
      #   summarise(`True Survival` = sum(TP),
      #             `False Survival` = sum(FP),
      #             `True Mortality` = sum(TN),
      #             `False Mortality` = sum(FN))%>%
      #   mutate(`True survival rate` = `True Survival`/(`True Survival`+`False Mortality`),
      #          `True mortality rate` = `True Mortality`/(`True Mortality`+`False Survival`))
      # 
    
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
      
      # ySurv = actuals == 1
      # yMort = actuals == 0
      
     
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
               model.type = "Species", 
               SPCD = SPCD.id
               )
   
   
   # out-of-sample:
   #preds.oos <-  #%>% select(-.chain, -.iteration, -.draw) %>% as.matrix()
   preds.oos.class <- y_hat_samps == 1
   
   # ySurv = actuals.oos == 1
   # yMort = actuals.oos == 0
   
   
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
            model.type = "Species", 
            SPCD = SPCD.id
     )
      
     # hist(confusion.oos_draws$`True survival rate`)
     # hist(confusion.oos_draws$`True mortality rate`)
     
     
    # combine the AUCs and confusion matrices together with draws;
   confusion.is_draws$AUC <- AUC.is.samples.df
   confusion.oos_draws$AUC <- AUC.oos.samples.df 
   
   AUC.confusion_draws <- rbind(confusion.is_draws, confusion.oos_draws)
   
   qs2::qs_save(AUC.confusion_draws, paste0(output.dir,"SPCD_stanoutput_cmdstan/AUC/AUC_draws_", model.name, ".qs"))
     
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
   
      # PLOTS: model fit -----
     
      gc()
      
   
 #   }
 }
 
 # run.species.models(i = 17, 
 #                    m = 1, 
 #                    nparallel = nparallel,
 #                    niter = niter, 
 #                    nwarmup = nwarmup , 
 #                    nchain = nchain, 
 #                    output.dir = output.dir, 
 #                    print.progress = 0)

 
# library(progressr)
# progressr::handlers(handler_txtprogressbar(char = cli::col_red(cli::symbol$heart)))
# options(cli.progress_handlers = "progressr")
# #y <- purrr::map(1:30, slow_sqrt, .progress = TRUE)
# # Define the rainbow heart format
# rainbow_hearts <-  function(x) {
#   # Map 7 colors using cli theme syntax
#   colors <- c("#FF3B30", "#FF9500", "#FFCC00", "#4CD964", "#5AC8FA", "#5856D6", "#AF52DE")
#   
#   # Return custom formatted cli progress bar
#   cli_progress_bar(
#     name = "Processing",
#     total = length(x),
#     format = paste0(
#       "{name} {cli::pb_bar} {cli::pb_current}/{cli::pb_total} | ",
#       "{col_red('\u2665')}{col_orange('\u2665')}{col_yellow('\u2665')}{col_green('\u2665')}{col_blue('\u2665')}{col_indigo('\u2665')}{col_violet('\u2665')}"
#     )
#   )
# }
# 
# # # Example usage with lapply
# # data_list <- 1:50
# # 
# # result <- function(data_list){
# #   lapply(cli_progress_along(data_list, format = rainbow_hearts), function(i) {
# #   # Simulate a computation
# #   Sys.sleep(1) 
# #   return(i * 2)
# # })
# # }
# 
# model_list <- 1:9
# 
# run.spp.17 <- function(model_list){
#   lapply(cli_progress_along(model_list, 
#                             format = rainbow_hearts), FUN = function(x){
#    run.species.models(i = 17, 
#                       m = x, 
#                       nparallel = nparallel,
#                       niter = niter, 
#                       nwarmup = nwarmup , 
#                       nchain = nchain, 
#                       output.dir = output.dir, 
#                       print.progress = 0)
#    }
#  )
# }
# run.spp.17(model_list)



 # run.species.models(i = 13, 
 #                    m = 2, 
 #                    nparallel = nparallel,
 #                    niter = niter, 
 #                    nwarmup = nwarmup , 
 #                    nchain = nchain, 
 #                    output.dir = output.dir, 
 #                    print.progress = 0)

for (species.num in 9:1){
  if(species.num == 9){
    lapply(9, FUN = function(x){
      run.species.models(i = species.num, 
                         m = x, 
                         nparallel = nparallel,
                         niter = niter, 
                         nwarmup = nwarmup , 
                         nchain = nchain, 
                         output.dir = output.dir, 
                         print.progress = 0)
    })
  }else{

  lapply(1:9, FUN = function(x){
    run.species.models(i = species.num, 
                       m = x, 
                       nparallel = nparallel,
                       niter = niter, 
                       nwarmup = nwarmup , 
                       nchain = nchain, 
                       output.dir = output.dir, 
                       print.progress = 0)
  })  
}
}  
  
  
# do the summaries on diagnostics, time, etc
#fit.1 <- readRDS(paste0(output.dir,"SPCD_stanoutput_cmdstan/", model.name, ".rds"))



# get convergence stats---
convergence.files <- list.files(path = paste0(output.dir,"SPCD_stanoutput_cmdstan/diagnostics/"), 
                               pattern = paste0("Rhats_ESS"), full.names = TRUE)


convergence_all <- do.call(rbind, lapply(convergence.files, FUN = function(x)
  read.csv(x) %>% 
    group_by(model.number, model.type, SPCD, remper.correction)%>%
    summarise(rhat.median = median(rhat, na.rm = TRUE), 
              rhat.ci.lo = quantile(rhat, 0.025, na.rm = TRUE),
              rhat.ci.hi = quantile(rhat, 0.975, na.rm = TRUE),
              ess_bulk.median = median(ess_bulk, na.rm =TRUE), 
              ess_bulk.ci.lo = quantile(ess_bulk, 0.025, na.rm = TRUE),
              ess_bulk.ci.hi = quantile(ess_bulk, 0.975, na.rm = TRUE),
              ess_tail.median = median(ess_tail, na.rm =TRUE), 
              ess_tail.ci.lo = quantile(ess_tail, 0.025, na.rm = TRUE),
              ess_tail.ci.hi = quantile(ess_tail, 0.975, na.rm = TRUE), .groups = "drop_last")
))





diagnostic.files <- list.files(path = paste0(output.dir,"SPCD_stanoutput_cmdstan/diagnostics/"), pattern = paste0("sample_diagnostics"), full.names = TRUE)
diags_all <- do.call(rbind, lapply(diagnostic.files, read.csv))
diag.summary<- diags_all %>%  #group_by(SPCD, model.type, model.number, chain_id)%>%
  mutate(chain_core_hours = total*(1/60)*(1/60)*(ncores/nchain)) %>%
  group_by(SPCD, model.type, model.number)%>%
  summarise(total_core_hours = sum(chain_core_hours), 
            total_iter = sum(niter), 
            total_warmup = sum(nwarmup), 
            
            total_divergent = sum(num_divergent),
            total_max_treedepth = sum(num_max_treedepth))

# check that there are no divergent transitions
diag.summary %>% filter(total_divergent > 0)  
diags_all %>% filter(num_divergent > 0)  

# only one model has a few samples that exceeded max tree depth
diag.summary %>% filter( total_max_treedepth > 0)  %>% select(SPCD, model.number, total_max_treedepth)

# combine and save all the species outputs.
diag_converg.df <- left_join(diag.summary, convergence_all)

hist(diag_converg.df$rhat.median)
hist(diag_converg.df$rhat.ci.hi)
hist(diag_converg.df$rhat.ci.lo)
hist(diag_converg.df$ess_bulk.median)
hist(diag_converg.df$ess_bulk.ci.lo)


# get loo results and compare---
# need to do for each species...!!
SPCD.id <- 97

# function to read in LOO output and summarise for each species:
LOO_summarise_SPCD <- function(SPCD.id){
    spp.loo.files <- paste0(output.dir,"SPCD_stanoutput_cmdstan/LOO/LOO_results_mort_model_", 1:9, 
                              "_SPCD_", SPCD.id, "_remper_correction_0.5_niter_1000_nchain_4.qs")
      
    loo_results_all <- lapply(spp.loo.files, qs_read)
    
    
    # check pareto-k estimates:
    pareto.k.checks <- do.call(rbind, lapply(loo_results_all, function(x){data.frame(good = sum(x$diagnostics$pareto_k <= 0.7), 
                                                   bad = sum(x$diagnostics$pareto_k > 0.7), 
                                                   total = length(x$diagnostics$pareto_k))}))%>%
      mutate(percent.bad = (bad/total)*100)%>%
      mutate(model = paste0("model", 1:9), 
             SPCD = SPCD.id)
      
      
    
    loo_compare.out <- loo::loo_compare(loo_results_all) # best fit based on loo elpd differences
    loo_comparisons <- loo_compare.out %>% data.frame()%>%left_join(., pareto.k.checks) %>%
      mutate(elpd_se_ratio = abs(elpd_diff)/se_diff)%>%
      mutate(model.number = substr(model, start = 6, stop = 6))
    
    
    # Get model weights---
    # get pointwise log predictive densities
    lpd_point <- do.call(cbind,lapply(loo_results_all, function(x){x$pointwise[,"elpd_loo"]}))
    pbma_wts <- pseudobma_weights(lpd_point, BB=FALSE)
    pbma_BB_wts <- pseudobma_weights(lpd_point) # default is BB=TRUE
    stacking_wts <- stacking_weights(lpd_point)
    mod.weights <- round(cbind(pbma_wts, pbma_BB_wts, stacking_wts),3)%>% data.frame() %>% 
      mutate(model = paste0("model", 1:9), 
                                                          SPCD = SPCD.id)
    loo_comparisons <- loo_comparisons %>% left_join(., mod.weights, by = c("model", "SPCD"))
    return(loo_comparisons)
}


LOO_ELPD.df <- do.call(rbind, lapply(nspp$SPCD,LOO_summarise_SPCD))
# make sure none of the models have >5-10 percent bad pareto k values
max(LOO_ELPD.df$percent.bad)


# Order Species common names by from most to least abundant
LOO_ELPD.df$COMMON_NAME <- FIESTA::ref_species[match(LOO_ELPD.df$SPCD, FIESTA::ref_species$SPCD),]$COMMON_NAME
LOO_ELPD.df$COMMON_NAME <- factor(LOO_ELPD.df$COMMON_NAME, levels = unique(nspp$COMMON_NAME))


# plot of all species ELPD diff +/- SE, with significance
LOO_ELPD.df %>%
  mutate(elpd_diff_sig = ifelse(elpd_se_ratio >= 2, "ELPD_diff >= 2 SE", 
                                ifelse(elpd_se_ratio <2, "ELPD_diff < 2 SE", "Best-fit")))%>%
  mutate(elpd_diff_sig = ifelse(is.na(elpd_diff_sig), "Best-fit ELPD", elpd_diff_sig))|>
  ggplot()+geom_pointrange(aes(x = model, y = elpd_diff, ymin = elpd_diff+se_diff, ymax = elpd_diff-se_diff, color = elpd_diff_sig))+
  facet_wrap(~COMMON_NAME, scales = "free_y")+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust =0.5))

# Plot of ELPD differences vs model weights (BMA + bootstrapping)
LOO_ELPD.df %>%
  mutate(elpd_diff_sig = ifelse(elpd_se_ratio >= 2, "ELPD_diff >= 2 SE", 
                                ifelse(elpd_se_ratio <2, "ELPD_diff < 2 SE", "Best-fit")))%>%
  mutate(elpd_diff_sig = ifelse(is.na(elpd_diff_sig), "Best-fit ELPD", elpd_diff_sig))|>
  ggplot()+
  geom_text(aes(x = elpd_diff, y = pbma_BB_wts, color = elpd_diff_sig, label = model.number))+
  facet_wrap(~COMMON_NAME, scales = "free_x")+theme_bw()+
  ylab("BMA weights")+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust =0.5))
ggsave(filename = paste0(output.dir, "SPCD_stanoutput_cmdstan/summary/Species_model_ELPDdiff_BMA_weights.png"))


LOO_ELPD.df %>%
  mutate(elpd_diff_sig = ifelse(elpd_se_ratio >= 2, "ELPD_diff >= 2 SE", 
                                ifelse(elpd_se_ratio <2, "ELPD_diff < 2 SE", "Best-fit")))%>%
  mutate(elpd_diff_sig = ifelse(is.na(elpd_diff_sig), "Best-fit ELPD", elpd_diff_sig))|>
  ggplot()+
  geom_text(aes(x = elpd_diff, y = stacking_wts, color = elpd_diff_sig, label = model.number))+
  facet_wrap(~COMMON_NAME, scales = "free_x")+theme_bw()+
  ylab("Stacking weights")+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust =0.5))


# stacked barplot of model weights by species
LOO_ELPD.df |>  ggplot()+
  geom_bar(aes(x = COMMON_NAME, y = pbma_BB_wts, fill = model), stat = "identity")+theme_bw()+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust =0.5))

LOO_ELPD.df |>  ggplot()+
  geom_bar(aes(x = COMMON_NAME, y = stacking_wts, fill = model), stat = "identity")+theme_bw()+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust =0.5))


# unstacked barplot of model weights by species
LOO_ELPD.df |>  ggplot()+
  geom_bar(aes(x = model, y = pbma_BB_wts), stat = "identity")+theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust =0.5))+
  facet_wrap(~COMMON_NAME)

# line plot of model weights by species
LOO_ELPD.df |>  ggplot()+
  geom_point(aes(x = model, y = pbma_BB_wts, color = COMMON_NAME))+
  geom_line(aes(x = model, y = pbma_BB_wts, color = COMMON_NAME, group = COMMON_NAME))+theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust =0.5))

LOO_ELPD.df |>  ggplot()+
  geom_tile(aes(x = model, y = COMMON_NAME, fill = pbma_BB_wts))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust =0.5))+
  scale_fill_distiller(palette = "YlGnBu", direction = 1, name = "BMA Weight")+
  ylab("Species")+xlab("Model")

ggsave(filename = paste0(output.dir, "SPCD_stanoutput_cmdstan/summary/Species_model_BMA_weights.png"))

# save LOO_ELPD.df 

# get auc results ---

AUC_summarise_SPCD <- function(SPCD.id){
  spp.AUC.files <- paste0(output.dir,"SPCD_stanoutput_cmdstan/AUC/AUC_draws_mort_model_", 1:9, 
                          "_SPCD_", SPCD.id, "_remper_correction_0.5_niter_1000_nchain_4.qs")
  
  AUC_results_all <- lapply(spp.AUC.files, qs_read)
 
  AUC_confusion_summary <- do.call(rbind, lapply(AUC_results_all, function(x){
   x |> group_by(model.number, type, model.type, SPCD)%>%
     summarise(AUC_median = median(AUC),
               AUC_ci.lo = quantile(AUC, 0.025),
               AUC_ci.hi = quantile(AUC, 0.975),
               
               True_surv_rate = median(`True survival rate`),
               True_surv_rate_ci.lo = quantile(`True survival rate`, 0.025),
               True_surv_rate_ci.hi = quantile(`True survival rate`, 0.975),
               
               True_mort_rate = median(`True mortality rate`),
               True_mort_rate_ci.lo = quantile(`True mortality rate`, 0.025),
               True_mort_rate_ci.hi = quantile(`True mortality rate`, 0.975), .groups = "drop_last")
 }))
  return(AUC_confusion_summary)
}


AUC.df <- do.call(rbind, lapply(nspp$SPCD, AUC_summarise_SPCD))
AUC.df$COMMON_NAME <- FIESTA::ref_species[match(AUC.df$SPCD, FIESTA::ref_species$SPCD),]$COMMON_NAME
AUC.df$COMMON_NAME <- factor(AUC.df$COMMON_NAME, levels = unique(nspp$COMMON_NAME))


AUC.df %>% filter(type == "in-sample")|> 
  ggplot()+geom_pointrange(aes(x = as.character(model.number), y = AUC_median, ymin = AUC_ci.lo, ymax = AUC_ci.hi))+
  facet_wrap(~COMMON_NAME)

AUC.df %>% filter(type == "out-of-sample")|> 
  ggplot()+geom_pointrange(aes(x = as.character(model.number), y = AUC_median, ymin = AUC_ci.lo, ymax = AUC_ci.hi))+
  facet_wrap(~COMMON_NAME, scales = "free_y")


AUC.df %>% filter(type == "out-of-sample") |> 
  ggplot()+geom_pointrange(aes(x = as.character(model.number), y = True_surv_rate, ymin = True_surv_rate_ci.lo, ymax = True_surv_rate_ci.hi))+
  facet_wrap(~COMMON_NAME, scales = "free_y")

AUC.df %>% filter(type == "out-of-sample") |> 
  ggplot()+geom_pointrange(aes(x = as.character(model.number), y = True_mort_rate, ymin = True_mort_rate_ci.lo, ymax = True_mort_rate_ci.hi))+
  facet_wrap(~COMMON_NAME, scales = "free_y")


AUC.df %>% filter(type == "in-sample") |> 
  ggplot()+geom_pointrange(aes(x = as.character(model.number), y = True_surv_rate, ymin = True_surv_rate_ci.lo, ymax = True_surv_rate_ci.hi))+
  facet_wrap(~COMMON_NAME, scales = "free_y")

AUC.df %>% filter(type == "in-sample") |> 
  ggplot()+geom_pointrange(aes(x = as.character(model.number), y = True_mort_rate, ymin = True_mort_rate_ci.lo, ymax = True_mort_rate_ci.hi))+
  facet_wrap(~COMMON_NAME, scales = "free_y")

# link up AUC scores with the ELPD and weights for comparisons

AUC.oos.df <- AUC.df %>% filter(type == "out-of-sample")
AUC.is.df <- AUC.df %>% filter(type == "out-of-sample")

OOS.AUC.ELPD.df <- left_join(LOO_ELPD.df%>% mutate(model.number = as.numeric(model.number)),AUC.oos.df)
IS.AUC.ELPD.df <- left_join(LOO_ELPD.df%>% mutate(model.number = as.numeric(model.number)),AUC.is.df)

OOS.AUC.ELPD.df <- OOS.AUC.ELPD.df %>% 
  mutate(elpd_diff_sig = ifelse(elpd_se_ratio >= 2, "ELPD_diff >= 2 SE", 
                                ifelse(elpd_se_ratio <2, "ELPD_diff < 2 SE", "Best-fit")))%>%
  mutate(elpd_diff_sig = ifelse(is.na(elpd_diff_sig), "Best-fit ELPD", elpd_diff_sig))

IS.AUC.ELPD.df <- IS.AUC.ELPD.df %>% 
  mutate(elpd_diff_sig = ifelse(elpd_se_ratio >= 2, "ELPD_diff >= 2 SE", 
                                ifelse(elpd_se_ratio <2, "ELPD_diff < 2 SE", "Best-fit")))%>%
  mutate(elpd_diff_sig = ifelse(is.na(elpd_diff_sig), "Best-fit ELPD", elpd_diff_sig))



ggplot(data = OOS.AUC.ELPD.df)+
  geom_pointrange(aes(x = model, y = AUC_median, ymin = AUC_ci.lo, ymax = AUC_ci.hi))+
  facet_wrap(~COMMON_NAME, scales = "free")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  ylab("AUC")+xlab("")


ggplot(data = OOS.AUC.ELPD.df)+
  geom_pointrange(aes(x = elpd_loo, y = AUC_median, ymin = AUC_ci.lo, ymax = AUC_ci.hi, color = model))+
  facet_wrap(~COMMON_NAME, scales = "free")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  ylab("AUC")+xlab("Expeced Log Predictive Density-LOO")

ggplot(data = OOS.AUC.ELPD.df)+
  geom_pointrange(aes(x = elpd_diff, y = AUC_median, ymin = AUC_ci.lo, ymax = AUC_ci.hi, color = model))+
  facet_wrap(~COMMON_NAME, scales = "free")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  ylab("AUC")+xlab("Expeced Log Predictive Density-LOO")


OOS.AUC.ELPD.df <- OOS.AUC.ELPD.df%>%left_join(.,
  OOS.AUC.ELPD.df %>% group_by(SPCD, COMMON_NAME)%>%
    mutate(max.AUC = max(AUC_median,na.rm = TRUE))%>%
    filter(AUC_median == max.AUC)%>% select(SPCD, COMMON_NAME, type, AUC_median, AUC_ci.lo, AUC_ci.hi)%>%
    rename("bf_AUC"= "AUC_median", 
           "bf_AUC.ci.lo"= "AUC_ci.lo", 
           "bf_AUC_ci.hi"= "AUC_ci.hi")
)%>%
  mutate(AUC_sig = ifelse(AUC_median == bf_AUC, "Best-fit AUC",
                          ifelse(AUC_ci.hi > bf_AUC.ci.lo, "overlapping CI", "non-overlapping CI")))


IS.AUC.ELPD.df <- IS.AUC.ELPD.df%>%left_join(.,
                                               IS.AUC.ELPD.df %>% group_by(SPCD, COMMON_NAME)%>%
                                                 mutate(max.AUC = max(AUC_median,na.rm = TRUE))%>%
                                                 filter(AUC_median == max.AUC)%>% select(SPCD, COMMON_NAME, type, AUC_median, AUC_ci.lo, AUC_ci.hi)%>%
                                                 rename("bf_AUC"= "AUC_median", 
                                                        "bf_AUC.ci.lo"= "AUC_ci.lo", 
                                                        "bf_AUC_ci.hi"= "AUC_ci.hi")
)%>%
  mutate(AUC_sig = ifelse(AUC_median == bf_AUC, "Best-fit AUC",
                          ifelse(AUC_ci.hi > bf_AUC.ci.lo, "overlapping CI", "non-overlapping CI")))

# How consistent is AUC with elpd and model weights?---
# out-of-sample AUC, colored by best fit elpd, shapes are AUC sig.
OOS.AUC.ELPD.df|>
  ggplot()+
  geom_pointrange(aes(x = model, y = AUC_median, ymin = AUC_ci.lo, ymax = AUC_ci.hi, color = elpd_diff_sig, shape = AUC_sig))+
  #facet_wrap(~COMMON_NAME, scales = "free")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  ylab("Out-of-Sample AUC")+
  facet_wrap(~COMMON_NAME, scales = "free_y")
ggsave(filename = paste0(output.dir, "SPCD_stanoutput_cmdstan/summary/Species_model_AUC_OOS.png"))


# in-sample AUC, colored by best fit elpd
IS.AUC.ELPD.df|>
  ggplot()+
  geom_pointrange(aes(x = model, y = AUC_median, ymin = AUC_ci.lo, ymax = AUC_ci.hi, color = elpd_diff_sig, shape = AUC_sig))+
  #facet_wrap(~COMMON_NAME, scales = "free")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  ylab("Out-of-Sample AUC")+
  facet_wrap(~COMMON_NAME, scales = "free_y")
ggsave(filename = paste0(output.dir, "SPCD_stanoutput_cmdstan/summary/Species_model_AUC_IS.png"))


OOS.AUC.ELPD.df|>
  ggplot()+
  geom_bar(aes(x = AUC_sig, fill = elpd_diff_sig ))+
 
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  ylab("Out-of-Sample AUC")+
  facet_wrap(~COMMON_NAME, scales = "free_y")
IS.AUC.ELPD.df %>% select(model, COMMON_NAME, AUC_sig)%>%
  rename("IS_AUC_sig"= "AUC_sig")%>%left_join(.,OOS.AUC.ELPD.df) %>% select(model, COMMON_NAME, IS_AUC_sig, AUC_sig, elpd_diff_sig)%>%
  pivot_longer( cols = c("AUC_sig","IS_AUC_sig", "elpd_diff_sig"))%>%
  mutate(Comparison.stat = ifelse(name %in% "AUC_sig", "out-of-sample AUC", 
                                  ifelse(name %in% "IS_AUC_sig","in-sample AUC","ELPD-diff")), 
         Cat.signficance = ifelse(value %in% c("Best-fit ELPD", "Best-fit AUC"), "Best-fit", 
                                  ifelse(value %in% c("overlapping CI", "ELPD_diff < 2 SE"), "Overlapping with Best-fit", "Non-overlapping")))|>
  ggplot()+
  geom_tile(aes(y = Comparison.stat, x = model, fill = Cat.signficance ))+
  
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  scale_fill_manual(values = c("Best-fit" = "darkred", 
                               "Overlapping with Best-fit" = "red", 
                               "Non-overlapping"= "lightgrey"
                               ), name = "Significance")+
  facet_wrap(~COMMON_NAME)+
  ylab("Model Comparison Statistic")

ggsave(filename = paste0(output.dir, "SPCD_stanoutput_cmdstan/summary/Species_model_AUC_ELPD_comparison_tile.png"))

