library(tidyverse)
library(cmdstanr)
library(posterior)
library(bayesplot)
color_scheme_set("brightblue")
#check_cmdstan_toolchain()
#cmdstan_path()


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
nchain <- 4
for(m in 1:9){ 
#for(i in 9:1){# run for each of the 17 species
#for(m in 1:9){ 
  
  model.number <- model.list[m]

  common.name <- nspp[1:17, ] %>% filter(SPCD %in% SPCD.df[i,]$SPCD) %>% dplyr::select(COMMON)
  SPCD.id <- SPCD.df[i,]$SPCD
    #for(m in 1:length(model.list)){  # run each of the 9 models
    
    
   # for (j in 1:length(remper.cor.vector)){ # for the growth only model explore the consequences of other assumptions about remeasurement period
      cat(paste("running stan mortality model ", model.number, " for SPCD", SPCD.df[i,]$SPCD, common.name$COMMON, " remper correction", remper.cor.vector[j]))
     
      
       fit.1 <- species.mod$sample(
        data = paste0("SPCD_standata_json/SPCD_",SPCD.id,"remper_correction_",
                      remper.cor.vector[j],"model_",model.number,".json"), # path to json data files
        seed = 123,
        chains = nchain,
        iter_warmup = 500,
        iter_sampling = niter,
        parallel_chains = nchain,
        adapt_delta = 0.95,
        init = 0.5,
        refresh = 50 # print update every 500 iters
        
      )
      
      # # get and save the samples:
      # fit_ssm_df <-  fit.1$draws(format = "df")
      # 
      # # get and save the model diagnositics
      # summary.params <-  fit.1$summary(variables = c("alpha_SPP", "u_beta"))
      # 
      # # get loo results and save the log_lik, ypredictions, and mmhat and mmrep
      # loo_results <-  fit.1$loo()
      # log_lik_samps <-  fit.1$draws(variables = c("log_lik"))
      # y_rep_samps <-  fit.1$draws(variables = c("y_rep"))
      # y_hat_samps <-  fit.1$draws(variables = c("y_hat"))
      # 
      # pSurv_rep_samps <-  fit.1$draws(variables = c("mMrep"))
      # pSurv_hat_samps <-  fit.1$draws(variables = c("mMhat"))
      # 
      # #summary.all <-  fit.1$summary()
      # 
      # set up to make plots of the stan outputs 
      model.name <- paste0("mort_model_",model.number,"_SPCD_", SPCD.id, 
                           "_remper_correction_", remper.cor.vector[j], "_niter_", niter, "_nchain_", nchain)
      remp.cor <- remper.cor.vector[j]
      remper.correction <- remper.cor.vector[j]
      
      output.dir = "C:/Users/KellyHeilman/Box/01. kelly.heilman Workspace/mortality/Eastern-Mortality/mortality_models/"
      fit.1$save_object(file = paste0(output.dir,"SPCD_stanoutput_cmdstan/", model.name, ".rds"))
      # source("R/speciesModels/SPCD_plot_stan.R")
      rm(fit.1)
#     }
#   }
}

# do the summaries on diagnostics, time, etc
fit.1 <- readRDS(paste0(output.dir,"SPCD_stanoutput_cmdstan/", model.name, ".rds"))
fit.1$sampler_diagnostics()
mod.data <- fromJSON(fit.1$data_file())
summary(mod.data$xM)
summary(mod.data$Remper)
mod.data$N
sum(mod.data$y == 0) # only 139 mortalities and estimating 78 parameters?


loo.results <- fit.1$loo()

# do the summaries on model parameter estimates
fit.1$summary(variables = c("lp__", "alpha_SPP","u_beta"))
# 
