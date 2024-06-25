SPCD_run_stan <- function(SPCD.id, model.no = 1, niter = 1000, nchains = 2, remper.correction = 0.5,  model.file = 'modelcode/mort_model3_SPCD.stan'){

  model.name <- paste0("mort_model", model.no, "_single_SPCD_", SPCD.id, "remper_", remper.correction)
  load(paste0("SPCD_standata_general_full/SPCD_",SPCD.id, "remper_correction_", remper.correction,"model_",model.no, ".Rdata")) # load the species code data
  mod.data$K <- ncol(mod.data$xM)
  # y == 0 is mortality and y == 1 is survival
  num_cores <-  parallel::detectCores()
  # null model:
  start.time <- Sys.time()
  fit.1 <- stan(file = model.file , 
                data = mod.data,
                iter = niter, 
                chains = nchains, 
                verbose=FALSE, 
                ##control =  list(max_treedepth = 15),#list(adapt_delta = 0.99, stepsize = 0.5, max_treedepth = 15),#, stepsize = 0.01, max_treedepth = 15),
                #sample_file = model.name, 
                #adapt_delta = 0.99, 
                pars =c("alpha_SPP", "u_beta", 
                        "y_rep", "mMrep",## in sample predictions
                        "y_hat", "mMhat", ## out of sample predictions
                        "log_lik")) #, "y_hat", 
  
 
 

end.time <- Sys.time()


# Calculate elapsed time in seconds
elapsed_time <- as.numeric(difftime(end.time, start.time, units = "secs"))

# Calculate core hours
core_hours <- elapsed_time * num_cores / 3600  # convert seconds to hours

time.diag <- data.frame(model = model.no, 
                        SPCD = SPCD.id, 
                        remper = remper.correction,
                        core.hours = core_hours, 
                        elapsed.time = elapsed_time, 
                        cores = num_cores)

write.csv(time.diag, paste0("SPCD_stanoutput_full/computational_resources/time_diag_SPCD_",SPCD.id, "_model_", model.no, "_remper_", remper.correction,".csv"))
# get the sampler diagnostics and save:
saveRDS(fit.1, paste0("SPCD_stanoutput_full/samples/model_",model.no,"_SPCD_",SPCD.id, "_remper_correction_", remper.correction, ".RDS"))

}

save_diagnostics <- function(stanfitobj = fit.1, nchains = 2, model.no = 1, remper.correction = 0.5){
  model.name <-paste0("model_", model.no)
    # get model diagnostics and save these to look at
    sampler_params <- get_sampler_params(stanfitobj, inc_warmup = FALSE)
    
    mean_accept_stat_by_chain <- sapply(sampler_params, function(x) mean(x[, "accept_stat__"]))
    sum_divergent_transitions_by_chain <- sapply(sampler_params, function(x) sum(x[, "divergent__"]))
    sampler_diag <- data.frame(model.name = rep(model.name, nchains),
                               chain = 1:nchains, 
                               accept = mean_accept_stat_by_chain, 
                               sum_divergent_transitions_by_chain = sum_divergent_transitions_by_chain, 
                               n.samples = nrow(sampler_params[[1]]))
    # we want few (no) divergent transitions so this is good
    sampler_diag
    # there are no divergent transistions in either chain, and acceptance rate > 0.89
    write.csv(sampler_diag, paste0("SPCD_stanoutput/sample_diagnostics_", model.name,"_remper_",remper.correction, "_species_", SPCD.id , ".csv"), row.names = FALSE)
    
    
    # get the convergence statistics of the model:
    fit_ssm_df <- as.data.frame(stanfitobj) # takes awhile to convert to df
    Rhats <- apply(fit_ssm_df, 2, Rhat)
    hist(Rhats)
    # most of the R hat values are below 1.01
    ESS_bulks <- apply(fit_ssm_df, 2, ess_bulk)
    hist(ESS_bulks)
    ESS_tails <- apply(fit_ssm_df, 2, ess_tail)
    hist(ESS_tails)
    
    convergence.stats <- as.data.frame(rbind(Rhats, ESS_bulks, ESS_tails))
    convergence.stats$Statistic <- c("Rhat", "ESS_bulk", "ESS_tail")
    
    write.csv(convergence.stats,paste0("SPCD_stanoutput/Rhats_diagnostics_", model.name,"_remper_",remper.correction,  "_species_", SPCD.id , ".csv"))
}
  #end.time <- Sys.time()
  #mod1REtime <-  end.time - start.time 
plot.stan.mort <- function(fit = fit.1, SPCD.id){
  
  names(fit.1) <- c("alpha_SPP", colnames(mod.data$xM), "lp__")
  par.names = c("alpha_SPP", colnames(mod.data$xM)) #,
  
  
  
  png(height = 12, width = 12, units = "in", res = 100, paste0("SPGRP_stanoutput/images/traceplots_mortality_", model.name,"_species_", SPCD.id , ".png"))
  #par(mfrow = c(5, 3))
  traceplot (fit.1, pars = par.names, nrow = 7, ncol = 6, inc_warmup = FALSE) 
  dev.off()
  
  png(height = 12, width = 12, units = "in", res = 100, paste0("SPCD_stanoutput/pairs_plot_mortality_", model.name, "_species_", SPCD.id , ".png"))
  pairs(fit.1, pars = "alpha", "u_beta")
  dev.off()
  
  
  species.table <- unique(train.data[,c("SPCD","SPP")])
  species.table$COMMON <- FIESTA::ref_species[match(species.table$SPCD, FIESTA::ref_species$SPCD),]$COMMON_NAME
  species.table$SPP <- as.character(species.table$SPP)
  
  # get the main effects by species:
  fit_ssm_df <- as.data.frame(fit.1) # takes awhile to convert to df
  #covariates = c(
  # "alpha", "beta_DIA","beta_si", "beta_growth")
  
  
  # note for model 0 there are no covariates
  cov.estimates <- fit_ssm_df %>% dplyr::select(-lp__, -alpha_SPP) 
  cov.m <- reshape2::melt(cov.estimates)
  
  betas.quant <- cov.m %>% group_by(variable) %>% summarise(median = quantile(value, 0.5, na.rm =TRUE),
                                                            ci.lo = quantile(value, 0.025, na.rm =TRUE),
                                                            ci.hi = quantile(value, 0.975, na.rm =TRUE))
  
  
  # clean up the naming structure of this
 
  betas.quant$SPCD <- SPCD.id
  betas.quant$Covariate <-  colnames(mod.data$xM)
  
  
  betas.quant <- left_join(betas.quant, species.table)
  
  
  ggplot(data = na.omit(betas.quant), aes(x = Covariate, y = median))+geom_point()+
    geom_errorbar(data = na.omit(betas.quant), aes(x = Covariate , ymin = ci.lo, ymax = ci.hi), width = 0.1)+
    geom_abline(aes(slope = 0, intercept = 0), color = "grey", linetype = "dashed")+facet_wrap(~COMMON)+theme_bw(base_size = 12)+
    theme( axis.text.x = element_text(angle = 45, hjust = 1))+ylab("Effect on mortality")+xlab("Parameter")
  ggsave(height = 5, width = 7, units = "in",paste0("SPCD_stanoutput/images/Estimated_effects_on_mortality_",model.name, "_species_", SPCD.id , ".png"))
}
