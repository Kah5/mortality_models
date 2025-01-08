load(paste0("SPCD_standata_general_full_standardized/SPCD_",SPCD.id, "remper_correction_", remper.correction,"model_",model.number, ".Rdata")) # load the species code data
# read in the model fit
#if(model.number < 4){

#fit.1 <- readRDS(url( paste0("https://data.cyverse.org/dav-anon/iplant/home/kellyheilman/analyses/Species_level_mortality-2024-10-09-21-52-40.1/SPCD_stanoutput_full_standardized/samples/model_",model.number,"_SPCD_",SPCD.id,"_remper_correction_",remper.cor.vector[j],".RDS")))

#}else{
# if(model.number %in% c(1,2,3,4,5)){
#   fit.1 <- readRDS( paste0(output.folder, "SPCD_stanoutput_full_cyverse_4_5/SPCD_stanoutput_full_standardized/samples/model_",model.number,"_SPCD_",SPCD.id, "_remper_correction_", remper.cor.vector[j], ".RDS"))
#   
# }else{
fit.1 <- readRDS( paste0(output.folder, "SPCD_stanoutput_full_standardized/samples/model_",model.number,"_SPCD_",SPCD.id, "_remper_correction_", remper.cor.vector[j], ".RDS"))
#}
# }
species.table <- unique(train.data[,c("SPCD","SPP")])
species.table$COMMON <- FIESTA::ref_species[match(species.table$SPCD, FIESTA::ref_species$SPCD),]$COMMON_NAME
species.table$SPP <- as.character(species.table$SPP)

# if (model.number >= 7){names(fit.1) <- c("alpha_SPP", colnames(mod.data$xM),
#                                         ## in sample predicted status
#                                         # paste0("yrep[",1:mod.data$N, "]"),
#                                         #paste0("psurv[",1:mod.data$N, "]"),
#                                         #paste0("psurv.annual[",1:mod.data$N, "]"),
#                                         # 
#                                         # ## out of  sample predicted status
#                                         # paste0("yhat[",1:mod.data$Nrep, "]"),
#                                         # ## out of sample predicted prob mor
#                                         # paste0("psurv.hat[",1:mod.data$Nrep, "]"),
#                                         # paste0("psurv.hat.annual[",1:mod.data$Nrep, "]"),
#                                         ## in sample predicted status
#                                         #paste0("log_lik[",1:mod.data$N, "]"),
#                                         "lp__")
# }else{

names(fit.1) <- c("alpha_SPP", colnames(mod.data$xM),
                  ## in sample predicted status
                  paste0("yrep[",1:mod.data$N, "]"),
                  paste0("psurv[",1:mod.data$N, "]"),
                  #paste0("psurv.annual[",1:mod.data$N, "]"),
                  # 
                  # ## out of  sample predicted status
                  paste0("yhat[",1:mod.data$Nrep, "]"),
                  # ## out of sample predicted prob mor
                  paste0("psurv.hat[",1:mod.data$Nrep, "]"),
                  # paste0("psurv.hat.annual[",1:mod.data$Nrep, "]"),
                  ## in sample predicted status
                  # paste0("log_lik[",1:mod.data$N, "]"),
                  "lp__")

#}
par.names = c("alpha_SPP", colnames(mod.data$xM)) #,

nvariables <- length(par.names)

fit_ssm_df <- as_draws_df(fit.1) # takes awhile to convert to df

# ##################################
# # in-sample pred vs obs plots
# ##################################
# # estimate the annual probability of survival and convert to the remper probability of survival for each tree
# # get alpha estimates
alpha.estimates <- fit_ssm_df %>% dplyr::select("alpha_SPP")

# note for model 0 there are no covariates
cov.estimates <- fit_ssm_df %>% select(colnames(mod.data$xM))
rm(fit.1)

# # extract the posterior samples of the stan model coefficients:
psurv.quant <- summarise_draws(fit_ssm_df %>% select(paste0("psurv[",1:mod.data$N, "]")), 
                               #variables = c(colnames(mod.data$xM)), 
                               median, ~quantile(.x, probs = c(0.025, 0.975)))%>% rename(`ci.lo` = `2.5%`, `ci.hi` = `97.5%`)
psurv.quant$S <- mod.data$y

p.surv.all <- fit_ssm_df %>% select(paste0("psurv[",1:mod.data$N, "]"))

# param.sample <- cbind(alpha.estimates, cov.estimates)
# 
# posterior.predict.annual <- function(x){
#   
#   
# # system.time( result.sample <- t(apply(param.sample, 1, function(x) x["alpha_SPP"] + x[paste0(colnames(mod.data$xM))]*mod.data$xM[1:10] )))
#   
#   # linear predictions
#   logit_p_annual = as.data.frame(alpha.estimates) %>% dplyr::select(alpha_SPP) +
#     rowSums(as.data.frame(cov.estimates) %>% dplyr::select(paste0(colnames(mod.data$xM)))*mod.data$xM[x,])
#   p.surv.annual <- exp(logit_p_annual)/ (1 + exp(logit_p_annual)) # inverse.logit to calculate annual survival probabilities
#   #p.surv.annual1 <- inv.logit(logit_p_annual$alpha_SPP)# inverse.logit to calculate annual survival probabilities
#   
#   # calcualte total surivvial probability
#   p.surv.remper <- (p.surv.annual$alpha_SPP)^(mod.data$Remper[x])
#   
#   # calculate mortality probabilities
#   p.mort.annual <- 1-p.surv.annual
#   p.mort.remper <- 1-p.surv.remper
#   
#   
#   p.mort.df <- data.frame(
#     annual.survival = p.surv.annual$alpha_SPP,
#     remper.survival = p.surv.remper,
#     annual.mortality = p.mort.annual,
#     remper.mortality = p.mort.remper,
#     S = mod.data$y[x],
#     remper = mod.data$Remper[x],
#     tree = x)
#  p.mort.summary.df <-  p.mort.df %>% group_by(tree) %>% summarise(Mob = mean(S),
#                                                 median = quantile(remper.survival, 0.5, na.rm =TRUE),
#                                                 ci.lo = quantile(remper.survival, 0.025, na.rm =TRUE),
#                                                 ci.hi = quantile(remper.survival, 0.975, na.rm =TRUE))
#  p.mort.summary.df
# }
# pred.list <- list()
# 
# for(s in 1:mod.data$N){
#   cat(paste(s, "/", mod.data$N, "\n"))
#   pred.list[[s]] <-  posterior.predict.annual(x = s)
# }
# # system.time(pred.list <- lapply(1:10, FUN = posterior.predict.annual))
# # 
# # system.time(pred.list <- lapply(1:mod.data$N, FUN = posterior.predict.annual))
# psurv.quant <- do.call(rbind, pred.list)

# psurv.quant <- pmort.rep.is %>% group_by(tree) %>% summarise(Mob = mean(S),
#                                                              median = quantile(remper.survival, 0.5, na.rm =TRUE),
#                                                              ci.lo = quantile(remper.survival, 0.025, na.rm =TRUE),
#                                                              ci.hi = quantile(remper.survival, 0.975, na.rm =TRUE))



# pmort.hat.oos.m <- reshape2::melt(pmort.hat.oos)
# colnames(pmort.hat.oos.m) <- c("sample", "treeid", "pmort")
# pmort.hat.oos.summary <- pmort.hat.oos.m %>% group_by(treeid) %>% summarise(median = quantile(pmort, 0.5, na.rm =TRUE),
#                                                                             ci.lo = quantile(pmort, 0.025, na.rm =TRUE),
#                                                                             ci.hi = quantile(pmort, 0.975, na.rm =TRUE))
#
# psurv.estimates <- fit_ssm_df %>% dplyr::select( paste0("psurv[",1:mod.data$N, "]"))
# psurv.m <- reshape2::melt(psurv.estimates)

# psurv.quant <- psurv.m %>% group_by(variable) %>% summarise(median = quantile(value, 0.5, na.rm =TRUE),
#                                                           ci.lo = quantile(value, 0.005, na.rm =TRUE),
#                                                           ci.hi = quantile(value, 0.975, na.rm =TRUE))

#hist(psurv.quant$median)
psurv.quant$Mobs <- as.character(mod.data$y)
psurv.quant$COMMON <- unique(species.table$COMMON)


ggplot(psurv.quant, aes(as.character(Mobs), y = median, fill =Mobs ))+geom_violin()+ylab("Median annual predicted probability of survival")+
  xlab("In-sample Observed Survival Status")+scale_fill_manual(values = c("0" = "#a6611a", "1"= "forestgreen"))+theme_bw()+theme(legend.position = "none")+facet_wrap(~COMMON)
ggsave(height = 4, width = 4, units = "in",paste0(output.folder, "SPCD_stanoutput_full_standardized/images/psurv_vs_obs_in_sample_violin_model_",model.number, "_species_", SPCD.id ,"_remper_corr_", remper.cor.vector[j], ".png"))

# # compare the observed status to the predicted status
yrep.estimates <- fit_ssm_df %>% dplyr::select( paste0("yrep[",1:mod.data$N, "]"))
yrep.m <- reshape2::melt(yrep.estimates)

yrep.quant <- yrep.m %>% group_by(variable) %>% summarise(
  median = quantile(value, 0.5, na.rm =TRUE),
  ci.lo = quantile(value, 0.005, na.rm =TRUE),
  ci.hi = quantile(value, 0.975, na.rm =TRUE))


yrep.quant$Mobs <- as.character(mod.data$y)
yrep.quant$COMMON <- unique(species.table$COMMON)

ggplot(yrep.quant, aes(as.character(Mobs), y = median, fill =Mobs ))+geom_violin()+ylab("Median predicted survival status")+
  xlab("Observed in-sample tree survival status")+scale_fill_manual(values =  c("0" = "#a6611a", "1"= "forestgreen"))+theme_bw()+theme(legend.position = "none")+facet_wrap(~COMMON)
ggsave(height = 4, width = 4, units = "in",paste0(output.folder, "SPCD_stanoutput_full_standardized/images/Yrepsurv_insample_median_vs_obs_violin_model_",model.number, "_species_", SPCD.id ,"_remper_corr_", remper.cor.vector[j], ".png"))

ggplot(yrep.quant, aes(as.character(Mobs), y = ci.hi, fill =Mobs ))+geom_violin()+ylab("97.5% quantile of predicted survival status")+
  xlab("Observed in-sample tree survival status")+scale_fill_manual(values =  c("0" = "#a6611a", "1"= "forestgreen"))+theme_bw()+theme(legend.position = "none")+facet_wrap(~COMMON)
ggsave(height = 4, width = 4, units = "in",paste0(output.folder, "SPCD_stanoutput_full_standardized/images/Yrepsurv_insample_95pct_vs_obs_violin_model_",model.number, "_species_", SPCD.id ,"_remper_corr_", remper.cor.vector[j], ".png"))

ggplot(yrep.quant, aes(as.character(Mobs), y = ci.lo, fill =Mobs ))+geom_violin()+ylab("2.5% quantile of predicted survival status")+
  xlab("Observed in-sample tree survival status")+scale_fill_manual(values =  c("0" = "#a6611a", "1"= "forestgreen"))+theme_bw()+theme(legend.position = "none")+facet_wrap(~COMMON)
ggsave(height = 4, width = 4, units = "in",paste0(output.folder, "SPCD_stanoutput_full_standardized/images/Yrepsurv_insample_2.5pct_vs_obs_violin_model_",model.number, "_species_", SPCD.id ,"_remper_corr_", remper.cor.vector[j], ".png"))

# ##################################
# # out of sample plots
# ##################################
#alpha.estimates <- fit_ssm_df %>% dplyr::select("alpha_SPP")

# posterior.predict.annual.oos <- function(x){
#   
#   # linear predictions
#   logit_p_annual = as.data.frame(alpha.estimates) %>% dplyr::select(alpha_SPP) +
#     rowSums(as.data.frame(cov.estimates) %>% dplyr::select(paste0(colnames(mod.data$xMrep)))*mod.data$xMrep[x,])
#   p.surv.annual <- exp(logit_p_annual)/ (1 + exp(logit_p_annual))  # inverse.logit to calculate annual survival probabilities
#   
#   # calcualte total surivvial probability
#   p.surv.remper <- (p.surv.annual$alpha_SPP)^(mod.data$Remper[x])
#   
#   # calculate mortality probabilities
#   p.mort.annual <- 1-p.surv.annual
#   p.mort.remper <- 1-p.surv.remper
#   
#   
#   p.mort.df <- data.frame(
#     annual.survival = p.surv.annual$alpha_SPP,
#     remper.survival = p.surv.remper,
#     annual.mortality = p.mort.annual,
#     remper.mortality = p.mort.remper,
#     S = mod.data$ytest[x],
#     remper = mod.data$Remper[x],
#     tree = x)
#   p.mort.summary.df <-  p.mort.df %>% group_by(tree) %>% summarise(Mob = mean(S),
#                                                                    median = quantile(remper.survival, 0.5, na.rm =TRUE),
#                                                                    ci.lo = quantile(remper.survival, 0.025, na.rm =TRUE),
#                                                                    ci.hi = quantile(remper.survival, 0.975, na.rm =TRUE))
#   p.mort.summary.df
# }
# 
# pred.list <- list()
# 
# for(s in 1:mod.data$Nrep){
#   cat(paste(s, "/", mod.data$Nrep, "\n"))
#   pred.list[[s]] <-  posterior.predict.annual.oos(x = s)
# }
# 
# #pred.list <- lapply(1:mod.data$Nrep, FUN = posterior.predict.annual.oos)
# psurv.hat.quant <- do.call(rbind, pred.list)
psurv.hat.all <- fit_ssm_df %>% select(paste0("psurv.hat[",1:mod.data$Nrep, "]"))
psurv.hat.quant <- summarise_draws(fit_ssm_df %>% select(paste0("psurv.hat[",1:mod.data$Nrep, "]")), 
                                   #variables = c(colnames(mod.data$xM)), 
                                   median, ~quantile(.x, probs = c(0.025, 0.975)))%>% rename(`ci.lo` = `2.5%`, `ci.hi` = `97.5%`)
psurv.hat.quant$Mobs <- as.character(mod.data$ytest)
psurv.hat.quant$COMMON <- unique(species.table$COMMON)


ggplot(psurv.hat.quant, aes(as.character(Mobs), y = median, fill = Mobs ))+geom_violin()+ylab("Median predicted probability of survival")+
  xlab("Out-of-sample Observed Tree Status")+scale_fill_manual(values = c("0" = "#a6611a", "1"= "forestgreen"))+theme_bw()+theme(legend.position = "none")+facet_wrap(~COMMON)
ggsave(height = 4, width = 4, units = "in",paste0(output.folder, "SPCD_stanoutput_full_standardized/images/psurv.hat_vs_obs_out_of_sample_violin_model_",model.number, "_species_", SPCD.id ,"_remper_corr_", remper.cor.vector[j], ".png"))
#
yhat.estimates <- fit_ssm_df %>% dplyr::select( paste0("yhat[",1:mod.data$Nrep, "]"))
yhat.m <- reshape2::melt(yhat.estimates)

yhat.quant <- yhat.m %>% group_by(variable) %>% summarise(median = quantile(value, 0.5, na.rm =TRUE),
                                                          ci.lo = quantile(value, 0.005, na.rm =TRUE),
                                                          ci.hi = quantile(value, 0.975, na.rm =TRUE))


yhat.quant$Mobs <- as.character(mod.data$ytest)
yhat.quant$COMMON <- unique(species.table$COMMON)

ggplot(yhat.quant, aes(as.character(Mobs), y = median, fill =Mobs ))+geom_violin()+ylab("Median predicted survival status")+
  xlab("Observed in-sample tree status")+scale_fill_manual(values =  c("0" = "#a6611a", "1"= "forestgreen"))+theme_bw()+theme(legend.position = "none")+facet_wrap(~COMMON)
ggsave(height = 4, width = 4, units = "in",paste0(output.folder, "SPCD_stanoutput_full_standardized/images/YhatMort_out_of_sample_median_vs_obs_violin_model_",model.number, "_species_", SPCD.id ,"_remper_corr_", remper.cor.vector[j], ".png"))

ggplot(yhat.quant, aes(as.character(Mobs), y = ci.hi, fill =Mobs ))+geom_violin()+ylab("97.5% quantile of predicted survival status")+
  xlab("Observed in-sample tree status")+scale_fill_manual(values =  c("0" = "#a6611a", "1"= "forestgreen"))+theme_bw()+theme(legend.position = "none")+facet_wrap(~COMMON)
ggsave(height = 4, width = 4, units = "in",paste0(output.folder, "SPCD_stanoutput_full_standardized/images/YhatMort_out_of_sample_95pct_vs_obs_violin_model_",model.number, "_species_", SPCD.id ,"_remper_corr_", remper.cor.vector[j], ".png"))


ggplot(yhat.quant, aes(as.character(Mobs), y = ci.hi, fill =Mobs ))+geom_violin()+ylab("2.5% quantile of predicted survival status")+
  xlab("Observed in-sample tree status")+scale_fill_manual(values =  c("0" = "#a6611a", "1"= "forestgreen"))+theme_bw()+theme(legend.position = "none")+facet_wrap(~COMMON)
ggsave(height = 4, width = 4, units = "in",paste0(output.folder, "SPCD_stanoutput_full_standardized/images/YhatMort_out_of_sample_2.5pct_vs_obs_violin_model_",model.number, "_species_", SPCD.id ,"_remper_corr_", remper.cor.vector[j], ".png"))

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
# save(log_lik_1, r_eff, loo_1, file = paste0(output.folder, "SPCD_stanoutput_full_standardized/LOO_model_",model.number,"remper.corr_",remper.correction, "_species_", SPCD.id ,"_remper_corr_", remper.cor.vector[j], ".Rdata"))
#
# ###------------------------------------
# psis <- loo_1$psis_object
# keep_obs <- sample(1:length(mod.data$y), 100)
# lw <- weights(psis)
#
#
# ppc1 <- bayesplot::ppc_loo_intervals(mod.data$y,
#                                      yrep = as.matrix(yrep.estimates),
#                                      psis_object = psis, subset = keep_obs, order = "median")
# # ppc0 <- rstantools::ppc_loo_pit_overlay(mod.data$y, yrep =  as.matrix(yrep.estimates), lw = lw)
# # ppc3 <- bayesplot::ppc_loo_pit_qq(mod.data$y, yrep =  as.matrix(yrep.estimates), lw = lw)
#
# ############################################################
# # get accuracy of prediction
# # Accuracy
# ext_fit <- rstan::extract(fit.1)
# accuracy.is <- mean(as.vector(yrep.quant$median) == mod.data$y)
# accuracy.oos <- mean(as.vector(yhat.quant$median) == mod.data$ytest)

# # AUC using mltools auc_roc function
# # for in sample data
actuals = mod.data$y
preds = as.vector(psurv.quant$median)
auc.is <-auc_roc(preds, actuals)

# get the range of responses for each sample:
auc.is.list <- list()

auc.is.list <- lapply(1:nrow(p.surv.all), FUN = function(x){
  auc_roc(preds, as.vector(as.numeric(p.surv.all[x,])))
})

AUC.is.df <- do.call(rbind, auc.is.list)
full.auc.quants <- quantile(AUC.is.df[,1], c(0.025, 0.5, 0.975) )



#
#
# ## for out of sampled data
actuals = mod.data$ytest
preds = as.vector(psurv.hat.quant$median)

auc.oos <- auc_roc(preds, actuals)

# get the range of responses for each sample:
auc.oos.list <- list()

auc.oos.list <- lapply(1:nrow(psurv.hat.all), FUN = function(x){
  auc_roc(preds, as.vector(as.numeric(psurv.hat.all[x,])))
})

AUC.oos.df <- do.call(rbind, auc.oos.list)
full.oos.auc.quants <- quantile(AUC.oos.df[,1], c(0.025, 0.5, 0.975) )

#
# # save in one model summary table:
#
model.assessment.df <- data.frame(SPCD = SPCD.id,
                                  model = model.number,
                                  remper.correction = remper.correction,
                                  
                                  # loo estimates
                                  #elpd_loo = loo_1$estimates[1,1],
                                  #p_loo = loo_1$estimates[2,1],
                                  #looic = loo_1$estimates[3,1],
                                  # in sample
                                  auc.insample = auc.is,
                                  auc.insample.median = full.auc.quants[2],
                                  auc.insample.lo = full.auc.quants[1], 
                                  auc.insample.hi = full.auc.quants[3],
                                  # out of sample
                                  auc.oosample = auc.oos, 
                                  auc.oosample.median = full.oos.auc.quants[2],
                                  auc.oosample.lo = full.oos.auc.quants[1], 
                                  auc.oosample.hi = full.oos.auc.quants[3]
                                  #accuracy.oos = accuracy.oos
)
write.csv(model.assessment.df , paste0(output.folder, "SPCD_stanoutput_full_standardized/Accuracy_df_model_",model.number, "_remper_0.5_species_", SPCD.id,"_remper_corr_", remper.cor.vector[j], ".csv" ), row.names = FALSE)

# save the predicted values in sample and out of sampel
saveRDS(psurv.quant, paste0(output.folder, "SPCD_stanoutput_full_standardized/predicted_mort/psurv_quant_", model.number, "remper_", remper.cor.vector[j], "_SPCD_", SPCD.id, ".rds"))
saveRDS(psurv.hat.quant, paste0(output.folder, "SPCD_stanoutput_full_standardized/predicted_mort/psurv_hat_quant_", model.number, "remper_", remper.cor.vector[j], "_SPCD_", SPCD.id, ".rds"))
