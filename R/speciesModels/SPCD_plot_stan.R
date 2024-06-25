library(loo)
#model.name <- paste0("mort_model3alphaRE_allspp_SPGRPCD_", SPGRPCD.id)
load(paste0("SPCD_standata_general_full/SPCD_",SPCD.id, "remper_correction_", remper.correction,"model_",model.number, ".Rdata")) # load the species code data
# read in the model fit
fit.1 <- readRDS( paste0("SPCD_stanoutput_full/samples/model_",model.number,"_SPCD_",SPCD.id, "_remper_correction_", remper.cor.vector[j], ".RDS"))


species.table <- unique(train.data[,c("SPCD","SPP")])
species.table$COMMON <- FIESTA::ref_species[match(species.table$SPCD, FIESTA::ref_species$SPCD),]$COMMON_NAME
species.table$SPP <- as.character(species.table$SPP)


names(fit.1) <- c("alpha_SPP", colnames(mod.data$xM),
                  ## in sample predicted status
                  paste0("yrep[",1:mod.data$N, "]"),
                  paste0("psurv[",1:mod.data$N, "]"),
                  ## out of  sample predicted status
                  paste0("yhat[",1:mod.data$Nrep, "]"),
                  ## out of sample predicted prob mor
                  paste0("psurv.hat[",1:mod.data$Nrep, "]"),
                  ## in sample predicted status
                  paste0("log_lik[",1:mod.data$N, "]"),
                  "lp__")
par.names = c("alpha_SPP", colnames(mod.data$xM)) #,

nvariables <- length(par.names)
nvariables
#if(nvariables < 10){
pdf( paste0("SPCD_stanoutput_full/images/traceplots_survival_model_",model.number,"_species_", SPCD.id ,"_remper_corr_", remper.cor.vector[j], ".pdf"))
#specify to save plots in 0x0 grid
par(mfrow = c(8,3))
for (p in 1:length(par.names)) {   
  print(traceplot (fit.1,pars = par.names[p], inc_warmup = FALSE))
}
dev.off()

#png(height = (nvariables/0)*3, width = (nvariables/0)*3, units = "in", res = 100, paste0("SPCD_stanoutput_full/images/pairs_plot_survival_", model.name, "_species_", SPCD.id , ".png"))
#pairs(fit.1, pars = par.names)
#dev.off()

species.table <- unique(train.data[,c("SPCD","SPP")])
species.table$COMMON <- FIESTA::ref_species[match(species.table$SPCD, FIESTA::ref_species$SPCD),]$COMMON_NAME
species.table$SPP <- as.character(species.table$SPP)

# get the main effects by species:
fit_ssm_df <- as.data.frame(fit.1) # takes awhile to convert to df
#covariates = c(
# "alpha", "beta_DIA","beta_si", "beta_growth")


# note for model 0 there are no covariates
cov.estimates <- fit_ssm_df %>% dplyr::select(colnames(mod.data$xM)) 
cov.m <- reshape2::melt(cov.estimates)

betas.quant <- cov.m %>% group_by(variable) %>% summarise(median = quantile(value, 0.5, na.rm =TRUE),
                                                          ci.lo = quantile(value, 0.005, na.rm =TRUE),
                                                          ci.hi = quantile(value, 0.975, na.rm =TRUE))


# clean up the naming structure of this

betas.quant$SPCD <- SPCD.id
betas.quant$Covariate <-  colnames(mod.data$xM)


betas.quant <- left_join(betas.quant, species.table)

# reorder by the value of the covariate
betas.quant <- betas.quant %>% arrange(by = median)
betas.quant$Covariate <- factor(betas.quant$Covariate, levels = betas.quant$Covariate)

# get overlapping zero to color the error bars
betas.quant$`significance` <- ifelse(betas.quant$ci.lo < 0 & betas.quant$ci.hi < 0, "significant", 
                                     ifelse(betas.quant$ci.lo > 0 & betas.quant$ci.hi > 0, "significant", "not overlapping zero"))

ggplot(data = na.omit(betas.quant), aes(x = Covariate, y = median, color = significance))+geom_point()+
  geom_errorbar(data = na.omit(betas.quant), aes(x = Covariate , ymin = ci.lo, ymax = ci.hi, color = significance), width = 0.1)+
  geom_abline(aes(slope = 0, intercept = 0), color = "grey", linetype = "dashed")+facet_wrap(~COMMON)+theme_bw(base_size = 10)+
  theme( axis.text.x = element_text(angle = 65, hjust = 1), panel.grid  = element_blank(), legend.position = "none")+ylab("Effect on survival")+xlab("Parameter")+
  scale_color_manual(values = c("not overlapping zero"="darkgrey", "significant"="black"))

ggsave(height = 5, width = 10, units = "in",paste0("SPCD_stanoutput_full/images/Estimated_effects_on_survival_model_",model.number, "_species_", SPCD.id ,"_remper_corr_", remper.cor.vector[j], ".png"))

##################################
# in-sample pred vs obs plots
##################################
psurv.estimates <- fit_ssm_df %>% dplyr::select( paste0("psurv[",1:mod.data$N, "]")) 
psurv.m <- reshape2::melt(psurv.estimates)

psurv.quant <- psurv.m %>% group_by(variable) %>% summarise(median = quantile(value, 0.5, na.rm =TRUE),
                                                          ci.lo = quantile(value, 0.005, na.rm =TRUE),
                                                          ci.hi = quantile(value, 0.975, na.rm =TRUE))

#hist(psurv.quant$median)
psurv.quant$Mobs <- as.character(mod.data$y)
psurv.quant$COMMON <- unique(species.table$COMMON)


ggplot(psurv.quant, aes(as.character(Mobs), y = median, fill =Mobs ))+geom_violin()+ylab("Median predicted probability of survival")+
  xlab("In-sample Observed Survival Status")+scale_fill_manual(values = c("0" = "#a6611a", "1"= "forestgreen"))+theme_bw()+theme(legend.position = "none")+facet_wrap(~COMMON)
ggsave(height = 4, width = 4, units = "in",paste0("SPCD_stanoutput_full/images/psurv_vs_obs_in_sample_violin_model_",model.number, "_species_", SPCD.id ,"_remper_corr_", remper.cor.vector[j], ".png"))

# compare the observed status to the predicted status
yrep.estimates <- fit_ssm_df %>% dplyr::select( paste0("yrep[",1:mod.data$N, "]")) 
yrep.m <- reshape2::melt(yrep.estimates)

yrep.quant <- yrep.m %>% group_by(variable) %>% summarise(median = quantile(value, 0.5, na.rm =TRUE),
                                                            ci.lo = quantile(value, 0.005, na.rm =TRUE),
                                                            ci.hi = quantile(value, 0.975, na.rm =TRUE))


yrep.quant$Mobs <- as.character(mod.data$y)
yrep.quant$COMMON <- unique(species.table$COMMON)

ggplot(yrep.quant, aes(as.character(Mobs), y = median, fill =Mobs ))+geom_violin()+ylab("Median predicted survival status")+
  xlab("Observed in-sample tree survival status")+scale_fill_manual(values =  c("0" = "#a6611a", "1"= "forestgreen"))+theme_bw()+theme(legend.position = "none")+facet_wrap(~COMMON)
ggsave(height = 4, width = 4, units = "in",paste0("SPCD_stanoutput_full/images/Yrepsurv_insample_median_vs_obs_violin_model_",model.number, "_species_", SPCD.id ,"_remper_corr_", remper.cor.vector[j], ".png"))

ggplot(yrep.quant, aes(as.character(Mobs), y = ci.hi, fill =Mobs ))+geom_violin()+ylab("97.5% quantile of predicted survival status")+
  xlab("Observed in-sample tree survival status")+scale_fill_manual(values =  c("0" = "#a6611a", "1"= "forestgreen"))+theme_bw()+theme(legend.position = "none")+facet_wrap(~COMMON)
ggsave(height = 4, width = 4, units = "in",paste0("SPCD_stanoutput_full/images/Yrepsurv_insample_95pct_vs_obs_violin_model_",model.number, "_species_", SPCD.id ,"_remper_corr_", remper.cor.vector[j], ".png"))

##################################
# out of sample plots
##################################
psurv.hat.estimates <- fit_ssm_df %>% dplyr::select( paste0("psurv.hat[",1:mod.data$Nrep, "]")) 
psurv.hat.m <- reshape2::melt(psurv.hat.estimates)

psurv.hat.quant <- psurv.hat.m %>% group_by(variable) %>% summarise(median = quantile(value, 0.5, na.rm =TRUE),
                                                            ci.lo = quantile(value, 0.005, na.rm =TRUE),
                                                            ci.hi = quantile(value, 0.975, na.rm =TRUE))

#hist(psurv.hat.quant$median)
psurv.hat.quant$Mobs <- as.character(mod.data$ytest)
psurv.hat.quant$COMMON <- unique(species.table$COMMON)


ggplot(psurv.hat.quant, aes(as.character(Mobs), y = median, fill =Mobs ))+geom_violin()+ylab("Median predicted probability of survival")+
  xlab("Out-of-sample Observed Tree Status")+scale_fill_manual(values = c("0" = "#a6611a", "1"= "forestgreen"))+theme_bw()+theme(legend.position = "none")+facet_wrap(~COMMON)
ggsave(height = 4, width = 4, units = "in",paste0("SPCD_stanoutput_full/images/psurv.hat_vs_obs_out_of_sample_violin_model_",model.number, "_species_", SPCD.id ,"_remper_corr_", remper.cor.vector[j], ".png"))

yhat.estimates <- fit_ssm_df %>% dplyr::select( paste0("yhat[",1:mod.data$Nrep, "]")) 
yhat.m <- reshape2::melt(yhat.estimates)

yhat.quant <- yhat.m %>% group_by(variable) %>% summarise(median = quantile(value, 0.5, na.rm =TRUE),
                                                          ci.lo = quantile(value, 0.005, na.rm =TRUE),
                                                          ci.hi = quantile(value, 0.975, na.rm =TRUE))


yhat.quant$Mobs <- as.character(mod.data$ytest)
yhat.quant$COMMON <- unique(species.table$COMMON)

ggplot(yhat.quant, aes(as.character(Mobs), y = median, fill =Mobs ))+geom_violin()+ylab("Median predicted survival status")+
  xlab("Observed in-sample tree status")+scale_fill_manual(values =  c("0" = "#a6611a", "1"= "forestgreen"))+theme_bw()+theme(legend.position = "none")+facet_wrap(~COMMON)
ggsave(height = 4, width = 4, units = "in",paste0("SPCD_stanoutput_full/images/YhatMort_out_of_sample_median_vs_obs_violin_model_",model.number, "_species_", SPCD.id ,"_remper_corr_", remper.cor.vector[j], ".png"))

ggplot(yhat.quant, aes(as.character(Mobs), y = ci.hi, fill =Mobs ))+geom_violin()+ylab("97.5% quantile of predicted survival status")+
  xlab("Observed in-sample tree status")+scale_fill_manual(values =  c("0" = "#a6611a", "1"= "forestgreen"))+theme_bw()+theme(legend.position = "none")+facet_wrap(~COMMON)
ggsave(height = 4, width = 4, units = "in",paste0("SPCD_stanoutput_full/images/YhatMort_out_of_sample_95pct_vs_obs_violin_model_",model.number, "_species_", SPCD.id ,"_remper_corr_", remper.cor.vector[j], ".png"))

###########
# Use the Loo package to compute PSIS-LOO and check diagnositcs
# this may take awhile..
# Extract pointwise log-likelihood
log_lik_1 <- extract_log_lik(fit.1, merge_chains = FALSE)

# provide relative effective sample sizes, to bettin estimate PSIS 
r_eff <- relative_eff(exp(log_lik_1))

# preferably use more than 0 cores (as many cores as possible)
# will use value of 'mc.cores' option if cores is not specified
loo_1 <- loo(log_lik_1, r_eff = r_eff, save_psis = TRUE)
print(loo_1)
# save the loo_1 object

save(log_lik_1, r_eff, loo_1, file = paste0("SPCD_stanoutput_full/LOO_model_",model.number,"remper.corr_",remper.correction, "_species_", SPCD.id ,"_remper_corr_", remper.cor.vector[j], ".Rdata"))

###------------------------------------
psis <- loo_1$psis_object
keep_obs <- sample(1:length(mod.data$y), 100)
lw <- weights(psis)


ppc1 <- bayesplot::ppc_loo_intervals(mod.data$y, 
                                     yrep = as.matrix(yrep.estimates), 
                                     psis_object = psis, subset = keep_obs, order = "median") 
# ppc0 <- rstantools::ppc_loo_pit_overlay(mod.data$y, yrep =  as.matrix(yrep.estimates), lw = lw)
# ppc3 <- bayesplot::ppc_loo_pit_qq(mod.data$y, yrep =  as.matrix(yrep.estimates), lw = lw)

############################################################
# get accuracy of prediction
# Accuracy
ext_fit <- rstan::extract(fit.1)
accuracy.is <- mean(as.vector(yrep.quant$median) == mod.data$y)
accuracy.oos <- mean(as.vector(yhat.quant$median) == mod.data$ytest)

# AUC using mltools auc_roc function
# for in sample data
actuals = mod.data$y
preds = as.vector(psurv.quant$median)
auc.is <-auc_roc(preds, actuals)  


## for out of sampled data
actuals = mod.data$ytest
preds = as.vector(psurv.hat.quant$median)

auc.oos <- auc_roc(preds, actuals)  

# save in one model summary table:

model.assessment.df <- data.frame(SPCD = SPCD.id,
                          model = model.number, 
                          remper.correction = remper.correction, 
                          
                          # loo estimates
                          elpd_loo = loo_1$estimates[1,1],
                          p_loo = loo_1$estimates[2,1],
                          looic = loo_1$estimates[3,1],
                          # in sample
                          auc.insample = auc.is,
                          accuracy.is = accuracy.is,
                          # out of sample 
                          auc.oosample = auc.oos,
                          accuracy.oos = accuracy.oos
                          )
write.csv(model.assessment.df , paste0("SPCD_stanoutput_full/Accuracy_df_model_",model.number, "_remper_0.5_species_", SPCD.id,"_remper_corr_", remper.cor.vector[j], ".csv" ), row.names = FALSE)


