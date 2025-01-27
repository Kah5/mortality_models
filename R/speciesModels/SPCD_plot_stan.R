library(loo)
#model.name <- paste0("mort_model3alphaRE_allspp_SPGRPCD_", SPGRPCD.id)
load(paste0("SPCD_standata_general_full_standardized/SPCD_",SPCD.id, "remper_correction_", remper.correction,"model_",model.number, ".Rdata")) # load the species code data
# read in the model fit
fit.1 <- readRDS( paste0(output.folder, "samples/model_",model.number,"_SPCD_",SPCD.id, "_remper_correction_", remper.cor.vector[j], ".RDS"))


species.table <- unique(train.data[,c("SPCD","SPP")])
species.table$COMMON <- FIESTA::ref_species[match(species.table$SPCD, FIESTA::ref_species$SPCD),]$COMMON_NAME
species.table$SPP <- as.character(species.table$SPP)


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
par.names = c("alpha_SPP", colnames(mod.data$xM)) #,

nvariables <- length(par.names)
nvariables
#if(nvariables < 10){
pdf( paste0(output.folder, "/images/traceplots_survival_model_",model.number,"_species_", SPCD.id ,"_remper_corr_", remper.cor.vector[j], ".pdf"))
#specify to save plots in 0x0 grid
par(mfrow = c(8,3))
for (p in 1:length(par.names)) {   
  print(traceplot (fit.1,pars = par.names[p], inc_warmup = FALSE))
}
dev.off()

png(height = (nvariables/2)*3, width = (nvariables/2)*3, units = "in", res = 100, paste0(output.folder, "images/pairs_plot_survival_", model.name, "_species_", SPCD.id , ".png"))
pairs(fit.1, pars = par.names)
dev.off()

species.table <- unique(train.data[,c("SPCD","SPP")])
species.table$COMMON <- FIESTA::ref_species[match(species.table$SPCD, FIESTA::ref_species$SPCD),]$COMMON_NAME
species.table$SPP <- as.character(species.table$SPP)

# get the main effects by species:
fit_ssm_df <- as_draws_df(fit.1) # takes awhile to convert to df
#covariates = c(
# "alpha", "beta_DIA","beta_si", "beta_growth")


# note for model 0 there are no covariates
betas.quant <- summarise_draws(fit_ssm_df %>% select(colnames(mod.data$xM)), 
                               #variables = c(colnames(mod.data$xM)), 
                               median, ~quantile(.x, probs = c(0.025, 0.975)))%>% rename(`ci.lo` = `2.5%`, `ci.hi` = `97.5%`)


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

ggsave(height = 5, width = 10, units = "in",paste0(output.folder, "/images/Estimated_effects_on_survival_model_",model.number, "_species_", SPCD.id ,"_remper_corr_", remper.cor.vector[j], ".png"))

rm(fit.1)

