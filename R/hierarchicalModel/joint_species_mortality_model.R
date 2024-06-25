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

options(mc.cores = parallel::detectCores())

# get the complete spcies list
cleaned.data <- readRDS( "data-store/data/iplant/home/kellyheilman/mort_data/cleaned.data.mortality.TRplots.RDS")
cleaned.data <- cleaned.data %>% dplyr::select(state, county, pltnum, cndtn, point, tree, PLOT.ID, cycle, spp, dbhcur, dbhold, damage, Species, SPCD,
                                               remper, LAT_FIADB, LONG_FIADB, elev, DIA_DIFF, annual.growth, M, relative.growth, si, physio:RD) %>% distinct()

nspp <- cleaned.data %>% group_by(SPCD) %>% summarise(n = n(), 
                                                      pct = n/nrow(cleaned.data)) %>% arrange (desc(`pct`))

nspp$cumulative.pct <- cumsum(nspp$pct)



# link up to the species table:
nspp$COMMON <- FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD),]$COMMON

#View(nspp)

nspp[1:17,]$COMMON


# read in the test data for all the species
spp.table <- data.frame(SPCD.id = nspp[1:17,]$SPCD, 
                        spp = 1:17, 
                        COMMON = nspp[1:17,]$COMMON)
model.no <- 6
SPCD.id <- spp.table[1,]$SPCD.id

xM.list <- xMrep.list <-y.list <- y.test.list <- nSPP.list <- nSPP.rep.list <- list()
for(i in 1:17){
  # SPCD.id
  SPCD.id <- spp.table[i,]$SPCD.id
  load(paste0("data-store/data/iplant/home/kellyheilman/SPCD_standata_general_full/SPCD_",SPCD.id, "remper_correction_0.5model_",model.no, ".Rdata")) # load the species code data
  #mod.data$K <- ncol(mod.data$xM)
  
  
  xM.list[[i]] <- mod.data$xM
  xMrep.list[[i]] <- mod.data$xMrep
  y.list[[i]] <- mod.data$y
  y.test.list[[i]] <- mod.data$ytest
  nSPP.list[[i]] <- rep(i, length(mod.data$y))
  nSPP.rep.list[[i]] <- rep(i, length(mod.data$ytest))
}




mod.data.full <-  list(xM = do.call(rbind, xM.list),
                       xMrep = do.call(rbind, xMrep.list),
                       y = unlist(y.list), 
                       ytest = unlist(y.test.list), 
                       SPP = unlist(nSPP.list), 
                       SPPrep = unlist( nSPP.rep.list))
mod.data.full$K <- ncol(mod.data.full$xM)
mod.data.full$N <- length(mod.data.full$y)
#mod.data.full$Nspp <- 3
mod.data.full$Nspp <- 17
mod.data.full$Nrep <- length(mod.data.full$ytest)
saveRDS (mod.data.full, paste0("SPCD_stanoutput_joint/all_SPCD_model_", model.no,".RDS"))
mod.data.full <- readRDS ( paste0("SPCD_stanoutput_joint/all_SPCD_model_", model.no,".RDS"))


saveRDS(spp.table, "SPCD_stanoutput_joint/spp.table.rds")
spp.table <- readRDS("SPCD_stanoutput_joint/spp.table.rds")
#species.table <- data.frame(spp = unlist(nSPP.list))
#species.table <- data.frame(spp = unlist(nSPP.list))

unique(mod.data.full$SPP)
mod.data.full$y 

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



# null model:

num_cores <-  parallel::detectCores()
# null model:
start.time <- Sys.time()
fit.1 <- stan(file = "mort_model_general_heiarchical.stan" , 
              data = mod.data.full,
              iter = 3000, 
              chains = 2, 
              verbose=FALSE, 
              ##control =  list(max_treedepth = 15),#list(adapt_delta = 0.99, stepsize = 0.5, max_treedepth = 15),#, stepsize = 0.01, max_treedepth = 15),
              #sample_file = model.name, 
              #adapt_delta = 0.99, 
              pars =c("alpha_SPP", "u_beta", # the species-specific params
                      "alpha", "mu_beta",
                      "y_rep", "mMrep",## in sample predictions
                      "y_hat", "mMhat", ## out of sample predictions
                      "log_lik")) #, "y_hat", 
end.time <- Sys.time()
# Calculate elapsed time in seconds
elapsed_time <- as.numeric(difftime(end.time, start.time, units = "secs"))

# Calculate core hours
core_hours <- elapsed_time * num_cores / 3600  # convert seconds to hours
remper.correction <- 0.5

time.diag <- data.frame(SPP = mod.data.full$SPP, 
                        y = mod.data.full$y) %>% group_by(SPP) %>% summarise(n()) %>% ungroup()%>%
  mutate(total = sum(`n()`)) %>% group_by(SPP) %>% mutate(weight = `n()`/total) %>% ungroup()%>%
  mutate(total.core.hours = core_hours, 
         total.elapsed.time = elapsed_time, 
         cores = num_cores, 
         model = model.no, 
         remper = 0.5) %>% group_by(SPP) %>%
  mutate(core.hours = weight*total.core.hours)%>%rename(`spp` = "SPP") %>%
  left_join(spp.table)


# time.diag <- data.frame(model = model.no, 
#                         SPCD = SPCD.id, 
#                         remper = 0.5,
#                         core.hours = core_hours, 
#                         elapsed.time = elapsed_time, 
#                         cores = num_cores)

write.csv(time.diag, paste0("SPCD_stanoutput_joint/joint_model_time_diag_SPCD_joint_model_", model.no, "_remper_", remper.correction,".csv"))

# get the sampler diagnostics and save:
# Time difference of 22.91259 mins for all 17 species but only 30% of the data

saveRDS(fit.1, paste0("SPCD_stanoutput_joint/samples/model_",model.no,"all_SPCD_remper_correction_0.5.RDS"))

joint.samples <- as_draws_df(fit.1)

alpha.p <- subset_draws(joint.samples, variable = "alpha", chain = 1:2, iteration = 500:1500)
alpha.spp <- subset_draws(joint.samples, variable = "alpha_SPP", chain = 1:2, iteration = 500:1500)
saveRDS(alpha.p, paste0("SPCD_stanoutput_joint/samples/alpha.p_model_",model.no,"_1000samples.rds"))
saveRDS(alpha.spp, paste0("SPCD_stanoutput_joint/samples/alpha.spp_model_",model.no,"_1000samples.rds"))

beta.p <- subset_draws(joint.samples, variable = "mu_beta", chain = 1:2, iteration = 500:1500)
bet0a.spp <- subset_draws(joint.samples, variable = "u_beta", chain = 1:2, iteration = 500:1500)
saveRDS(beta.p, paste0("SPCD_stanoutput_joint/samples/beta_model_",model.no,"_1000samples.rds"))
saveRDS(bet0a.spp, paste0("SPCD_stanoutput_joint/samples/u_betas_model_",model.no,"_1000samples.rds"))

yrep <- subset_draws(joint.samples, variable = "y_rep", chain = 1:2, iteration = 500:1500)
yhat <- subset_draws(joint.samples, variable = "y_hat", chain = 1:2, iteration = 500:1500)
saveRDS(yrep, paste0("SPCD_stanoutput_joint/samples/yrep_model_",model.no,"_1000samples.rds"))
saveRDS(yhat, paste0("SPCD_stanoutput_joint/samples/yhat_model_",model.no,"_1000samples.rds"))

yhat <- readRDS(paste0("SPCD_stanoutput_joint/samples/yhat_model_",model.no,"_1000samples.rds"))
yhat <- readRDS(paste0("SPCD_stanoutput_joint/samples/yrep_model_",model.no,"_1000samples.rds"))

mMrep <- subset_draws(joint.samples, variable = "mMrep", chain = 1:2, iteration = 500:1500)
mMhat <- subset_draws(joint.samples, variable = "mMhat", chain = 1:2, iteration = 500:1500)
saveRDS(mMrep, paste0("SPCD_stanoutput_joint/samples/mMrep_model_",model.no,"_1000samples.rds"))
saveRDS(mMhat, paste0("SPCD_stanoutput_joint/samples/mMhat_model_",model.no,"_1000samples.rds"))

log_lik <- subset_draws(joint.samples, variable = "log_lik", chain = 1:2, iteration = 500:1500)
saveRDS(log_lik, paste0("SPCD_stanoutput_joint/samples/log_lik_model_",model.no,"_1000samples.rds"))




growth.params = c(paste0("u_beta[", 1:17, ",1]"), "mu_beta[1]")
dia.params = c(paste0("u_beta[", 1:17, ",2]"), "mu_beta[2]")
RD.params = c(paste0("u_beta[", 1:17, ",3]"), "mu_beta[3]")
ba.params = c(paste0("u_beta[", 1:17, ",4]"), "mu_beta[4]")

traceplot(fit.1, pars = "alpha_SPP")
traceplot(fit.1, pars = growth.params)
traceplot(fit.1, pars = dia.params)
traceplot(fit.1, pars = RD.params)
traceplot(fit.1, pars = ba.params)
#traceplot(fit.1, pars = "u_beta[6,1]")
traceplot(fit.1, pars = "mu_beta")

beta.names <- data.frame(parameter = colnames(mod.data.full$xM), 
                         param.no = 1:length(colnames(mod.data.full$xM)))

nvariables <- length(names(fit.1))

nvariables
#if(nvariables < 10){
pdf( paste0("SPCD_stanoutput/images/traceplots_mortality_model_6_all.species.pdf"))
#specify to save plots in 0x0 grid
par(mfrow = c(8,3))
for (p in 1:nvariables) {   
  print(traceplot (fit.1,pars = names(fit.1)[p], inc_warmup = FALSE))
}
dev.off()

# names(fit.1) <- c(paste0("alpha_SPP_", 1:17),
#                   colnames(mod.data$xM),
#                   # ## in sample predicted status
#                   paste0("yrep[",1:mod.data$N, "]"),
#                   paste0("pmort[",1:mod.data$N, "]"),
#                   ## out of  sample predicted status
#                   paste0("yhat[",1:mod.data$Nrep, "]"),
#                   ## out of sample predicted prob mor
#                   paste0("pmort.hat[",1:mod.data$Nrep, "]"),
#                   ## in sample predicted status
#                   paste0("log_lik[",1:mod.data$N, "]"),
#                   "lp__")
#---------------------------------------------------------------------------------------
# look at outputs and plot traces
bet0a.spp <- readRDS( paste0("SPCD_stanoutput_joint/u_betas_model_",model.no,"_1000samples.rds"))

betas.quant <- summarise_draws(bet0a.spp, median, ~quantile(.x, probs = c(0.025, 0.975)))%>% rename(`ci.lo` = `2.5%`, `ci.hi` = `97.5%`)

spp.ids <- rep(1:17, 48)
param.ids <- rep(1:48, each = 17)
# # note for model 0 there are no covariates
# cov.estimates <- fit_ssm_df %>% dplyr::select(paste0("u_beta[", spp.ids, ",", param.ids, "]")) 
# cov.m <- reshape2::melt(cov.estimates)
# 
# betas.quant <- cov.m %>% group_by(variable) %>% summarise(median = quantile(value, 0.5, na.rm =TRUE),
#                                                           ci.lo = quantile(value, 0.025, na.rm =TRUE),
#                                                           ci.hi = quantile(value, 0.975, na.rm =TRUE))


# clean up the naming structure of this
betas.quant$spp <- rep(1:17, 48)
betas.quant$param.no <- rep(1:48, each = 17)

betas.quant <- left_join(betas.quant, beta.names)

betas.quant <- left_join(betas.quant, spp.table)

# reorder by the value of the covariate
betas.quant <- betas.quant %>% arrange(by = median) 
# betas.quant$parameter <- factor(betas.quant$parameter, levels = betas.quant$parameter)

# get overlapping zero to color the error bars
betas.quant$`significance` <- ifelse(betas.quant$ci.lo < 0 & betas.quant$ci.hi < 0, "significant", 
                                     ifelse(betas.quant$ci.lo > 0 & betas.quant$ci.hi > 0, "significant", "not overlapping zero"))

ggplot(data = na.omit(betas.quant), aes(x = parameter, y = median, color = significance))+geom_point()+
  geom_errorbar(data = na.omit(betas.quant), aes(x = parameter , ymin = ci.lo, ymax = ci.hi, color = significance), width = 0.1)+
  geom_abline(aes(slope = 0, intercept = 0), color = "grey", linetype = "dashed")+facet_wrap(~COMMON)+theme_bw(base_size = 10)+
  theme( axis.text.x = element_text(angle = 65, hjust = 1), panel.grid  = element_blank(), legend.position = "none")+ylab("Effect on mortality")+xlab("Parameter")+
  scale_color_manual(values = c("not overlapping zero"="darkgrey", "significant"="black"))



ggplot(data = na.omit(betas.quant), aes(x = COMMON, y = median, color = significance))+geom_point()+
  geom_errorbar(data = na.omit(betas.quant), aes(x = COMMON , ymin = ci.lo, ymax = ci.hi, color = significance), width = 0.1)+
  geom_abline(aes(slope = 0, intercept = 0), color = "grey", linetype = "dashed")+facet_wrap(~parameter)+theme_bw(base_size = 10)+
  theme( axis.text.x = element_text(angle = 65, hjust = 1), panel.grid  = element_blank(), legend.position = "none")+ylab("Effect on mortality")+xlab("Parameter")+
  scale_color_manual(values = c("not overlapping zero"="darkgrey", "significant"="black"))
ggsave(height = 10, width = 10, units = "in",paste0("SPCD_stanoutput_joint/images/Estimated_effects_on_mortality_model_model6_all_species_betas.png"))

# get the population estimates

mubetas.estimates <- readRDS( paste0("SPCD_stanoutput_joint/beta_model_6_1000samples.RDS"))

mubetas.quant <- summarise_draws(mubetas.estimates, median, ~quantile(.x, probs = c(0.025, 0.975)))%>% rename(`ci.lo` = `2.5%`, `ci.hi` = `97.5%`)

# mubetas.estimates <- fit_ssm_df %>% dplyr::select(paste0("mu_beta[",1:48,"]")) 
# mub.m <- reshape2::melt(mubetas.estimates)
# 
# 
# mubetas.quant <- mub.m %>% group_by(variable) %>% summarise(median = quantile(value, 0.5, na.rm =TRUE),
#                                                             ci.lo = quantile(value, 0.005, na.rm =TRUE),
#                                                             ci.hi = quantile(value, 0.975, na.rm =TRUE))
mubetas.quant$spp <- 18
mubetas.quant$param.no <- 1:48

mubetas.quant <- left_join(mubetas.quant, beta.names)
main.table <- data.frame(SPCD.id = "1000", 
                         COMMON = "population", 
                         spp = 18)

mubetas.quant <- left_join(mubetas.quant, main.table)


mubetas.quant$`significance` <- ifelse(mubetas.quant$ci.lo < 0 & mubetas.quant$ci.hi < 0, "significant", 
                                       ifelse(mubetas.quant$ci.lo > 0 & mubetas.quant$ci.hi > 0, "significant", "not overlapping zero"))

mubetas.quant <- mubetas.quant %>% arrange(by = median)
mubetas.quant$parameter <- factor(mubetas.quant$parameter, levels = mubetas.quant$parameter)

ggplot(data = na.omit(mubetas.quant), aes(x = parameter, y = median, color = significance))+geom_point()+
  geom_errorbar(data = na.omit(mubetas.quant), aes(x = parameter , ymin = ci.lo, ymax = ci.hi, color = significance), width = 0.1)+
  geom_abline(aes(slope = 0, intercept = 0), color = "grey", linetype = "dashed")+facet_wrap(~COMMON)+theme_bw(base_size = 10)+
  theme( axis.text.x = element_text(angle = 65, hjust = 1), panel.grid  = element_blank(), legend.position = "none")+ylab("Effect on mortality")+xlab("Parameter")+
  scale_color_manual(values = c("not overlapping zero"="darkgrey", "significant"="black"))
ggsave(height = 5, width = 10, units = "in",paste0("SPCD_stanoutput_joint/images/Estimated_effects_on_mortality_model_model6_all_species_population_betas.png"))



### combine the betas together
all.joint.betas <- rbind(mubetas.quant, betas.quant) 
# reorder factors for better plotting
all.joint.betas$COMMON <- factor(all.joint.betas$COMMON, levels = c("population", unique(betas.quant$COMMON)[order(unique(betas.quant$COMMON))]))
all.joint.betas$parameter <- factor(all.joint.betas$parameter, levels = beta.names$parameter)

ggplot(data = na.omit(all.joint.betas), aes(x = COMMON, y = median, color = significance, shape = COMMON %in% "population"))+geom_point()+
  geom_errorbar(data = na.omit(all.joint.betas), aes(x = COMMON , ymin = ci.lo, ymax = ci.hi, color = significance), width = 0.1)+
  geom_abline(aes(slope = 0, intercept = 0), color = "grey", linetype = "dashed")+facet_wrap(~parameter)+theme_bw(base_size = 10)+
  theme( axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.grid  = element_blank(), legend.position = "none")+ylab("Effect on survival")+xlab("Parameter")+
  scale_color_manual(values = c("not overlapping zero"="darkgrey", "significant"="black"))+
scale_shape_manual(values = c("TRUE" = 15, "FALSE" = 19))

ggsave(height = 10, width = 11, units = "in",paste0("SPCD_stanoutput_joint/images/Estimated_betas_all_model6.png"))



# get the species alpha estimates
alphas.estimates <- readRDS( paste0("SPCD_stanoutput_joint/alpha.spp_model_",model.no,"_1000samples.rds"))

alphas.quant <- summarise_draws(alphas.estimates, median, ~quantile(.x, probs = c(0.025, 0.975)))%>% rename(`ci.lo` = `2.5%`, `ci.hi` = `97.5%`)

# alphas.estimates <- fit_ssm_df %>% dplyr::select(paste0("alpha_SPP[",1:17,"]")) 
# mub.m <- reshape2::melt(alphas.estimates)
# 
# 
# alphas.quant <- mub.m %>% group_by(variable) %>% summarise(median = quantile(value, 0.5, na.rm =TRUE),
#                                                            ci.lo = quantile(value, 0.005, na.rm =TRUE),
#                                                            ci.hi = quantile(value, 0.975, na.rm =TRUE))
 alphas.quant$spp <- 1:17


alphas.quant <- left_join(alphas.quant, spp.table)


alphas.quant$`significance` <- ifelse(alphas.quant$ci.lo < 0 & alphas.quant$ci.hi < 0, "significant", 
                                      ifelse(alphas.quant$ci.lo > 0 & alphas.quant$ci.hi > 0, "significant", "not overlapping zero"))

#alphas.quant <- alphas.quant %>% arrange(by = median)
#alphas.quant$parameter <- factor(alphas.quant$parameter, levels = alphas.quant$parameter)

ggplot(data = na.omit(alphas.quant), aes(x = COMMON, y = median, color = significance))+geom_point()+
  geom_errorbar(data = na.omit(alphas.quant), aes(x = COMMON , ymin = ci.lo, ymax = ci.hi, color = significance), width = 0.1)+
  geom_abline(aes(slope = 0, intercept = 0), color = "grey", linetype = "dashed")+theme_bw(base_size = 10)+
  theme( axis.text.x = element_text(angle = 65, hjust = 1), panel.grid  = element_blank(), legend.position = "none")+ylab("alpha estimate")+xlab("SPECIES")+
  scale_color_manual(values = c("not overlapping zero"="darkgrey", "significant"="black"))
ggsave(height = 5, width = 10, units = "in",paste0("SPCD_stanoutput_joint/images/Estimated_alpha_SPP_model6_all_species_alphas.png"))

# combine the population and the species level betas

## get population level alpha estimates:
alphas.pop.estimates <- readRDS( paste0("SPCD_stanoutput_joint/alpha.p_model_",model.no,"_1000samples.rds"))

alphas.pop.quant <- summarise_draws(alphas.pop.estimates, median, ~quantile(.x, probs = c(0.025, 0.975)))%>% rename(`ci.lo` = `2.5%`, `ci.hi` = `97.5%`)

# alphas.estimates <- fit_ssm_df %>% dplyr::select(paste0("alpha_SPP[",1:17,"]")) 
# mub.m <- reshape2::melt(alphas.estimates)
# 
# 
# alphas.quant <- mub.m %>% group_by(variable) %>% summarise(median = quantile(value, 0.5, na.rm =TRUE),
#                                                            ci.lo = quantile(value, 0.005, na.rm =TRUE),
#                                                            ci.hi = quantile(value, 0.975, na.rm =TRUE))
alphas.pop.quant$spp <- 18


alphas.pop.quant <- left_join(alphas.pop.quant, main.table)


alphas.pop.quant$`significance` <- ifelse(alphas.pop.quant$ci.lo < 0 & alphas.pop.quant$ci.hi < 0, "significant", 
                                      ifelse(alphas.pop.quant$ci.lo > 0 & alphas.pop.quant$ci.hi > 0, "significant", "not overlapping zero"))

#alphas.quant <- alphas.quant %>% arrange(by = median)
#alphas.quant$parameter <- factor(alphas.quant$parameter, levels = alphas.quant$parameter)

ggplot(data = na.omit(alphas.pop.quant), aes(x = COMMON, y = median, color = significance))+geom_point()+
  geom_errorbar(data = na.omit(alphas.pop.quant), aes(x = COMMON , ymin = ci.lo, ymax = ci.hi, color = significance), width = 0.1)+
  geom_abline(aes(slope = 0, intercept = 0), color = "grey", linetype = "dashed")+theme_bw(base_size = 10)+
  theme( axis.text.x = element_text(angle = 65, hjust = 1), panel.grid  = element_blank(), legend.position = "none")+ylab("alpha estimate")+xlab("SPECIES")+
  scale_color_manual(values = c("not overlapping zero"="darkgrey", "significant"="black"))


### combine the alphas together
all.joint.alphas <- rbind(alphas.pop.quant, alphas.quant) 
all.joint.alphas$COMMON <- factor(all.joint.alphas$COMMON, levels = c("population", alphas.quant$COMMON[order(unique(alphas.quant$COMMON))]))

ggplot(data = na.omit(all.joint.alphas), aes(x = COMMON, y = median, color = significance, shape = COMMON %in% "population"))+geom_point()+
  geom_errorbar(data = na.omit(all.joint.alphas), aes(x = COMMON , ymin = ci.lo, ymax = ci.hi, color = significance), width = 0.1)+
  geom_abline(aes(slope = 0, intercept = 0), color = "grey", linetype = "dashed")+theme_bw(base_size = 10)+
  theme( axis.text.x = element_text(angle = 65, hjust = 1), panel.grid  = element_blank(), legend.position = "none")+ylab("alpha estimate")+xlab("SPECIES")+
  scale_color_manual(values = c("not overlapping zero"="darkgrey", "significant"="black"))+
  scale_shape_manual(values = c("TRUE" = 15, "FALSE" = 19))

ggsave(height = 5, width = 6, units = "in",paste0("SPCD_stanoutput_joint/images/Estimated_alpha_all_model6.png"))

## combine alphas and betas together to make the figure for the paper

all.joint.alphas <- all.joint.alphas %>% mutate(parameter = "alpha") %>% select(variable, median, ci.lo, ci.hi, spp, parameter, SPCD.id, COMMON, significance)
all.joint.betas <- all.joint.betas %>% select(-param.no)
all.params <- rbind(all.joint.betas, all.joint.alphas)

ggplot(data = na.omit(all.params), aes(x = COMMON, y = median, color = significance, shape = COMMON %in% "population"))+geom_point()+
  geom_errorbar(data = na.omit(all.params), aes(x = COMMON , ymin = ci.lo, ymax = ci.hi, color = significance), width = 0.1)+
  geom_abline(aes(slope = 0, intercept = 0), color = "grey", linetype = "dashed")+facet_wrap(~parameter)+theme_bw(base_size = 10)+
  theme( axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.grid  = element_blank(), legend.position = "none")+ylab("Effect on survival")+xlab("Parameter")+
  scale_color_manual(values = c("not overlapping zero"="darkgrey", "significant"="black"))+
  scale_shape_manual(values = c("TRUE" = 15, "FALSE" = 19))

ggsave(height = 10, width = 11, units = "in",dpi = 350,paste0("SPCD_stanoutput_joint/images/Estimated_parameters_model6_joint.png"))

# compare the observed status to the predicted status
yrep.estimates <- readRDS( paste0("SPCD_stanoutput_joint/yrep_model_",model.no,"_1000samples.rds"))

yrep.quant <- summarise_draws(yrep.estimates, median, ~quantile(.x, probs = c(0.025, 0.975)))%>% rename(`ci.lo` = `2.5%`, `ci.hi` = `97.5%`)

yrep.quant$Mobs <- as.character(mod.data.full$y)
yrep.quant$spp <- mod.data.full$SPP
yrep.quant <- left_join(yrep.quant, spp.table)

ggplot(yrep.quant, aes(as.character(Mobs), y = median, fill =Mobs ))+geom_violin()+ylab("Median predicted survival status")+
  xlab("Observed in-sample tree survival status")+scale_fill_manual(values =  c("0" = "#a6611a", "1"= "forestgreen"))+theme_bw()+theme(legend.position = "none")+facet_wrap(~COMMON)
ggsave(height = 8, width = 8, units = "in",paste0("SPCD_stanoutput_joint/images/Yrepsurv_insample_median_vs_obs_violin_model_",model.no, "_species_joint_model_remper_corr_", remper.correction, ".png"))

ggplot(yrep.quant, aes(as.character(Mobs), y = ci.hi, fill =Mobs ))+geom_violin()+ylab("97.5% quantile of predicted survival status")+
  xlab("Observed in-sample tree survival status")+scale_fill_manual(values =  c("0" = "#a6611a", "1"= "forestgreen"))+theme_bw()+theme(legend.position = "none")+facet_wrap(~COMMON)
ggsave(height = 8, width = 8, units = "in",paste0("SPCD_stanoutput_full/images/Yrepsurv_insample_95pct_vs_obs_violin_model_",model.no, "_species_joint_model_remper_corr_", remper.correction, ".png"))


# compare the observed status to the predicted status
yhat.estimates <- readRDS( paste0("SPCD_stanoutput_joint/yhat_model_",model.no,"_1000samples.rds"))

yhat.quant <- summarise_draws(yhat.estimates, median, ~quantile(.x, probs = c(0.025, 0.975)))%>% rename(`ci.lo` = `2.5%`, `ci.hi` = `97.5%`)

yhat.quant$Mobs <- as.character(mod.data.full$ytest)
yhat.quant$spp <- mod.data.full$SPPrep
yhat.quant <- left_join(yhat.quant, spp.table)

ggplot(yhat.quant, aes(as.character(Mobs), y = median, fill =Mobs ))+geom_violin()+ylab("Median predicted survival status")+
  xlab("Observed out-of-sample tree survival status")+scale_fill_manual(values =  c("0" = "#a6611a", "1"= "forestgreen"))+theme_bw()+theme(legend.position = "none")+facet_wrap(~COMMON)
ggsave(height = 8, width = 8, units = "in",paste0("SPCD_stanoutput_joint/images/Yhatsurv_oossample_median_vs_obs_violin_model_",model.no, "_species_joint_model_remper_corr_", remper.correction, ".png"))

ggplot(yhat.quant, aes(as.character(Mobs), y = ci.hi, fill =Mobs ))+geom_violin()+ylab("97.5% quantile of predicted survival status")+
  xlab("Observed out-of-sample tree survival status")+scale_fill_manual(values =  c("0" = "#a6611a", "1"= "forestgreen"))+theme_bw()+theme(legend.position = "none")+facet_wrap(~COMMON)
ggsave(height = 8, width = 8, units = "in",paste0("SPCD_stanoutput_joint/images/Yhatsurv_oossample_95pct_vs_obs_violin_model_",model.no, "_species_joint_model_remper_corr_", remper.correction, ".png"))


# survival probability for in-sample data
psurv.estimates <- readRDS( paste0("SPCD_stanoutput_joint/mMrep_model_",model.no,"_1000samples.rds"))

psurv.quant <- summarise_draws(psurv.estimates, median, ~quantile(.x, probs = c(0.025, 0.975)))%>% rename(`ci.lo` = `2.5%`, `ci.hi` = `97.5%`)

#hist(psurv.quant$median)

psurv.quant$Mobs <- as.character(mod.data.full$y)
psurv.quant$spp <- mod.data.full$SPP
psurv.quant <- left_join(psurv.quant, spp.table)

ll.train.pmort <- psurv.quant
saveRDS(ll.train.pmort, "SPCD_stanoutput_joint/ll.train.pmort.RDS")


# survival probability for held out-sample data
psurv.hat.estimates <- readRDS( paste0("SPCD_stanoutput_joint/mMhat_model_",model.no,"_1000samples.rds"))

psurv.hat.quant <- summarise_draws(psurv.hat.estimates, median, ~quantile(.x, probs = c(0.025, 0.975)))%>% rename(`ci.lo` = `2.5%`, `ci.hi` = `97.5%`)


psurv.hat.quant$Mobs <- as.character(mod.data.full$ytest)
psurv.hat.quant$spp <- mod.data.full$SPPrep
psurv.hat.quant <- left_join(psurv.hat.quant, spp.table)


ll.test.pmort <- psurv.hat.quant

saveRDS(ll.test.pmort, "SPCD_stanoutput_joint/ll.test.pmort.RDS")

############################################################
# get accuracy of prediction
# Accuracy
#ext_fit <- rstan::extract(fit.1)
accuracy.is <- mean(as.vector(yrep.quant$median) == mod.data.full$y)
accuracy.oos <- mean(as.vector(yhat.quant$median) == mod.data.full$ytest)

# AUC using mltools auc_roc function
# for in sample data
actuals = mod.data.full$y
preds = as.vector(psurv.quant$median)
auc.is <-auc_roc(preds, actuals)  


## for out of sampled data
actuals = mod.data.full$ytest
preds = as.vector(psurv.hat.quant$median)

auc.oos <- auc_roc(preds, actuals)  

# now get the species-level auc and accuracies
oos.all <- data.frame(actuals = mod.data.full$ytest,
                      preds = as.vector(psurv.hat.quant$median), 
                      ypreds = as.vector(yhat.quant$median),
                      spp = mod.data.full$SPPrep)
oos.all <- left_join(oos.all, spp.table)


auc.oos.spp <- oos.all %>% group_by(COMMON, SPCD.id, spp) %>% summarise(auc.oosample = auc_roc(preds, actuals), 
                                                                        accuracy.oos =  mean(as.vector(ypreds) == actuals))

oos.population <- data.frame(COMMON= "population", 
                            SPCD.id = 1000, 
                            spp = 18, 
                            auc.oosample = auc.oos, 
                            accuracy.oos = accuracy.oos)
auc.oos.spp <- rbind(auc.oos.spp, oos.population)


# get the species level in sample summaries
is.all <- data.frame(actuals = mod.data.full$y,
                      preds = as.vector(psurv.quant$median),
                      ypreds = as.vector(yrep.quant$median),
                      spp = mod.data.full$SPP)
is.all <- left_join(is.all, spp.table)


auc.is.spp <- is.all %>% group_by(COMMON, SPCD.id, spp) %>% summarise(auc.insample = auc_roc(preds, actuals),
                                                                      accuracy.is =  mean(as.vector(ypreds) == actuals))
is.population <- data.frame(COMMON= "population", 
                            SPCD.id = 1000, 
                            spp = 18, 
                            auc.insample = auc.is, 
                            accuracy.is = accuracy.is)
auc.is.spp <- rbind(auc.is.spp, is.population)
auc.df <- left_join(auc.is.spp, auc.oos.spp)


# save in one model summary table:
model.assessment.df <- auc.df %>% rename(`SPCD`="SPCD.id") %>%
  mutate(model = model.no, 
         remper.correction = remper.correction, 
         elpd_loo =NA,   p_loo =NA,    looic = NA)%>% select(SPCD, model, remper.correction, elpd_loo, 
                                                             p_loo, looic, auc.insample, accuracy.is, 
                                                          auc.oosample, accuracy.oos)

write.csv(model.assessment.df , paste0("SPCD_stanoutput_full/Accuracy_df_model_",model.no, "_remper_0.5_species_joint_model_remper_corr_",remper.correction, ".csv" ), row.names = FALSE)
write.csv(model.assessment.df, "SPCD_stanoutput_joint/Accuracy_df_model_6_remper_0.5_species_joint_model_remper_corr_0.5.csv", row.names = FALSE)

# make a big plot of the predictions

library(maps)
library(mapdata)
states <- map_data("state")
#9=CT, 25=MA, 33=NH, 23=ME, 50=VT, 44=RI, 42=PA, 39=OH, 54=WV
state_sub <- filter(states, region %in% c("connecticut","maine","new hampshire","vermont","new york", "new jersey",
                                          "rhode island","pennsylvania","ohio","west virginia", "massachusetts", "virginia", "delaware", 
                                          "north carolina", "kentucky", "tennessee", "michigan", "indiana", "district of columbia", "south carolina", 
                                          "georgia", "maryland"))

canada <- map_data("worldHires", "Canada")

# plot distribution 
ll.test.pmort  <- ll.test.pmort  %>%mutate(`p(mort)` = 1- median) %>%  mutate(Mort.quantiles = cut(`p(mort)`, 
                                                                                                   breaks = c(0,0.01,0.05, 0.1, 0.2, 0.3, 0.40, 0.50, 1), 
                                                                                                   include.lowest=TRUE))

ggplot() +
  geom_polygon(data = canada, 
               aes(x=long, y=lat, group = group), 
               color = "black", fill = "white") +
  geom_polygon(data = state_sub, 
               aes(x=long, y=lat, group = group), 
               color = "black", fill = "white") +
  
  geom_point(data = ll.test.pmort, aes( x = LONG_FIADB, y = LAT_FIADB, color = Mort.quantiles), size = 0.5)+theme_bw()+facet_wrap(~COMMON)+
  scale_color_manual(values = c(
    "[0,0.01]" = "#d0d1e6", 
    "(0.01,0.05]" = "yellow",
    "(0.05,0.1]" = "yellow2",
    "(0.1,0.2]" = "#fee090", 
    "(0.22,0.3]" = "#fec44f", 
    "(0.3,0.4]" = "#fc8d59", 
    "(0.4,0.5]" =  "#f46d43", 
    "(0.5,1]" = "#d73027" 
    
    
  ))+
  coord_sf(xlim = c(-85, -68), ylim = c(35, 48))+theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'lightblue'))

ggsave(height = 10, width = 18, units = "in", "SPCD_stanoutput_joint/images/species_pmort_hat_distribution_maps.png")

# plot distribution 

# plot distribution 
ll.train.pmort  <- ll.train.pmort  %>% mutate(Mort.quantiles = cut(`p(mort)`, 
                                                                   breaks = c(0, 0.1, 0.2, 0.3, 0.40, 0.50, 1), 
                                                                   include.lowest=TRUE))

ggplot() +
  geom_polygon(data = canada, 
               aes(x=long, y=lat, group = group), 
               color = "black", fill = "white") +
  geom_polygon(data = state_sub, 
               aes(x=long, y=lat, group = group), 
               color = "black", fill = "white") +
  
  geom_point(data = ll.train.pmort, aes( x = LONG_FIADB, y = LAT_FIADB, color = Mort.quantiles), size = 0.5)+theme_bw()+facet_wrap(~COMMON)+
  scale_color_manual(values = c(
    "[0,0.1]" = "#d0d1e6", 
    "(0.1,0.2]" = "#fee090", 
    "(0.22,0.3]" = "#fec44f", 
    "(0.3,0.4]" = "#fc8d59", 
    "(0.4,0.5]" =  "#f46d43", 
    "(0.5,1]" = "#d73027" 
    
    
  ))+
  coord_sf(xlim = c(-85, -68), ylim = c(35, 48))+theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'lightblue'))




ggsave(height = 10, width = 18, units = "in", "SPCD_stanoutput/images/species_pmort_insample_distribution_maps.png")

# plot distribution 

ggplot() +
  geom_polygon(data = canada, 
               aes(x=long, y=lat, group = group), 
               color = "black", fill = "white") +
  geom_polygon(data = state_sub, 
               aes(x=long, y=lat, group = group), 
               color = "black", fill = "white") +
  
  geom_point(data = ll.train.pmort %>% filter(spp %in% spp.table[i,]$SPCD.id), aes( x = LONG_FIADB, y = LAT_FIADB, color = `p(mort)`), size = 0.5)+theme_bw()+facet_wrap(~COMMON)+
  scale_color_viridis()+
  coord_sf(xlim = c(-85, -68), ylim = c(35, 48))+theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'lightblue'))

ggsave(height = 5, width = 9, units = "in", paste0("SPCD_stanoutput/images/species_",spp.table[i,]$SPCD.id,"_pmort_rep_distribution_maps.png"))




ggplot(data = ll.test.pmort, aes( x = LONG_FIADB, y = LAT_FIADB, color = median))+
  geom_point()+facet_wrap(~COMMON)

ggplot(data = ll.train.pmort, aes( x = LONG_FIADB, y = LAT_FIADB, color = median))+
  geom_point()+facet_wrap(~COMMON)


##########################################################
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

save(log_lik_1, r_eff, loo_1, file = paste0("SPCD_stanoutput_joint/LOO_model_",model.no,"remper.corr_", "_species_", SPCD.id ,"_remper_corr_",remper.correction, ".Rdata"))

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

