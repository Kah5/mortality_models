library(rstan)
library(MASS)
#library(here)
library(tidyverse)
#library(gt)
library(FIESTA)
library(dplyr)
library(posterior)
library(mltools)
output.folder <- "C:/Users/KellyHeilman/Box/01. kelly.heilman Workspace/mortality/Eastern-Mortality/mortality_models/"

# Use posterior estimates to make predictions for annual probability of mortality:---

# get the complete species list
nspp <- data.frame(SPCD = c(316, 318, 833, 832, 261, 531, 802, 129, 762,  12, 541,  97, 621, 400, 371, 241, 375))
nspp$Species <- paste(FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD),]$GENUS, FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD),]$SPECIES)
nspp$COMMON <- paste(FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD),]$GENUS, FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD),]$SPECIES)


# read in the test data for all the species
spp.table <- data.frame(SPCD.id = nspp[1:17,]$SPCD, 
                        spp = 1:17, 
                        COMMON = nspp[1:17,]$COMMON)
model.no <- 6

set.seed(22)

# load the posteriors:---

# read in alphas by species
alpha_species <- readRDS(paste0(output.folder,"/SPCD_stanoutput_joint_v3/samples/alpha.spp_model_6_5000samples.rds"))

# get random sample of each parameter:
#rand.draws <- sample(1:max(alpha_species$.iteration), 1000)
#alpha_species <- alpha_species #%>% posterior::subset_draws(iteration = rand.draws)


# read in betas by species
beta_species <- readRDS(paste0(output.folder,"/SPCD_stanoutput_joint_v3/samples/u_betas_model_6_5000samples.rds"))#%>% 
  #posterior::subset_draws(iteration = rand.draws)


# population betas, alphas:
alpha_population <- readRDS(paste0(output.folder,"/SPCD_stanoutput_joint_v3/samples/alpha.p_model_6_5000samples.rds"))#%>%
  #posterior::subset_draws(iteration = rand.draws)
beta_population <- readRDS(paste0(output.folder,"/SPCD_stanoutput_joint_v3/samples/beta_model_6_5000samples.rds"))#%>%
  #posterior::subset_draws(iteration = rand.draws)



# load the data that we want to predict survival for:---
mod.data.full <-
  readRDS (
    paste0(
      output.folder,
      "SPCD_stanoutput_joint_v3/all_SPCD_model_",
      model.no,
      ".RDS"
    )
  )

plot.data.train <- readRDS (
  
  paste0(
    output.folder,
    "/SPCD_stanoutput_joint_v3/train_plot_data_SPCD_model_",
    model.no,
    ".RDS"
  )
)
# --- COMPUTE TREE-LEVEL POSTERIOR SURVIVAL PROBABILITIES----

 # function to generate posterior predictions for a species
 generate_survival_preds <- function(slope, # posterior samples of betas
                                     Xmat_tree, #a vector of the xM covariate values for that tree
                                     alpha_re, # species or population alphas
                                     remper,  # tree remper
                                     observed){ # observed tree status
   
   # calculate the logit
   logit.p.annual = alpha_re + slope %*% Xmat_tree 
   # calcuate the annual survival probability 
   pSannual = 1/(1 + exp(-logit.p.annual))
   
   # calculate the remper survival probability
   mMhat = as.data.frame(pSannual^remper)
   colnames(mMhat) <- "mMhat"
   
   # sample from binomial distribution to get predicted status over the remper
   y_hat = apply(mMhat, MARGIN=1, FUN = function(x){rbinom(1, 1, x)})
   
   
   cbind(y_hat, mMhat, pSannual, logit.p.annual)
 }
 
 
# Generate predictions of tree survival probability for all 17 of the species:--
# hat = in sample 
# rep = out of sample
# calculate AUC estimates with uncertainty

for(sp_i in 17:1){ 
      
    cat(paste("Generating in-sample predictions for tree survival for ", spp.table[match(sp_i, spp.table$spp),]$COMMON, "\n"))
      
    # get species indices
    species.index <- mod.data.full$SPP %in% sp_i # in sample
    species.index.rep <- mod.data.full$SPPrep %in% sp_i # out of sample
    
    # get species betas and alphas   
    slope_spp <- beta_species[, paste0("u_beta[", sp_i, ",", 1:33, "]")]
    alpha_spp <- alpha_species[paste0("alpha_SPP[", sp_i, "]")]
    
    
    # predictor variables for the species:
    xM_spp <- mod.data.full$xM[species.index,]
    xM_spp.rep <- mod.data.full$xMrep[species.index.rep,]
    
    # tree level remper variables
    remper_spp <- mod.data.full$Remper[species.index]
    remper_spp.rep <- mod.data.full$Remperoos[species.index.rep]
    
    
    
    # for in sample predictions:-----
    # create arrays for all the trees of this species  
    yhat_spp <- mMhat_spp <- pSannual_spp <- AUC_spp <- array(data = NA, dim = c(nrow(alpha_spp), nrow(xM_spp)))
    dimnames(yhat_spp) <- list(c(1:nrow(alpha_spp)),c(paste0("yhat[", 1:nrow(xM_spp),"]")))
    
    
    for(j in 1:nrow(xM_spp)){
        tree_posteriors <- generate_survival_preds(slope = as.matrix(slope_spp), # posterior samples of betas
                                Xmat_tree = as.vector(xM_spp[j,]), #a vector of the xM covariate values for that tree
                                alpha_re = alpha_spp, # species or population alphas
                                remper = remper_spp[j])
        colnames(tree_posteriors) <- c("y_hat", "mMhat", "pSannual", "logit.p.annual")
          
        yhat_spp[,j] <- tree_posteriors$y_hat
        mMhat_spp[,j] <- tree_posteriors$mMhat
        pSannual_spp[,j] <- tree_posteriors$pSannual
      
    }
    
    # save the posterior estimates for yhat and mMhat of this species:
    saveRDS(yhat_spp, paste0(output.folder, "SPCD_stanoutput_joint_v3/samples/Yhat_spp_",sp_i, ".RDS" ))
    saveRDS(mMhat_spp, paste0(output.folder, "SPCD_stanoutput_joint_v3/samples/mMhat_spp_",sp_i, ".RDS" ))
    saveRDS(pSannual_spp, paste0(output.folder, "SPCD_stanoutput_joint_v3/samples/pSannual_spp_",sp_i, ".RDS" ))
    
    # save the tree indices for each species
    species.tree.indices <- data.frame(SPP = mod.data.full$SPP[species.index], 
                                       y = mod.data.full$y[species.index],
                                       Remper = mod.data.full$Remper[species.index],
                                       tree.id = 1:length(mod.data.full$y[species.index]))
    saveRDS(species.tree.indices, paste0(output.folder, "SPCD_stanoutput_joint_v3/samples/species_indices_",sp_i, ".RDS" ))
    
    ############################################################
    # AUC using mltools auc_roc function
    # for in sample data
    actuals = mod.data.full$y[species.index]
    
    
    # get the range of responses for each sample:
    auc.is.list <- list()
    #p.surv.hat <- mMhat_spp #%>% select(paste0("mMhat[", 1:length(actuals), "]"))
    
    # calculate the AUC for each sample:
    AUC.is.species <- do.call(rbind, 
                        lapply(1:nrow(mMhat_spp), FUN =function(x){
      overall.auc <- auc_roc( as.vector(as.numeric(mMhat_spp[x, ])), mod.data.full$y[species.index])
      overall.auc
             })) %>% as.data.frame() %>% mutate(spp = sp_i) %>% mutate( 
                            SPCD = spp.table[match(sp_i, spp.table$spp),]$SPCD.id, 
                            AUC_type = "in-sample")
    colnames(AUC.is.species)[1] <- c("AUC")
    
    saveRDS(AUC.is.species, paste0(output.folder, "SPCD_stanoutput_joint_v3/samples/AUC_is_spp_",sp_i, ".RDS" ))
    
    
    
    
    # for the out of sample predictions:---
    cat(paste("Generating out-of-sample predictions for tree survival for ", spp.table[match(sp_i, spp.table$spp),]$COMMON, "\n"))
    
    # create arrays for all the trees of this species  
    yrep_spp <- mMrep_spp <- pSannual_spp_rep <- array(data = NA, dim = c(nrow(alpha_spp), nrow(xM_spp.rep)))
    dimnames(yrep_spp) <- list(c(1:nrow(alpha_spp)),c(paste0("yrep[", 1:nrow(xM_spp.rep),"]")))
    
    
    for(j in 1:nrow(xM_spp.rep)){
      tree_posteriors <- generate_survival_preds(slope = as.matrix(slope_spp), # posterior samples of betas
                                                 Xmat_tree = as.vector(xM_spp.rep[j,]), #a vector of the xM covariate values for trep tree
                                                 alpha_re = alpha_spp, # species or population alphas
                                                 remper = remper_spp.rep[j])
      colnames(tree_posteriors) <- c("y_rep", "mMrep", "pSannual_rep", "logit.p.annual")
      
      yrep_spp[,j] <- tree_posteriors$y_rep
      mMrep_spp[,j] <- tree_posteriors$mMrep
      pSannual_spp_rep[,j] <- tree_posteriors$pSannual_rep
      
    }
    
    # save the posterior estimates for yrep and mMrep of this species:
    saveRDS(yrep_spp, paste0(output.folder, "SPCD_stanoutput_joint_v3/samples/Yrep_spp_",sp_i, ".RDS" ))
    saveRDS(mMrep_spp, paste0(output.folder, "SPCD_stanoutput_joint_v3/samples/mMrep_spp_",sp_i, ".RDS" ))
    saveRDS(pSannual_spp_rep, paste0(output.folder, "SPCD_stanoutput_joint_v3/samples/pSannual_rep_spp_",sp_i, ".RDS" ))
    
    # save the tree indices for each species
    species.tree.indices <- data.frame(SPP = mod.data.full$SPPrep[species.index.rep], 
                                       y = mod.data.full$ytest[species.index.rep],
                                       Remper = mod.data.full$Remperoos[species.index.rep],
                                       tree.id = 1:length(mod.data.full$ytest[species.index.rep]))
    saveRDS(species.tree.indices, paste0(output.folder, "SPCD_stanoutput_joint_v3/samples/species_indices_",sp_i, "_rep.RDS" ))
    
    ############################################################
    # AUC using mltools auc_roc function
    # for in sample data
    actuals = mod.data.full$ytest[species.index.rep]
    
    
    # get the range of responses for each sample:
    auc.oos.list <- list()
    
    
    # calculate the AUC for each sample:
    AUC.oos.species <- do.call(rbind, 
                              lapply(1:nrow(mMrep_spp), FUN =function(x){
                                overall.auc <- auc_roc( as.vector(as.numeric(mMrep_spp[x, ])), mod.data.full$ytest[species.index.rep])
                                overall.auc
                              })) %>% as.data.frame() %>% mutate(spp = sp_i) %>% mutate( 
                                SPCD = spp.table[match(sp_i, spp.table$spp),]$SPCD.id, 
                                AUC_type = "out-of-sample")
    colnames(AUC.oos.species)[1] <- c("AUC")
    
    saveRDS(AUC.oos.species, paste0(output.folder, "SPCD_stanoutput_joint_v3/samples/AUC_oos_spp_",sp_i, ".RDS" ))

}

###############################################################
# Combine all species AUC scores for the posterior predictions
###############################################################


###############################################################
# Plot remper mortality estimated probabilities vs regional mortality rates
###############################################################

########################################################################################
# Generating posterior predictions from population estimates for the rest of the species
########################################################################################
