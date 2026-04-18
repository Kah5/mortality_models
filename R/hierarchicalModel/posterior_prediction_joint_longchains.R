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
nspp$COMMON <- FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD),]$COMMON_NAME


# read in the test data for all the species
spp.table <- data.frame(SPCD.id = nspp[1:17,]$SPCD, 
                        spp = 1:17, 
                        COMMON = nspp[1:17,]$COMMON)

# set the species order using the factors:
SP.TRAITS <- read.csv("data/NinemetsSpeciesTraits.csv") %>% filter(COMMON_NAME %in% unique(nspp[1:17,]$COMMON))
# order the trait db by softwood-hardwood, then shade tolerance, then name (this puts all the oaks together b/c hickory and red oak have the same tolerance values)
SP.TRAITS <- SP.TRAITS %>% group_by(SFTWD_HRDWD) %>% arrange(desc(SFTWD_HRDWD), desc(ShadeTol), desc(COMMON_NAME))

SP.TRAITS$Color <- c(# softwoods
  "#b2df8a",
  "#003c30", 
  "#b2182b", 
  "#fee090", 
  "#33a02c",
  
  
  # sugar  maples
  "#a6cee3",
  "#1f78b4",
  
  # red maple
  "#e31a1c",
  # yellow birch
  "#fdbf6f",
  # oaks
  "#cab2d6",
  "#8073ac",
  "#6a3d9a",
  
  # hickory
  "#7f3b08",
  # white ash
  "#bababa",
  # black cherry
  "#4d4d4d",
  # yellow poplar
  "#ff7f00",
  "#fccde5" # paper birch
  
  
)

SP.TRAITS$`Shade Tolerance`  <- ifelse(SP.TRAITS$ShadeTol >=4, "High", 
                                       ifelse(SP.TRAITS$ShadeTol <=2.5, "Low", "Moderate"))

# set up custom colors for species
sppColors <- SP.TRAITS$Color 
names(sppColors) <- unique(SP.TRAITS$COMMON_NAME) 

species_fill <- scale_fill_manual(values = sppColors)
species_color <- scale_color_manual(values = sppColors)

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
    "/SPCD_stanoutput_joint_v3/train_data_all_SPCD_model_",
    model.no,
    ".RDS"
  )
)

plot.data.test <- readRDS (
  
  paste0(
    output.folder,
    "/SPCD_stanoutput_joint_v3/test_data_all_SPCD_model_",
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
# read in the in-sample AUC scores:
auc.is.files <- list.files(path = paste0(output.folder, "SPCD_stanoutput_joint_v3/samples/"), pattern = 
"AUC_is_spp_", full.names = TRUE)


AUC.is.df <-
  do.call(rbind, lapply(auc.is.files, readRDS)) %>% group_by(spp, AUC_type) %>%
  summarise(
    median = median(AUC),
    auc.ci.lo = quantile(AUC, 0.025),
    auc.ci.hi = quantile(AUC, 0.975)
  ) %>% left_join(., spp.table) 

AUC.is.df$COMMON <- FIESTA::ref_species[match(AUC.is.df$SPCD.id, ref_species$SPCD),]$COMMON_NAME

# save as a combined AUC in-sample
saveRDS(
  AUC.is.df,
  paste0(
    output.folder,
    "SPCD_stanoutput_joint_v3/AUC_is_with_uncertainty.rds"
  )
)

# read in the out of sample species AUC scores:
auc.oos.files <- list.files(path = paste0(output.folder, "SPCD_stanoutput_joint_v3/samples/"), pattern = 
                             "AUC_oos_spp_", full.names = TRUE)


AUC.oos.df <-
  do.call(rbind, lapply(auc.oos.files, readRDS)) %>% group_by(spp, AUC_type) %>%
  summarise(
    median = median(AUC),
    auc.ci.lo = quantile(AUC, 0.025),
    auc.ci.hi = quantile(AUC, 0.975)
  ) %>% left_join(., spp.table) 

AUC.oos.df$COMMON <- FIESTA::ref_species[match(AUC.oos.df$SPCD.id, ref_species$SPCD),]$COMMON_NAME

# save as a combined AUC out of sample
saveRDS(
  AUC.oos.df,
  paste0(
    output.folder,
    "SPCD_stanoutput_joint_v3/AUC_oos_with_uncertainty.rds"
  )
)



###############################################################
# Plot remper mortality estimated probabilities vs regional mortality rates
###############################################################
# for each species, read in the in-sample and out-of-sample annual probability of mortality:
# read in tree_remeas to get the volfac for each represented tree:
TREE.remeas <- readRDS( "data/unfiltered_TREE.remeas.rds")

all.remeas <- TREE.remeas %>%
  # do the filtering section
  filter( exprem > 0 & # if exprem == 0, these could be modeled plots?
            dbhold >= 5 & # need an initial dbh greater than 5
            ! remper == 0 & # if remper is listed as zero, filter out
            DIA_DIFF >= 0 & # filter diameter differences >= 0
            # !status == 3 & # keep the cut trees for this
            SPCD %in% nspp[1:17,]$SPCD & # filter species in the top 17 of all species
            !is.na(status) & # filter out trees with no status
            !is.na(elev)) 

i <- 14

# just take the last 500 draws from each of 3 chains
draws2keep <- c(2001:2500, 4501:5000, 7001:7500)

for(i in 13:1){
  
cat("\n",paste0("estimating weighted county and region means for ", nspp[i,]$COMMON))
  
  # read in the in-sample and out of sample
  pSannual_s <- readRDS(paste0(output.folder, "SPCD_stanoutput_joint_v3/samples/pSannual_spp_",i,".RDS"))[draws2keep,]%>%
     reshape2::melt(.) %>% 
    rename("tree.id"  = "Var2", 
           "sample" = "Var1", 
           "p1year" = "value")%>%
    group_by(tree.id)%>%
    mutate(spp = i) %>%
    mutate(data.type = "in-sample")
  
  # read the out of sample data
  pSannual_rep <- readRDS(paste0(output.folder, "SPCD_stanoutput_joint_v3/samples/pSannual_rep_spp_",i,".RDS"))[draws2keep,]%>%
    reshape2::melt(.) %>% 
    rename("tree.id"  = "Var2", 
           "sample" = "Var1", 
           "p1year" = "value")%>%
    group_by(tree.id)%>%
    mutate(spp = i) %>%
    mutate(data.type = "in-sample")
  
  # get tree level volfac infromation
  spp.state.id.is <- plot.data.train %>% filter(SPCD %in% spp.table[i,]$SPCD.id)%>%
    mutate(tree.id = 1:length(county))%>%
    dplyr::select(PLOT.ID, SPCD, tree.id, state, county, pltnum, cndtn, point, tree, cycle, dbhcur, dbhold, remper)%>%
    mutate(data.type = "in-sample") %>% 
    left_join(., all.remeas %>% dplyr::select(state, county, pltnum, point, tree, status, SPCD, dbhcur, dbhold, volfac) %>% distinct())
  
  spp.state.id.oos <- plot.data.test %>% filter(SPCD %in% spp.table[i,]$SPCD.id)%>%
    mutate(tree.id = 1:length(county))%>%
    dplyr::select(PLOT.ID, SPCD, tree.id, state, county, pltnum, cndtn, point, tree, cycle, dbhcur, dbhold, remper)%>%
    mutate(data.type = "out-of-sample")%>% 
    left_join(., all.remeas %>% dplyr::select(state, county, pltnum, point, tree, status, SPCD, dbhcur, dbhold, volfac) %>% distinct())
  
spp.state.ids.all  <- rbind(spp.state.id.is, spp.state.id.oos)
  
  # calculate county specific weights for all the trees
weights.volfac <- spp.state.ids.all %>% group_by(state, county) %>%
  mutate(total.volfac.spp.county = sum(volfac, na.rm =TRUE), 
         ntrees = n())%>%
  ungroup()%>%
  group_by(PLOT.ID, tree.id, SPCD, state, county, pltnum, point, tree, data.type)%>%
  mutate(weight_co_volfac = volfac/total.volfac.spp.county)%>%
  dplyr::select(tree, PLOT.ID, tree.id, state, county, data.type, volfac, status, weight_co_volfac, total.volfac.spp.county,ntrees)
  

pannual_all <- rbind(pSannual_rep, pSannual_s) %>%
  mutate(p_mort_1year = 1-p1year, 
         p_surv_10year = p1year^10, 
         p_mort_10year = 1-(p1year^10))%>% left_join(., weights.volfac) %>%
  
  # get weights based on volfac and counties
  mutate(p_mort_1year_volfac = p_mort_1year*volfac, 
         p_surv_10year_volfac = p_surv_10year*volfac, 
         p_mort_10year_volfac = p_mort_10year*volfac)

pMort_county_weighted_samples <- pannual_all %>% 
  ungroup() %>% 
  group_by(SPCD, sample, state, county) %>%
  
  # get the weighted average for each state-county-mcmc sample
  summarise(pmort_co_1year = sum(p_mort_1year_volfac, na.rm =TRUE)/sum(volfac, na.rm =TRUE),
            pmort_co_10year = sum(p_mort_10year_volfac, na.rm =TRUE)/sum(volfac, na.rm =TRUE)) #

pMort_county_weighted <- pMort_county_weighted_samples %>% 
  ungroup()%>%
  group_by(SPCD, state, county)%>%
  summarise(pmort_weighted_1 = median(pmort_co_1year, na.rm =TRUE), 
            pmort_weighted_1.ci.lo = quantile(pmort_co_1year, 0.025, na.rm =TRUE), 
            pmort_weighted_1.ci.hi = quantile(pmort_co_1year, 0.975, na.rm =TRUE), 
            pmort_weighted_10 = median(pmort_co_10year, na.rm =TRUE), 
            pmort_weighted_10.ci.lo = quantile(pmort_co_10year, 0.025, na.rm =TRUE), 
            pmort_weighted_10.ci.hi = quantile(pmort_co_10year, 0.975, na.rm =TRUE))



pMort_region_weighted_samples <- pannual_all %>% 
  ungroup() %>% 
  group_by(SPCD, sample) %>%
  
  # get the weighted average for each state-county-mcmc sample
  summarise(pmort_co_1year = sum(p_mort_1year_volfac, na.rm =TRUE)/sum(volfac, na.rm =TRUE),
            pmort_co_10year = sum(p_mort_10year_volfac, na.rm =TRUE)/sum(volfac, na.rm =TRUE)) #

# save the samples weighted by volface of each species in the county
saveRDS(
  pMort_county_weighted_samples,
  paste0(
    output.folder,
    "SPCD_stanoutput_joint_v3/predicted_mort/Mort_weighted_county_mortality_samps_",nspp[i,]$SPCD,".rds"
  )
)

saveRDS(
  pMort_county_weighted,
  paste0(
    output.folder,
    "SPCD_stanoutput_joint_v3/predicted_mort/Mort_weighted_county_mortality_averages_",nspp[i,]$SPCD,".rds"
  )
)

saveRDS(
  pMort_region_weighted_samples,
  paste0(
    output.folder,
    "SPCD_stanoutput_joint_v3/predicted_mort/Mort_weighted_region_mortality_samps_",nspp[i,]$SPCD,".rds"
  )
)

rm(pMort_county_weighted, pMort_county_weighted_samples, pMort_region_weighted_samples, pSannual_s, pSannual_rep, pannual_all)

}

# do the summaries another way

for(i in 17:1){
    
    # read in the in-sample and out of sample
    pSannual_s <- readRDS(paste0(output.folder, "SPCD_stanoutput_joint_v3/samples/pSannual_spp_",i,".RDS"))#%>%
    
   
    
    # get the tree annual survival probabilities:
    pSannualyear_tree <- pSannual_s %>% reshape2::melt(.) %>% 
      rename("tree.id"  = "Var2", 
             "sample" = "Var1", 
             "p1year" = "value")%>%
      group_by(tree.id)%>%
      summarise(p1year.med = median(p1year), 
                p1year.ci.lo = quantile(p1year, 0.025), 
                p1year.ci.hi = quantile(p1year, 0.975), 
                p1year.sd = sd(p1year))%>%
      mutate(spp = i) 
    
    # convert to 10 year survival probabilities:
    pS_10year <- pSannual_s^10
    
    # get the tree and species-level summaries
    pS10year_tree <- pS_10year %>% reshape2::melt(.) %>% 
      rename("tree.id"  = "Var2", 
             "sample" = "Var1", 
             "p10year" = "value")%>%
      group_by(tree.id)%>%
      summarise(p10year.med = median(p10year), 
                p10year.ci.lo = quantile(p10year, 0.025), 
                p10year.ci.hi = quantile(p10year, 0.975), 
                p10year.sd = sd(p10year))%>%
      mutate(spp = i) %>%
     # left_join(., pSannualyear_tree)%>%
      left_join(., spp.table) %>%
      group_by(tree.id)%>%
      # use pS_10median to get a survival for each tree over a ten year interval:
      mutate(survival.draw = rbinom(1,1, prob = p10year.med))%>%
      mutate(data.type = "in-sample") %>%
      rename("SPCD" = "SPCD.id")
    
    spp.state.id.is <- plot.data.train %>% filter(SPCD %in% unique(pS10year_tree$SPCD))%>%
      mutate(tree.id = 1:length(county))%>%
      select(PLOT.ID, SPCD, tree.id, state, county, pltnum, cndtn, point, tree, cycle, dbhcur, dbhold, remper)%>%
      mutate(data.type = "in-sample")
    
    pS10year_tree <- left_join(pS10year_tree, spp.state.id.is)
    
    # read the out of sample data
    pSannual_rep <- readRDS(paste0(output.folder, "SPCD_stanoutput_joint_v3/samples/pSannual_rep_spp_",i,".RDS"))#%>%
    
    # convert to 10 year survival probabilities:
    pS_10year_rep <- pSannual_rep^10
    
    # get the tree and species-level summaries
    pS10year_tree_rep <- pS_10year_rep %>% reshape2::melt(.) %>% 
      rename("tree.id"  = "Var2", 
             "sample" = "Var1", 
             "p10year" = "value")%>%
      group_by(tree.id)%>%
      summarise(p10year.med = median(p10year), 
                p10year.ci.lo = quantile(p10year, 0.025), 
                p10year.ci.hi = quantile(p10year, 0.975), 
                p10year.sd = sd(p10year))%>%
      mutate(spp = i) %>%
      left_join(., spp.table) %>%
      group_by(tree.id)%>%
      # use pS_10median to get a survival for each tree over a ten year interval:
      mutate(survival.draw = rbinom(1,1, prob = p10year.med)) %>%
      mutate(data.type = "out-of-sample")%>%
      rename("SPCD" = "SPCD.id")
    
    spp.state.id.oos <- plot.data.test %>% filter(SPCD %in% unique(pS10year_tree$SPCD))%>%
      mutate(tree.id = 1:length(county))%>%
      select(PLOT.ID, SPCD, tree.id, state, county, pltnum, cndtn, point, tree, cycle, dbhcur, dbhold, remper)%>%
      mutate(data.type = "out-of-sample")
    
    pS10year_tree_rep <- left_join(pS10year_tree_rep, spp.state.id.oos)
    
    # combine all the trees
    pS10year_all <- rbind(pS10year_tree, pS10year_tree_rep) %>% ungroup() %>%
      mutate(mortality.draw = 1-survival.draw)
    
    # calculate a whole region summary: note that these are not scaled by volfac
    regional.spp.summary <-  pS10year_all %>% 
      group_by(spp, SPCD, COMMON)%>% 
      summarise(ntree = n(), 
                nmort = sum(mortality.draw),
                pct.mort.10year = ((sum(mortality.draw)/n())*100))%>%
      mutate(pct.mort.annual = pct.mort.10year/10)
    
    
    
    # calculate a state-level summaries: note that these are not scaled by volfac
    state.spp.summary <-  pS10year_all %>% ungroup() %>%
      
      mutate(mortality.draw = 1-survival.draw) %>% 
      group_by(spp, state, SPCD, COMMON)%>% 
      summarise(ntree = n(), 
                nmort = sum(mortality.draw),
                pct.mort.10year = ((sum(mortality.draw)/n())*100))%>%
      mutate(pct.mort.annual = pct.mort.10year/10)
    
    # save these combined summaries:
    saveRDS(
      state.spp.summary,
      paste0(
        output.folder,
        "SPCD_stanoutput_joint_v3/predicted_mort/Mort_10yr/state_10yr_mortality_",nspp[i,]$SPCD,".rds"
      )
    )
    
    saveRDS(
      regional.spp.summary,
      paste0(
        output.folder,
        "SPCD_stanoutput_joint_v3/predicted_mort/Mort_10yr/regional_10yr_mortality_",nspp[i,]$SPCD,".rds"
      )
    )
    
    saveRDS(
      pS10year_all,
      paste0(
        output.folder,
        "SPCD_stanoutput_joint_v3/predicted_mort/Mort_10yr/tree_10yr_mortality_",nspp[i,]$SPCD,".rds"
      )
    )
}

# read in predicted mortality summaries: note that these are not scaled by volfac
all.state.species.10year <- do.call(rbind, lapply(list.files(path = paste0(
  output.folder,
  "SPCD_stanoutput_joint_v3/predicted_mort/Mort_10yr/"), 
  pattern = "state_10yr_mortality_", 
  full.names = TRUE), 
  readRDS))

all.region.species.10year <- do.call(rbind, lapply(list.files(path = paste0(
  output.folder,
  "SPCD_stanoutput_joint_v3/predicted_mort/Mort_10yr/"), 
  pattern = "regional_10yr_mortality_", 
  full.names = TRUE), 
  readRDS))


# refactor the species so they are plotted in our order

ggplot()+
  geom_bar(data = all.region.species.10year, aes(x = COMMON, y = pct.mort.annual, fill = COMMON), stat = "identity")+
  theme_bw(base_size = 16)+
  species_fill+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        #panel.grid = element_blank(), 
        legend.position = "none", 
        axis.title.x = element_blank())+
  ylab("Predicted annual mortality rate")
# %>% filter(n_plots_SPCD > 50)
ggplot(data = all.state.species.10year )+
  geom_bar(aes(x = COMMON, y = pct.mort.annual, fill = COMMON), stat = "summary", fun = median,  color = "black")+
  theme_bw(base_size = 16)+
  species_fill+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        axis.title.x = element_blank(), 
        legend.position = "none", 
        panel.grid = element_blank())+
  geom_point(aes(x = COMMON, y = pct.mort.annual, group = state), alpha = 0.9) +
  #geom_text(aes(x = Species, y = species_volfac_mort, label = region),position = position_jitter(width = 0.1), alpha = 0.6)
  ggrepel::geom_text_repel(aes(x = COMMON, y = pct.mort.annual, label = state), color = "black", size = 2.5, segment.color = "grey", max.overlaps = 12) +  
  ylab("Mortality Rate (% per year)")




########################################################################################
# Generating posterior predictions from population estimates for the rest of the species
########################################################################################
