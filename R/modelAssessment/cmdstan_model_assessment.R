library(tidyverse)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(qs2)
library(jsonlite)
library(pROC)
library(FIESTA)
color_scheme_set("brightblue")
# script to read in cmdstan model outputs, generate model assessement and comparisons

output.dir = "C:/Users/KellyHeilman/Box/01. kelly.heilman Workspace/mortality/Eastern-Mortality/mortality_models/"

spp.table <- read.csv(file = paste0(output.dir, "/data/Hierarchical_obs_model_7.csv"))
spp.table
nspp <- spp.table
# nspp <- data.frame(SPCD = c(316, 318, 833, 832, 261, 531, 802, 129, 762,  12, 541,  97, 621, 400, 371, 241, 375))
# nspp$Species <- paste(FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD),]$GENUS, FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD),]$SPECIES)
# nspp$COMMON_NAME <- FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD),]$COMMON_NAME
# 


options(mc.cores = parallel::detectCores())
SPCD.df <- data.frame(spp.table[,c("SPCD", "SPP")])
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


unique(convergence_all$model.type)


ggplot(convergence_all)+geom_histogram(aes(x = rhat.median))+
  facet_wrap(~model.type)

ggplot(convergence_all)+geom_histogram(aes(x = ess_bulk.median))+
  facet_wrap(~model.type)

saveRDS(convergence_all, paste0(output.dir, "SPCD_stanoutput_cmdstan/summary/convergence_stats_all_models.RDS"))


# get the species model diagnostic files with core hours
diagnostic.files <- list.files(path = paste0(output.dir,"SPCD_stanoutput_cmdstan/diagnostics/"), 
                               pattern = paste0("sample_diagnostics_mort_model"), full.names = TRUE)

diags_all <- do.call(rbind, lapply(diagnostic.files, read.csv))
diag.summary <- diags_all %>%  #group_by(SPCD, model.type, model.number, chain_id)%>%
  mutate(chain_core_hours = total*(1/60)*(1/60)*(ncores/nchain), 
         chain_core_hours_sampling = sampling*(1/60)*(1/60)*(ncores/nchain), 
         chain_core_hours_warmup = warmup*(1/60)*(1/60)*(ncores/nchain)) %>%
  group_by(SPCD, model.type, model.number)%>%
  summarise(total_core_hours = sum(chain_core_hours), 
            sampling_core_hours = sum(chain_core_hours_sampling),
            warmup_core_hours = sum(chain_core_hours_warmup),
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


# check hierarchical models diagnostic files
hier.diagnostic.files <- list.files(path = paste0(output.dir,"SPCD_stanoutput_cmdstan/diagnostics/"), pattern = paste0("sample_diagnostics_hierarchical_mort_model"), full.names = TRUE)

diag.test <-read.csv(hier.diagnostic.files[1])

diags_heir <- do.call(rbind, lapply(hier.diagnostic.files, read.csv))%>%
  rename("model.number"= "model.no")

# get the number of individuals per species and weight by species--
hier.data <- fromJSON(paste0("SPCD_standata_json/hierarchical_data_model_",1,".json"))
# sample.trees <- data.frame(SPP = hier.data$SPP) %>% 
#   left_join(.,data.frame(SPP = 1:17, SPCD = nspp$SPCD)) %>%
#   group_by(SPCD) %>% 
#   summarise(ntrees = n()) %>% ungroup()%>%
#   mutate(total_trees = sum(ntrees)) %>% 
#   mutate(weight_species = ntrees/total_trees)

sample.trees <- spp.table %>%
  rename("ntrees" = "n")%>%
  group_by(SPCD, COMMON) %>% 
  ungroup()%>%
  mutate(total_trees = sum(ntrees)) %>% 
  mutate(weight_species = ntrees/total_trees)

diagnostic_heir <- expand.grid(model.number = unique(diags_heir$model.number), 
                               chain_id = unique(diags_heir$chain_id), 
                               SPCD = sample.trees$SPCD)%>%
  left_join(.,sample.trees) %>% 
  mutate(model.type = "Hierarchical")%>% 
  left_join(., diags_heir) %>%
  mutate(niter = 1000, 
         nwarmup = 1000) 



diag_heir.summary <- diagnostic_heir %>% 
  group_by(SPCD, model.type, model.number, chain_id)%>%
  mutate(chain_core_hours = total*weight_species*(1/60)*(1/60)*(ncores/nchain), 
         chain_core_hours_sampling = sampling*weight_species*(1/60)*(1/60)*(ncores/nchain), 
         chain_core_hours_warmup = warmup*weight_species*(1/60)*(1/60)*(ncores/nchain)) %>%
  group_by(SPCD, model.type, model.number)%>%
  summarise(total_core_hours = sum(chain_core_hours), 
            sampling_core_hours = sum(chain_core_hours_sampling),
            warmup_core_hours = sum(chain_core_hours_warmup),
            total_iter = sum(niter), 
            total_warmup = sum(nwarmup), 
            
            total_divergent = sum(n_divergent),
            total_max_treedepth = sum(tot.max.treedepth))



# check that there are no divergent transitions
diag_heir.summary %>% filter(total_divergent > 0)  
diags_heir %>% filter(n_divergent > 0)  

# no hierarchical models have samples that exceeded max tree depth
diag_heir.summary %>% filter( total_max_treedepth > 0)  %>% 
  select(SPCD, model.number, total_max_treedepth)

# combine and save all the species outputs.
diag_heir.converg.df <- left_join(diag_heir.summary, convergence_all)


hier.convergence.df <- convergence_all %>% filter(model.type %in% "Hierarchical")
hist(hier.convergence.df$rhat.median)
hist(hier.convergence.df$rhat.ci.hi)
hist(hier.convergence.df$rhat.ci.lo)
hist(hier.convergence.df$ess_bulk.median)
hist(hier.convergence.df$ess_bulk.ci.lo)

all_diagnostics <- rbind(diag_heir.converg.df, diag_converg.df)
all_diagnostics$COMMON <- ref_species[match(all_diagnostics$SPCD, ref_species$SPCD),]$COMMON_NAME

saveRDS(all_diagnostics, paste0(output.dir, "SPCD_stanoutput_cmdstan/summary/sampler_diagnostics_all_models.RDS"))

# total core hours
all_diagnostics |>
  ggplot()+geom_bar(aes(x = model.number, y = total_core_hours, fill = model.type), stat = "identity", position = "dodge")+
  facet_wrap(~SPCD, scales = "free_y")

# sampling core hours--same # of samples
all_diagnostics |>
  ggplot()+geom_bar(aes(x = as.character(model.number), y = sampling_core_hours, fill = model.type), stat = "identity", position = "dodge")+
  facet_wrap(~COMMON, scales = "free_y")+
  theme_bw(base_size = 12)+ylab("Core Hours for Sampling (4000 draws)")+xlab("Model")+
  theme(panel.background = element_blank())

ggsave(filename = paste0(output.dir, "SPCD_stanoutput_cmdstan/summary/Core_hours_all_models.png"), 
       height = 7, width = 10)

# standardized to 500 warmup samples 
all_diagnostics |>
  ggplot()+geom_bar(aes(x = model.number, y = (warmup_core_hours/total_warmup)*2000, fill = model.type), stat = "identity", position = "dodge")+
  facet_wrap(~SPCD, scales = "free_y")

# get loo results and compare---
# need to do for each species...!!
SPCD.id <- 261

# function to read in LOO output and summarise for each species:
# udpate hardcoded file numbers when hierarchical models are run!

LOO_summarise_SPCD <- function(SPCD.id){
  
    cat(paste0("Doing LOO comparison for ", SPCD.id))
  
   # all file names for species models and hierarchical models
   all.species.files <- list.files(path = paste0(output.dir,"SPCD_stanoutput_cmdstan/LOO/"), 
               pattern ="LOO_results_mort_model", full.names = T)
   
   all.Hierarchical.files <- list.files(path = paste0(output.dir,"SPCD_stanoutput_cmdstan/LOO/"), 
                                        pattern ="LOO_results_hierarchical_mort_model", full.names = T)
     
   # get only the species we are looking at
   spp.loo.files <- grep( paste0("_", SPCD.id, "(_|\\.)"),  
                          all.species.files, 
                          value = TRUE, 
                          perl = TRUE)
   Hierarchical.loo.files <-   grep( paste0("_", SPCD.id, "(_|\\.)"),  
                                      all.Hierarchical.files, 
                                      value = TRUE, 
                                      perl = TRUE)
   

   
   # get the model numbers for each species and hierarchical models
   # we were using length to name before, but this will ensure we don't accidently rename a model incorrectly
   spp.model_numbers <- as.numeric(str_extract(spp.loo.files, "(?<=LOO_results_mort_model_)\\d+"))
   hier.model_numbers <- as.numeric(str_extract( Hierarchical.loo.files, "(?<=LOO_results_hierarchical_mort_model_)\\d+"))
  
   all.mod.numbers <- 1:9
   
   # add some warnings to list the models that need to be run
   if(length(spp.loo.files) < 9){
     cat(paste0("SPCD ", SPCD.id, " is missing Species-level model ", all.mod.numbers[!1:9 %in% spp.model_numbers], "\n"))
   }  
   
   if(length(Hierarchical.loo.files) < 9){
     cat(paste0("SPCD ", SPCD.id, " is missing Hierarchical model ", all.mod.numbers[!1:9 %in% hier.model_numbers], "\n"))
   }  
   
    loo.files.all <- c(spp.loo.files, Hierarchical.loo.files)
    loo_results_all <- lapply(loo.files.all, qs_read)
    
  
    model.full.names <-  c(paste0("Species_model_",  spp.model_numbers), 
      paste0("Hierarchical_model_", hier.model_numbers))
    
    # add names to the list:
    names(loo_results_all) <- model.full.names
   
    # check pareto-k estimates:
    pareto.k.checks <- do.call(rbind, lapply(loo_results_all, function(x){data.frame(good = sum(x$diagnostics$pareto_k <= 0.7),
                                                                                     bad = sum(x$diagnostics$pareto_k > 0.7),
                                                                                     total = length(x$diagnostics$pareto_k))}))%>%
      mutate(percent.bad = (bad/total)*100)%>%
      mutate(model.name = c(paste0("model ", spp.model_numbers), paste0("model ", hier.model_numbers)),
             model =model.full.names,
             model.type = c(rep("Species", length(spp.loo.files)), rep("Hierarchical", length( Hierarchical.loo.files))),
             SPCD = SPCD.id)
  
   # mod.data.7 <- fromJSON(paste0("SPCD_standata_json/SPCD_", SPCD.id,"remper_correction_0.5model_7.json"))
   
    
    mod.data <- fromJSON(paste0("SPCD_standata_json/hierarchical_data_model_7.json"))
    
    
    cat(paste("generating posterior predictions for SPCD", SPCD.id, "species number", s, "\n"))
    
    spec.idx     <- mod.data$SPP    %in% s
    spec_rep.idx <- mod.data$SPPrep %in% s
    
    spp.mod.data <- mod.data
    
    # set up species-specific data for gen quants
    spp.mod.data$SPP    <- mod.data$SPP[spec.idx]
    spp.mod.data$xM     <- as.matrix(mod.data$xM[spec.idx, ])
    spp.mod.data$Remper <- mod.data$Remper[spec.idx]
    spp.mod.data$N      <- length(spp.mod.data$SPP)
    spp.mod.data$y      <- mod.data$y[spec.idx]
    
    spp.mod.data$SPPrep    <- mod.data$SPPrep[spec_rep.idx]
    spp.mod.data$xMrep     <- as.matrix(mod.data$xMrep[spec_rep.idx, ])
    spp.mod.data$Remperoos <- mod.data$Remperoos[spec_rep.idx]
    spp.mod.data$Nrep      <- length(spp.mod.data$SPPrep)
    spp.mod.data$ytest     <- mod.data$ytest[spec_rep.idx]
    
    #spp.mod.data$N
    
    loo_compare.out <- loo::loo_compare(loo_results_all) # best fit based on loo elpd differences
    loo_comparisons <- loo_compare.out %>% data.frame()%>%
      left_join(., pareto.k.checks) %>%
      mutate(elpd_se_ratio = abs(elpd_diff)/se_diff)%>%
      mutate(model.number = substr(model.name, start = 7, stop = 7))
  
  
    # Get model weights---
    # get pointwise log predictive densities
    lpd_point <- do.call(cbind,lapply(loo_results_all, function(x){x$pointwise[,"elpd_loo"]}))
    pbma_wts <- loo::pseudobma_weights(lpd_point, BB=FALSE)
    pbma_BB_wts <- loo::pseudobma_weights(lpd_point) # default is BB=TRUE
    stacking_wts <- loo::stacking_weights(lpd_point)
    mod.weights <- round(cbind(pbma_wts, pbma_BB_wts, stacking_wts),3)%>% data.frame() %>%
      mutate(model = model.full.names,
             SPCD = SPCD.id)
    loo_comparisons_weighted <- loo_comparisons %>% left_join(., mod.weights, by = c("model", "SPCD"))
    return(loo_comparisons_weighted)
}

LOO_ELPD_list <- list()
for(s in 1:length(nspp$SPCD)){
  LOO_ELPD_list[[s]] <- LOO_summarise_SPCD(SPCD.id = nspp$SPCD[s])
}

LOO_ELPD.df <- do.call(rbind, LOO_ELPD_list)



# make sure none of the models have >5-10 percent bad pareto k values
max(LOO_ELPD.df$percent.bad)


# Order Species common names by from most to least abundant
LOO_ELPD.df$COMMON_NAME <- FIESTA::ref_species[match(LOO_ELPD.df$SPCD, FIESTA::ref_species$SPCD),]$COMMON_NAME
LOO_ELPD.df$COMMON_NAME <- factor(LOO_ELPD.df$COMMON_NAME, levels = unique(nspp$COMMON))
LOO_ELPD.df$model.full <- paste(LOO_ELPD.df$model.name, LOO_ELPD.df$model.type)

# save LOO_ELPD.df 
saveRDS(LOO_ELPD.df, paste0(output.dir, "SPCD_stanoutput_cmdstan/summary/All_LOO_comparisons.rds"))
LOO_ELPD.df <- readRDS( paste0(output.dir, "SPCD_stanoutput_cmdstan/summary/All_LOO_comparisons.rds"))


# plot of all species ELPD diff +/- SE, with significance
LOO_ELPD.df %>%
  mutate(elpd_diff_sig = ifelse(elpd_se_ratio >= 2, "ELPD_diff >= 2 SE", 
                                ifelse(elpd_se_ratio <2, "ELPD_diff < 2 SE", "Best-fit")))%>%
  mutate(elpd_diff_sig = ifelse(is.na(elpd_diff_sig), "Best-fit ELPD", elpd_diff_sig))|>
    ggplot()+geom_pointrange(aes(x = model.full, y = elpd_diff, ymin = elpd_diff+se_diff, ymax = elpd_diff-se_diff, color = elpd_diff_sig))+
      facet_wrap(~COMMON_NAME, scales = "free_y")+
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust =0.5))

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
  geom_bar(aes(x = COMMON_NAME, y = pbma_BB_wts, fill = model.full), stat = "identity")+theme_bw()+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust =0.5))

LOO_ELPD.df |>  ggplot()+
  geom_bar(aes(x = COMMON_NAME, y = stacking_wts, fill = model.full), stat = "identity")+theme_bw()+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust =0.5))


# unstacked barplot of model weights by species
LOO_ELPD.df |>  ggplot()+
  geom_bar(aes(x = model.full, y = pbma_BB_wts, fill = model.type), stat = "identity", position = "dodge")+theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust =0.5))+
  facet_wrap(~COMMON_NAME)


LOO_ELPD.df |>  ggplot()+
  geom_tile(aes(x = model.full, y = COMMON_NAME, fill = pbma_BB_wts))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust =0.5))+
  scale_fill_distiller(palette = "YlGnBu", direction = 1, name = "BMA Weight")+
  ylab("Species")+xlab("Model")

ggsave(filename = paste0(output.dir, "SPCD_stanoutput_cmdstan/summary/Species_model_BMA_weights.png"))

# plot up ELP by model with uncertainty
LOO_ELPD.df$model.type <- factor(LOO_ELPD.df$model.type, levels = c( "Species", "Hierarchical"))
LOO_ELPD.df |>  ggplot()+
  geom_tile(aes(x = model.full, y = COMMON_NAME, fill = stacking_wts))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust =0.5))+
  scale_fill_distiller(palette = "YlGnBu", direction = 1, name = "Stacking Weight")+
  ylab("Species")+xlab("Model")+facet_wrap(~model.type, scales = "free_x")

ggsave(filename = paste0(output.dir, "SPCD_stanoutput_cmdstan/summary/Species_model_stacking_weights.png"))



LOO_ELPD.df |>  ggplot()+
  geom_pointrange(aes(x = model.name, y = elpd_loo, ymin = elpd_loo - se_elpd_loo, ymax = elpd_loo+ se_elpd_loo, color = model.type), position = position_dodge(width = 1))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust =0.5))+xlab("Model")+
  facet_wrap(~COMMON_NAME, scales = "free_y")


LOO_ELPD.df %>% mutate(elpd_diff_sig = ifelse(elpd_se_ratio >= 2, "ELPD_diff >= 2 SE", 
                                          ifelse(elpd_se_ratio <2, "ELPD_diff < 2 SE", "Best-fit")))%>%
  mutate(elpd_diff_sig = ifelse(is.na(elpd_diff_sig), "Best-fit ELPD", elpd_diff_sig)) |>  ggplot()+
  geom_pointrange(aes(x = model.name, y = elpd_diff, ymin = elpd_diff - se_diff, ymax = elpd_diff + se_diff, shape = model.type, color = elpd_diff_sig), position = position_dodge(width = 1))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust =0.5))+xlab("Model")+
  facet_wrap(~COMMON_NAME, scales = "free_y") +
  ylab("ELPD Difference (+/- SE difference)")+
  scale_color_manual(values = c("ELPD_diff >= 2 SE" = "lightgrey", 
                                "ELPD_diff < 2 SE" = "black", 
                                "Best-fit ELPD" = "red"))
ggsave(filename = paste0(output.dir, "SPCD_stanoutput_cmdstan/summary/ELPD_differences_all_models.png"), 
       height = 7, width = 10)


LOO_ELPD.df %>% mutate(elpd_diff_sig = ifelse(elpd_se_ratio >= 4, "ELPD_diff >= 4 SE", 
                                              ifelse(elpd_se_ratio <4, "ELPD_diff < 4 SE", "Best-fit")))%>%
  mutate(elpd_diff_sig = ifelse(is.na(elpd_diff_sig), "Best-fit ELPD", elpd_diff_sig)) |>  ggplot()+
  geom_pointrange(aes(x = model.name, y = elpd_diff, ymin = elpd_diff - se_diff, ymax = elpd_diff + se_diff, shape = model.type, color = elpd_diff_sig), position = position_dodge(width = 1))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust =0.5))+xlab("Model")+
  facet_wrap(~COMMON_NAME, scales = "free_y") +
  ylab("ELPD Difference (+/- SE difference)")+
  scale_color_manual(values = c("ELPD_diff >= 4 SE" = "lightgrey", 
                                "ELPD_diff < 4 SE" = "black", 
                                "Best-fit ELPD" = "red"))
ggsave(filename = paste0(output.dir, "SPCD_stanoutput_cmdstan/summary/ELPD_differences_all_models_4SE.png"), 
       height = 7, width = 10)

LOO_ELPD.df %>% mutate(elpd_diff_sig = ifelse(elpd_se_ratio >= 2, "ELPD_diff >= 2 SE", 
                                              ifelse(elpd_se_ratio <2, "ELPD_diff < 2 SE", "Best-fit")))%>%
  mutate(elpd_diff_sig = ifelse(is.na(elpd_diff_sig), "Best-fit ELPD", elpd_diff_sig)) |>  ggplot()+
  geom_point(aes(x = model.name, y = stacking_wts, shape = model.type, color = elpd_diff_sig), position = position_dodge(width = 1))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust =0.5))+xlab("Model")+
  facet_wrap(~COMMON_NAME, scales = "free_y") +
  ylab("Stacking weights")+
  scale_color_manual(values = c("ELPD_diff >= 2 SE" = "lightgrey", 
                                "ELPD_diff < 2 SE" = "black", 
                                "Best-fit ELPD" = "red"))
ggsave(filename = paste0(output.dir, "SPCD_stanoutput_cmdstan/summary/Stacking_weights_elpd_all_models.png"), 
       height = 7, width = 10)

# get auc results ---

AUC_summarise_SPCD <- function(SPCD.id){
  
  # all file names for species models and hierarchical models
  all.species.files <- list.files(path = paste0(output.dir,"SPCD_stanoutput_cmdstan/AUC/"), 
                                  pattern ="AUC_draws_mort_model_", full.names = T)
  
  all.Hierarchical.files <- list.files(path = paste0(output.dir,"SPCD_stanoutput_cmdstan/AUC/"), 
                                       pattern ="_hierarchical_mort_model_", full.names = T)
  
  # get only the species we are looking for...the paste0("_", SPCD.id, "(_|\\.)") ensures we dont read in 129 for spcd == 12
  spp.AUC.files <- grep( paste0("_", SPCD.id, "(_|\\.)"),  
                         all.species.files, 
                         value = TRUE, 
                         perl = TRUE)
  hierarchical.AUC.files <-   grep( paste0("_", SPCD.id, "(_|\\.)"),  
                                    all.Hierarchical.files, 
                                    value = TRUE, 
                                    perl = TRUE)
 
  all.AUC.files <- c(spp.AUC.files, hierarchical.AUC.files)
  
  AUC_results_all <- lapply(all.AUC.files, qs_read)
  
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
AUC.df$COMMON_NAME <- factor(AUC.df$COMMON_NAME, levels = unique(nspp$COMMON))


AUC.df %>% filter(type == "in-sample")|> 
  ggplot()+geom_pointrange(aes(x = as.character(model.number), y = AUC_median, ymin = AUC_ci.lo, ymax = AUC_ci.hi, shape = model.type))+
  facet_wrap(~COMMON_NAME)

AUC.df %>% filter(type == "out-of-sample")|> 
  ggplot()+geom_pointrange(aes(x = as.character(model.number), y = AUC_median, ymin = AUC_ci.lo, ymax = AUC_ci.hi, shape = model.type))+
  facet_wrap(~COMMON_NAME, scales = "free_y")

AUC.df %>% filter(type == "in-sample")|> 
  ggplot()+geom_pointrange(aes(x = as.character(model.number), y = AUC_median, ymin = AUC_ci.lo, ymax = AUC_ci.hi, shape = model.type, color = model.type), position = position_dodge(width = 1))+
  facet_wrap(~COMMON_NAME)+theme_bw()+
  ylab("In-sample AUC score")+
  xlab("Model Number")

ggsave(filename = paste0(output.dir, "SPCD_stanoutput_cmdstan/summary/AUC_is_all_models.png"), 
       height = 7, width = 10)

AUC.df %>% filter(type == "out-of-sample")|> 
  ggplot()+geom_pointrange(aes(x = as.character(model.number), y = AUC_median, ymin = AUC_ci.lo, ymax = AUC_ci.hi, shape = model.type, color = model.type), position = position_dodge(width = 1))+
  facet_wrap(~COMMON_NAME)+theme_bw()+
  ylab("Held-out AUC score")+
  xlab("Model Number")

ggsave(filename = paste0(output.dir, "SPCD_stanoutput_cmdstan/summary/AUC_oos_all_models.png"), 
       height = 7, width = 10)

AUC.df %>% filter(type == "out-of-sample") |> 
  ggplot()+geom_pointrange(aes(x = as.character(model.number), y = True_surv_rate, ymin = True_surv_rate_ci.lo, ymax = True_surv_rate_ci.hi,shape = model.type, color = model.type), position = position_dodge(width = 1))+
  facet_wrap(~COMMON_NAME, scales = "free_y")

AUC.df %>% filter(type == "out-of-sample") |> 
  ggplot()+geom_pointrange(aes(x = as.character(model.number), y = True_mort_rate, ymin = True_mort_rate_ci.lo, ymax = True_mort_rate_ci.hi,shape = model.type, color = model.type), position = position_dodge(width = 1))+
  facet_wrap(~COMMON_NAME, scales = "free_y")


AUC.df %>% filter(type == "in-sample") |> 
  ggplot()+geom_pointrange(aes(x = as.character(model.number), y = True_surv_rate, ymin = True_surv_rate_ci.lo, ymax = True_surv_rate_ci.hi))+
  facet_wrap(~COMMON_NAME, scales = "free_y")

AUC.df %>% filter(type == "in-sample") |> 
  ggplot()+geom_pointrange(aes(x = as.character(model.number), y = True_mort_rate, ymin = True_mort_rate_ci.lo, ymax = True_mort_rate_ci.hi))+
  facet_wrap(~COMMON_NAME, scales = "free_y")

# link up AUC scores with the ELPD and weights for comparisons

AUC.oos.df <- AUC.df %>% filter(type == "out-of-sample")
AUC.is.df <- AUC.df %>% filter(type == "in-sample")

AUC.oos.df |>
  ggplot()+geom_bar(aes( x= model.number, fill = SPCD))

OOS.AUC.ELPD.df <- left_join(LOO_ELPD.df %>% mutate(model.number = as.numeric(model.number)),AUC.oos.df)
IS.AUC.ELPD.df <- left_join(LOO_ELPD.df %>% mutate(model.number = as.numeric(model.number)),AUC.is.df)

OOS.AUC.ELPD.df <- OOS.AUC.ELPD.df %>% 
  mutate(elpd_diff_sig = ifelse(elpd_se_ratio >= 2, "ELPD_diff >= 2 SE", 
                                ifelse(elpd_se_ratio <2, "ELPD_diff < 2 SE", "Best-fit")))%>%
  mutate(elpd_diff_sig = ifelse(is.na(elpd_diff_sig), "Best-fit ELPD", elpd_diff_sig))

IS.AUC.ELPD.df <- IS.AUC.ELPD.df %>% 
  mutate(elpd_diff_sig = ifelse(elpd_se_ratio >= 2, "ELPD_diff >= 2 SE", 
                                ifelse(elpd_se_ratio <2, "ELPD_diff < 2 SE", "Best-fit")))%>%
  mutate(elpd_diff_sig = ifelse(is.na(elpd_diff_sig), "Best-fit ELPD", elpd_diff_sig))






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
                                                 filter(AUC_median == max.AUC)%>% 
                                                 select(SPCD, COMMON_NAME, type, AUC_median, AUC_ci.lo, AUC_ci.hi)%>%
                                                 rename("bf_AUC"= "AUC_median", 
                                                        "bf_AUC_ci.lo"= "AUC_ci.lo", 
                                                        "bf_AUC_ci.hi"= "AUC_ci.hi"))%>%

  mutate(AUC_sig = ifelse(AUC_median == bf_AUC, "Best-fit AUC",
                          ifelse(AUC_ci.lo <= bf_AUC_ci.hi & AUC_ci.hi >= bf_AUC_ci.lo,
                                 "overlapping CI", "non-overlapping CI")))


IS.AUC.ELPD.df <- IS.AUC.ELPD.df%>%left_join(.,
                                             IS.AUC.ELPD.df %>% group_by(SPCD, COMMON_NAME)%>%
                                               mutate(max.AUC = max(AUC_median,na.rm = TRUE))%>%
                                               filter(AUC_median == max.AUC)%>% 
                                               select(SPCD, COMMON_NAME, type, AUC_median, AUC_ci.lo, AUC_ci.hi)%>%
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
  geom_pointrange(aes(x = model.full, y = AUC_median, ymin = AUC_ci.lo, ymax = AUC_ci.hi, color = elpd_diff_sig, shape = AUC_sig))+
  #facet_wrap(~COMMON_NAME, scales = "free")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  ylab("Out-of-Sample AUC")+
  facet_wrap(~COMMON_NAME, scales = "free_y")
ggsave(filename = paste0(output.dir, "SPCD_stanoutput_cmdstan/summary/Species_model_AUC_OOS_elpd.png"), 
       height = 7, width = 10)


OOS.AUC.ELPD.df |>
  ggplot()+
  geom_bar(aes(x = model.name, fill =elpd_diff_sig, group = model.type), stat = "identity")+
  #facet_wrap(~COMMON_NAME, scales = "free")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  ylab("Out-of-Sample AUC")+
  facet_wrap(~COMMON_NAME, scales = "free_y")

# in-sample AUC, colored by best fit elpd
IS.AUC.ELPD.df|>
  ggplot()+
  geom_pointrange(aes(x = model.full, y = AUC_median, ymin = AUC_ci.lo, ymax = AUC_ci.hi, color = elpd_diff_sig, shape = AUC_sig))+
  #facet_wrap(~COMMON_NAME, scales = "free")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  ylab("Out-of-Sample AUC")+
  facet_wrap(~COMMON_NAME, scales = "free_y")
ggsave(filename = paste0(output.dir, "SPCD_stanoutput_cmdstan/summary/Species_model_AUC_IS_elpd.png"), 
       height = 7, width = 10)




OOS.AUC.ELPD.df|>
  ggplot()+
  geom_bar(aes(x = AUC_sig, fill = elpd_diff_sig ))+
  
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  facet_wrap(~COMMON_NAME, scales = "free_y")

IS.AUC.ELPD.df %>% select(model, model.name, model.type, COMMON_NAME, AUC_sig)%>%
  rename("IS_AUC_sig"= "AUC_sig")%>%
  left_join(.,OOS.AUC.ELPD.df) %>% 
  select(model, model.name, model.type, COMMON_NAME, IS_AUC_sig, AUC_sig, elpd_diff_sig)%>%
  pivot_longer( cols = c("AUC_sig","IS_AUC_sig", "elpd_diff_sig"))%>%
  mutate(Comparison.stat = ifelse(name %in% "AUC_sig", "out-of-sample AUC", 
                                  ifelse(name %in% "IS_AUC_sig","in-sample AUC","ELPD-diff")), 
         Cat.signficance = ifelse(value %in% c("Best-fit ELPD", "Best-fit AUC"), "Best-fit", 
                                  ifelse(value %in% c("overlapping CI", "ELPD_diff < 2 SE"), "Overlapping with Best-fit", "Non-overlapping")))|>
  ggplot()+
  geom_tile(aes(y = Comparison.stat, x = model.name, fill = Cat.signficance ), color = "white")+
  
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  scale_fill_manual(values = c("Best-fit" = "darkred", 
                               "Overlapping with Best-fit" = "red", 
                               "Non-overlapping"= "lightgrey"
  ), name = "Significance")+
  facet_grid(vars(COMMON_NAME),vars(model.type), space = "free_y")+
  ylab("Model Comparison Statistic")+
  theme(strip.text.y = element_text(angle = 0))

ggsave(filename = paste0(output.dir, "SPCD_stanoutput_cmdstan/summary/Species_model_AUC_ELPD_comparison_tile.png"), 
       width = 7, height = 15)

OOS.AUC.ELPD.df |>
  group_by(model.name, model.type) %>% filter(type %in% "out-of-sample")%>% 
  summarise(n_ELPD_bf = sum(elpd_diff_sig %in% "Best-fit ELPD"), 
            n_AUC_bf = sum(AUC_sig %in% "Best-fit AUC"))

##########################################################################
# Compare beta estimates -----------
# for each species, load in the betas and alphas. Highlight best fit model based on ELPD_difference and AUC scores

betas.heir <- list.files(paste0(output.dir, "SPCD_stanoutput_cmdstan/betas/"), pattern = "hierarchical_mort_model_", full.names = T)

beta_sumary_df <- do.call(rbind, lapply(1:9, function(i){
  

  ubetas <- qs2::qs_read(betas.heir[i]) %>% subset_draws(., "u_beta")%>%  summarise_draws() %>% mutate(model.number = i) %>% 
    mutate(param = "u_beta")
  
  mu_beta <- qs2::qs_read(betas.heir[i]) %>% subset_draws(., "mu_beta")%>%  summarise_draws() %>% mutate(model.number = i) %>% 
    mutate(param = "mu_beta")
  
  alpha_spp <- qs2::qs_read(betas.heir[i]) %>% subset_draws(., "alpha_SPP")%>%  summarise_draws() %>% mutate(model.number = i) %>% 
    mutate(param = "alpha_SPP")
  
  mu_alpha = qs2::qs_read(betas.heir[i]) %>% subset_draws(., "alpha_SPP")%>%  summarise_draws() %>% mutate(model.number = i) %>% 
    mutate(param = "mu_alpha" )
  rbind(ubetas, mu_beta, alpha_spp, mu_alpha)
})
)

u_betas <- beta_sumary_df %>% filter(param %in% "u_beta")

u_beta_names <- data.frame(variable = unique(u_betas$variable), 
                           SPP = rep(1:17, 78),
                           cov.number = rep(1:78, each = 17))


u_betas_alpha.quant <- u_betas %>% rename("ci.lo"="q5", "ci.hi" = "q95")%>% left_join(., u_beta_names)%>% left_join(., spp.table)

u_betas_alpha.quant$Covariate <-  u_betas_alpha.quant$variable
# get overlapping zero to color the error bars
u_betas_alpha.quant$`significance` <- ifelse(u_betas_alpha.quant$ci.lo < 0 & u_betas_alpha.quant$ci.hi < 0, "significant", 
                                             ifelse(u_betas_alpha.quant$ci.lo > 0 & u_betas_alpha.quant$ci.hi > 0, "significant", "overlapping zero"))


# 
# u_betas_alpha.quant%>% group_by(COMMON, SPCD, cov.number) %>% 
#   summarise(n_pos = sum(median > 0 & significance %in% "significant"), 
#             n_neg = sum(median < 0 & significance %in% "significant")) %>% View()


u_betas_alpha.quant %>% filter(COMMON %in% "sugar maple") |>
  ggplot()+geom_point(aes(x = model.number, y = median, color = model.number, group = COMMON))+
  facet_wrap(~cov.number)

ggplot(data = u_betas_alpha.quant, aes(x = Covariate, y = median, color = model.number))+geom_point()+
  geom_errorbar(data = u_betas_alpha.quant, aes(x = Covariate , ymin = ci.lo, ymax = ci.hi, color = model.number), width = 0.1)+
  geom_abline(aes(slope = 0, intercept = 0), color = "grey", linetype = "dashed")+theme_bw(base_size = 10)+
  theme( axis.text.x = element_text(angle = 65, hjust = 1), panel.grid  = element_blank(), legend.position = "none")+
  ylab("Effect on survival")+xlab("Parameter")+
 # scale_color_manual(values = c("overlapping zero"="darkgrey", "significant"="black"))+
  ggtitle(paste0("Species Model Posterior Estimates"))+
  facet_wrap(~COMMON)

ggsave(height = 5, width = 10, units = "in",
       paste0(output.dir, "/images/Estimated_effects_on_survival_",model.name,".png"))


ggplot(beta_sumary_df %>% filter(param %in% u_beta))+geom_pointrange(aes(x = model.number, y = median, ymin = q5, ymax = q95))+
  facet_wrap(~variable)

# get the original, unscaled datasets

SPCD.id <- 97
load(file = paste0("SPCD_standata_general_full_standardized_v3/SPCD_",SPCD.id,"remper_correction_0.5model_9.Rdata"))
train.data

colnames(train.data)


##########################################################################
# Scaling tree-level survival predictions to plot-level survival/mortality estimates ----
library(qs2)
post_draws_dir <- paste0(output.dir, "SPCD_stanoutput_cmdstan/predicted_mort/")

# step 1: lookup to deal with misnaming of species pSurv-hat vs rep
y_type_lookup <- data.frame(
  actual_y_type = c("y_rep", "y_hat", "pSurv_rep", "pSurv_hat"),
  heir_y_type = c("y_rep", "y_hat", "pSurv_rep", "pSurv_hat"),
  
  SPP_y_type = c("y_hat", "y_rep", "pSurv_hat", "pSurv_rep"),# the mixed up labels for species model
  data_type = c("out-of-sample", "in-sample", "out-of-sample", "in-sample"),
  parameter = c("Survival Assignment", "Survival Assignment", "annual survival probability", "annual survival probability")
)


# function to read the species paths
species_path <- function(y_type, spcd, k) {
  # change the naming of the species files to read in the correct file:
  y_type_real <- y_type_lookup[match(y_type, y_type_lookup$actual_y_type),]$SPP_y_type
  
  # get the filename
  paste0(post_draws_dir, y_type_real, "_samps_mort_model_", k, "_SPCD_", spcd, "_remper_correction_0.5_niter_1000_nchain_4.qs")
}

# function to read the hierarchical paths
hier_path     <- function(y_type, spcd, k) {
  # just to be consistent with species paths
  y_type_real <- y_type_lookup[match(y_type, y_type_lookup$actual_y_type),]$heir_y_type
  paste0(post_draws_dir, y_type_real, "_samps_SPCD_", spcd, "_hierarchical_mort_model_",k,"_niter_1000_nchain_4.qs")
}



# get posterior matrix reads and outputs these paths
get_posterior_matrix_renamed <- function(structure,  spcd, y_type, k, spp_index = NULL) {
  
  if (structure == "species_only") {
    full <- qs2::qs_read(species_path( y_type, spcd, k))
    colnames(full) <- paste0(y_type, "_",1:length(colnames(full))) # rename species column names
    full
    
  } else {
    full <- qs2::qs_read(hier_path(y_type, spcd, k))        
    colnames(full) <- paste0(y_type, "_",1:length(colnames(full))) # rename species column names
    full
  }
}



SPCD.id <- 531
k <- 7


# for each species, read in the in-sample and out-of-sample annual probabilities:
# read in tree_remeas to get the volfac for each represented tree:
TREE.remeas <- readRDS( "data/unfiltered_TREE.remeas.rds")

all.remeas <- TREE.remeas %>%
  # do the filtering section
  filter( exprem > 0 & # if exprem == 0, these could be modeled plots?
            dbhold >= 5 & # need an initial dbh greater than 5
            ! remper == 0 & # if remper is listed as zero, filter out
            DIA_DIFF > 0 & # filter diameter differences >= 0
            # !status == 3 & # keep the cut trees for this
            SPCD %in% nspp[1:17,]$SPCD & # filter species in the top 17 of all species
            !is.na(status) & # filter out trees with no status
            !is.na(elev)) 

PLOT <- read_delim(paste0(output.dir,"data/formatted_older_matching_plts_PLOT.txt"))

plot_expansions <- PLOT %>% select(PLOT.ID, state, county, pltnum, cndtn, cycle,expacr) %>% distinct()


calculate_state_county_rates <- function(k, SPCD.id, model.type){
        
        cat("\n",paste0("estimating county and state mortality for: ", SPCD.id,",",model.type, " model ", k))
          
        # read in the in-sample and out of sample probability of survival pSurv over remper for each model k
        load(file = paste0("SPCD_standata_general_full_standardized_v3/SPCD_",SPCD.id,"remper_correction_0.5model_9.Rdata"))
        
        pSurv_hat_samps <- get_posterior_matrix_renamed(structure = ifelse(model.type %in% "Species","species_only","Hierarchical"),  
                             spcd = SPCD.id, 
                             y_type = "pSurv_hat",
                             k = k)
        
        pSurv_rep_samps <- get_posterior_matrix_renamed(structure = ifelse(model.type %in% "Species","species_only","Hierarchical"),   
                                                spcd = SPCD.id, 
                                                y_type = "pSurv_rep",
                                                k = k)
        
        actuals     <- mod.data$y
        actuals.oos <- mod.data$ytest
      
        # convert remper probabilities to annualized probabilities of survival:
        Remper_matrix <- matrix(mod.data$Remper, nrow = nrow(pSurv_hat_samps), ncol = ncol(pSurv_hat_samps), byrow = TRUE)
        
        pSannual_hat <- pSurv_hat_samps^(1/Remper_matrix) 
        
        Remperoos_matrix <- matrix(mod.data$Remperoos, nrow = nrow(pSurv_rep_samps), ncol = ncol(pSurv_rep_samps), byrow = TRUE)
        pSannual_rep <- pSurv_rep_samps^(1/Remperoos_matrix) 
        
        qs2::qs_save(pSannual_rep, paste0(output.dir, "SPCD_stanoutput_cmdstan/predicted_mort/pSurv_annual_rep_samps_mort_model_",k,"_SPCD_", SPCD.id, "_remper_correction_0.5_niter_1000_nchain_4.qs"))
        qs2::qs_save(pSannual_hat, paste0(output.dir, "SPCD_stanoutput_cmdstan/predicted_mort/pSurv_annual_hat_samp_mort_model_",k,"_SPCD_", SPCD.id, "_remper_correction_0.5_niter_1000_nchain_4.qs"))
        
          # convert pSannuals to posterior expected deaths over the remper scale
          volfac_hat_matrix <- matrix(train.data$volfac, nrow = nrow(pSurv_hat_samps), ncol = ncol(pSurv_hat_samps), byrow = TRUE)
          volfac_rep_matrix <- matrix(test.data$volfac, nrow = nrow(pSurv_rep_samps), ncol = ncol(pSurv_rep_samps), byrow = TRUE)
          
          # calculate the posterior expected # of trees dead based on tree-level volfac
          Post_Emort_remp_volfac_rep <- (1-pSurv_rep_samps)*volfac_rep_matrix
          Post_Emort_remp_volfac_hat <- (1-pSurv_hat_samps)*volfac_hat_matrix
          
          # calculate matrices of observed exposure for these predictions (#trees/acre * n years) tree-years
          Exposure_mat_rep <-  volfac_rep_matrix*Remperoos_matrix
          Exposure_mat_hat <-  volfac_hat_matrix*Remper_matrix
          
          
          # aggregate by states and counties:----------
          
          state_county_lookup_train <- train.data %>% mutate(tree.id = 1:length(train.data$state))%>% 
            left_join(.,plot_expansions)%>%
            select(tree.id, remper, volfac, S, pltnum, state, county, expacr) %>% 
            mutate(data.type = "in-sample", 
                   ST_CTY = paste0(state, "_", county))
          
          state_county_lookup_test <- test.data %>% mutate(tree.id = 1:length(test.data$state))%>% 
            left_join(.,plot_expansions)%>%
            select(tree.id, remper, volfac, S, pltnum, state, county, expacr) %>% 
            mutate(data.type = "out-of-sample", 
                   ST_CTY = paste0(state, "_", county))
          
          
          # get expacr matrices:
          Expacr_rep_matrix <- matrix(state_county_lookup_test$expacr, 
                                      nrow = nrow(pSurv_rep_samps), ncol = ncol(pSurv_rep_samps), byrow = TRUE)
          
          Expacr_hat_matrix  <- matrix(state_county_lookup_train$expacr, 
                                      nrow = nrow(pSurv_hat_samps), ncol = ncol(pSurv_hat_samps), byrow = TRUE)
          
          
          
          # weight posterior mortality rates by expacr for county and state scales
          # calculate the posterior expected # of trees dead based on tree-level volfac
          Post_Emort_remp_volfac_expacre_rep <- Expacr_rep_matrix*((1-pSurv_rep_samps)*volfac_rep_matrix)
          Post_Emort_remp_volfac_expacre_hat <- Expacr_hat_matrix*(1-pSurv_hat_samps)*volfac_hat_matrix
          
          # calculate matrices of observed exposure for these predictions (#trees/acre * n years) tree-years
          Exposure_expacre_mat_rep <-   Expacr_rep_matrix*(volfac_rep_matrix*Remperoos_matrix)
          Exposure_expacre_mat_hat <-   Expacr_hat_matrix*(volfac_hat_matrix*Remper_matrix)
          
   
  # overall estimates of observed and predicted mortality rates       
  
          
          # calculate tree-level mortality rates
          st_Pmort_rep <- Post_Emort_remp_volfac_rep
          st_Expos_rep <- Exposure_mat_rep
          st_obs_train <- state_county_lookup_test%>% 
            mutate(Exposure_i = volfac*remper)%>%
            mutate(Deaths_i = volfac*(1-S))%>% 
            mutate(Exposure_EXPN_i = Exposure_i*expacr, 
                   Deaths_EXPN_i = Deaths_i*expacr)
          
          st_Pmort_hat <- Post_Emort_remp_volfac_hat
          st_Expos_hat <- Exposure_mat_hat
          st_obs_test <- state_county_lookup_train %>% 
            mutate(Exposure_i = volfac*remper)%>%
            mutate(Deaths_i = volfac*(1-S)) %>% 
            mutate(Exposure_EXPN_i = Exposure_i*expacr, 
                   Deaths_EXPN_i = Deaths_i*expacr)
          
          
          # do the same with expacre weighted values:
          st_Pmort_rep_EXPN <- Post_Emort_remp_volfac_expacre_rep
          st_Expos_rep_EXPN <- Exposure_expacre_mat_rep
          
          
          st_Pmort_hat_EXPN <- Post_Emort_remp_volfac_expacre_hat
          st_Expos_hat_EXPN <- Exposure_expacre_mat_hat
          
          # # tree-level mortality rates--
          # tree_E_M_expn_i <- cbind(st_Pmort_hat_EXPN, st_Pmort_rep_EXPN)
          # tree_Expos_expn_i <- cbind(st_Expos_hat_EXPN, st_Expos_rep_EXPN)
          # tree_Pred_mort_expn <- tree_E_M_expn_i/tree_Expos_expn_i
          # 
          # 
          # tree_obs_Expos_expn_i = c(st_obs_train$Exposure_EXPN_i, st_obs_test$Exposure_EXPN_i)# calculating the "exposure for each tree" remper* volfac, in tree/acre-years
          # tree_obs_Deaths_expn_i= c(st_obs_train$Deaths_EXPN_i, st_obs_test$Deaths_EXPN_i) # calculate the number of deaths per tree observed over remper
          # 
          # tree_obs_mort_expn <- tree_obs_Deaths_expn_i/tree_obs_Expos_expn_i
          # 
          # tree_mort_rate_summary <- data.frame(Obs_mort_expn = tree_obs_mort_expn, 
          #                                      Pred_mort_expn.mean = colMeans(tree_Pred_mort_expn))
          # 
          # 
          
          # Regional mortality rates:
          OVERALL_mort_rate <- data.frame(E_mort = rowSums(cbind(st_Pmort_hat, st_Pmort_rep)), 
                                     Exposure = rowSums(cbind(st_Expos_hat, st_Expos_rep)),
                                     E_mort_expn = rowSums(cbind(st_Pmort_hat_EXPN, st_Pmort_rep_EXPN)), 
                                     Exposure_expn = rowSums(cbind(st_Expos_hat_EXPN, st_Expos_rep_EXPN)), 
                                     
                                     
                                     Exposure_obs = sum(st_obs_test$Exposure_i, st_obs_train$Exposure_i),# calculating the "exposure for each tree" remper* volfac, in tree/acre-years
                                     Deaths_obs = sum(st_obs_test$Deaths_i, st_obs_train$Deaths_i),
                                     
                                     Exposure_expn_obs = sum(st_obs_test$Exposure_EXPN_i, st_obs_train$Exposure_EXPN_i),# calculating the "exposure for each tree" remper* volfac, in tree/acre-years
                                     Deaths_expn_obs = sum(st_obs_test$Deaths_EXPN_i, st_obs_train$Deaths_EXPN_i),#, # calculate the number of deaths per tree observed over remper
                                     total_obs = length(c(st_obs_test$Deaths_i, st_obs_train$Deaths_i)),
                                     
                                     SPCD = SPCD.id, 
                                     model.number = k, 
                                     model.type = model.type) %>%
            mutate(Obs_mort_rate =  Deaths_obs/Exposure_obs,
                   Pred_mort_rate = E_mort/Exposure, 
                   Obs_mort_rate_expn = Deaths_expn_obs/Exposure_expn_obs, 
                   Pred_mort_rate_expn = E_mort_expn/Exposure_expn)
  
          
  
  
  # calculate state-level estimates of observed and predicted mortality rates
          
          # unique states
          state_ids <- unique(c(state_county_lookup_train$state, state_county_lookup_test$state))
          
          state_mort_rate_list <- lapply(state_ids, FUN = function(state_cd){
            
            # get index for the focal state
            st_index_train <- state_county_lookup_train$state == state_cd
            st_index_test <- state_county_lookup_test$state == state_cd
            
            # calculate tree level mortality rates
            st_Pmort_rep <- Post_Emort_remp_volfac_rep[,st_index_test]
            st_Expos_rep <- Exposure_mat_rep[,st_index_test]
            st_obs_train <- state_county_lookup_test[st_index_test,]%>% 
              mutate(Exposure_i = volfac*remper)%>%
              mutate(Deaths_i = volfac*(1-S))%>% 
              mutate(Exposure_EXPN_i = Exposure_i*expacr, 
                     Deaths_EXPN_i = Deaths_i*expacr)
            
            st_Pmort_hat <- Post_Emort_remp_volfac_hat[,st_index_train]
            st_Expos_hat <- Exposure_mat_hat[,st_index_train]
            st_obs_test <- state_county_lookup_train[st_index_train,] %>% 
              mutate(Exposure_i = volfac*remper)%>%
              mutate(Deaths_i = volfac*(1-S)) %>% 
              mutate(Exposure_EXPN_i = Exposure_i*expacr, 
                     Deaths_EXPN_i = Deaths_i*expacr)
            
            
            # do the same with expacre weighted values:
            st_Pmort_rep_EXPN <- Post_Emort_remp_volfac_expacre_rep[,st_index_test]
            st_Expos_rep_EXPN <- Exposure_expacre_mat_rep[,st_index_test]
           
            
            st_Pmort_hat_EXPN <- Post_Emort_remp_volfac_expacre_hat[,st_index_train]
            st_Expos_hat_EXPN <- Exposure_expacre_mat_hat[,st_index_train]
         
            
            # combine together:
            st_mort_rate <- data.frame(E_mort = rowSums(cbind(st_Pmort_hat, st_Pmort_rep)), 
                                       Exposure = rowSums(cbind(st_Expos_hat, st_Expos_rep)),
                                       E_mort_expn = rowSums(cbind(st_Pmort_hat_EXPN, st_Pmort_rep_EXPN)), 
                                       Exposure_expn = rowSums(cbind(st_Expos_hat_EXPN, st_Expos_rep_EXPN)), 
                                       
                                       state = state_cd,
                                       Exposure_obs = sum(st_obs_test$Exposure_i, st_obs_train$Exposure_i),# calculating the "exposure for each tree" remper* volfac, in tree/acre-years
                                       Deaths_obs = sum(st_obs_test$Deaths_i, st_obs_train$Deaths_i),
                                       
                                       Exposure_expn_obs = sum(st_obs_test$Exposure_EXPN_i, st_obs_train$Exposure_EXPN_i),# calculating the "exposure for each tree" remper* volfac, in tree/acre-years
                                       Deaths_expn_obs = sum(st_obs_test$Deaths_EXPN_i, st_obs_train$Deaths_EXPN_i),#, # calculate the number of deaths per tree observed over remper
                                       total_obs = length(c(st_obs_test$Deaths_i, st_obs_train$Deaths_i)),
                                       
                                       SPCD = SPCD.id, 
                                       model.number = k, 
                                       model.type = model.type) %>%
              mutate(Obs_mort_rate =  Deaths_obs/Exposure_obs,
                     Pred_mort_rate = E_mort/Exposure, 
                     Obs_mort_rate_expn = Deaths_expn_obs/Exposure_expn_obs, 
                     Pred_mort_rate_expn = E_mort_expn/Exposure_expn)
            return(st_mort_rate)
          })
          
          state_mort_df <- do.call(rbind,state_mort_rate_list) 
  
          # get county level predicted mortality summaries:
          STCTY_ids <- unique(c(state_county_lookup_train$ST_CTY, state_county_lookup_test$ST_CTY))
         # ST_CT_id <- STCTY_ids[1]
          
          county_mort_rate_list <- lapply(STCTY_ids, FUN = function(ST_CT_id){
            
            # get index for the focal state
            st_index_train <- state_county_lookup_train$ST_CTY == ST_CT_id
            st_index_test <- state_county_lookup_test$ST_CTY == ST_CT_id
            
            # calculate tree level expected mortality and exposure 
            st_Pmort_rep <- Post_Emort_remp_volfac_rep[,st_index_test]
            st_Expos_rep <- Exposure_mat_rep[,st_index_test]
            st_obs_train <- state_county_lookup_test[st_index_test,]%>% 
              mutate(Exposure_i = volfac*remper)%>%
              mutate(Deaths_i = volfac*(1-S))%>% 
              mutate(Exposure_EXPN_i = Exposure_i*expacr, 
                     Deaths_EXPN_i = Deaths_i*expacr)
            
            st_Pmort_hat <- Post_Emort_remp_volfac_hat[,st_index_train]
            st_Expos_hat <- Exposure_mat_hat[,st_index_train]
            st_obs_test <- state_county_lookup_train[st_index_train,] %>% 
              mutate(Exposure_i = volfac*remper)%>%
              mutate(Deaths_i = volfac*(1-S)) %>% 
              mutate(Exposure_EXPN_i = Exposure_i*expacr, 
                     Deaths_EXPN_i = Deaths_i*expacr)
            
            
            # do the same with expacre weighted values:
            st_Pmort_rep_EXPN <- Post_Emort_remp_volfac_expacre_rep[,st_index_test]
            st_Expos_rep_EXPN <- Exposure_expacre_mat_rep[,st_index_test]
            
            
            st_Pmort_hat_EXPN <- Post_Emort_remp_volfac_expacre_hat[,st_index_train]
            st_Expos_hat_EXPN <- Exposure_expacre_mat_hat[,st_index_train]
            
            
            
            # combine together:
            st_mort_rate <- data.frame(E_mort = rowSums(cbind(st_Pmort_hat, st_Pmort_rep)), 
                                       Exposure = rowSums(cbind(st_Expos_hat, st_Expos_rep)), 
                                       E_mort_expn = rowSums(cbind(st_Pmort_hat_EXPN, st_Pmort_rep_EXPN)), 
                                       Exposure_expn = rowSums(cbind(st_Expos_hat_EXPN, st_Expos_rep_EXPN)), 
                                       
                                       
                                       ST_CTY = ST_CT_id,
                                       
                                       Exposure_obs = sum(st_obs_test$Exposure_i, st_obs_train$Exposure_i),# calculating the "exposure for each tree" remper* volfac, in tree/acre-years
                                       Deaths_obs = sum(st_obs_test$Deaths_i, st_obs_train$Deaths_i),#, # calculate the number of deaths per tree observed over remper
                                       Exposure_expn_obs = sum(st_obs_test$Exposure_EXPN_i, st_obs_train$Exposure_EXPN_i),# calculating the "exposure for each tree" remper* volfac, in tree/acre-years
                                       Deaths_expn_obs = sum(st_obs_test$Deaths_EXPN_i, st_obs_train$Deaths_EXPN_i),#, # calculate the number of deaths per tree observed over remper
                                       
                                       total_obs = length(c(st_obs_test$Deaths_i, st_obs_train$Deaths_i)),
                                       SPCD = SPCD.id, 
                                       model.number = k, 
                                       model.type = model.type) %>%
              
              mutate(Obs_mort_rate =  Deaths_obs/Exposure_obs,
                     Pred_mort_rate = E_mort/Exposure,
                     Obs_mort_rate_expn = Deaths_expn_obs/Exposure_expn_obs, 
                     Pred_mort_rate_expn = E_mort_expn/Exposure_expn)
            return(st_mort_rate)
          })
          
          county_mort_df <- do.call(rbind, county_mort_rate_list) 
          
          county_mort_summary <- county_mort_df %>%  
            group_by(SPCD, ST_CTY, model.number, model.type)%>%
            summarise(obs_M_median = median(Obs_mort_rate, na.rm =TRUE), 
                    n_obs = median(total_obs, na.rm =TRUE),
                    pred_M_median = median(Pred_mort_rate, na.rm =TRUE), 
                    pred_M_5.ci.lo = quantile(Pred_mort_rate, 0.05, na.rm =TRUE), 
                    pred_M_95.ci.hi = quantile(Pred_mort_rate, 0.95, na.rm =TRUE),
                    pred_M_25.ci.lo = quantile(Pred_mort_rate, 0.25, na.rm =TRUE), 
                    pred_M_75.ci.hi = quantile(Pred_mort_rate, 0.75, na.rm =TRUE), 
                    
                    obs_M_expn_median = median(Obs_mort_rate_expn, na.rm =TRUE), 
                  
                    pred_M_expn_median = median(Pred_mort_rate_expn, na.rm =TRUE), 
                    pred_M_expn_5.ci.lo = quantile(Pred_mort_rate_expn, 0.05, na.rm =TRUE), 
                    pred_M_expn_95.ci.hi = quantile(Pred_mort_rate_expn, 0.95, na.rm =TRUE),
                    pred_M_expn_25.ci.lo = quantile(Pred_mort_rate_expn, 0.25, na.rm =TRUE), 
                    pred_M_expn_75.ci.hi = quantile(Pred_mort_rate_expn, 0.75, na.rm =TRUE), .groups = "drop")
        
     
        # save the county-level outputs for each species:
        # save the samples weighted by volface of each species in the county
        
         qs_save(county_mort_df, paste0(
            output.dir,
            "SPCD_stanoutput_cmdstan/predicted_mort/ST_CTY_mort_rate_samps_",SPCD.id,"_",model.type, "_", k,".qs"
          ))
         qs_save(state_mort_df, paste0(
           output.dir,
           "SPCD_stanoutput_cmdstan/predicted_mort/State_mort_rate_samps_",SPCD.id,"_",model.type, "_", k,".qs"
         ))
         
         qs_save(OVERALL_mort_rate, paste0(
           output.dir,
           "SPCD_stanoutput_cmdstan/predicted_mort/Regional_mort_rate_samps_",SPCD.id,"_",model.type, "_", k,".qs"
         ))
         
         return(county_mort_summary)
}


# get all county-level estimates for each species:
species.ests <- hierarchical.ests <- list()
for(s in length(spp.table$SPCD):1){
  
  species.ests[[s]] <- do.call(rbind, lapply(1:9, FUN  = function(i){calculate_state_county_rates(i, spp.table$SPCD[s], "Species")}))
  hierarchical.ests[[s]] <- do.call(rbind, lapply(1:9, FUN  = function(i){calculate_state_county_rates(i, spp.table$SPCD[s], "Hierarchical")}))

}


county_mort_summary%>% filter(n_obs > 25)|>
  ggplot()+geom_point(aes(x = obs_M_expn_median, y = pred_M_expn_median))+
  geom_errorbar(aes(x = obs_M_expn_median, ymin = pred_M_expn_5.ci.lo, ymax = pred_M_expn_95.ci.hi))+
  geom_errorbar(aes(x = obs_M_expn_median, ymin = pred_M_expn_25.ci.lo, ymax = pred_M_expn_75.ci.hi))
  
