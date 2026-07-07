library(tidyverse)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(qs2)
library(jsonlite)
library(pROC)
color_scheme_set("brightblue")
# script to read in cmdstan model outputs, generate model assessement and comparisons

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
# LOO_summarise_SPCD <- function(SPCD.id){
#   spp.loo.files <- paste0(output.dir,"SPCD_stanoutput_cmdstan/LOO/LOO_results_mort_model_", 1:9, 
#                           "_SPCD_", SPCD.id, "_remper_correction_0.5_niter_1000_nchain_4.qs")
#   
#   loo_results_all <- lapply(spp.loo.files, qs_read)
#   
#   
#   # check pareto-k estimates:
#   pareto.k.checks <- do.call(rbind, lapply(loo_results_all, function(x){data.frame(good = sum(x$diagnostics$pareto_k <= 0.7), 
#                                                                                    bad = sum(x$diagnostics$pareto_k > 0.7), 
#                                                                                    total = length(x$diagnostics$pareto_k))}))%>%
#     mutate(percent.bad = (bad/total)*100)%>%
#     mutate(model = paste0("model", 1:9), 
#            SPCD = SPCD.id)
#   
#   
#   
#   loo_compare.out <- loo::loo_compare(loo_results_all) # best fit based on loo elpd differences
#   loo_comparisons <- loo_compare.out %>% data.frame()%>%left_join(., pareto.k.checks) %>%
#     mutate(elpd_se_ratio = abs(elpd_diff)/se_diff)%>%
#     mutate(model.number = substr(model, start = 6, stop = 6))
#   
#   
#   # Get model weights---
#   # get pointwise log predictive densities
#   lpd_point <- do.call(cbind,lapply(loo_results_all, function(x){x$pointwise[,"elpd_loo"]}))
#   pbma_wts <- pseudobma_weights(lpd_point, BB=FALSE)
#   pbma_BB_wts <- pseudobma_weights(lpd_point) # default is BB=TRUE
#   stacking_wts <- stacking_weights(lpd_point)
#   mod.weights <- round(cbind(pbma_wts, pbma_BB_wts, stacking_wts),3)%>% data.frame() %>% 
#     mutate(model = paste0("model", 1:9), 
#            SPCD = SPCD.id)
#   loo_comparisons <- loo_comparisons %>% left_join(., mod.weights, by = c("model", "SPCD"))
#   return(loo_comparisons)
# }
# 

LOO_summarise_SPCD <- function(SPCD.id){

  cat(paste0("Doing LOO comparison for ", SPCD.id))
  spp.loo.files <- paste0(output.dir,"SPCD_stanoutput_cmdstan/LOO/LOO_results_mort_model_", 1:9,
                          "_SPCD_", SPCD.id, "_remper_correction_0.5_niter_1000_nchain_4.qs")
  Hierarchical.loo.files <- paste0(output.dir,"SPCD_stanoutput_cmdstan/LOO/LOO_results_hierarchical_mort_model_",1:6,"_niter_1000_nchain_4_SPCD_",
                                   SPCD.id, ".qs")
  loo.files.all <- c(spp.loo.files, Hierarchical.loo.files)
  loo_results_all <- lapply(loo.files.all, qs_read)


  # check pareto-k estimates:
  pareto.k.checks <- do.call(rbind, lapply(loo_results_all, function(x){data.frame(good = sum(x$diagnostics$pareto_k <= 0.7),
                                                                                   bad = sum(x$diagnostics$pareto_k > 0.7),
                                                                                   total = length(x$diagnostics$pareto_k))}))%>%
    mutate(percent.bad = (bad/total)*100)%>%
    mutate(model.name = c(paste0("model ", 1:9),paste0("model ", 1:6)),
           model = c(paste0("model", 1:15)),
           model.type = c(rep("Species", 9), rep("Hierarchical", 6)),
           SPCD = SPCD.id)



  loo_compare.out <- loo::loo_compare(loo_results_all) # best fit based on loo elpd differences
  loo_comparisons <- loo_compare.out %>% data.frame()%>%left_join(., pareto.k.checks) %>%
    mutate(elpd_se_ratio = abs(elpd_diff)/se_diff)%>%
    mutate(model.number = substr(model, start = 6, stop = 6))


  # Get model weights---
  # get pointwise log predictive densities
  lpd_point <- do.call(cbind,lapply(loo_results_all, function(x){x$pointwise[,"elpd_loo"]}))
  pbma_wts <- loo::pseudobma_weights(lpd_point, BB=FALSE)
  pbma_BB_wts <- loo::pseudobma_weights(lpd_point) # default is BB=TRUE
  stacking_wts <- loo::stacking_weights(lpd_point)
  mod.weights <- round(cbind(pbma_wts, pbma_BB_wts, stacking_wts),3)%>% data.frame() %>%
    mutate(model = paste0("model", 1:15),
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
  hierarchical.AUC.files <- paste0(output.dir,"SPCD_stanoutput_cmdstan/AUC/AUC_draws_hierarchical_mort_model_",1:6,"_niter_1000_nchain_4_SPCD_",
                                   SPCD.id, ".qs")
  all.AUC.fils <- c(spp.AUC.files, hierarchical.AUC.files)
  
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


##########################################################################
# TODO: Using model weights (stacking or BMA) to average posterior predictive distributions:
# do this for parameter inference on alphas and betas (?)

# outline of the code:
# for each species, pull up the model stacking weights

spcd_loo_weights <- LOO_ELPD.df %>% filter(SPCD == 316)

get_weighted_species_mod_betas <- function(SPCD.id){
    # loo weights
    spcd_loo_weights <- LOO_ELPD.df %>% filter(SPCD == SPCD.id)
    
    # beta coefficients
    # 1. weight posterior betas by their model stacking weight
    # read in all u_betas_for each species:
    spp_beta.files <- paste0(output.dir,"SPCD_stanoutput_cmdstan/betas/u_beta_alpha_samps_mort_model_", 1:9, 
                                                     "_SPCD_", SPCD.id, "_remper_correction_0.5_niter_1000_nchain_4.qs")
    
    
    spp_betas <- lapply(spp_beta.files, FUN = qs2::qs_read)
    
    
    spp.stacking.weights <- spcd_loo_weights %>% arrange(model.number) %>% select(model.number, stacking_wts)
    
    spp_betas_weighted <- lapply(1:9, FUN = function(x) {
      (spp_betas[[x]]*spp.stacking.weights[x,]$stacking_wts) %>% as_draws_df() %>% reshape2::melt(id.vars = c(".chain", ".iteration", ".draw"))})
    
    spp_betas_weighted_all <- do.call(rbind, spp_betas_weighted) %>% group_by(variable, .chain, .iteration, .draw)%>%
      summarise(weighted_draw = sum(value, na.rm =TRUE))
    
    spp_betas_weighted_summary <- spp_betas_weighted_all %>% ungroup()%>%
      group_by(variable)%>%
      summarise(weighted_median = median(weighted_draw, na.rm =TRUE), 
                weighted_ci_lo.2.5 = quantile(weighted_draw, 0.025, na.rm =TRUE), 
                weighted_ci_hi.97.5 = quantile(weighted_draw, 0.975, na.rm = TRUE), 
                weighted_ci_lo.5 = quantile(weighted_draw, 0.05, na.rm =TRUE), 
                weighted_ci_hi.95 = quantile(weighted_draw, 0.95, na.rm =TRUE))%>%
      mutate(SPCD = SPCD.id)
    
    return(spp_betas_weighted_summary)
}

# get for all of the species
species_mod_weighted_betas <- lapply(unique(nspp$SPCD), get_weighted_species_mod_betas)
spcd_mod_wt_betas <- do.call(rbind, species_mod_weighted_betas)


# 2. compare across species
spcd_mod_wt_betas %>%  ungroup()%>% mutate(significant = ifelse(weighted_ci_lo.2.5 < 0 & weighted_ci_hi.97.5 < 0 , "significant", 
                                                                   ifelse(weighted_ci_lo.2.5 > 0 & weighted_ci_hi.97.5 > 0, "significant", "n.s.")))%>%
  filter(!variable %in% "alpha_SPP")|>
  ggplot()+geom_pointrange(aes(x = variable, y = weighted_median, ymin = weighted_ci_lo.2.5, ymax = weighted_ci_hi.97.5, color = significant))+
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1))+facet_wrap(~SPCD, scales = "free_y")


# model averaged marginal effects---
# 1. for each model & species, run generated quantities over a prediction grid 
#       - annual mortality predictions (p_annual) in response to covariate
#       - 10 year mortality probability predictions in response to covariate

# 2. multiply draws by stacking weights 
# p_annual_1 * w_1 + p_annual_2 * w_2 + .....p_annual_9 + w_9

# 3. Plot up weighted marginal effects by species....



