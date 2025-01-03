###########################################################################################
# Model Validation and comparison plots with both the joint model and species-level models
###########################################################################################
library(loo)
library(ggplot2)
library(tidyverse)
library(FIESTA)
options(mc.cores = parallel::detectCores())
##### SINGLE SPECIES MODELS
#read in all the accuracy dataframes from the single species models
output.folder <- "C:/Users/KellyHeilman/Box/01. kelly.heilman Workspace/mortality/Eastern-Mortality/mortality_models/"


#paste0(output.folder, "SPCD_stanoutput_full/Accuracy_df_model_",,SPCD.df$SPCD
accuracy.files <- list.files(path = paste0(output.folder, "SPCD_stanoutput_full_standardized/computational_resources/"), pattern = c("Accuracy_df_model_"))

accuracy.files.full <- paste0(output.folder, "SPCD_stanoutput_full_standardized/computational_resources/", accuracy.files)

accuracy.list <- lapply(accuracy.files.full, read.csv)
accuracy.df <- do.call(rbind, accuracy.list) %>% filter(remper.correction %in% 0.5)
accuracy.df$Species <- FIESTA::ref_species[match(accuracy.df$SPCD, FIESTA::ref_species$SPCD),]$COMMON
accuracy.df$Model.name <- paste0("model ", accuracy.df$model)
accuracy.df$Size_effect <- "Linear"

AUC.singlespecies <- accuracy.df %>% dplyr::select(SPCD, auc.oosample, auc.insample, Species, Model.name)

is.auc <- ggplot(AUC.singlespecies, aes(x = Model.name, y = auc.insample, shape = Model.name %in% "model 6"))+geom_point()+
 # geom_hline(data = AUC.all %>% filter(Model.name %in% "hierarchical"), aes(yintercept =  auc.insample), linetype = "dashed", color = "red")+
  facet_wrap(~Species, scales =  "free_y")+
  scale_color_manual( values = c( "black" , 
                                  "red" ))+
  theme_bw()+theme(axis.text.x = element_text(angle = 45,  hjust=1), legend.position = "none") +
  xlab("")+ylab("In Sample AUC")

ggsave(paste0("model_summary_full/All_species_models_all9models_compare-auc-insample.png"), 
       is.auc,
       width = 10, height = 6)


oos.auc <- ggplot(AUC.singlespecies, aes(x = Model.name, y = auc.oosample))+geom_point()+facet_wrap(~Species, scales =  "free_y")+
  theme_bw()+theme(axis.text.x = element_text(angle = 45,  hjust=1)) +
  xlab("")+ylab("Out of Sample AUC")
ggsave(paste0("model_summary_full/All_species_models_all9models_compare-auc-outofsample.png"), 
       oos.auc,
       width = 10, height = 6)



AUC.summary <- AUC.singlespecies %>% group_by(Model.name) %>% 
  summarise(auc.oos.median = median(auc.oosample), 
            auc.is.median = median(auc.insample), 
            auc.oos.mean = mean(auc.oosample), 
            auc.is.mean = mean(auc.insample), 
            auc.oos.total = sum(auc.oosample), 
            auc.is.total = sum(auc.insample))


ggplot(AUC.summary, aes(x = Model.name, y = auc.is.mean, shape = Model.name %in% "model 6"))+geom_point()+
  geom_hline(data = AUC.summary %>% filter(Model.name %in% "hierarchical"), aes(yintercept =  auc.is.mean), linetype = "dashed", color = "red")+
  #facet_wrap(~Species, scales =  "free_y")+
  scale_color_manual( values = c( "black" , 
                                  "red" ))+
  theme_bw()+theme(axis.text.x = element_text(angle = 45,  hjust=1), legend.position = "none") +
  xlab("")+ylab("Mean In Sample AUC")


ggplot(AUC.summary, aes(x = Model.name, y = auc.oos.mean, shape = Model.name %in% "model 6"))+geom_point()+
  geom_hline(data = AUC.summary %>% filter(Model.name %in% "hierarchical"), aes(yintercept =  auc.oos.mean), linetype = "dashed", color = "red")+
  #facet_wrap(~Species, scales =  "free_y")+
  scale_color_manual( values = c( "black" , 
                                  "red" ))+
  theme_bw()+theme(axis.text.x = element_text(angle = 45,  hjust=1), legend.position = "none") +
  xlab("")+ylab("out of Sample AUC mean")
ggsave("model_summary_full/AUC_oos_mean_species_model.png")

ggplot(AUC.summary, aes(x = Model.name, y = auc.oos.median, shape = Model.name %in% "model 6"))+geom_point()+
  geom_hline(data = AUC.summary %>% filter(Model.name %in% "hierarchical"), aes(yintercept =  auc.oos.median), linetype = "dashed", color = "red")+
  #facet_wrap(~Species, scales =  "free_y")+
  scale_color_manual( values = c( "black" , 
                                  "red" ))+
  theme_bw()+theme(axis.text.x = element_text(angle = 45,  hjust=1), legend.position = "none") +
  xlab("")+ylab("out of Sample AUC median")
ggsave("model_summary_full/AUC_oos_median_species_model.png")

ggplot(AUC.summary, aes(x = Model.name, y = auc.oos.total, shape = Model.name %in% "model 6"))+geom_point()+
  geom_hline(data = AUC.summary %>% filter(Model.name %in% "hierarchical"), aes(yintercept =  auc.oos.median), linetype = "dashed", color = "red")+
  #facet_wrap(~Species, scales =  "free_y")+
  scale_color_manual( values = c( "black" , 
                                  "red" ))+
  theme_bw()+theme(axis.text.x = element_text(angle = 45,  hjust=1), legend.position = "none") +
  xlab("")+ylab("out of Sample AUC total")
ggsave("model_summary_full/AUC_oos_total_species_model.png")

# read in the comptational efficiency results
compute.files <- list.files(path = paste0(output.folder,"SPCD_stanoutput_full_standardized/computational_resources/"), pattern = "time_diag_SPCD_")

compute.files.full <- paste0(output.folder,"SPCD_stanoutput_full_standardized/computational_resources/", compute.files)
compute.list <- lapply(compute.files.full, read.csv)

compute.df <- do.call(rbind, compute.list) %>% filter(remper %in% 0.5)
compute.df$Species <- FIESTA::ref_species[match(compute.df$SPCD, FIESTA::ref_species$SPCD),]$COMMON
compute.df$Model.name <- paste0("model ", compute.df$model)



# get the total number of in-sample observations
all.full.data <- list.files(paste0(output.folder, "SPCD_standata_general_full/"), pattern = "remper_correction_0.5model_1.Rdata")
Nobs <- list()

for(i in 1:length(all.full.data)){
  load(paste0(output.folder,"SPCD_standata_general_full/", all.full.data[[i]]))
   
  Nobs[[i]] <- data.frame(SPCD = unique(test.data$SPCD),
            N = mod.data$N, 
             Nrep = mod.data$Nrep)
}

Nobs <- do.call(rbind, Nobs)

# chec these values
model.complexity <- data.frame(model = 1:9, 
                               nbetas = c(1, 2, 5, 13, 18, 51, 103, 145, 159))
compute.df <- left_join(compute.df, model.complexity)
compute.df <- left_join(compute.df, Nobs)

compute.df <- compute.df %>% mutate(Nparams = N*3 + Nrep*2 + nbetas + 1)

ggplot()+geom_text(data = compute.df, aes(x =  nbetas, y = core.hours, label = model ))+
  facet_wrap(~Species, scales = "free_y")+
  theme_bw()+theme(panel.grid = element_blank())+ylab("Core Hours")+xlab("Number of Fixed Effects")
ggsave(height = 7, width = 9,"model_summary_full/Core_hours_fixed_effects_total_species_model.png")

# combine core hours with accuracy
compute.df <- left_join(compute.df, accuracy.df)

ggplot()+geom_text(data = compute.df, aes( x = core.hours, y =  auc.oosample, label = model ))+
  facet_wrap(~Species, scales = "free")+
  theme_bw()+theme(panel.grid = element_blank())+ylab("Out-of-sample AUC")+xlab("Core Hours")
ggsave(height = 7, width = 9,"model_summary_full/Core_hours_oos_accuracy_total_species_model.png")

## basis function accuracy 
# read in the DBH spline model computational and accuracy metrics
accuracy.basis.files <- list.files(path = "SPCD_stanoutput_full/", pattern = "Accuracy_df_basis")

accuracy.basis.files.full <- paste0("SPCD_stanoutput_full/", accuracy.basis.files)
accuracy.basis.list <- lapply(accuracy.basis.files.full, read.csv)
accuracy.basis.df <- do.call(rbind, accuracy.basis.list) %>% filter(remper.correction %in% 0.5)
accuracy.basis.df$Species <- FIESTA::ref_species[match(accuracy.basis.df$SPCD, FIESTA::ref_species$SPCD),]$COMMON
accuracy.basis.df$Model.name <- paste0("model ", accuracy.basis.df$model)
accuracy.basis.df$Size_effect <- "Spline"

accuracy.both.df <- rbind(accuracy.df, accuracy.basis.df)

is.auc <- ggplot()+
  geom_point(data = accuracy.both.df, aes(x = Model.name, y = auc.insample,  color = Size_effect))+
  #geom_point(data = accuracy.df, aes(x = Model.name, y = auc.insample, shape = Model.name %in% "model 6"), color = "red")+
  # geom_hline(data = AUC.all %>% filter(Model.name %in% "hierarchical"), aes(yintercept =  auc.insample), linetype = "dashed", color = "red")+
  facet_wrap(~Species, scales =  "free_y")+
  scale_color_manual( values = c( "black" , 
                                  "red" ))+
  theme_bw()+theme(axis.text.x = element_text(angle = 45,  hjust=1), legend.position = "bottom") +
  xlab("")+ylab("In Sample AUC")

ggsave(paste0("model_summary_full/All_species_models_6_basis_dbh_models_compare-auc-insample.png"), 
       is.auc,
       width = 10, height = 6)


oos.auc <- ggplot()+
  geom_point(data = accuracy.both.df, aes(x = Model.name, y = auc.oosample, color =Size_effect))+
 # geom_point(data = accuracy.df, aes(x = Model.name, auc.oosample), color = "red")+facet_wrap(~Species, scales =  "free_y")+
  theme_bw()+theme(axis.text.x = element_text(angle = 45,  hjust=1), legend.position = "bottom") +
  xlab("")+ylab("Out of Sample AUC")
ggsave(paste0("model_summary_full/All_species_models_6_basis_dbh_models_compare-auc-outofsample.png"), 
       oos.auc,
       width = 10, height = 6)

# read in the computational hours for the basis function models
# read in the comptational efficiency results
compute.basis.files <- list.files(path = "SPCD_stanoutput_full/computational_resources/", pattern = "time_diag_BASIS_SPCD_")

compute.basis.files.full <- paste0("SPCD_stanoutput_full/computational_resources/", compute.basis.files)
compute.basis.list <- lapply(compute.basis.files.full, read.csv)

compute.basis.df <- do.call(rbind, compute.basis.list) %>% filter(remper %in% 0.5)
compute.basis.df$Species <- FIESTA::ref_species[match(compute.basis.df$SPCD, FIESTA::ref_species$SPCD),]$COMMON
compute.basis.df$Model.name <- paste0("model ", compute.basis.df$model)

compute.basis.df <- left_join(compute.basis.df, model.complexity)
compute.basis.df <- left_join(compute.basis.df, Nobs)

compute.basis.df <- compute.basis.df %>% mutate(nbetas = nbetas + 8)%>%  mutate(Nparams = N*3 + Nrep*2 + nbetas + 1)

compute.basis.df <- left_join(compute.basis.df, accuracy.basis.df)

ggplot()+geom_text(data = compute.basis.df, aes(x =  nbetas, y = core.hours, label = model ))+
  facet_wrap(~Species, scales = "free_y")+
  theme_bw()+theme(panel.grid = element_blank())+ylab("Core Hours")+xlab("Number of Fixed Effects")
ggsave(height = 7, width = 9,"model_summary_full/Core_hours_fixed_effects_total_species__basis_model.png")


ggplot()+geom_text(data = compute.basis.df, aes( x = core.hours, y =  auc.oosample, label = model ))+
  facet_wrap(~Species, scales = "free")+
  theme_bw()+theme(panel.grid = element_blank())+ylab("Out-of-sample AUC")+xlab("Core Hours")
ggsave(height = 7, width = 9,"model_summary_full/Core_hours_oos_accuracy_total_species_basis_model.png")

# is there any computational difference between basis function models and non basis function models?

compute.all.df <- rbind(compute.basis.df,compute.df)
ggplot()+geom_text(data = compute.all.df, aes( x = core.hours, y =  auc.oosample, label = model , color = Size_effect))+
  scale_color_manual( values = c( "black" , 
                                  "red" ))+
  facet_wrap(~Species, scales = "free")+
  theme_bw()+theme(panel.grid = element_blank())+ylab("Out-of-sample AUC")+xlab("Core Hours")
ggsave(height = 7, width = 9,"model_summary_full/Core_hours_oos_accuracy_total_species_basis_and_linear_model.png")


# Plot up % increase in accuracy relative to the last model fit vs core hour cost difference
# for both size effect types
compute.all.diff <- compute.all.df %>% group_by(SPCD, Species, Size_effect) %>%
  mutate(accuracy.diff = auc.oosample - lag(auc.oosample), 
         core.hr.diff = core.hours - lag(core.hours)) %>%
  mutate(pct.accuracy.diff = accuracy.diff/lag(auc.oosample), 
         pct.core.hr.diff = core.hr.diff/lag(core.hours))

ggplot()+geom_text(data = compute.all.diff, aes(x = pct.accuracy.diff, y = pct.core.hr.diff, label = model, color= Size_effect))+
  facet_wrap(~Species, scales = "free")+
  theme_bw()+theme(panel.grid = element_blank())+ylab("% change in core hours")+xlab("% change in accuracy")
ggsave(height = 7, width = 9,"model_summary_full/Change_in_Core_hours_change_oos_accuracy_total_species_basis_and_linear_model.png")


ggplot()+geom_point(data = compute.all.diff, aes(x = as.character(model), y = pct.core.hr.diff, color= Size_effect))+
  facet_wrap(~Species, scales = "free")+
  theme_bw()+theme(panel.grid = element_blank())+ylab("% change in core hours")+xlab("% change in accuracy")
#ggsave(height = 7, width = 9,"model_summary_full/Change_in_Core_hours_change_oos_accuracy_total_species_basis_and_linear_model.png")
compute.all.diff <- compute.all.diff %>% mutate(ratio.diff = pct.accuracy.diff/pct.core.hr.diff)
ggplot()+geom_point(data = compute.all.diff, aes(x = as.character(model), y = ratio.diff, color= Size_effect))+
  facet_wrap(~Species, scales = "free")+
  theme_bw()+theme(panel.grid = element_blank())+ylab("% change in core hours")+xlab("% change in accuracy")

# read in the AUC values from the joint model:
AUC.joint <- read.csv(paste0(output.folder,"SPCD_stanoutput_joint_v2/Accuracy_df_model_6_remper_0.5_species_joint_model_remper_corr_0.5.csv"))
AUC.joint <- AUC.joint %>% rename(Species = COMMON) %>% dplyr::select(SPCD, auc.oosample, auc.insample, Species) %>%
  mutate(Model.name = "hierarchical")
AUC.all <- rbind(AUC.singlespecies, AUC.joint)
AUC.all$Model <- ifelse(AUC.all$Model.name %in% "hierarchical", "hierarchical", "species model")
AUC.all$Model.name <- factor(AUC.all$Model.name, levels = unique(AUC.all$Model.name))
AUC.all <- AUC.all %>% filter(! Species %in% "population")
is.auc <- ggplot(AUC.all , aes(x = Model.name, y = auc.insample, color = Model, shape = Model.name %in% "model 6"))+geom_point()+
  geom_hline(data = AUC.all %>% filter(Model.name %in% "hierarchical"), aes(yintercept =  auc.insample), linetype = "dashed", color = "red")+
  facet_wrap(~Species, scales =  "free_y")+
  scale_color_manual( values = c("species model" = "black" , 
                                "hierarchical"="red" ))+
  theme_bw()+theme(axis.text.x = element_text(angle = 45,  hjust=1), legend.position = "none") +
  xlab("")+ylab("In Sample AUC")

ggsave(paste0("model_summary_full/All_species_models_all9models_compare-auc-insample_plus_joint.png"), 
       is.auc,
       width = 10, height = 6)




oos.auc <- ggplot(AUC.all, aes(x = Model.name, y = auc.oosample, color = Model, shape = Model.name %in% "model 6"))+geom_point()+
  geom_hline(data = AUC.all %>% filter(Model.name %in% "hierarchical"), aes(yintercept =  auc.oosample), linetype = "dashed", color = "red")+
  facet_wrap(~Species, scales =  "free_y")+
  scale_color_manual( values = c("species model" = "black" , 
                                 "hierarchical"="red" ))+
  theme_bw()+theme(axis.text.x = element_text(angle = 45,  hjust=1), legend.position = "none") +
  xlab("")+ylab("Out of Sample AUC")
ggsave(paste0("model_summary_full/All_species_models_all9models_compare-auc-outofsample_plus_joint.png"), 
       oos.auc,
       width = 10, height = 6)

# make figure with the total AUC, mean AUC, and median AUC for oos at the beginning
AUC.summary.m <- AUC.summary %>% select(Model.name, auc.oos.median, auc.oos.mean, auc.oos.total)%>%
  rename(`Median AUC` = auc.oos.median, 
         `Mean AUC` = auc.oos.mean,
         `Total AUC` = auc.oos.total)%>%
  reshape2::melt(., id.vars = "Model.name")

median.auc.p <- ggplot(data = AUC.summary.m %>% filter(variable %in% "Median AUC"), aes(x = Model.name, y = value, shape = Model.name %in% "model 6"))+geom_point()+
  theme_bw()+theme(axis.text.x = element_text(angle = 45,  hjust=1), legend.position = "none")+xlab("")+
  ylab("Out of Sample AUC")+facet_wrap(~variable)+coord_flip()

mean.auc.p <- ggplot(data = AUC.summary.m %>% filter(variable %in% "Mean AUC"), aes(x = Model.name, y = value, shape = Model.name %in% "model 6"))+geom_point()+
  theme_bw()+theme(axis.text.x = element_text(angle = 45,  hjust=1), legend.position = "none")+xlab("")+
  ylab("Out of Sample AUC")+facet_wrap(~variable)+coord_flip()

total.auc.p <- ggplot(data =AUC.summary.m %>% filter(variable %in% "Total AUC"), aes(x = Model.name, y = value, shape = Model.name %in% "model 6"))+geom_point()+
  theme_bw()+theme(axis.text.x = element_text(angle = 45,  hjust=1), legend.position = "none")+xlab("")+
  ylab("Out of Sample AUC")+facet_wrap(~variable)+coord_flip()


# Then individually plot the joint model and single model outputs:
AUC.best.fit.df <- AUC.all %>% filter(!Model.name %in% "hierarchical")%>% group_by(Species) %>% mutate(max.model = max(auc.oosample))%>%
  mutate(best.fit = ifelse(auc.oosample == max.model, TRUE, FALSE))


AUC.all <- left_join(AUC.all, AUC.best.fit.df)
AUC.all <- AUC.all %>% mutate(Model = ifelse(is.na(best.fit), "hierarchical", 
                                  ifelse(best.fit == TRUE, "best fit species", "non best fit model")))

beech.AUC <- ggplot(AUC.all %>% filter(Species %in% "American beech"), aes(x = Model.name, y = auc.oosample, color = Model, shape = Model.name %in% "model 6"))+geom_point()+facet_wrap(~Species, scales =  "free_y")+
  theme_bw()+
  xlab("")+ylab("Out of Sample AUC")+
  geom_hline(data = AUC.all %>% filter(Species %in% "American beech" & Model.name %in% "hierarchical"), aes(yintercept =  auc.oosample), linetype = "dashed", color = "red")+
  scale_color_manual( values = c( "black", 
                                  "red",
                                  "darkgrey" ))+theme(axis.text.x = element_text(angle = 45,  hjust=1), legend.position = "none") +coord_flip()

balsam.AUC <-  ggplot(AUC.all %>% filter(Species %in% "balsam fir"), aes(x = Model.name, y = auc.oosample, color = Model, shape = Model.name %in% "model 6"))+geom_point()+facet_wrap(~Species, scales =  "free_y")+
  theme_bw()+ xlab("")+ylab("Out of Sample AUC")+
  geom_hline(data = AUC.all %>% filter(Species %in% "balsam fir" & Model.name %in% "hierarchical"), aes(yintercept =  auc.oosample), linetype = "dashed", color = "red")+
  scale_color_manual(  values = c( "black", 
                                   "red",
                                   "darkgrey" ))+theme(axis.text.x = element_text(angle = 45,  hjust=1), legend.position = "none") +coord_flip()


blackcherry.AUC <-  ggplot(AUC.all %>% filter(Species %in% "black cherry"), aes(x = Model.name, y = auc.oosample, color = Model, shape = Model.name %in% "model 6"))+geom_point()+facet_wrap(~Species, scales =  "free_y")+
  theme_bw()+ xlab("")+ylab("Out of Sample AUC")+
  geom_hline(data = AUC.all %>% filter(Species %in% "black cherry" & Model.name %in% "hierarchical"), aes(yintercept =  auc.oosample), linetype = "dashed", color = "red")+
  scale_color_manual( values = c( "black", 
                                  "red",
                                  "darkgrey" ))+theme(axis.text.x = element_text(angle = 45,  hjust=1), legend.position = "none") +coord_flip()


blackoak.AUC <-  ggplot(AUC.all %>% filter(Species %in% "black oak"), aes(x = Model.name, y = auc.oosample, color = Model, shape = Model.name %in% "model 6"))+geom_point()+facet_wrap(~Species, scales =  "free_y")+
  theme_bw()+ xlab("")+ylab("Out of Sample AUC")+
  geom_hline(data = AUC.all %>% filter(Species %in% "black oak" & Model.name %in% "hierarchical"), aes(yintercept =  auc.oosample), linetype = "dashed", color = "red")+
  scale_color_manual(  values = c( "black", 
                                   "red",
                                   "darkgrey" ))+theme(axis.text.x = element_text(angle = 45,  hjust=1), legend.position = "none") +coord_flip()


chestnutoak.AUC <-  ggplot(AUC.all %>% filter(Species %in% "chestnut oak"), aes(x = Model.name, y = auc.oosample, color = Model, shape = Model.name %in% "model 6"))+geom_point()+facet_wrap(~Species, scales =  "free_y")+
  theme_bw()+ xlab("")+ylab("Out of Sample AUC")+
  geom_hline(data = AUC.all %>% filter(Species %in% "chestnut oak" & Model.name %in% "hierarchical"), aes(yintercept =  auc.oosample), linetype = "dashed", color = "red")+
  scale_color_manual(  values = c( "black", 
                                   "red",
                                   "darkgrey" ))+theme(axis.text.x = element_text(angle = 45,  hjust=1), legend.position = "none") +coord_flip()


hemlock.AUC <- ggplot(AUC.all %>% filter(Species %in% "eastern hemlock"), aes(x = Model.name, y = auc.oosample, color = Model, shape = Model.name %in% "model 6"))+geom_point()+facet_wrap(~Species, scales =  "free_y")+
  theme_bw()+ xlab("")+ylab("Out of Sample AUC")+
  geom_hline(data = AUC.all %>% filter(Species %in% "eastern hemlock" & Model.name %in% "hierarchical"), aes(yintercept =  auc.oosample), linetype = "dashed", color = "red")+
  scale_color_manual(  values = c( "black", 
                                   "red",
                                   "darkgrey" ))+theme(axis.text.x = element_text(angle = 45,  hjust=1), legend.position = "none") +coord_flip()

easternWP.AUC <-  ggplot(AUC.all %>% filter(Species %in% "eastern white pine"), aes(x = Model.name, y = auc.oosample, color = Model, shape = Model.name %in% "model 6"))+geom_point()+facet_wrap(~Species, scales =  "free_y")+
  theme_bw()+ xlab("")+ylab("Out of Sample AUC")+
  geom_hline(data = AUC.all %>% filter(Species %in% "eastern white pine" & Model.name %in% "hierarchical"), aes(yintercept =  auc.oosample), linetype = "dashed", color = "red")+
  scale_color_manual(  values = c( "black", 
                                   "red",
                                   "darkgrey" ))+theme(axis.text.x = element_text(angle = 45,  hjust=1), legend.position = "none") +coord_flip()


hickory.AUC <-  ggplot(AUC.all %>% filter(Species %in% "hickory spp."), aes(x = Model.name, y = auc.oosample, color = Model, shape = Model.name %in% "model 6"))+geom_point()+facet_wrap(~Species, scales =  "free_y")+
  theme_bw()+ xlab("")+ylab("Out of Sample AUC")+
  geom_hline(data = AUC.all %>% filter(Species %in% "hickory spp." & Model.name %in% "hierarchical"), aes(yintercept =  auc.oosample), linetype = "dashed", color = "red")+
  scale_color_manual(  values = c( "black", 
                                   "red",
                                   "darkgrey" ))+theme(axis.text.x = element_text(angle = 45,  hjust=1), legend.position = "none") +coord_flip()


northernRedOak.AUC <-  ggplot(AUC.all %>% filter(Species %in% "northern red oak"), aes(x = Model.name, y = auc.oosample, color = Model, shape = Model.name %in% "model 6"))+geom_point()+facet_wrap(~Species, scales =  "free_y")+
  theme_bw()+ xlab("")+ylab("Out of Sample AUC")+
  geom_hline(data = AUC.all %>% filter(Species %in% "northern red oak" & Model.name %in% "hierarchical"), aes(yintercept =  auc.oosample), linetype = "dashed", color = "red")+
  scale_color_manual( values = c( "black", 
                                  "red",
                                  "darkgrey" ))+theme(axis.text.x = element_text(angle = 45,  hjust=1), legend.position = "none") +coord_flip()


NWhiteCedar.AUC <-  ggplot(AUC.all %>% filter(Species %in% "northern white-cedar"), aes(x = Model.name, y = auc.oosample, color = Model, shape = Model.name %in% "model 6"))+geom_point()+facet_wrap(~Species, scales =  "free_y")+
  theme_bw()+ xlab("")+ylab("Out of Sample AUC")+
  geom_hline(data = AUC.all %>% filter(Species %in% "northern white-cedar" & Model.name %in% "hierarchical"), aes(yintercept =  auc.oosample), linetype = "dashed", color = "red")+
  scale_color_manual(  values = c( "black", 
                                   "red",
                                   "darkgrey" ))+theme(axis.text.x = element_text(angle = 45,  hjust=1), legend.position = "none") +coord_flip()


RedMaple.AUC <-  ggplot(AUC.all %>% filter(Species %in% "red maple"), aes(x = Model.name, y = auc.oosample, color = Model, shape = Model.name %in% "model 6"))+geom_point()+facet_wrap(~Species, scales =  "free_y")+
  theme_bw()+ xlab("")+ylab("Out of Sample AUC")+
  geom_hline(data = AUC.all %>% filter(Species %in% "red maple" & Model.name %in% "hierarchical"), aes(yintercept =  auc.oosample), linetype = "dashed", color = "red")+
  scale_color_manual( values = c( "black", 
                                  "red",
                                  "darkgrey" ))+theme(axis.text.x = element_text(angle = 45,  hjust=1), legend.position = "none") +coord_flip()


RedSpruce.AUC <-  ggplot(AUC.all %>% filter(Species %in% "red spruce"), aes(x = Model.name, y = auc.oosample, color = Model, shape = Model.name %in% "model 6"))+geom_point()+facet_wrap(~Species, scales =  "free_y")+
  theme_bw()+ xlab("")+ylab("Out of Sample AUC")+
  geom_hline(data = AUC.all %>% filter(Species %in% "red spruce" & Model.name %in% "hierarchical"), aes(yintercept =  auc.oosample), linetype = "dashed", color = "red")+
  scale_color_manual(  values = c( "black", 
                                   "red",
                                   "darkgrey" ))+theme(axis.text.x = element_text(angle = 45,  hjust=1), legend.position = "none") +coord_flip()


SugarMaple.AUC <-  ggplot(AUC.all %>% filter(Species %in% "sugar maple"), aes(x = Model.name, y = auc.oosample, color = Model, shape = Model.name %in% "model 6"))+geom_point()+facet_wrap(~Species, scales =  "free_y")+
  theme_bw()+ xlab("")+ylab("Out of Sample AUC")+
  geom_hline(data = AUC.all %>% filter(Species %in% "sugar maple" & Model.name %in% "hierarchical"), aes(yintercept =  auc.oosample), linetype = "dashed", color = "red")+
  scale_color_manual(  values = c( "black", 
                                   "red",
                                   "darkgrey" ))+theme(axis.text.x = element_text(angle = 45,  hjust=1), legend.position = "none") +coord_flip()


WhiteAsh.AUC <-  ggplot(AUC.all %>% filter(Species %in% "white ash"), aes(x = Model.name, y = auc.oosample, color = Model, shape = Model.name %in% "model 6"))+geom_point()+facet_wrap(~Species, scales =  "free_y")+
  theme_bw()+ xlab("")+ylab("Out of Sample AUC")+
  geom_hline(data = AUC.all %>% filter(Species %in% "white ash" & Model.name %in% "hierarchical"), aes(yintercept =  auc.oosample), linetype = "dashed", color = "red")+
  scale_color_manual(  values = c( "black", 
                                   "red",
                                   "darkgrey" ))+theme(axis.text.x = element_text(angle = 45,  hjust=1), legend.position = "none") +coord_flip()

WhiteOak.AUC <-  ggplot(AUC.all %>% filter(Species %in% "white oak"), aes(x = Model.name, y = auc.oosample, color = Model, shape = Model.name %in% "model 6"))+geom_point()+facet_wrap(~Species, scales =  "free_y")+
  theme_bw()+ xlab("")+ylab("Out of Sample AUC")+
  geom_hline(data = AUC.all %>% filter(Species %in% "white oak" & Model.name %in% "hierarchical"), aes(yintercept =  auc.oosample), linetype = "dashed", color = "red")+
  scale_color_manual(  values = c( "black", 
                                   "red",
                                   "darkgrey" ))+theme(axis.text.x = element_text(angle = 45,  hjust=1), legend.position = "none") +coord_flip()


YellowPoplar.AUC <-  ggplot(AUC.all %>% filter(Species %in% "yellow-poplar"), aes(x = Model.name, y = auc.oosample, color = Model, shape = Model.name %in% "model 6"))+geom_point()+facet_wrap(~Species, scales =  "free_y")+
  theme_bw()+ xlab("")+ylab("Out of Sample AUC")+
  geom_hline(data = AUC.all %>% filter(Species %in% "yellow-poplar" & Model.name %in% "hierarchical"), aes(yintercept =  auc.oosample), linetype = "dashed", color = "red")+
  scale_color_manual(  values = c( "black", 
                                   "red",
                                   "darkgrey" ))+theme(axis.text.x = element_text(angle = 45,  hjust=1), legend.position = "none")+coord_flip() 


YellowBirch.AUC <-  ggplot(AUC.all %>% filter(Species %in% "yellow birch"), aes(x = Model.name, y = auc.oosample, color = Model, shape = Model.name %in% "model 6"))+geom_point()+facet_wrap(~Species, scales =  "free_y")+
  theme_bw()+ xlab("")+ylab("Out of Sample AUC")+
  geom_hline(data = AUC.all %>% filter(Species %in% "yellow birch" & Model.name %in% "hierarchical"), aes(yintercept =  auc.oosample), linetype = "dashed", color = "red")+
  scale_color_manual(  values = c( "black", 
                                   "red",
                                   "darkgrey" ))+theme(axis.text.x = element_text(angle = 45,  hjust=1), legend.position = "none") +coord_flip()


png(height = 10, width = 11, units = "in", res = 300, "model_summary_full/all_AUC_sumary_with_joint_summary.png")
cowplot::plot_grid(median.auc.p, mean.auc.p, total.auc.p, beech.AUC, balsam.AUC, 
          blackcherry.AUC, blackoak.AUC, chestnutoak.AUC, hemlock.AUC, easternWP.AUC, 
          hickory.AUC, northernRedOak.AUC, NWhiteCedar.AUC, RedMaple.AUC, RedSpruce.AUC, 
          SugarMaple.AUC,  WhiteAsh.AUC, WhiteOak.AUC, YellowPoplar.AUC, YellowBirch.AUC,
          align = "hv", ncol = 5, labels = "AUTO")
dev.off()

#--------------------------------------------------------------------------------------------------
# make the same computational efficiency figures, but now with the joint-level model as a comparison
# need to redo with basis function update
# read in the computational files for the joint model
joint.compute <- read.csv(paste0(output.folder, "SPCD_stanoutput_joint_v2/joint_model_time_diag_SPCD_joint_model_6_remper_0.5.csv"))
joint.accuracy <- read.csv(paste0(output.folder, "SPCD_stanoutput_joint_v2/Accuracy_df_model_6_remper_0.5_species_joint_model_remper_corr_0.5.csv"))

joint.all.df <- left_join(joint.accuracy, joint.compute)
#colnames(compute.all.df)
colnames(joint.all.df)

joint.all.df  <- joint.all.df %>% rename(`Species` = "COMMON") %>% 
  select(Species, SPCD, model, remper.correction, auc.insample, auc.oosample, accuracy.oos, core.hours) %>%
  mutate(Size_effect = "Linear") %>% 
  mutate(variety = "hierarchical")

compute.all.df <- compute.df %>% 
  select(Species, SPCD, model, remper.correction, auc.insample, auc.oosample, accuracy.oos, core.hours, Size_effect) %>%

    mutate(variety = ifelse(Size_effect == "Linear", "Species", "Species basis"))

joint.species.compute <- rbind(joint.all.df, compute.all.df) %>% filter(!Species %in% "population")
joint.species.compute <- joint.species.compute %>% mutate(variety = ifelse(variety %in% "hierarchical", "Hierarchical", variety))
ggplot()+geom_text(data = joint.species.compute, aes( x = core.hours, y =  auc.oosample, label = model , color = variety))+
  scale_color_manual( values = c("red", "black", "darkgrey"))+ # "#1b9e77",
                                   # "#d95f02",
                                    #"#7570b3"))+
  facet_wrap(~Species, scales = "free")+
  theme_bw()+theme(panel.grid = element_blank(), legend.title = element_blank())+ylab("Out-of-sample AUC")+xlab("Core Hours")
ggsave(height = 7, width = 9,"model_summary_full/Core_hours_oos_accuracy_total_species_basis_and_linear_and_joint_model.png")


# Plot up % increase in accuracy relative to the last model fit vs core hour cost difference
# for both size effect types
joint.species.compute.diff <- joint.species.compute %>% group_by(SPCD, Species, Size_effect, variety) %>%
  mutate(accuracy.diff = auc.oosample - lag(auc.oosample), 
         core.hr.diff = core.hours - lag(core.hours)) %>%
  mutate(pct.accuracy.diff = accuracy.diff/lag(auc.oosample), 
         pct.core.hr.diff = core.hr.diff/lag(core.hours))

# for model 6, lets compare model computational efficiency:
model6.summary <- joint.species.compute %>% filter(model == 6) %>% group_by(variety) %>% 
  summarise(median.auc.oos = median(auc.oosample), 
         median.auc.ins = median(auc.insample), 
         total.core.hours = sum(core.hours))
model6.summary$variety <- factor(model6.summary$variety, levels = c("Species", "Species basis", "hierarchical"))

auc.plot.6 <- ggplot(model6.summary, aes(x = variety, y = median.auc.oos, fill = variety))+geom_bar(stat = "identity")+
  scale_fill_manual( values = c( "black", "darkgrey", "red"))+
  ylab("Median out-of-sample AUC \n (model 6)")+xlab("Model type")+theme_bw(base_size = 14) + theme(legend.position = "none")

core.hours.plot.6 <- ggplot(model6.summary, aes(x = variety, y = total.core.hours, fill = variety))+geom_bar(stat = "identity")+
  scale_fill_manual( values = c( "black", "darkgrey", "red"))+
  ylab("Total Core Hours")+xlab("Model type")+theme_bw(base_size = 14) + theme( legend.position = "none")

png(height = 4, width = 8, units = "in", res = 300, "model_summary_full/model_6_computational_efficiency_time.png")
cowplot::plot_grid(auc.plot.6, core.hours.plot.6, labels = "AUTO")
dev.off()  
#--------------------------------------------------------------------------------------------------

# compare model estimates to species-level betas
# get the complete spcies list
cleaned.data <- readRDS( "data/cleaned.data.mortality.TRplots.RDS")
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


joint.betas <- readRDS("SPCD_stanoutput/samples/model_6_u_beta_last1000samples.RDS")
joint.betas.summary <- joint.betas %>% select(-.iteration, -.chain, -.draw)%>% reshape2::melt() %>%
                                group_by(variable)%>% summarise(median = quantile(value, 0.5), 
                                                                ci.lo = quantile(value, 0.025), 
                                                                ci.hi = quantile(value, 0.975), 
                                                                log.odds.median = quantile(exp(value), 0.5), 
                                                                log.odds.ci.lo = quantile(exp(value), 0.025), 
                                                                log.odds.ci.hi = quantile(exp(value), 0.975))

# read in the data used to fit
mod.data.full <- readRDS ( paste0("full_stan_data/three_SPCD_model_6.RDS"))
param.id = rep(1:48, each = 17)
spp.id = rep(1:17, 48)
beta.names.df <- data.frame(variable = unique(joint.betas.summary$variable), 
                            spp = spp.id, 
                            Param.id = param.id, 
                            Parameter.name = rep(colnames(mod.data.full$xM), each = 17))
joint.betas.summary <- joint.betas.summary %>% left_join(.,beta.names.df) %>% left_join(., spp.table)
joint.betas.summary$significance <- ifelse(joint.betas.summary$ci.lo > 0 & joint.betas.summary$ci.hi > 0, "significant", 
                                           ifelse(joint.betas.summary$ci.lo < 0 & joint.betas.summary$ci.hi < 0, "significant", "not significant"))

ggplot()+geom_point(data = joint.betas.summary %>% filter(significance %in% "significant"), aes(x = COMMON, y = log.odds.median))+
  geom_errorbar(data = joint.betas.summary %>% filter(significance %in% "significant"), aes(x = COMMON, ymin = log.odds.ci.lo, ymax = log.odds.ci.hi))+
  facet_wrap(~Parameter.name, scales = "free_y")+theme(axis.text.x = element_text(hjust = 1, angle = 45))


ggplot()+geom_point(data = joint.betas.summary  %>% filter(significance %in% "significant"), aes(x = Parameter.name, y = log.odds.median, color = log.odds.median >=1))+
  geom_errorbar(data = joint.betas.summary  %>% filter(significance %in% "significant"), aes(x = Parameter.name, ymin = log.odds.ci.lo, ymax = log.odds.ci.hi, color = log.odds.median >=1))+
  geom_hline(yintercept = 1, color = "red", linetype = "dashed")+facet_wrap(~COMMON, scales = "free")+theme(axis.text.x = element_text(hjust = 1, angle = 45), legend.position = "none")
ggsave(filename = "SPCD_stanoutput/images/joint_model_log_odds_betas_significant.png", 
       height = 10, width =12, units = "in")

