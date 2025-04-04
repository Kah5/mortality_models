###########################################################################################
# Model Validation and comparison plots with both the joint model and species-level models
###########################################################################################
#library(loo)
library(ggplot2)
library(tidyverse)
library(FIESTA)
options(mc.cores = parallel::detectCores())
##### SINGLE SPECIES MODELS
#read in all the accuracy dataframes from the single species models
output.folder <- "/home/rstudio"


#paste0(output.folder, "SPCD_stanoutput_full/Accuracy_df_model_",,SPCD.df$SPCD
accuracy.files <- list.files(path = paste0(output.folder, "/SPCD_stanoutput_full_standardized_v3/"), pattern = c("Accuracy_df_model_"))

accuracy.files.full <- paste0(output.folder, "/SPCD_stanoutput_full_standardized_v3/", accuracy.files)

accuracy.list <- lapply(accuracy.files.full, read.csv)
accuracy.df <- do.call(rbind, accuracy.list) %>% filter(remper.correction %in% 0.5)
accuracy.df$Species <- FIESTA::ref_species[match(accuracy.df$SPCD, FIESTA::ref_species$SPCD),]$COMMON
accuracy.df$Model.name <- paste0("model ", accuracy.df$model)
accuracy.df$Size_effect <- "Linear"

AUC.singlespecies <- accuracy.df #%>% dplyr::select(SPCD, auc.oosample, auc.insample, Species, Model.name)

is.auc <- ggplot(AUC.singlespecies, aes(x = Model.name, y = auc.insample.median, shape = Model.name %in% "model 6"))+geom_point()+
  geom_errorbar(data = AUC.singlespecies, aes(x = Model.name, ymin = auc.insample.lo, ymax = auc.insample.hi))+
  # geom_hline(data = AUC.all %>% filter(Model.name %in% "hierarchical"), aes(yintercept =  auc.insample), linetype = "dashed", color = "red")+
  facet_wrap(~Species, scales =  "free_y")+
  scale_color_manual( values = c( "black" , 
                                  "red" ))+
  theme_bw()+theme(axis.text.x = element_text(angle = 45,  hjust=1), legend.position = "none") +
  xlab("")+ylab("In Sample AUC")

ggsave(paste0("model_summary_full/All_species_models_all9models_compare-auc-insample.png"), 
       is.auc,
       width = 10, height = 6)


oos.auc <- ggplot(AUC.singlespecies, aes(x = Model.name, y = auc.oosample.median))+geom_point()+facet_wrap(~Species, scales =  "free_y")+
  geom_errorbar(data = AUC.singlespecies, aes(x = Model.name, ymin = auc.oosample.lo, ymax = auc.oosample.hi))+
  theme_bw()+theme(axis.text.x = element_text(angle = 45,  hjust=1)) +
  xlab("")+ylab("Out of Sample AUC")
ggsave(paste0("model_summary_full/All_species_models_all9models_compare-auc-outofsample.png"), 
       oos.auc,
       width = 10, height = 6)



AUC.summary <- AUC.singlespecies %>% group_by(Model.name) %>% 
  summarise(auc.oos.median = median(auc.oosample.median), 
            auc.is.median = median(auc.insample.median), 
            auc.oos.mean = mean(auc.oosample.median), 
            auc.is.mean = mean(auc.insample.median), 
            auc.oos.total = sum(auc.oosample.median), 
            auc.is.total = sum(auc.insample.median))


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
compute.files <- list.files(path = paste0(output.folder,"/SPCD_stanoutput_full_standardized_v3/computational_resources/"), pattern = "time_diag_SPCD_")

compute.files.full <- paste0(output.folder,"/SPCD_stanoutput_full_standardized_v3/computational_resources/", compute.files)
compute.list <- lapply(compute.files.full, read.csv)

compute.df <- do.call(rbind, compute.list) %>% filter(remper %in% 0.5)
compute.df$Species <- FIESTA::ref_species[match(compute.df$SPCD, FIESTA::ref_species$SPCD),]$COMMON
compute.df$Model.name <- paste0("model ", compute.df$model)



# get the total number of in-sample observations
all.full.data <- list.files(paste0(output.folder, "/SPCD_standata_general_full_standardized_v3/"), pattern = "remper_correction_0.5model_1.Rdata")
Nobs <- list()

for(i in 1:length(all.full.data)){
  load(paste0(output.folder,"/SPCD_standata_general_full_standardized_v3/", all.full.data[[i]]))
  
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


# copy the data-store files
system(paste("cp -r",  "model_summary_full/",
             "data-store/data/output/"))

