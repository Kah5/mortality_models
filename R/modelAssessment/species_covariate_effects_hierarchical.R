##############################################################################
# Script to take the posterior beta estimates and generate a figure of effects
##############################################################################
library(rstan)
library(tidyverse)
library(posterior)
################################################################################
# Read in mortality data for 17 species
################################################################################
cleaned.data <- readRDS( "data/cleaned.data.mortality.TRplots.RDS")
unique(cleaned.data$SPCD)

# get the top species
nspp <- cleaned.data %>% group_by(SPCD) %>% summarise(n = n(), 
                                                      pct = n/nrow(cleaned.data)) %>% arrange (desc(`pct`))

nspp$cumulative.pct <- cumsum(nspp$pct)



# link up to the species table:
nspp$COMMON <- FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD),]$COMMON

nspp[1:17,]$COMMON



SPCD.id <- 318
model.no <- 6

  output.folder = "C:/Users/KellyHeilman/Box/01. kelly.heilman Workspace/mortality/Eastern-Mortality/mortality_models/"
  
  # Extract posterior samples for the model and species
  fit <- readRDS(paste0(output.folder,"/SPCD_stanoutput_joint_v2/u_betas_model_",model.no,"_1000samples.rds"))
  fit_ssm_df <- as_draws_df(fit) # takes awhile to convert to df
  
  # get all the covariates using posterior package
  betas.quant <- subset_draws(fit_ssm_df, variable = "u_beta") %>% summarise_draws(median, ~quantile(., probs = c(0.025, 0.975))) %>%
    rename(`ci.lo` = "2.5%", `ci.hi` = "97.5%") %>% 
    mutate(remper.cor = 0.5)
  # relabel u_betas to meaningful species ids names
  betas.quant$spp <- rep(1:17, 51)
  
  alpha.fit <- readRDS(paste0(output.folder,"/SPCD_stanoutput_joint_v2/alpha.spp_model_",model.no,"_1000samples.rds"))
  alpha_df <- as_draws_df(alpha.fit) # takes awhile to convert to df
  
  # get all the covariates using posterior package
  betas.quant <- subset_draws(fit_ssm_df, variable = "u_beta") %>% summarise_draws(median, ~quantile(., probs = c(0.025, 0.975))) %>%
    rename(`ci.lo` = "2.5%", `ci.hi` = "97.5%") %>% 
    mutate(remper.cor = 0.5)
  # relabel u_betas to meaningful species ids names
  betas.quant$spp <- rep(1:17, 51)
  betas.quant$cov <- rep(1:51, each = 17)
  load(paste0("SPCD_standata_general_full_standardized/SPCD_",SPCD.id, "remper_correction_0.5model_",model.no, ".Rdata")) # load the species code data
  mod.data$K <- ncol(mod.data$xM)
  
  covariate_names <- c(colnames(mod.data$xM))  # Replace with your covariate names
  betas.quant$Covariate <- rep(covariate_names, each = 17)
  betas.quant$Species <- rep(nspp[1:17,]$COMMON, 51)
  
  # make a figure  
  ggplot()+geom_point(data = betas.quant, aes(x = Species, y = median ))+
    geom_errorbar(data= betas.quant, aes(x = Species, ymin = ci.lo, ymax = ci.hi))+
    facet_wrap(~Covariate, scales = "free_y")+
      theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90))

  #bet0a.spp <- readRDS( paste0("SPCD_stanoutput_joint/u_betas_model_",model.no,"_1000samples.rds"))
  
  
  # reorder by the value of the covariate
  #betas.quant <- betas.quant %>% arrange(by = median) 
  # betas.quant$parameter <- factor(betas.quant$parameter, levels = betas.quant$parameter)
  
  # get overlapping zero to color the error bars
  betas.quant$`significance` <- ifelse(betas.quant$ci.lo < 0 & betas.quant$ci.hi < 0, "significant", 
                                       ifelse(betas.quant$ci.lo > 0 & betas.quant$ci.hi > 0, "significant", "not overlapping zero"))
  
 
  
  betas.quant$Covariate <- factor(betas.quant$Covariate, levels = unique(betas.quant$Covariate))
  ggplot(data = na.omit(betas.quant), aes(x = Species, y = median, color = significance))+geom_point()+
    geom_errorbar(data = na.omit(betas.quant), aes(x = Species , ymin = ci.lo, ymax = ci.hi, color = significance), width = 0.1)+
    geom_abline(aes(slope = 0, intercept = 0), color = "grey", linetype = "dashed")+
    facet_wrap(~Covariate, scales= "free_y")+
    theme_bw(base_size = 10)+
    theme( axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.grid  = element_blank(), legend.position = "none")+ylab("Effect on mortality")+xlab("Parameter")+
    scale_color_manual(values = c("not overlapping zero"="darkgrey", "significant"="black"))
  ggsave(height = 10, width = 15, units = "in",paste0(output.folder,"SPCD_stanoutput_joint_v2/images/Estimated_effects_on_mortality_model_model6_all_species_betas.png"))
  
  # get the population estimates
  
  mubetas.estimates <- readRDS( paste0(output.folder, "SPCD_stanoutput_joint_v2/beta_model_6_1000samples.RDS"))
  
  mubetas.quant <- summarise_draws(mubetas.estimates, median, ~quantile(.x, probs = c(0.025, 0.975)))%>% rename(`ci.lo` = `2.5%`, `ci.hi` = `97.5%`)
  
  mubetas.quant$remper.cor <- 0.5                                                        
  mubetas.quant$spp <- 18
  mubetas.quant$cov <- 1:51
  beta.names <- data.frame(cov = 1:51, 
                           Covariate = unique(betas.quant$Covariate))
  
  mubetas.quant <- left_join(mubetas.quant, beta.names)
  main.table <- data.frame( 
                           Species = "population", 
                           spp = 18, 
                           SPCD = 1000,
                           COMMON = 1000)
  
  mubetas.quant <- left_join(mubetas.quant, main.table)
  
  
  mubetas.quant$`significance` <- ifelse(mubetas.quant$ci.lo < 0 & mubetas.quant$ci.hi < 0, "significant", 
                                         ifelse(mubetas.quant$ci.lo > 0 & mubetas.quant$ci.hi > 0, "significant", "not overlapping zero"))
  
  mubetas.quant <- mubetas.quant %>% arrange(by = median)
  mubetas.quant$Covariate <- factor(mubetas.quant$Covariate, levels = mubetas.quant$Covariate)
  
  ggplot(data = na.omit(mubetas.quant), aes(x = Covariate, y = median, color = significance))+geom_point()+
    geom_errorbar(data = na.omit(mubetas.quant), aes(x = Covariate , ymin = ci.lo, ymax = ci.hi, color = significance), width = 0.1)+
    geom_abline(aes(slope = 0, intercept = 0), color = "grey", linetype = "dashed")+facet_wrap(~Species)+theme_bw(base_size = 10)+
    theme( axis.text.x = element_text(angle = 65, hjust = 1), panel.grid  = element_blank(), legend.position = "none")+ylab("Effect on mortality")+xlab("Covariate")+
    scale_color_manual(values = c("not overlapping zero"="darkgrey", "significant"="black"))
  ggsave(height = 5, width = 10, units = "in",paste0(output.folder, "SPCD_stanoutput_joint_v2/images/Estimated_effects_on_mortality_model_model6_all_species_population_betas.png"))
  
  
  
  ### combine the betas together
  all.joint.betas <- rbind(mubetas.quant %>% select(colnames(betas.quant)), betas.quant) 
  # reorder factors for better plotting
  all.joint.betas$Species <- factor(all.joint.betas$Species, levels = c("population", unique(betas.quant$Species)[order(unique(betas.quant$Species))]))
  all.joint.betas$Covariate <- factor(all.joint.betas$Covariate, levels = beta.names$Covariate)
  
  ggplot(data = na.omit(all.joint.betas), aes(x = Species, y = median, color = significance, shape = Species %in% "population"))+geom_point()+
    geom_errorbar(data = na.omit(all.joint.betas), aes(x = Species , ymin = ci.lo, ymax = ci.hi, color = significance), width = 0.1)+
    geom_abline(aes(slope = 0, intercept = 0), color = "grey", linetype = "dashed")+facet_wrap(~Covariate, scales = "free_y")+theme_bw(base_size = 10)+
    theme( axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.grid  = element_blank(), legend.position = "none")+ylab("Effect on survival")+xlab("Covariate")+
    scale_color_manual(values = c("not overlapping zero"="darkgrey", "significant"="black"))+
    scale_shape_manual(values = c("TRUE" = 15, "FALSE" = 19))
  
  ggsave(height = 10, width = 11, units = "in", paste0(output.folder,"SPCD_stanoutput_joint_v2/images/Estimated_betas_all_model6.png"))
  
  ggplot(data = na.omit(all.joint.betas) %>% filter(! Species %in% "population"), aes(x = Species, y = median, color = significance, shape = Species %in% "population"))+geom_point()+
    geom_errorbar(data = na.omit(all.joint.betas) %>% filter(! Species %in% "population"), aes(x = Species , ymin = ci.lo, ymax = ci.hi, color = significance), width = 0.1)+
    geom_abline(aes(slope = 0, intercept = 0), color = "grey", linetype = "dashed")+facet_wrap(~Covariate, scales = "free_y")+theme_bw(base_size = 10)+
    theme( axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.grid  = element_blank(), legend.position = "none")+ylab("Effect on survival")+xlab("Covariate")+
    scale_color_manual(values = c("not overlapping zero"="darkgrey", "significant"="black"))+
    scale_shape_manual(values = c("TRUE" = 15, "FALSE" = 19))
  
  ggsave(height = 10, width = 11, units = "in", paste0(output.folder,"SPCD_stanoutput_joint_v2/images/Estimated_betas_all_model6_no_population.png"))
  
  
  # get the main effects
  growth.diam <- c(unique(betas.quant$Covariate)[1:2])
  competition <- c(unique(betas.quant$Covariate)[3:7])
  climate <- c(unique(betas.quant$Covariate)[8:13])
  site.vars  <- c(unique(betas.quant$Covariate)[14:18])
  main.effects <- c(growth.diam, competition, climate, site.vars)
  
  ggplot(data = na.omit(all.joint.betas) %>% filter(! Species %in% "population" & Covariate %in% main.effects), aes(x = Species, y = median, color = significance, shape = Species %in% "population"))+geom_point()+
    geom_errorbar(data = na.omit(all.joint.betas) %>% filter(! Species %in% "population"& Covariate %in% main.effects), aes(x = Species , ymin = ci.lo, ymax = ci.hi, color = significance), width = 0.1)+
    geom_abline(aes(slope = 0, intercept = 0), color = "grey", linetype = "dashed")+
    facet_wrap(~Covariate, scales = "free_y", ncol = 6)+theme_bw(base_size = 10)+
    theme( axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.grid  = element_blank(), legend.position = "none")+ylab("Effect on survival")+xlab("Covariate")+
    scale_color_manual(values = c("not overlapping zero"="darkgrey", "significant"="black"))+
    scale_shape_manual(values = c("TRUE" = 15, "FALSE" = 19))
  
  ggsave(height = 6, width = 11, units = "in", paste0(output.folder,"SPCD_stanoutput_joint_v2/images/Estimated_betas_main_effects_model6_no_population.png"))
  
  # reorganize to have broader categories:

  class.definitions <- data.frame(Covariate = unique(betas.quant$Covariate), 
             Class = c("Growth, Diameter, Disturbance", "Growth, Diameter, Disturbance", 
                       "Competition", "Site Conditions", "Competition", "Competition", "Growth, Diameter, Disturbance", 
                       "Climate Normals", "Climate Normals", "Climate Normals", 
                      "Climate Anomalies", "Climate Anomalies", "Climate Anomalies", 
                      "Site Position", "Site Position", "Site Position", "Site Conditions", "Site Conditions", 
                      rep("Growth Interactions", 17),  rep("Diameter Interactions", 16)
                      ))
 all.joint.betas <-  left_join(all.joint.betas, class.definitions)
  ggplot(data = na.omit(all.joint.betas) %>% filter(! Species %in% "population" & Covariate %in% main.effects), aes(x =Covariate, y = median, group = Species, color = Species, shape = Species %in% "population"))+geom_point(position= position_dodge(width = 1))+
    #geom_errorbar(data = na.omit(all.joint.betas) %>% filter(! Species %in% "population"& Covariate %in% main.effects), aes(x = Species , ymin = ci.lo, ymax = ci.hi, group = Covariate, color = significance), width = 0.1)+
    geom_abline(aes(slope = 0, intercept = 0), color = "grey", linetype = "dashed")+
    facet_wrap(~Class, scales = "free", ncol = 6)+theme_bw(base_size = 10)+
    theme( axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.grid  = element_blank(), legend.position = "none")+ylab("Effect on survival")+xlab("Covariate")+
    #scale_color_manual(values = c("not overlapping zero"="darkgrey", "significant"="black"))+
    scale_shape_manual(values = c("TRUE" = 15, "FALSE" = 19))+coord_flip()
  
  #ggsave(height = 6, width = 11, units = "in", paste0(output.folder,"SPCD_stanoutput_joint_v2/images/Estimated_betas_main_effects_model6_no_population.png"))
  
 growth.diam <-  ggplot(data = na.omit(all.joint.betas) %>% filter(! Species %in% "population" & Class %in% "Growth, Diameter, Disturbance"), aes(y =median , x = Species, group = Covariate, color = significance, shape = Species %in% "population"))+geom_point()+
    geom_errorbar(data = na.omit(all.joint.betas) %>% filter(! Species %in% "population"& Class %in% "Growth, Diameter, Disturbance"), aes(x = Species , ymin = ci.lo, ymax = ci.hi, group = Covariate, color = significance), width = 0.1)+
    geom_abline(aes(slope = 0, intercept = 0), color = "grey", linetype = "dashed")+
    facet_grid(rows = vars(Covariate), cols = vars(Class),, scales = "free_y")+
    theme_bw(base_size = 10)+
    theme( axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.grid  = element_blank(), legend.position = "none")+
   ylab("Effect on survival")+xlab(" ")+
    scale_color_manual(values = c("not overlapping zero"="darkgrey", "significant"="black"))+
    scale_shape_manual(values = c("TRUE" = 15, "FALSE" = 19))+coord_flip()
  
  
  competition <- ggplot(data = na.omit(all.joint.betas) %>% filter(! Species %in% "population" & Class %in% "Competition"), aes(y =median , x = Species, group = Covariate, color = significance, shape = Species %in% "population"))+geom_point()+
    geom_errorbar(data = na.omit(all.joint.betas) %>% filter(! Species %in% "population"& Class %in% "Competition"), aes(x = Species , ymin = ci.lo, ymax = ci.hi, group = Covariate, color = significance), width = 0.1)+
    geom_abline(aes(slope = 0, intercept = 0), color = "grey", linetype = "dashed")+
    facet_grid(rows = vars(Covariate), cols = vars(Class),, scales = "free_y")+
    theme_bw(base_size = 10)+
    theme( axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.grid  = element_blank(), legend.position = "none")+
    ylab("Effect on survival")+xlab(" ")+
    scale_color_manual(values = c("not overlapping zero"="darkgrey", "significant"="black"))+
    scale_shape_manual(values = c("TRUE" = 15, "FALSE" = 19))+coord_flip()
  
  
  normals <- ggplot(data = na.omit(all.joint.betas) %>% filter(! Species %in% "population" & Class %in% "Climate Normals"), aes(y =median , x = Species, group = Covariate, color = significance, shape = Species %in% "population"))+geom_point()+
    geom_errorbar(data = na.omit(all.joint.betas) %>% filter(! Species %in% "population"& Class %in% "Climate Normals"), aes(x = Species , ymin = ci.lo, ymax = ci.hi, group = Covariate, color = significance), width = 0.1)+
    geom_abline(aes(slope = 0, intercept = 0), color = "grey", linetype = "dashed")+
    facet_grid(rows = vars(Covariate), cols = vars(Class),, scales = "free_y")+
    theme_bw(base_size = 10)+
    theme( axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.grid  = element_blank(), legend.position = "none")+
    ylab("Effect on survival")+xlab(" ")+
    scale_color_manual(values = c("not overlapping zero"="darkgrey", "significant"="black"))+
    scale_shape_manual(values = c("TRUE" = 15, "FALSE" = 19))+coord_flip()
  
  
  anomalies <- ggplot(data = na.omit(all.joint.betas) %>% filter(! Species %in% "population" & Class %in% "Climate Anomalies"), aes(y =median , x = Species, group = Covariate, color = significance, shape = Species %in% "population"))+geom_point()+
    geom_errorbar(data = na.omit(all.joint.betas) %>% filter(! Species %in% "population"& Class %in% "Climate Anomalies"), aes(x = Species , ymin = ci.lo, ymax = ci.hi, group = Covariate, color = significance), width = 0.1)+
    geom_abline(aes(slope = 0, intercept = 0), color = "grey", linetype = "dashed")+
    facet_grid(rows = vars(Covariate), cols = vars(Class),, scales = "free_y")+
    theme_bw(base_size = 10)+
    theme( axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.grid  = element_blank(), legend.position = "none")+
    ylab("Effect on survival")+xlab(" ")+
    scale_color_manual(values = c("not overlapping zero"="darkgrey", "significant"="black"))+
    scale_shape_manual(values = c("TRUE" = 15, "FALSE" = 19))+coord_flip()
  
  positions <- ggplot(data = na.omit(all.joint.betas) %>% filter(! Species %in% "population" & Class %in% "Site Position"), aes(y =median , x = Species, group = Covariate, color = significance, shape = Species %in% "population"))+geom_point()+
    geom_errorbar(data = na.omit(all.joint.betas) %>% filter(! Species %in% "population"& Class %in% "Site Position"), aes(x = Species , ymin = ci.lo, ymax = ci.hi, group = Covariate, color = significance), width = 0.1)+
    geom_abline(aes(slope = 0, intercept = 0), color = "grey", linetype = "dashed")+
    facet_grid(rows = vars(Covariate), cols = vars(Class),, scales = "free_y")+
    theme_bw(base_size = 10)+
    theme( axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.grid  = element_blank(), legend.position = "none")+
    ylab("Effect on survival")+xlab(" ")+
    scale_color_manual(values = c("not overlapping zero"="darkgrey", "significant"="black"))+
    scale_shape_manual(values = c("TRUE" = 15, "FALSE" = 19))+coord_flip()
  
  site.cond <- ggplot(data = na.omit(all.joint.betas) %>% filter(! Species %in% "population" & Class %in% "Site Conditions"), aes(y =median , x = Species, group = Covariate, color = significance, shape = Species %in% "population"))+geom_point()+
    geom_errorbar(data = na.omit(all.joint.betas) %>% filter(! Species %in% "population"& Class %in% "Site Conditions"), aes(x = Species , ymin = ci.lo, ymax = ci.hi, group = Covariate, color = significance), width = 0.1)+
    geom_abline(aes(slope = 0, intercept = 0), color = "grey", linetype = "dashed")+
    facet_grid(rows = vars(Covariate), cols = vars(Class), scales = "free_y")+
    theme_bw(base_size = 10)+
    theme( axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
           panel.grid  = element_blank(), legend.position = "none")+ylab("Effect on survival")+xlab("")+
    scale_color_manual(values = c("not overlapping zero"="darkgrey", "significant"="black"))+
    scale_shape_manual(values = c("TRUE" = 15, "FALSE" = 19))+coord_flip()
  
 png(height = 12, width = 12, units = "in",res = 300,  paste0(output.folder, "SPCD_stanoutput_joint_v2/images/full_main_effects_summary.png")) 
 cowplot::plot_grid(growth.diam, normals, anomalies, 
           competition, site.cond, positions, align = "hv")
dev.off()

# do the same but reorder in terms of species type/forest type
betas.quant$E_SPGRPCD <- ref_species[match(betas.quant$Species, ref_species$COMMON_NAME), ]$E_SPGRPCD
betas.quant$MAJOR_SPGRPCD <- ref_species[match(betas.quant$Species, ref_species$COMMON_NAME), ]$MAJOR_SPGRPCD
betas.quant$GENUS <- ref_species[match(betas.quant$Species, ref_species$COMMON_NAME), ]$GENUS

#betas.quant$Jenkins_SPGRPCD <- ref_species[match(betas.quant$Species, ref_codes$COMMON_NAME), ]$JENKINS_SPGRPCD
unique(betas.quant[, c("MAJOR_SPGRPCD", "Species")])
SPECIES.groups <- unique(betas.quant[,c("Species", "GENUS", "MAJOR_SPGRPCD")]) %>% arrange(MAJOR_SPGRPCD)
SPECIES.groups$Species

all.joint.betas$Species <- factor(all.joint.betas$Species, levels = c(
  "balsam fir", "red spruce", "northern white-cedar", "eastern hemlock", "eastern white pine", 
  "yellow birch", "American beech", "sugar maple", "red maple",
  "northern red oak","chestnut oak", "black oak", "white oak", 
  "hickory spp.","white ash", "yellow-poplar"))



growth.diam <-  ggplot(data = na.omit(all.joint.betas) %>% filter(! Species %in% "population" & Class %in% "Growth, Diameter, Disturbance"), aes(y =median , x = Species, group = Covariate, color = significance, shape = Species %in% "population"))+geom_point()+
  geom_errorbar(data = na.omit(all.joint.betas) %>% filter(! Species %in% "population"& Class %in% "Growth, Diameter, Disturbance"), aes(x = Species , ymin = ci.lo, ymax = ci.hi, group = Covariate, color = significance), width = 0.1)+
  geom_abline(aes(slope = 0, intercept = 0), color = "grey", linetype = "dashed")+
  facet_grid(rows = vars(Covariate), cols = vars(Class),, scales = "free_y")+
  theme_bw(base_size = 10)+
  theme( axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.grid  = element_blank(), legend.position = "none")+
  ylab("Effect on survival")+xlab(" ")+
  scale_color_manual(values = c("not overlapping zero"="darkgrey", "significant"="black"))+
  scale_shape_manual(values = c("TRUE" = 15, "FALSE" = 19))+coord_flip()


competition <- ggplot(data = na.omit(all.joint.betas) %>% filter(! Species %in% "population" & Class %in% "Competition"), aes(y =median , x = Species, group = Covariate, color = significance, shape = Species %in% "population"))+geom_point()+
  geom_errorbar(data = na.omit(all.joint.betas) %>% filter(! Species %in% "population"& Class %in% "Competition"), aes(x = Species , ymin = ci.lo, ymax = ci.hi, group = Covariate, color = significance), width = 0.1)+
  geom_abline(aes(slope = 0, intercept = 0), color = "grey", linetype = "dashed")+
  facet_grid(rows = vars(Covariate), cols = vars(Class),, scales = "free_y")+
  theme_bw(base_size = 10)+
  theme( axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.grid  = element_blank(), legend.position = "none")+
  ylab("Effect on survival")+xlab(" ")+
  scale_color_manual(values = c("not overlapping zero"="darkgrey", "significant"="black"))+
  scale_shape_manual(values = c("TRUE" = 15, "FALSE" = 19))+coord_flip()


normals <- ggplot(data = na.omit(all.joint.betas) %>% filter(! Species %in% "population" & Class %in% "Climate Normals"), aes(y =median , x = Species, group = Covariate, color = significance, shape = Species %in% "population"))+geom_point()+
  geom_errorbar(data = na.omit(all.joint.betas) %>% filter(! Species %in% "population"& Class %in% "Climate Normals"), aes(x = Species , ymin = ci.lo, ymax = ci.hi, group = Covariate, color = significance), width = 0.1)+
  geom_abline(aes(slope = 0, intercept = 0), color = "grey", linetype = "dashed")+
  facet_grid(rows = vars(Covariate), cols = vars(Class),, scales = "free_y")+
  theme_bw(base_size = 10)+
  theme( axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.grid  = element_blank(), legend.position = "none")+
  ylab("Effect on survival")+xlab(" ")+
  scale_color_manual(values = c("not overlapping zero"="darkgrey", "significant"="black"))+
  scale_shape_manual(values = c("TRUE" = 15, "FALSE" = 19))+coord_flip()


anomalies <- ggplot(data = na.omit(all.joint.betas) %>% filter(! Species %in% "population" & Class %in% "Climate Anomalies"), aes(y =median , x = Species, group = Covariate, color = significance, shape = Species %in% "population"))+geom_point()+
  geom_errorbar(data = na.omit(all.joint.betas) %>% filter(! Species %in% "population"& Class %in% "Climate Anomalies"), aes(x = Species , ymin = ci.lo, ymax = ci.hi, group = Covariate, color = significance), width = 0.1)+
  geom_abline(aes(slope = 0, intercept = 0), color = "grey", linetype = "dashed")+
  facet_grid(rows = vars(Covariate), cols = vars(Class),, scales = "free_y")+
  theme_bw(base_size = 10)+
  theme( axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.grid  = element_blank(), legend.position = "none")+
  ylab("Effect on survival")+xlab(" ")+
  scale_color_manual(values = c("not overlapping zero"="darkgrey", "significant"="black"))+
  scale_shape_manual(values = c("TRUE" = 15, "FALSE" = 19))+coord_flip()

positions <- ggplot(data = na.omit(all.joint.betas) %>% filter(! Species %in% "population" & Class %in% "Site Position"), aes(y =median , x = Species, group = Covariate, color = significance, shape = Species %in% "population"))+geom_point()+
  geom_errorbar(data = na.omit(all.joint.betas) %>% filter(! Species %in% "population"& Class %in% "Site Position"), aes(x = Species , ymin = ci.lo, ymax = ci.hi, group = Covariate, color = significance), width = 0.1)+
  geom_abline(aes(slope = 0, intercept = 0), color = "grey", linetype = "dashed")+
  facet_grid(rows = vars(Covariate), cols = vars(Class),, scales = "free_y")+
  theme_bw(base_size = 10)+
  theme( axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.grid  = element_blank(), legend.position = "none")+
  ylab("Effect on survival")+xlab(" ")+
  scale_color_manual(values = c("not overlapping zero"="darkgrey", "significant"="black"))+
  scale_shape_manual(values = c("TRUE" = 15, "FALSE" = 19))+coord_flip()

site.cond <- ggplot(data = na.omit(all.joint.betas) %>% filter(! Species %in% "population" & Class %in% "Site Conditions"), aes(y =median , x = Species, group = Covariate, color = significance, shape = Species %in% "population"))+geom_point()+
  geom_errorbar(data = na.omit(all.joint.betas) %>% filter(! Species %in% "population"& Class %in% "Site Conditions"), aes(x = Species , ymin = ci.lo, ymax = ci.hi, group = Covariate, color = significance), width = 0.1)+
  geom_abline(aes(slope = 0, intercept = 0), color = "grey", linetype = "dashed")+
  facet_grid(rows = vars(Covariate), cols = vars(Class), scales = "free_y")+
  theme_bw(base_size = 10)+
  theme( axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
         panel.grid  = element_blank(), legend.position = "none")+ylab("Effect on survival")+xlab("")+
  scale_color_manual(values = c("not overlapping zero"="darkgrey", "significant"="black"))+
  scale_shape_manual(values = c("TRUE" = 15, "FALSE" = 19))+coord_flip()

png(height = 12, width = 12, units = "in",res = 300,  paste0(output.folder, "SPCD_stanoutput_joint_v2/images/full_main_effects_summary_reorganized.png")) 
cowplot::plot_grid(growth.diam, normals, anomalies, 
                   competition, site.cond, positions, align = "hv")
dev.off()

# get the species alpha estimates
  alphas.estimates <- readRDS( paste0(output.folder,"SPCD_stanoutput_joint_v2/alpha.spp_model_",model.no,"_1000samples.rds"))
  
  alphas.quant <- summarise_draws(alphas.estimates, median, ~quantile(.x, probs = c(0.025, 0.975)))%>% rename(`ci.lo` = `2.5%`, `ci.hi` = `97.5%`)
  
  # alphas.estimates <- fit_ssm_df %>% dplyr::select(paste0("alpha_SPP[",1:17,"]")) 
  # mub.m <- reshape2::melt(alphas.estimates)
  # 
  # 
  # alphas.quant <- mub.m %>% group_by(variable) %>% summarise(median = quantile(value, 0.5, na.rm =TRUE),
  #                                                            ci.lo = quantile(value, 0.005, na.rm =TRUE),
  #                                                            ci.hi = quantile(value, 0.975, na.rm =TRUE))
  alphas.quant$spp <- 1:17
  spp.table <- data.frame(SPCD = nspp[1:17,]$SPCD, 
                          COMMON = nspp[1:17,]$COMMON, 
                          spp = 1:17)
  
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
  ggsave(height = 5, width = 10, units = "in",paste0(output.folder,"SPCD_stanoutput_joint_v2/images/Estimated_alpha_SPP_model6_all_species_alphas.png"))
  
  # combine the population and the species level betas
  
  ## get population level alpha estimates:
  alphas.pop.estimates <- readRDS( paste0(output.folder,"SPCD_stanoutput_joint_v2/alpha.p_model_",model.no,"_1000samples.rds"))
  
  alphas.pop.quant <- summarise_draws(alphas.pop.estimates, median, ~quantile(.x, probs = c(0.025, 0.975)))%>% rename(`ci.lo` = `2.5%`, `ci.hi` = `97.5%`)
  
  # alphas.estimates <- fit_ssm_df %>% dplyr::select(paste0("alpha_SPP[",1:17,"]")) 
  # mub.m <- reshape2::melt(alphas.estimates)
  # 
  # 
  # alphas.quant <- mub.m %>% group_by(variable) %>% summarise(median = quantile(value, 0.5, na.rm =TRUE),
  #                                                            ci.lo = quantile(value, 0.005, na.rm =TRUE),
  #                                                            ci.hi = quantile(value, 0.975, na.rm =TRUE))
  alphas.pop.quant$spp <- 18
  
  # main.table$COMMON = "population"
  # main.table$SPCD = 1000
  alphas.pop.quant <- left_join(alphas.pop.quant, main.table %>% select(-Species))
  
  
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
  
  ggsave(height = 5, width = 6, units = "in",paste0(output.folder,"SPCD_stanoutput_joint_v2/images/Estimated_alpha_all_model6.png"))
  
  ## combine alphas and betas together to make the figure for the paper
  # 
  # all.joint.alphas <- all.joint.alphas %>% mutate(parameter = "alpha") %>% select(variable, median, ci.lo, ci.hi, spp, parameter, SPCD.id, COMMON, significance)
  # all.joint.betas <- all.joint.betas #%>% select(-param.no)
  # all.params <- rbind(all.joint.betas, all.joint.alphas)
  # 
  # ggplot(data = na.omit(all.params), aes(x = COMMON, y = median, color = significance, shape = COMMON %in% "population"))+geom_point()+
  #   geom_errorbar(data = na.omit(all.params), aes(x = COMMON , ymin = ci.lo, ymax = ci.hi, color = significance), width = 0.1)+
  #   geom_abline(aes(slope = 0, intercept = 0), color = "grey", linetype = "dashed")+facet_wrap(~parameter)+theme_bw(base_size = 10)+
  #   theme( axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.grid  = element_blank(), legend.position = "none")+ylab("Effect on survival")+xlab("Parameter")+
  #   scale_color_manual(values = c("not overlapping zero"="darkgrey", "significant"="black"))+
  #   scale_shape_manual(values = c("TRUE" = 15, "FALSE" = 19))
  # 
  # ggsave(height = 10, width = 11, units = "in",dpi = 350,paste0("SPCD_stanoutput_joint/images/Estimated_parameters_model6_joint.png"))
  # 
  #########################################################################################
  # Estimate conditional effects for each species and covariate
  #########################################################################################
  
  
  
  for(i in 17:1){
    cat(i)
    
    SPCD.id <- nspp[i,]$SPCD
    beta.species.names <- betas.quant %>% filter(spp %in% i) %>% select(variable)
    # select only the betas and intercepts for the species of interest:
    beta <- subset_draws(fit_ssm_df, variable = beta.species.names$variable) %>% select(-.chain, -.iteration, -.draw) 
    beta_0 <- data.frame(subset_draws(alpha_df, variable = paste0("alpha_SPP[",i,"]"))) %>% select(-.chain, -.iteration, -.draw)  # Intercept
  # read in the species data and covariates to get the min and max and get ranges
  load(paste0("SPCD_standata_general_full_standardized/SPCD_",SPCD.id, "remper_correction_0.5model_",model.no, ".Rdata")) # load the species code data
  mod.data$K <- ncol(mod.data$xM)
  
  covariate_names <- c(colnames(mod.data$xM))  # Replace with your covariate names
  
  var.mins <- as.vector(apply(data.frame(mod.data$xM),2 , function(x)quantile(x, 0.025)))
  var.maxes <- as.vector(apply(data.frame(mod.data$xM),2 , function(x) quantile(x, 0.975)))
  
  covariate_ranges_df <- data.frame(
    covariate  = covariate_names, 
    mins = var.mins, 
    maxes = var.maxes
  )
  
  covariate_ranges_df[1,]
  cov.list <- list()
  for(j in 1:length (covariate_ranges_df[,2])){
    cov.list[[j]] <- seq(covariate_ranges_df[j,2], covariate_ranges_df[j,3], length.out = 25)
  }
  names(cov.list) <- c(colnames(mod.data$xM))
  
  # Define a baseline for other covariates (mean or median values)
  #var.medians <- as.vector(apply(data.frame(mod.data$xM),2 , median))
  
  covariate_ranges_df$medians <- 0#as.vector(apply(data.frame(mod.data$xM),2 , median))
  
  
  # set up a function calculate probabilities
  inv_logit_fxn <- function(x) {
    1 / (1 + exp(-x))
  }
  cov_name <- "annual.growth.scaled"
  
  
  # get all the data and make prediction over the range of each individual covariate
  plot_data <- map_dfr(covariate_names, function(cov_name) {
    # Range of the focal covariate
    cov_range <- cov.list[[cov_name]]
    
    # Create a matrix for covariates
    covariate_matrix <- matrix(covariate_ranges_df$medians, nrow = length(cov_range), ncol = length(covariate_names), byrow = TRUE)
    colnames(covariate_matrix) <- covariate_names
    covariate_matrix[, cov_name] <- cov_range  # Vary only the focal covariate
    
    # Calculate probabilities for each sample
    #probabilities <- apply(beta, MARGIN = 1 , function(b) {
    probabilities <- list()
    for(j in 1:length(covariate_matrix[,1])){
      linear_predictor <- beta_0[,1] + rowSums( as.matrix(beta)%*% covariate_matrix[j,] )
      probabilities[[j]] <-  as.vector(inv_logit_fxn(linear_predictor))
    }
    
    
    
    # Summarize probabilities (mean and 90% credible interval)
    prob_summary <- lapply(probabilities, function(p) {
      c(mean = median(p), ci.lo = quantile(p, 0.1), ci.hi = quantile(p, 0.9))
    })
    
    prob_summary.df <- do.call(rbind, prob_summary)
    
    # Combine with covariate range for plotting
    data.frame(
      Covariate = cov_name,
      Value = cov_range,
      mean = prob_summary.df[, "mean"],
      ci.lo = prob_summary.df[, "ci.lo.10%"],
      ci.hi = prob_summary.df[, "ci.hi.90%"]
    )
  })
  
  # make a plot of covariate responses for each species
  ggplot(plot_data, aes(x = Value, y = 1-mean, color = Covariate)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Covariate), alpha = 0.2, color = NA) +
    facet_wrap(~ Covariate, scales = "free") +
    labs(
      x = "Covariate Value",
      y = "Annual Probability of Mortality",
      title = "Effect of Covariates on Probability of Mortality"
    ) +
    theme_minimal() +
    theme(legend.position = "none")
  ggsave(height = 10, width = 12, paste0(output.folder, "images/predicted_mortality/model_6_all_marginal_SPCD_", SPCD.id, "_joint.png"))
  
  # save for this species:
  plot_data$SPCD <- SPCD.id
  saveRDS(plot_data, paste0(output.folder, "SPCD_stanoutput_joint_v2/predicted_mort/model_6_all_marginal_SPCD_", SPCD.id, "_joint.RDS") )
}


################################################################################
# combine all of the species conditional responses together
################################################################################
marginal_response <- list()
for (i in 1:17){
  SPCD.id <- nspp[i,]$SPCD
  marginal_response[[i]] <- readRDS( paste0(output.folder, "SPCD_stanoutput_joint_v2/predicted_mort/model_6_all_marginal_SPCD_", SPCD.id, "_joint.RDS") )
  
}
marginal_response_df <- do.call(rbind, marginal_response)

# get species common names
marginal_response_df$Species <- FIESTA::ref_species[match(marginal_response_df$SPCD, FIESTA::ref_species$SPCD),]$COMMON

ggplot(marginal_response_df, aes(x = Value, y = 1-mean, color = Species)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
  facet_wrap(~ Covariate, scales = "free") +
  labs(
    x = "Covariate Value",
    y = "Annual Probability of Mortality",
    title = "Effect of Covariates on Probability of Mortality"
  ) +
  theme_bw(base_size = 12)+theme(panel.grid = element_blank()) 
ggsave(height = 12, width = 12, paste0(output.folder, "images/predicted_mortality/model_6_species_level_models_all_species_marginal_effects_joint.png"))


## split up by variable type:

marginal_response_df <- left_join(marginal_response_df, class.definitions)
# ## main effects:
main.effects <- c("Growth, Diameter, Disturbance", "Competition" ,"Site Conditions",              
 "Climate Normals","Climate Anomalies", "Site Position")
# generate figures for all of these separately 
growth.marginal <- ggplot(marginal_response_df %>% filter( Class %in% "Growth, Diameter, Disturbance"), aes(x = Value, y = 1-mean, color = Species)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
  facet_grid(cols = vars(Class), rows = vars(Covariate), scales = "free") +
  labs(
    x = "Covariate Value",
    y = "Annual Probability of Mortality",
    #title = "Effect of Covariates on Probability of Mortality"
  ) +
  theme_bw(base_size = 14)+theme(panel.grid = element_blank())

Competition.marginal <- ggplot(marginal_response_df %>% filter( Class %in% "Competition"), aes(x = Value, y = 1-mean, color = Species)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
  facet_grid(cols = vars(Class), rows = vars(Covariate), scales = "free") +
  labs(
    x = "Covariate Value",
    y = "Annual Probability of Mortality",
    #title = "Effect of Covariates on Probability of Mortality"
  ) +
  theme_bw(base_size = 14)+theme(panel.grid = element_blank())

site.cond.marginal <- ggplot(marginal_response_df %>% filter( Class %in% "Site Conditions"), aes(x = Value, y = 1-mean, color = Species)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
  facet_grid(cols = vars(Class), rows = vars(Covariate), scales = "free") +
  labs(
    x = "Covariate Value",
    y = "Annual Probability of Mortality",
    #title = "Effect of Covariates on Probability of Mortality"
  ) +
  theme_bw(base_size = 14)+theme(panel.grid = element_blank())

normals.marginal <- ggplot(marginal_response_df %>% filter( Class %in% "Climate Normals"), aes(x = Value, y = 1-mean, color = Species)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
  facet_grid(cols = vars(Class), rows = vars(Covariate), scales = "free") +
  labs(
    x = "Covariate Value",
    y = "Annual Probability of Mortality",
    #title = "Effect of Covariates on Probability of Mortality"
  ) +
  theme_bw(base_size = 14)+theme(panel.grid = element_blank())


anomalies.marginal <- ggplot(marginal_response_df %>% filter( Class %in% "Climate Anomalies"), aes(x = Value, y = 1-mean, color = Species)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
  facet_grid(cols = vars(Class), rows = vars(Covariate), scales = "free") +
  labs(
    x = "Covariate Value",
    y = "Annual Probability of Mortality",
    #title = "Effect of Covariates on Probability of Mortality"
  ) +
  theme_bw(base_size = 14)+theme(panel.grid = element_blank())

anomalies.marginal <- ggplot(marginal_response_df %>% filter( Class %in% "Climate Anomalies"), aes(x = Value, y = 1-mean, color = Species)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
  facet_grid(cols = vars(Class), rows = vars(Covariate), scales = "free") +
  labs(
    x = "Covariate Value",
    y = "Annual Probability of Mortality",
    #title = "Effect of Covariates on Probability of Mortality"
  ) +
  theme_bw(base_size = 14)+theme(panel.grid = element_blank())

site.position.marginal <- ggplot(marginal_response_df %>% filter( Class %in% "Site Position"), aes(x = Value, y = 1-mean, color = Species)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
  facet_grid(cols = vars(Class), rows = vars(Covariate), scales = "free") +
  labs(
    x = "Covariate Value",
    y = "Annual Probability of Mortality",
    #title = "Effect of Covariates on Probability of Mortality"
  ) +
  theme_bw(base_size = 14)+theme(panel.grid = element_blank())

species.legend <- cowplot::get_legend(growth.marginal)

cowplot::plot_grid(cowplot::plot_grid(growth.marginal + theme(legend.position = "none"), 
          Competition.marginal + theme(legend.position = "none"), 
          site.cond.marginal + theme(legend.position = "none"), 
          normals.marginal + theme(legend.position = "none"), 
          anomalies.marginal + theme(legend.position = "none"), 
          site.position.marginal + theme(legend.position = "none"), align = "hv"), species.legend, rel_widths = c(1,0.2))

ggsave(height = 12, width = 12, 
       paste0(output.folder, "images/predicted_mortality/model_6_species_level_model_marginal_main_effects_joint.png"))

# # balsam fir really domninates the graph so take it out to view
# ggplot(marginal_response_df %>% filter( Covariate %in% c(growth.diam, competition, climate, site.vars  ) & !Species %in% "balsam fir"), aes(x = Value, y = 1-mean, color = Species)) +
#   geom_line(size = 1) +
#   geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
#   facet_wrap(~ Covariate, scales = "free") +
#   labs(
#     x = "Covariate Value",
#     y = "Annual Probability of Mortality",
#     title = "Effect of Covariates on Probability of Mortality"
#   ) +
#   theme_bw(base_size = 12)+theme(panel.grid = element_blank())
# ggsave(height = 8, width = 12, 
#        paste0(output.folder, "images/predicted_mortality/model_6_species_level_model_marginal_main_effects_no_balsam_fir_joint.png"))
# 
# 
# ggplot(marginal_response_df %>% filter( Covariate %in% growth.diam ), aes(x = Value, y = 1-mean, color = Species)) +
#   geom_line(size = 1) +
#   geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
#   facet_wrap(~ Covariate, scales = "free") +
#   labs(
#     x = "Covariate Value",
#     y = "Annual Probability of Mortality",
#     title = "Effect of Covariates on Probability of Mortality"
#   ) +
#   theme_bw(base_size = 12)+theme(panel.grid = element_blank())
# ggsave(height = 5, width = 12, 
#        paste0(output.folder, "images/predicted_mortality/model_6_species_level_model_marginal_growth_diam_joint.png"))
# 
# ggplot(marginal_response_df %>% filter( Covariate %in% competition ), aes(x = Value, y = 1-mean, color = Species)) +
#   geom_line(size = 1) +
#   geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
#   facet_wrap(~ Covariate, scales = "free", ncol = 5) +
#   labs(
#     x = "Covariate Value",
#     y = "Annual Probability of Mortality",
#     title = "Effect of Covariates on Probability of Mortality"
#   ) +
#   theme_bw(base_size = 12)+theme(panel.grid = element_blank())
# ggsave(height = 5, width = 12, 
#        paste0(output.folder, "images/predicted_mortality/model_6_species_level_model_marginal_competition_joint.png"))
# 
# 
# ggplot(marginal_response_df %>% filter( Covariate %in% climate  ), aes(x = Value, y = 1-mean, color = Species)) +
#   geom_line(size = 1) +
#   geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
#   facet_wrap(~ Covariate, scales = "free", ncol = 3) +
#   labs(
#     x = "Covariate Value",
#     y = "Annual Probability of Mortality",
#     title = "Effect of Covariates on Probability of Mortality"
#   ) +
#   theme_bw(base_size = 12)+theme(panel.grid = element_blank())
# ggsave(height = 8, width = 12, 
#        paste0(output.folder, "images/predicted_mortality/model_6_species_level_model_marginal_climate_vars_joint.png"))
# 
# 
# ggplot(marginal_response_df %>% filter( Covariate %in% site.vars ), aes(x = Value, y = 1-mean, color = Species)) +
#   geom_line(size = 1) +
#   geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
#   facet_wrap(~ Covariate, scales = "free", ncol = 5) +
#   labs(
#     x = "Covariate Value",
#     y = "Annual Probability of Mortality",
#     title = "Effect of Covariates on Probability of Mortality"
#   ) +
#   theme_bw(base_size = 12)+theme(panel.grid = element_blank())
# ggsave(height = 5, width = 15, 
#        paste0(output.folder, "images/predicted_mortality/model_6_species_level_model_marginal_site_vars_joint.png"))
# 
# ### interaction terms
# ggplot(marginal_response_df %>% filter( Covariate %in% growth.competition.int & !Species %in% "balsam fir"), aes(x = Value, y = 1-mean, color = Species)) +
#   geom_line(size = 1) +
#   geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
#   facet_wrap(~ Covariate, scales = "free", ncol = 5) +
#   labs(
#     x = "Covariate Value",
#     y = "Annual Probability of Mortality",
#     title = "Effect of Covariates on Probability of Mortality"
#   ) +
#   theme_bw(base_size = 12)+theme(panel.grid = element_blank())
# ggsave(height = 5, width = 12, 
#        paste0(output.folder, "images/predicted_mortality/model_6_species_level_model_marginal_growth_competition_joint_noBF.png"))
# 
# 
# ggplot(marginal_response_df %>% filter( Covariate %in% growth.climate.int& !Species %in% "balsam fir" ), aes(x = Value, y = 1-mean, color = Species)) +
#   geom_line(size = 1) +
#   geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
#   facet_wrap(~ Covariate, scales = "free", ncol = 3) +
#   labs(
#     x = "Covariate Value",
#     y = "Annual Probability of Mortality",
#     title = "Effect of Covariates on Probability of Mortality"
#   ) +
#   theme_bw(base_size = 12)+theme(panel.grid = element_blank())
# ggsave(height = 8, width = 12, 
#        paste0(output.folder, "images/predicted_mortality/model_6_species_level_model_marginal_growth_cliamte_joint_noBF.png"))
# 
# ggplot(marginal_response_df %>% filter( Covariate %in% growth.site.int& !Species %in% "balsam fir" ), aes(x = Value, y = 1-mean, color = Species)) +
#   geom_line(size = 1) +
#   geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
#   facet_wrap(~ Covariate, scales = "free", ncol = 5) +
#   labs(
#     x = "Covariate Value",
#     y = "Annual Probability of Mortality",
#     title = "Effect of Covariates on Probability of Mortality"
#   ) +
#   theme_bw(base_size = 12)+theme(panel.grid = element_blank())
# ggsave(height = 5, width = 12, 
#        paste0(output.folder, "images/predicted_mortality/model_6_species_level_model_marginal_growth_site_joint_noBF.png"))
# 
# 
# # diameter interaction terms
# ggplot(marginal_response_df %>% filter( Covariate %in% diameter.competition.int & !Species %in% "balsam fir"), aes(x = Value, y = 1-mean, color = Species)) +
#   geom_line(size = 1) +
#   geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
#   facet_wrap(~ Covariate, scales = "free", ncol = 5) +
#   labs(
#     x = "Covariate Value",
#     y = "Annual Probability of Mortality",
#     title = "Effect of Covariates on Probability of Mortality"
#   ) +
#   theme_bw(base_size = 12)+theme(panel.grid = element_blank())
# ggsave(height = 5, width = 12, 
#        paste0(output.folder, "images/predicted_mortality/model_6_species_level_model_marginal_diameter_competition_joint_noBF.png"))
# 
# 
# ggplot(marginal_response_df %>% filter( Covariate %in% diameter.climate.int & !Species %in% "balsam fir"), aes(x = Value, y = 1-mean, color = Species)) +
#   geom_line(size = 1) +
#   geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
#   facet_wrap(~ Covariate, scales = "free", ncol = 3) +
#   labs(
#     x = "Covariate Value",
#     y = "Annual Probability of Mortality",
#     title = "Effect of Covariates on Probability of Mortality"
#   ) +
#   theme_bw(base_size = 12)+theme(panel.grid = element_blank())
# ggsave(height = 8, width = 12, 
#        paste0(output.folder, "images/predicted_mortality/model_6_species_level_model_marginal_diameter_climate_joint_noBF.png"))
# 
# ggplot(marginal_response_df %>% filter( Covariate %in% diameter.site.int & !Species %in% "balsam fir" ), aes(x = Value, y = 1-mean, color = Species)) +
#   geom_line(size = 1) +
#   geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA, color = NA) +
#   facet_wrap(~ Covariate, scales = "free", ncol = 5) +
#   labs(
#     x = "Covariate Value",
#     y = "Annual Probability of Mortality",
#     title = "Effect of Covariates on Probability of Mortality"
#   ) +
#   theme_bw(base_size = 12)+theme(panel.grid = element_blank())
# ggsave(height = 5, width = 12, 
#        paste0(output.folder, "images/predicted_mortality/model_6_species_level_model_marginal_diameter_site_joint_noBF.png"))


