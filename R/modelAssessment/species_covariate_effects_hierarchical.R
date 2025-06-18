##############################################################################
# Script to take the posterior beta estimates and generate a figure of effects
##############################################################################
library(rstan)
library(tidyverse)
library(posterior)
################################################################################
# Read in mortality data for 17 species
################################################################################
# cleaned.data <- readRDS( "data/cleaned.data.mortality.TRplots.RDS")
# unique(cleaned.data$SPCD)
# 
# # get the top species
# nspp <- cleaned.data %>% group_by(SPCD) %>% summarise(n = n(), 
#                                                       pct = n/nrow(cleaned.data)) %>% arrange (desc(`pct`))
# 
# nspp$cumulative.pct <- cumsum(nspp$pct)
# 
# 
# # get the complete species list
nspp <- data.frame(SPCD = c(316, 318, 833, 832, 261, 531, 802, 129, 762,  12, 541,  97, 621, 400, 371, 241, 375))
nspp$Species <- paste(FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD),]$GENUS, FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD),]$SPECIES)
# link up to the species table:
nspp$COMMON <- FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD),]$COMMON


# read in the test data for all the species
spp.table <- data.frame(SPCD.id = nspp[1:17,]$SPCD, 
                        spp = 1:17, 
                        COMMON = nspp[1:17,]$COMMON)
model.no <- 6

nspp[1:17,]$COMMON



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

# set up custom lines for shade tolerances
sppLinetype <- c("solid", "dashed", "dotted")
names(sppLinetype) <- unique(SP.TRAITS$`Shade Tolerance`)

species_linetype <- scale_linetype_manual(values = sppLinetype)
SP.TRAITS.shade <- SP.TRAITS %>% select(`Shade Tolerance`, COMMON_NAME) %>% 
  mutate(line.t = ifelse(`Shade Tolerance`%in% "High", "solid", 
                         ifelse(`Shade Tolerance` %in% "Moderate", "dashed", "dotted")))
sppLinetype.species <- c(SP.TRAITS.shade$line.t)
names(sppLinetype.species) <- SP.TRAITS$COMMON_NAME
named_species_linetype <- scale_linetype_manual(values = sppLinetype.species)

SPCD.id <- 318
model.no <- 6
input.folder = "/home/rstudio/data-store/data/output/"
input.folder = "/home/rstudio/data-store/data/iplant/home/kellyheilman/analyses/Mortality_hierarchical_models-2025-04-07-18-56-47.7/"
output.folder = "/home/rstudio/SPCD_stanoutput_joint_v3/"

m <- 6
#for(m in 1:9){

cat(paste0("generating figures for model ",m ))
model.no <- m



# get the number of variables
mod.data <- readRDS(paste0(input.folder, "SPCD_stanoutput_joint_v3_model_", model.no, "/all_SPCD_model_",model.no,".RDS"))
ncovar <- length(colnames(mod.data$xM))

full.model <- data.frame(Covariates = colnames(mod.data$xM), 
                         id = 1:length(colnames(mod.data$xM)))
# Extract posterior samples for the model and species
# fit <- readRDS(paste0(output.folder,"/SPCD_stanoutput_joint_v2/u_betas_model_",model.no,"_1000samples.rds"))
# fit_ssm_df <- as_draws_df(fit) # takes awhile to convert to df
# 
# # get all the covariates using posterior package
# betas.quant <- subset_draws(fit_ssm_df, variable = "u_beta") %>% summarise_draws(median, ~quantile(., probs = c(0.025, 0.975))) %>%
#   rename(`ci.lo` = "2.5%", `ci.hi` = "97.5%") %>% 
#   mutate(remper.cor = 0.5)
# # relabel u_betas to meaningful species ids names
# betas.quant$spp <- rep(1:17, ncovar)
alpha.fit <- readRDS(paste0(input.folder, "SPCD_stanoutput_joint_v3_model_", model.no, "/alpha.spp_model_", model.no, "_1000samples.rds"))
alpha_df <- as_draws_df(alpha.fit) # takes awhile to convert to df



# get all the covariates using posterior package
betas.df <- readRDS(paste0(input.folder, "SPCD_stanoutput_joint_v3_model_", model.no, "/u_betas_model_", model.no, "_1000samples.rds"))
betas.quant <- betas.df %>% summarise_draws(median, ~quantile(., probs = c(0.025, 0.975))) %>%
  rename(`ci.lo` = "2.5%", `ci.hi` = "97.5%") %>%
  mutate(remper.cor = 0.5)
# relabel u_betas to meaningful species ids names
betas.quant$spp <- rep(1:17, ncovar)
betas.quant$cov <- rep(1:ncovar, each = 17)


covariate_names <- c(colnames(mod.data$xM))  # Replace with your covariate names
betas.quant$Covariate <- rep(covariate_names, each = 17)
betas.quant$Species <- rep(nspp[1:17,]$COMMON, ncovar)

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
# order species by hardwood softwood, then shade tolence
betas.quant$Species <- factor(betas.quant$Species, levels = SP.TRAITS$COMMON_NAME)


ggplot(data = na.omit(betas.quant), aes(x = Species, y = median, color = significance))+geom_point()+
  geom_errorbar(data = na.omit(betas.quant), aes(x = Species , ymin = ci.lo, ymax = ci.hi, color = significance), width = 0.1)+
  geom_abline(aes(slope = 0, intercept = 0), color = "grey", linetype = "dashed")+
  facet_wrap(~Covariate, scales= "free_y")+
  theme_bw(base_size = 14)+
  theme( axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.grid  = element_blank(), legend.position = "none")+
  ylab("Effect on Survival")+xlab("Parameter")+
  scale_color_manual(values = c("not overlapping zero"="darkgrey", "significant"="black"))
ggsave(height = 10, width = 15,dpi = 350, units = "in",paste0(output.folder,"SPCD_stanoutput_joint_v3/images/Estimated_effects_on_mortality_model_model_",model.no,"_all_species_betas.png"))

# get the population estimates for the betas

mubetas.estimates <- readRDS( paste0(input.folder, "SPCD_stanoutput_joint_v3_model_", model.no, "/beta_model_", model.no, "_1000samples.rds"))

mubetas.quant <- summarise_draws(mubetas.estimates, median, ~quantile(.x, probs = c(0.025, 0.975)))%>% rename(`ci.lo` = `2.5%`, `ci.hi` = `97.5%`)

mubetas.quant$remper.cor <- 0.5                                                        
mubetas.quant$spp <- 18
mubetas.quant$cov <- 1:ncovar
beta.names <- data.frame(cov = 1:ncovar, 
                         Covariate = unique(betas.quant$Covariate))

mubetas.quant <- left_join(mubetas.quant, beta.names)
main.table <- data.frame( 
  Species = "population", 
  spp = 18, 
  SPCD = 1000,
  COMMON = "population")

mubetas.quant <- left_join(mubetas.quant, main.table)


mubetas.quant$`significance` <- ifelse(mubetas.quant$ci.lo < 0 & mubetas.quant$ci.hi < 0, "significant", 
                                       ifelse(mubetas.quant$ci.lo > 0 & mubetas.quant$ci.hi > 0, "significant", "not overlapping zero"))

mubetas.quant <- mubetas.quant %>% arrange(by = median)
mubetas.quant$Covariate <- factor(mubetas.quant$Covariate, levels = mubetas.quant$Covariate)

ggplot(data = na.omit(mubetas.quant), aes(x = Covariate, y = median, color = significance))+geom_point()+
  geom_errorbar(data = na.omit(mubetas.quant), aes(x = Covariate , ymin = ci.lo, ymax = ci.hi, color = significance), width = 0.1)+
  geom_abline(aes(slope = 0, intercept = 0), color = "grey", linetype = "dashed")+facet_wrap(~Species)+theme_bw(base_size = 14)+
  theme( axis.text.x = element_text(angle = 65, hjust = 1), panel.grid  = element_blank(), legend.position = "none")+ylab("Effect on mortality")+xlab("Covariate")+
  scale_color_manual(values = c("not overlapping zero"="darkgrey", "significant"="black"))
ggsave(height = 5, width = 10, units = "in",paste0(output.folder, "SPCD_stanoutput_joint_v3/images/Estimated_effects_on_mortality_model_model_", model.no, "_all_species_population_betas.png"))



### combine the betas together
all.joint.betas <- rbind(mubetas.quant %>% select(colnames(betas.quant)), betas.quant) 
# reorder factors for better plotting
all.joint.betas$Species <- factor(all.joint.betas$Species, levels = c("population", SP.TRAITS$COMMON_NAME))
all.joint.betas$Covariate <- factor(all.joint.betas$Covariate, levels = unique(beta.names$Covariate))

all.joint.betas$Covariate <- factor(all.joint.betas$Covariate, levels = unique(betas.quant$Covariate))


# read in file with pretty covariate names:
Covariate.types.df <- read.csv("data/model_covariate_types_v2.csv")
# # join to new variable names and make sure they are in order
all.joint.betas <- all.joint.betas %>% left_join(., Covariate.types.df)

# }

all.joint.betas$Predictor <- factor( all.joint.betas$Predictor, 
                                     levels = Covariate.types.df$Predictor)



ggplot(data = na.omit(all.joint.betas), aes(x = Species, y = median, color = significance, shape = Species %in% "population"))+geom_point()+
  geom_errorbar(data = na.omit(all.joint.betas), aes(x = Species , ymin = ci.lo, ymax = ci.hi, color = significance), width = 0.1)+
  geom_abline(aes(slope = 0, intercept = 0), color = "grey", linetype = "dashed")+facet_wrap(~Predictor, scales = "free_y")+
  theme_bw(base_size = 12)+
  theme( axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.grid  = element_blank(), legend.position = "none")+ylab("Effect on survival")+xlab("Species")+
  scale_color_manual(values = c("not overlapping zero"="darkgrey", "significant"="black"))+
  scale_shape_manual(values = c("TRUE" = 15, "FALSE" = 19))

ggsave(height = 22, width = 12, dpi = 350, units = "in", paste0(output.folder,"SPCD_stanoutput_joint_v3/images/Estimated_betas_all_model_",model.no,".png"))

ggplot(data = na.omit(all.joint.betas) %>% filter(! Species %in% "population"), aes(x = Species, y = median, color = significance, shape = Species %in% "population"))+geom_point()+
  geom_errorbar(data = na.omit(all.joint.betas) %>% filter(! Species %in% "population"), aes(x = Species , ymin = ci.lo, ymax = ci.hi, color = significance), width = 0.1)+
  geom_abline(aes(slope = 0, intercept = 0), color = "grey", linetype = "dashed")+
  facet_wrap(~Predictor, scales = "free_y")+
  theme_bw(base_size = 12)+
  theme( axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.grid  = element_blank(), legend.position = "none")+ylab("Effect on survival")+xlab("Species")+
  scale_color_manual(values = c("not overlapping zero"="darkgrey", "significant"="black"))+
  scale_shape_manual(values = c("TRUE" = 15, "FALSE" = 19))

ggsave(height = 12, width = 12, dpi = 350, units = "in", paste0(output.folder,"SPCD_stanoutput_joint_v3/images/Estimated_betas_all_model_",model.no,"_no_population.png"))

###################################################################################

#colnames(Covariate.types.df) <- c("Covariate", "Predictor", )
# if(model.no == 6){

# # get the main effects
growth.diam <- c("Growth x Size", "Size", "Change in Size")
# competition <- c(unique(betas.quant$Covariate)[3:7])
# climate <- c(unique(betas.quant$Covariate)[8:13])
# site.vars  <- c(unique(betas.quant$Covariate)[14:18])
# main.effects <- c(growth.diam, competition, climate, site.vars)

ggplot(data = all.joint.betas %>% filter(! Species %in% "population" & pred.type %in% "Main Effects"), aes(x = Species, y = median, color = significance, shape = Species %in% "population"))+geom_point()+
  geom_errorbar(data = all.joint.betas %>% filter(! Species %in% "population"  & pred.type %in% "Main Effects"), aes(x = Species , ymin = ci.lo, ymax = ci.hi, color = significance), width = 0.1)+
  geom_abline(aes(slope = 0, intercept = 0), color = "grey", linetype = "dashed")+
  facet_wrap(~Covariate, scales = "free_y", ncol = 6)+theme_bw(base_size = 14)+
  theme( axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.grid  = element_blank(), legend.position = "none")+ylab("Effect on survival")+xlab("Covariate")+
  scale_color_manual(values = c("not overlapping zero"="darkgrey", "significant"="black"))+
  scale_shape_manual(values = c("TRUE" = 15, "FALSE" = 19))

ggsave(height = 6, width = 11, dpi = 350, units = "in", paste0(output.folder,"SPCD_stanoutput_joint_v3/images/Estimated_betas_main_effects_model_",model.no,"_no_population.png"))

# reorganize to have broader categories:

# predictor.class.definitions <- data.frame(Covariate = unique(betas.quant$Covariate), 
#            predictor.class = c("Growth, Diameter, Disturbance", "Growth, Diameter, Disturbance", 
#                      "Competition", "Site Conditions", "Competition", "Competition", "Growth, Diameter, Disturbance", 
#                      "Climate Normals", "Climate Normals", "Climate Normals", 
#                     "Climate Anomalies", "Climate Anomalies", "Climate Anomalies", 
#                     "Site Position", "Site Position", "Site Position", "Site Conditions", "Site Conditions", 
#                     rep("Growth Interactions", 17),  rep("Diameter Interactions", 16)
#                     ))
#all.joint.betas <-  left_join(all.joint.betas, predictor.class.definitions)

all.joint.betas <- all.joint.betas %>% select(-Pred.1, -Pred.2)
ggplot(data = all.joint.betas %>% filter(! Species %in% "population" & pred.type %in% "Main Effects"), aes(x =Covariate, y = median, group = Species, color = Species, shape = Species %in% "population"))+geom_point(position= position_dodge(width = 1))+
  #geom_errorbar(data = na.omit(all.joint.betas) %>% filter(! Species %in% "population"& Covariate %in% main.effects), aes(x = Species , ymin = ci.lo, ymax = ci.hi, group = Covariate, color = significance), width = 0.1)+
  geom_abline(aes(slope = 0, intercept = 0), color = "grey", linetype = "dashed")+
  facet_wrap(~predictor.class, scales = "free", ncol = 6)+theme_bw(base_size = 14)+
  theme( axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.grid  = element_blank(), legend.position = "none")+ylab("Effect on survival")+xlab("Covariate")+
  #scale_color_manual(values = c("not overlapping zero"="darkgrey", "significant"="black"))+
  scale_shape_manual(values = c("TRUE" = 15, "FALSE" = 19))+coord_flip()

#ggsave(height = 6, width = 11, units = "in", paste0(output.folder,"SPCD_stanoutput_joint_v2/images/Estimated_betas_main_effects_model_", model.no, "_no_population.png"))

growth.diam <-  ggplot(data = all.joint.betas %>% filter(! Species %in% "population" & predictor.class %in% "Growth & Size"), aes(y =median , x = Species, group = Covariate, color = significance, shape = Species %in% "population"))+geom_point()+
  geom_errorbar(data = all.joint.betas %>% filter(! Species %in% "population"& predictor.class %in% "Growth & Size"), aes(x = Species , ymin = ci.lo, ymax = ci.hi, group = Covariate, color = significance), width = 0.1)+
  geom_abline(aes(slope = 0, intercept = 0), color = "grey", linetype = "dashed")+
  facet_grid(rows = vars(Predictor), scales = "free_y")+
  theme_bw(base_size = 14)+
  theme( axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.grid  = element_blank(), legend.position = "none")+
  ylab("Effect on survival")+xlab(" ")+
  scale_color_manual(values = c("not overlapping zero"="darkgrey", "significant"="black"))+
  scale_shape_manual(values = c("TRUE" = 15, "FALSE" = 19))+coord_flip()


competition <- ggplot(data = na.omit(all.joint.betas) %>% filter(! Species %in% "population" & predictor.class %in% "Competition"), aes(y =median , x = Species, group = Covariate, color = significance, shape = Species %in% "population"))+geom_point()+
  geom_errorbar(data = na.omit(all.joint.betas) %>% filter(! Species %in% "population"& predictor.class %in% "Competition"), aes(x = Species , ymin = ci.lo, ymax = ci.hi, group = Covariate, color = significance), width = 0.1)+
  geom_abline(aes(slope = 0, intercept = 0), color = "grey", linetype = "dashed")+
  facet_wrap(~Predictor, scales = "free_y")+
  theme_bw(base_size = 14)+
  theme( axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.grid  = element_blank(), legend.position = "none")+
  ylab("Effect on survival")+xlab(" ")+
  scale_color_manual(values = c("not overlapping zero"="darkgrey", "significant"="black"))+
  scale_shape_manual(values = c("TRUE" = 15, "FALSE" = 19))+coord_flip()


normals <- ggplot(data = na.omit(all.joint.betas) %>% filter(! Species %in% "population" & predictor.class %in% "Climate"), aes(y =median , x = Species, group = Covariate, color = significance, shape = Species %in% "population"))+geom_point()+
  geom_errorbar(data = na.omit(all.joint.betas) %>% filter(! Species %in% "population"& predictor.class %in% "Climate"), aes(x = Species , ymin = ci.lo, ymax = ci.hi, group = Covariate, color = significance), width = 0.1)+
  geom_abline(aes(slope = 0, intercept = 0), color = "grey", linetype = "dashed")+
  facet_grid(rows = vars(Predictor), cols = vars(predictor.class),, scales = "free_y")+
  theme_bw(base_size = 14)+
  theme( axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.grid  = element_blank(), legend.position = "none")+
  ylab("Effect on survival")+xlab(" ")+
  scale_color_manual(values = c("not overlapping zero"="darkgrey", "significant"="black"))+
  scale_shape_manual(values = c("TRUE" = 15, "FALSE" = 19))+coord_flip()


# anomalies <- ggplot(data = na.omit(all.joint.betas) %>% filter(! Species %in% "population" & predictor.class %in% "Climate Anomalies"), aes(y =median , x = Species, group = Covariate, color = significance, shape = Species %in% "population"))+geom_point()+
#   geom_errorbar(data = na.omit(all.joint.betas) %>% filter(! Species %in% "population"& predictor.class %in% "Climate Anomalies"), aes(x = Species , ymin = ci.lo, ymax = ci.hi, group = Covariate, color = significance), width = 0.1)+
#   geom_abline(aes(slope = 0, intercept = 0), color = "grey", linetype = "dashed")+
#   facet_grid(rows = vars(Predictor), cols = vars(predictor.class),, scales = "free_y")+
#   theme_bw(base_size = 14)+
#   theme( axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.grid  = element_blank(), legend.position = "none")+
#   ylab("Effect on survival")+xlab(" ")+
#   scale_color_manual(values = c("not overlapping zero"="darkgrey", "significant"="black"))+
#   scale_shape_manual(values = c("TRUE" = 15, "FALSE" = 19))+coord_flip()
# 
# positions <- ggplot(data = na.omit(all.joint.betas) %>% filter(! Species %in% "population" & predictor.class %in% "Site Position"), aes(y =median , x = Species, group = Covariate, color = significance, shape = Species %in% "population"))+geom_point()+
#   geom_errorbar(data = na.omit(all.joint.betas) %>% filter(! Species %in% "population"& predictor.class %in% "Site Position"), aes(x = Species , ymin = ci.lo, ymax = ci.hi, group = Covariate, color = significance), width = 0.1)+
#   geom_abline(aes(slope = 0, intercept = 0), color = "grey", linetype = "dashed")+
#   facet_grid(rows = vars(Predictor), cols = vars(predictor.class),, scales = "free_y")+
#   theme_bw(base_size = 14)+
#   theme( axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.grid  = element_blank(), legend.position = "none")+
#   ylab("Effect on survival")+xlab(" ")+
#   scale_color_manual(values = c("not overlapping zero"="darkgrey", "significant"="black"))+
#   scale_shape_manual(values = c("TRUE" = 15, "FALSE" = 19))+coord_flip()
# 
site.cond <- ggplot(data = na.omit(all.joint.betas) %>% filter(! Species %in% "population" & predictor.class %in% "Site Conditions"), aes(y =median , x = Species, group = Covariate, color = significance, shape = Species %in% "population"))+geom_point()+
  geom_errorbar(data = na.omit(all.joint.betas) %>% filter(! Species %in% "population"& predictor.class %in% "Site Conditions"), aes(x = Species , ymin = ci.lo, ymax = ci.hi, group = Covariate, color = significance), width = 0.1)+
  geom_abline(aes(slope = 0, intercept = 0), color = "grey", linetype = "dashed")+
  facet_grid(rows = vars(Predictor), cols = vars(predictor.class), scales = "free_y")+
  theme_bw(base_size = 14)+
  theme( axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
         panel.grid  = element_blank(), legend.position = "none")+ylab("Effect on survival")+xlab("")+
  scale_color_manual(values = c("not overlapping zero"="darkgrey", "significant"="black"))+
  scale_shape_manual(values = c("TRUE" = 15, "FALSE" = 19))+coord_flip()

if(model.no >= 2){  
  png(height = 5, width =5, units = "in",res = 350,  paste0(output.folder, "SPCD_stanoutput_joint_v3/images/full_main_effects_summary_model_",model.no,".png")) 
  cowplot::plot_grid(growth.diam)
  dev.off()
}

if(model.no ==3){  
  png(height = 5, width =5, units = "in",res = 350,  paste0(output.folder, "SPCD_stanoutput_joint_v3/images/full_main_effects_summary_model_",model.no,".png")) 
  cowplot::plot_grid(growth.diam, competition, align = "hv")
  dev.off()
}

if(model.no ==4){  
  png(height = 5, width =5, units = "in",res = 300,  paste0(output.folder, "SPCD_stanoutput_joint_v3/images/full_main_effects_summary_model_",model.no,".png")) 
  cowplot::plot_grid(growth.diam, competition,normals, align = "hv")
  dev.off()
}

if(model.no ==5){  
  png(height = 5, width =5, units = "in",res = 300,  paste0(output.folder, "SPCD_stanoutput_joint_v3/images/full_main_effects_summary_model_",model.no,".png")) 
  cowplot::plot_grid(growth.diam, competition,normals,site.cond, align = "hv")
  dev.off()
}

if(model.no >=5){  
  png(height = 15, width = 12, units = "in",res = 350,  paste0(output.folder, "SPCD_stanoutput_joint_v3/images/full_main_effects_summary_model_",model.no,".png")) 
  cowplot::plot_grid(growth.diam, normals, #anomalies, 
                     competition, site.cond, #positions, 
                     align = "hv")
  dev.off()
  
  png(height = 8, width = 15, units = "in",res = 350,  paste0(output.folder, "SPCD_stanoutput_joint_v3/images/full_main_effects_summary_model_",model.no,"_horizontal.png")) 
  cowplot::plot_grid(growth.diam, normals, #anomalies, 
                     competition, site.cond, #positions, 
                     ncol = 4, 
                     align = "hv")
  dev.off()
}
library(FIESTA)
# do the same but reorder in terms of species types
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



growth.diam <-  ggplot(data = na.omit(all.joint.betas) %>% filter(! Species %in% "population" & predictor.class %in% "Growth & Size"), aes(y =median , x = Species, group = Covariate, color = significance, shape = Species %in% "population"))+geom_point()+
  geom_errorbar(data = na.omit(all.joint.betas) %>% filter(! Species %in% "population"& predictor.class %in% "Growth & Size"), aes(x = Species , ymin = ci.lo, ymax = ci.hi, group = Covariate, color = significance), width = 0.1)+
  geom_abline(aes(slope = 0, intercept = 0), color = "grey", linetype = "dashed")+
  facet_grid(rows = vars(Predictor), cols = vars(predictor.class),, scales = "free_y")+
  theme_bw(base_size = 14)+
  theme( axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.grid  = element_blank(), legend.position = "none")+
  ylab("Effect on survival")+xlab(" ")+
  scale_color_manual(values = c("not overlapping zero"="darkgrey", "significant"="black"))+
  scale_shape_manual(values = c("TRUE" = 15, "FALSE" = 19))+coord_flip()


competition <- ggplot(data = na.omit(all.joint.betas) %>% filter(! Species %in% "population" & predictor.class %in% "Competition"), aes(y =median , x = Species, group = Covariate, color = significance, shape = Species %in% "population"))+geom_point()+
  geom_errorbar(data = na.omit(all.joint.betas) %>% filter(! Species %in% "population"& predictor.class %in% "Competition"), aes(x = Species , ymin = ci.lo, ymax = ci.hi, group = Covariate, color = significance), width = 0.1)+
  geom_abline(aes(slope = 0, intercept = 0), color = "grey", linetype = "dashed")+
  facet_grid(rows = vars(Predictor), cols = vars(predictor.class),, scales = "free_y")+
  theme_bw(base_size = 14)+
  theme( axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.grid  = element_blank(), legend.position = "none")+
  ylab("Effect on survival")+xlab(" ")+
  scale_color_manual(values = c("not overlapping zero"="darkgrey", "significant"="black"))+
  scale_shape_manual(values = c("TRUE" = 15, "FALSE" = 19))+coord_flip()


normals <- ggplot(data = na.omit(all.joint.betas) %>% filter(! Species %in% "population" & predictor.class %in% "Climate"), aes(y =median , x = Species, group = Covariate, color = significance, shape = Species %in% "population"))+geom_point()+
  geom_errorbar(data = na.omit(all.joint.betas) %>% filter(! Species %in% "population"& predictor.class %in% "Climate"), aes(x = Species , ymin = ci.lo, ymax = ci.hi, group = Covariate, color = significance), width = 0.1)+
  geom_abline(aes(slope = 0, intercept = 0), color = "grey", linetype = "dashed")+
  facet_grid(rows = vars(Predictor), cols = vars(predictor.class),, scales = "free_y")+
  theme_bw(base_size = 14)+
  theme( axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.grid  = element_blank(), legend.position = "none")+
  ylab("Effect on survival")+xlab(" ")+
  scale_color_manual(values = c("not overlapping zero"="darkgrey", "significant"="black"))+
  scale_shape_manual(values = c("TRUE" = 15, "FALSE" = 19))+coord_flip()


# anomalies <- ggplot(data = na.omit(all.joint.betas) %>% filter(! Species %in% "population" & predictor.class %in% "Climate Anomalies"), aes(y =median , x = Species, group = Covariate, color = significance, shape = Species %in% "population"))+geom_point()+
#   geom_errorbar(data = na.omit(all.joint.betas) %>% filter(! Species %in% "population"& predictor.class %in% "Climate Anomalies"), aes(x = Species , ymin = ci.lo, ymax = ci.hi, group = Covariate, color = significance), width = 0.1)+
#   geom_abline(aes(slope = 0, intercept = 0), color = "grey", linetype = "dashed")+
#   facet_grid(rows = vars(Predictor), cols = vars(predictor.class),, scales = "free_y")+
#   theme_bw(base_size = 14)+
#   theme( axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.grid  = element_blank(), legend.position = "none")+
#   ylab("Effect on survival")+xlab(" ")+
#   scale_color_manual(values = c("not overlapping zero"="darkgrey", "significant"="black"))+
#   scale_shape_manual(values = c("TRUE" = 15, "FALSE" = 19))+coord_flip()
# 
# positions <- ggplot(data = na.omit(all.joint.betas) %>% filter(! Species %in% "population" & predictor.class %in% "Site Position"), aes(y =median , x = Species, group = Covariate, color = significance, shape = Species %in% "population"))+geom_point()+
#   geom_errorbar(data = na.omit(all.joint.betas) %>% filter(! Species %in% "population"& predictor.class %in% "Site Position"), aes(x = Species , ymin = ci.lo, ymax = ci.hi, group = Covariate, color = significance), width = 0.1)+
#   geom_abline(aes(slope = 0, intercept = 0), color = "grey", linetype = "dashed")+
#   facet_grid(rows = vars(Predictor), cols = vars(predictor.class),, scales = "free_y")+
#   theme_bw(base_size = 14)+
#   theme( axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.grid  = element_blank(), legend.position = "none")+
#   ylab("Effect on survival")+xlab(" ")+
#   scale_color_manual(values = c("not overlapping zero"="darkgrey", "significant"="black"))+
#   scale_shape_manual(values = c("TRUE" = 15, "FALSE" = 19))+coord_flip()
# 
site.cond <- ggplot(data = na.omit(all.joint.betas) %>% filter(! Species %in% "population" & predictor.class %in% "Site Conditions"), aes(y =median , x = Species, group = Covariate, color = significance, shape = Species %in% "population"))+geom_point()+
  geom_errorbar(data = na.omit(all.joint.betas) %>% filter(! Species %in% "population"& predictor.class %in% "Site Conditions"), aes(x = Species , ymin = ci.lo, ymax = ci.hi, group = Covariate, color = significance), width = 0.1)+
  geom_abline(aes(slope = 0, intercept = 0), color = "grey", linetype = "dashed")+
  facet_grid(rows = vars(Predictor), cols = vars(predictor.class), scales = "free_y")+
  theme_bw(base_size = 14)+
  theme( axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
         panel.grid  = element_blank(), legend.position = "none")+ylab("Effect on survival")+xlab("")+
  scale_color_manual(values = c("not overlapping zero"="darkgrey", "significant"="black"))+
  scale_shape_manual(values = c("TRUE" = 15, "FALSE" = 19))+coord_flip()


if(model.no >=5){  
  png(height = 14, width = 12, units = "in",res = 300,  paste0(output.folder, "SPCD_stanoutput_joint_v3/images/full_main_effects_summary_reorganized_model_",model.no,".png")) 
  cowplot::plot_grid(growth.diam, normals, #anomalies, 
                     competition, site.cond, #positions, 
                     align = "hv")
  dev.off()
}


# get the species alpha estimates
alphas.estimates <- readRDS( paste0(input.folder,"SPCD_stanoutput_joint_v3_model_",model.no,"/alpha.spp_model_",model.no,"_1000samples.rds"))

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
  geom_abline(aes(slope = 0, intercept = 0), color = "grey", linetype = "dashed")+theme_bw(base_size = 14)+
  theme( axis.text.x = element_text(angle = 65, hjust = 1), panel.grid  = element_blank(), legend.position = "none")+ylab("alpha estimate")+xlab("SPECIES")+
  scale_color_manual(values = c("not overlapping zero"="darkgrey", "significant"="black"))
ggsave(height = 5, width = 10, units = "in",paste0(output.folder,"SPCD_stanoutput_joint_v3/images/Estimated_alpha_SPP_model_",model.no,"_all_species_alphas.png"))

# combine the population and the species level betas

## get population level alpha estimates:
alphas.pop.estimates <- readRDS( paste0(input.folder,"SPCD_stanoutput_joint_v3_model_",model.no,"/alpha.p_model_",model.no,"_1000samples.rds"))

alphas.pop.quant <- summarise_draws(alphas.pop.estimates, median, ~quantile(.x, probs = c(0.025, 0.975)))%>% rename(`ci.lo` = `2.5%`, `ci.hi` = `97.5%`)

alphas.pop.quant$spp <- 18

#main.table$COMMON = "population"
# main.table$SPCD = 1000
alphas.pop.quant <- left_join(alphas.pop.quant, main.table %>% select(-Species))


alphas.pop.quant$`significance` <- ifelse(alphas.pop.quant$ci.lo < 0 & alphas.pop.quant$ci.hi < 0, "significant", 
                                          ifelse(alphas.pop.quant$ci.lo > 0 & alphas.pop.quant$ci.hi > 0, "significant", "not overlapping zero"))

#alphas.quant <- alphas.quant %>% arrange(by = median)
#alphas.quant$parameter <- factor(alphas.quant$parameter, levels = alphas.quant$parameter)

ggplot(data = na.omit(alphas.pop.quant), aes(x = COMMON, y = median, color = significance))+geom_point()+
  geom_errorbar(data = na.omit(alphas.pop.quant), aes(x = COMMON , ymin = ci.lo, ymax = ci.hi, color = significance), width = 0.1)+
  geom_abline(aes(slope = 0, intercept = 0), color = "grey", linetype = "dashed")+theme_bw(base_size = 14)+
  theme( axis.text.x = element_text(angle = 65, hjust = 1), panel.grid  = element_blank(), legend.position = "none")+ylab("alpha estimate")+xlab("SPECIES")+
  scale_color_manual(values = c("not overlapping zero"="darkgrey", "significant"="black"))


### combine the alphas together
all.joint.alphas <- rbind(alphas.pop.quant, alphas.quant) 
all.joint.alphas$COMMON <- factor(all.joint.alphas$COMMON, levels = c("population", alphas.quant$COMMON[order(unique(alphas.quant$COMMON))]))

ggplot(data = na.omit(all.joint.alphas), aes(x = COMMON, y = median, color = significance, shape = COMMON %in% "population"))+geom_point()+
  geom_errorbar(data = na.omit(all.joint.alphas), aes(x = COMMON , ymin = ci.lo, ymax = ci.hi, color = significance), width = 0.1)+
  geom_abline(aes(slope = 0, intercept = 0), color = "grey", linetype = "dashed")+theme_bw(base_size = 14)+
  theme( axis.text.x = element_text(angle = 65, hjust = 1), panel.grid  = element_blank(), legend.position = "none")+ylab("alpha estimate")+xlab("SPECIES")+
  scale_color_manual(values = c("not overlapping zero"="darkgrey", "significant"="black"))+
  scale_shape_manual(values = c("TRUE" = 15, "FALSE" = 19))

ggsave(height = 5, width = 6, dpi = 350, units = "in",paste0(output.folder,"SPCD_stanoutput_joint_v3/images/Estimated_alpha_all_model_",model.no,".png"))

## combine alphas and betas together to make the figure for the paper
# 
# all.joint.alphas <- all.joint.alphas %>% mutate(parameter = "alpha") %>% select(variable, median, ci.lo, ci.hi, spp, parameter, SPCD.id, COMMON, significance)
# all.joint.betas <- all.joint.betas #%>% select(-param.no)
# all.params <- rbind(all.joint.betas, all.joint.alphas)
# 
# ggplot(data = na.omit(all.params), aes(x = COMMON, y = median, color = significance, shape = COMMON %in% "population"))+geom_point()+
#   geom_errorbar(data = na.omit(all.params), aes(x = COMMON , ymin = ci.lo, ymax = ci.hi, color = significance), width = 0.1)+
#   geom_abline(aes(slope = 0, intercept = 0), color = "grey", linetype = "dashed")+facet_wrap(~parameter)+theme_bw(base_size = 14)+
#   theme( axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.grid  = element_blank(), legend.position = "none")+ylab("Effect on survival")+xlab("Parameter")+
#   scale_color_manual(values = c("not overlapping zero"="darkgrey", "significant"="black"))+
#   scale_shape_manual(values = c("TRUE" = 15, "FALSE" = 19))
# 
# ggsave(height = 10, width = 11, units = "in",dpi = 350,paste0("SPCD_stanoutput_joint/images/Estimated_parameters_model_", model.no, "_joint.png"))
# 
#########################################################################################
# Estimate conditional effects for each species and covariate
#########################################################################################
# read in the model samples
# fit <- readRDS(paste0(output.folder,"/SPCD_stanoutput_joint_v2/u_betas_model_",model.no,"_1000samples.rds"))
# fit_ssm_df <- as_draws_df(fit) # takes awhile to convert to df
# 


for(i in 17:1){
  cat(i)
  
  SPCD.id <- nspp[i,]$SPCD
  beta.species.names <- betas.quant %>% filter(spp %in% i) %>% select(variable)
  # select only the betas and intercepts for the species of interest:
  beta <- subset_draws(betas.df, variable = beta.species.names$variable) %>% select(-.chain, -.iteration, -.draw) 
  beta_0 <- data.frame(subset_draws(alpha_df, variable = paste0("alpha_SPP[",i,"]"))) %>% select(-.chain, -.iteration, -.draw)  # Intercept
  # read in the species data and covariates to get the min and max and get ranges
  mod.data <-
    readRDS (
      paste0(
        input.folder,
        "SPCD_stanoutput_joint_v3_model_",model.no,"/all_SPCD_model_",
        model.no,
        ".RDS"
      )
    )
  
  
  
  covariate_names <- c(colnames(mod.data$xM))  # Replace with your covariate names
  
  var.mins <- as.vector(apply(data.frame(mod.data$xM),2 , function(x)quantile(x, 0.025)))
  var.maxes <- as.vector(apply(data.frame(mod.data$xM),2 , function(x) quantile(x, 0.975)))
  var.75s <- as.vector(apply(data.frame(mod.data$xM),2 , function(x) quantile(x, 0.75)))
  var.25s <- as.vector(apply(data.frame(mod.data$xM),2 , function(x) quantile(x, 0.25)))
  
  covariate_ranges_df <- data.frame(
    covariate  = covariate_names, 
    mins = var.mins, 
    maxes = var.maxes, 
    qq25 = var.25s, 
    qq75 = var.75s
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
  
  plot_data$Covariate <- factor(plot_data$Covariate, levels = unique(plot_data$Covariate))
  
  # make a plot of covariate responses for each species
  ggplot(plot_data, aes(x = Value, y = 1-mean, color = Covariate)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Covariate), alpha = 0.2, color = NA) +
    facet_wrap(~ Covariate, scales = "free") +
    labs(
      x = "Covariate Value",
      y = "Annual Probability of Mortality",
      title = paste0("Effect of Covariates on Probability of Mortality for ",nspp[i,]$Species)
    ) +
    theme_bw(base_size = 12)+theme(panel.grid = element_blank(), legend.position = "none") 
  
  ggsave(height = 10, width = 12, paste0(output.folder, "SPCD_stanoutput_joint_v3/predicted_mort/model_",model.no,"_all_marginal_SPCD_", SPCD.id, "_joint.png"))
  
  # save for this species:
  plot_data$SPCD <- SPCD.id
  saveRDS(plot_data, paste0(output.folder, "SPCD_stanoutput_joint_v3/predicted_mort/model_",model.no,"_all_marginal_SPCD_", SPCD.id, ".RDS") )
  
  # make interaction term plots:
  
  interaction.ranges <- left_join(covariate_ranges_df, Covariate.types.df %>% rename("covariate" = "Covariate"))#%>%
  #filter(!is.na(Pred.2)) 
  dia.ranges <- interaction.ranges %>% filter(covariate %in% "DIA_scaled")
  diadiff.ranges <- interaction.ranges %>% filter(covariate %in% "DIA_DIFF_scaled")
  
  main.effectsi <- interaction.ranges %>% filter(pred.type %in% "Main Effects") %>%
    select(Pred.1, mins:medians)%>%
    rename("p1.mins" = "mins", "p1.maxes" = "maxes", "p1.medians" = "medians")%>%
    select(Pred.1, p1.mins, p1.maxes, p1.medians)
  
  interactions <- interaction.ranges %>% filter(!is.na(Pred.2))%>%
    mutate(interaction.qq25 = ifelse(Pred.2 %in% "DIA_DIFF_scaled", diadiff.ranges$qq25, 
                                     dia.ranges$qq25), 
           interaction.qq75 = ifelse(Pred.2 %in% "DIA_DIFF_scaled", diadiff.ranges$qq75, 
                                     dia.ranges$qq75), 
           interaction.qq50 = ifelse(Pred.2 %in% "DIA_DIFF_scaled", mean(mod.data$xM[,"DIA_DIFF_scaled"]), 
                                     mean(mod.data$xM[,"DIA_scaled"])))%>% 
    left_join(.,main.effectsi)
  
  # for each interaction term, generate 3 prediction lines:
  interaction.list <- list()
  for(j in 1:length (interactions[,2])){
    
    p1.seq.median <-  seq(interactions[j,"p1.mins"], interactions[j,"p1.maxes"], length.out = 25)
    p1.seq.hi <-  seq(interactions[j,"p1.mins"], interactions[j,"p1.maxes"], length.out = 25)
    p1.seq.lo <-  seq(interactions[j,"p1.mins"], interactions[j,"p1.maxes"], length.out = 25)
    
    interaction.list[[j]] <- data.frame(p1.value = c(p1.seq.median, p1.seq.hi, p1.seq.lo), 
                                        p2.value = c(rep(interactions[j, "interaction.qq50"], 25), rep(interactions[j, "interaction.qq75"], 25), rep(interactions[j, "interaction.qq25"], 25)),
                                        p2.rank = c(rep("median", 25), rep("high", 25), rep("low", 25)),
                                        covariate = rep(interactions[j,"covariate"], 25*3), 
                                        Pred.1 = rep(interactions[j,"Pred.1"], 25*3), 
                                        Pred.2 = rep(interactions[j,"Pred.2"], 25*3))
  }
  names(interaction.list) <- c(interactions$covariate)
  
  interaction.list[[j]]
  
  interaction_effects <- function(j){
    
    int.list <- interaction.list[[j]]
    cov_name <- unique(int.list$covariate)
    Pred.1.name <- unique(int.list$Pred.1)
    Pred.2.name <- unique(int.list$Pred.2)
    
    covariate_matrix_int <- matrix(covariate_ranges_df$medians, 
                                   nrow = length(int.list$p1.value), 
                                   ncol = length(covariate_names), byrow = TRUE)
    colnames( covariate_matrix_int) <- covariate_names
    covariate_matrix_int[, cov_name] <- int.list$p1.value * int.list$p2.value
    covariate_matrix_int[, Pred.1.name] <- int.list$p1.value
    covariate_matrix_int[, Pred.2.name] <- int.list$p2.value
    #covariate_matrix[, Pred.1.name] <- int.list$ 
    probabilities <- list()
    for(d in 1:length(covariate_matrix_int[,1])){
      linear_predictor <- beta_0[,1] + rowSums( as.matrix(beta) %*%  covariate_matrix_int[d,] )
      probabilities[[d]] <-  as.vector(inv_logit_fxn(linear_predictor))
    }
    
    # Summarize probabilities (mean and 90% credible interval)
    prob_summary <- lapply(probabilities, function(p) {
      c(mean = median(p), ci.lo = quantile(p, 0.1), ci.hi = quantile(p, 0.9))
    })
    
    
    
    prob_summary.df <- do.call(rbind, prob_summary)
    prob_summary_covariate <- cbind(int.list, prob_summary.df)
    # Combine with covariate range for plotting
    prob_summary_covariate
  }
  
  allmarginal.interaction.spp <- do.call(rbind, lapply(1:length(interaction.list), interaction_effects))
  ggplot(data = allmarginal.interaction.spp %>% filter(covariate %in% "Ndep.scaled_DIA.int"))+
    geom_line(aes(x = p1.value, y = mean, color = p2.rank))+
    geom_ribbon(aes(x = p1.value, ymin = `ci.lo.10%`, ymax = `ci.hi.90%`, fill = p2.rank), alpha = 0.5)
  allmarginal.interaction.spp$SPCD <- SPCD.id
  saveRDS(allmarginal.interaction.spp, paste0(output.folder, "SPCD_stanoutput_joint_v3/predicted_mort/model_",model.no,"_interaction_marginal_SPCD_", SPCD.id, ".RDS") )
  
  
}


################################################################################
# combine all of the species conditional responses together
################################################################################

# for the main effects:
marginal_response <- interaction_response <- list()
for (i in 1:17){
  SPCD.id <- nspp[i,]$SPCD
  marginal_response[[i]] <- readRDS( paste0(output.folder, "SPCD_stanoutput_joint_v3/predicted_mort/model_",model.no,"_all_marginal_SPCD_", SPCD.id, ".RDS") )
  # read in the interaction effects:
  interaction_response[[i]] <-  readRDS(paste0(output.folder, "SPCD_stanoutput_joint_v3/predicted_mort/model_",model.no,"_interaction_marginal_SPCD_", SPCD.id, ".RDS") )
  
}
marginal_response_df <- do.call(rbind, marginal_response)

# get species common names
marginal_response_df$Species <- FIESTA::ref_species[match(marginal_response_df$SPCD, FIESTA::ref_species$SPCD),]$COMMON
marginal_response_df$Species <- factor(marginal_response_df$Species, levels = SP.TRAITS$COMMON_NAME)
marginal_response_df$Predictor <- Covariate.types.df[match(marginal_response_df$Covariate, Covariate.types.df$Covariate),]$Predictor
marginal_response_df$Predictor <- factor(marginal_response_df$Predictor, levels =  unique(marginal_response_df$Predictor))
marginal_response_df$Colors <- SP.TRAITS[match(marginal_response_df$Species, SP.TRAITS$COMMON_NAME),]$Color
marginal_response_df$`Shade Tolerance` <- SP.TRAITS[match(marginal_response_df$Species, SP.TRAITS$COMMON_NAME),]$`Shade Tolerance`


interaction_response_df <- do.call(rbind, interaction_response)

# get species common names
interaction_response_df$Species <- FIESTA::ref_species[match(interaction_response_df$SPCD, FIESTA::ref_species$SPCD),]$COMMON
interaction_response_df$Species <- factor(interaction_response_df$Species, levels = SP.TRAITS$COMMON_NAME)
interaction_response_df$Predictor <- Covariate.types.df[match(interaction_response_df$covariate, Covariate.types.df$Covariate),]$Predictor
interaction_response_df$Predictor <- factor(interaction_response_df$Predictor, levels =  unique(interaction_response_df$Predictor))
interaction_response_df$Colors <- SP.TRAITS[match(interaction_response_df$Species, SP.TRAITS$COMMON_NAME),]$Color
interaction_response_df$`Shade Tolerance` <- SP.TRAITS[match(interaction_response_df$Species, SP.TRAITS$COMMON_NAME),]$`Shade Tolerance`
interaction_response_df <- interaction_response_df %>% rename("ci.hi.90"= "ci.hi.90%", 
                                                              "ci.lo.10" = "ci.lo.10%")
class(interaction_response_df$ci.lo.10)
ggplot(data = interaction_response_df %>% filter(Predictor %in% "Diameter x DIA_DIFF")) +
  geom_ribbon(
    aes(x = p1.value, ymin = 1-ci.hi.90, ymax = 1-ci.lo.10, fill = p2.rank), alpha = 0.2, color = NA) +
  
  geom_line(aes(x = p1.value, y = 1-mean, color = p2.rank, group = p2.rank),size = 1) +
  facet_wrap(~ Species, scales = "free") +
  labs(
    x = "Covariate Value",
    y = "Annual Probability of Mortality",
    title = "Effect of Covariates on Probability of Mortality"
  ) +
  #species_fill + species_color + species_linetype+
  theme_bw(base_size = 14)+theme(panel.grid = element_blank(), legend.position = "bottom") 




ggplot(marginal_response_df, aes(x = Value, y = 1-mean, color = Species, linetype = `Shade Tolerance`)) +
  geom_line(size = 1) +
  facet_wrap(~ Predictor, scales = "free") +
  labs(
    x = "Covariate Value",
    y = "Annual Probability of Mortality",
    title = "Effect of Covariates on Probability of Mortality"
  ) +
  species_fill + species_color + species_linetype+
  theme_bw(base_size = 14)+theme(panel.grid = element_blank(), legend.position = "bottom") 
ggsave(height = 12, width = 15, dpi = 350,
       paste0(output.folder, "SPCD_stanoutput_joint_v3/predicted_mort/model_",model.no,"_species_level_models_all_species_marginal_effects_joint.png"))

# make the same plot but without BF
ggplot(marginal_response_df %>% filter(!Species %in% "balsam fir"), aes(x = Value, y = 1-mean, color = Species, linetype = `Shade Tolerance`)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
  facet_wrap(~ Predictor, scales = "free") +
  labs(
    x = "Covariate Value",
    y = "Annual Probability of Mortality",
    title = "Effect of Covariates on Probability of Mortality"
  ) +
  species_fill + species_color + species_linetype +
  theme_bw(base_size = 14)+theme(panel.grid = element_blank(), legend.position = "bottom") 
ggsave(height = 12, width = 15, dpi = 350,
       paste0(output.folder, "SPCD_stanoutput_joint_v3/predicted_mort/model_",model.no,"_species_level_models_all_species_marginal_effects_joint_noBF.png"))




## split up by variable type:

marginal_response_df <- marginal_response_df %>% left_join(., Covariate.types.df)
marginal_response_df$Predictor <- factor(marginal_response_df$Predictor, levels = unique(marginal_response_df$Predictor))


ggplot(marginal_response_df, aes(x = Value, y = 1-mean, color = Species, linetype = `Shade Tolerance`)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
  facet_wrap(~ Predictor, scales = "free") +
  labs(
    x = "Scaled Covariate Value",
    y = "Annual Probability of Mortality",
    title = "Effect of Covariates on Probability of Mortality"
  ) +
  species_fill + species_color + species_linetype +
  theme_bw(base_size = 14)+theme(panel.grid = element_blank(), legend.position = "bottom") 
ggsave(height = 12, width = 14.5, dpi = 350,
       paste0(output.folder, "SPCD_stanoutput_joint_v3/predicted_mort/summary/model_",model.no,"_species_level_models_all_species_marginal_effects.png"))

# same but without Balsam fir
ggplot(marginal_response_df %>% filter(!Species %in% "balsam fir"), aes(x = Value, y = 1-mean, color = Species, linetype = `Shade Tolerance`)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
  facet_wrap(~ Predictor, scales = "free") +
  labs(
    x = "Scaled Covariate Value",
    y = "Annual Probability of Mortality",
    title = "Effect of Covariates on Probability of Mortality"
  ) +
  species_fill + species_color + species_linetype +
  theme_bw(base_size = 14)+theme(panel.grid = element_blank(), legend.position = "bottom") 
ggsave(height = 12, width = 14.5, dpi = 350,
       paste0(output.folder, "SPCD_stanoutput_joint_v3/predicted_mort/summary/model_",model.no,"_species_level_models_all_species_marginal_effects_noBF.png"))




# generate figures for all of these separately 
ggplot(marginal_response_df %>% 
         filter( predictor.class %in% c("Growth & Size", "Competition", "Climate", "Site Conditions") & 
                   !Predictor %in% "Diameter x DIA_DIFF"), aes(x = Value, y = 1-mean, color = Species, linetype = `Shade Tolerance`)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
  facet_wrap(~ Predictor, scales = "free") +
  labs(
    x = "Scaled Covariate Value",
    y = "Annual Probability of Mortality",
    title = "Effect of Covariates on Probability of Mortality"
  ) +
  species_fill + species_color + species_linetype +
  theme_bw(base_size = 14)+theme(panel.grid = element_blank(), legend.position = "bottom")
ggsave(height = 10, width = 14, dpi = 350,
       paste0(output.folder,  "SPCD_stanoutput_joint_v3/predicted_mort/summary/model_",model.no,"_species_level_model_marginal_main_effects.png"))

# balsam fir really domninates the graph so take it out to view
ggplot(marginal_response_df %>% 
         filter( predictor.class %in% c("Growth & Size", "Competition", "Climate", "Site Conditions") & 
                   !Predictor %in% "Diameter x DIA_DIFF" & !Species %in% "balsam fir"), aes(x = Value, y = 1-mean, color = Species, linetype = `Shade Tolerance`)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
  facet_wrap(~ Predictor, scales = "free") +
  labs(
    x = "Scaled Covariate Value",
    y = "Annual Probability of Mortality",
    title = "Effect of Covariates on Probability of Mortality"
  ) +
  species_fill + species_color + species_linetype +
  theme_bw(base_size = 14)+theme(panel.grid = element_blank(), legend.position = "bottom")
ggsave(height = 10, width = 12.5, dpi = 350,
       paste0(output.folder, "SPCD_stanoutput_joint_v3/predicted_mort/summary/model_", model.no, "_species_level_model_marginal_main_effects_no_balsam_fir.png"))


ggplot(marginal_response_df %>% filter( predictor.class %in% "Growth & Size" ), aes(x = Value, y = 1-mean, color = Species, linetype = `Shade Tolerance`)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
  facet_wrap(~ Predictor, scales = "free") +
  labs(
    x = "Scaled Covariate Value",
    y = "Annual Probability of Mortality",
    title = "Growth & Size effects on Probability of Mortality"
  ) +
  species_fill + species_color + species_linetype +
  theme_bw(base_size = 14)+theme(panel.grid = element_blank())
ggsave(height = 6, width = 12, dpi = 350,
       paste0(output.folder, "SPCD_stanoutput_joint_v3/predicted_mort/summary/model_", model.no, "_species_level_model_marginal_growth_diam.png"))

if(model.no >=6){
  ggplot(marginal_response_df %>% filter( predictor.class %in% "Competition" ), aes(x = Value, y = 1-mean, color = Species, linetype = `Shade Tolerance`)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
    facet_wrap(~ Predictor, scales = "free", ncol = 5) +
    labs(
      x = "Scaled Covariate Value",
      y = "Annual Probability of Mortality",
      title = "Competition effects on Probability of Mortality"
    ) +
    species_fill + species_color + species_linetype +
    theme_bw(base_size = 14)+theme(panel.grid = element_blank())
  ggsave(height = 6, width = 12, dpi = 350,
         paste0(output.folder, "SPCD_stanoutput_joint_v3/predicted_mort/summary/model_", model.no, "_species_level_model_marginal_competition.png"))
  
  
  ggplot(marginal_response_df %>% filter(  predictor.class %in% "Climate" ), aes(x = Value, y = 1-mean, color = Species, linetype = `Shade Tolerance`)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
    facet_wrap(~ Predictor, scales = "free", ncol = 2) +
    labs(
      x = "Scaled Covariate Value",
      y = "Annual Probability of Mortality",
      title = "Climate effects on Probability of Mortality"
    ) +
    species_fill + species_color + species_linetype +
    theme_bw(base_size = 14)+theme(panel.grid = element_blank())
  ggsave(height = 8, width = 12, dpi = 350,
         paste0(output.folder, "SPCD_stanoutput_joint_v3/predicted_mort/summary/model_", model.no, "_species_level_model_marginal_climate_vars.png"))
  
  
  ggplot(marginal_response_df %>% filter( predictor.class %in% "Site Conditions" ), aes(x = Value, y = 1-mean, color = Species, linetype = `Shade Tolerance`)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
    facet_wrap(~ Predictor, scales = "free", ncol = 5) +
    labs(
      x = "Scaled Covariate Value",
      y = "Annual Probability of Mortality",
      title = "Site Condition effects on Probability of Mortality"
    ) +
    species_fill + species_color + species_linetype +
    theme_bw(base_size = 14)+theme(panel.grid = element_blank())
  ggsave(height = 6, width = 12, dpi = 350,
         paste0(output.folder, "SPCD_stanoutput_joint_v3/predicted_mort/summary/model_", model.no, "_species_level_model_marginal_site_vars.png"))
  
  ### interaction terms
  ggplot(marginal_response_df %>% filter( predictor.class %in% "Competiton x G & S" & !Species %in% "balsam fir"), aes(x = Value, y = 1-mean, color = Species, linetype = `Shade Tolerance`)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
    facet_wrap(~ Predictor, scales = "free", ncol = 3) +
    labs(
      x = "Covariate Value",
      y = "Annual Probability of Mortality" #,
      # title = "Competition x G & S  on Probability of Mortality"
    ) +
    species_fill + species_color + species_linetype +
    theme_bw(base_size = 14)+theme(panel.grid = element_blank())
  ggsave(height = 6, width = 12, dpi = 350,
         paste0(output.folder, "SPCD_stanoutput_joint_v3/predicted_mort/summary/model_", model.no, "_species_level_model_marginal_growth_competition_noBF.png"))
  
  
  ggplot(marginal_response_df %>% filter( predictor.class %in% "Climate x G & S" & !Species %in% "balsam fir" ), aes(x = Value, y = 1-mean, color = Species, linetype = `Shade Tolerance`)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
    facet_wrap(~ Predictor, scales = "free", ncol = 4) +
    labs(
      x = "Covariate Value",
      y = "Annual Probability of Mortality"#,
      #title = "Effect of Covariates on Probability of Mortality"
    ) +
    species_fill + species_color + species_linetype +
    theme_bw(base_size = 14)+theme(panel.grid = element_blank())
  ggsave(height = 8, width = 14, 
         paste0(output.folder, "SPCD_stanoutput_joint_v3/predicted_mort/summary/model_", model.no, "_species_level_model_marginal_growth_cliamte_noBF.png"))
  
  ggplot(marginal_response_df %>% filter(  predictor.class %in% "Site x G & S" & !Species %in% "balsam fir" ), aes(x = Value, y = 1-mean, color = Species, linetype = `Shade Tolerance`)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
    facet_wrap(~ Predictor, scales = "free", ncol = 3) +
    labs(
      x = "Covariate Value",
      y = "Annual Probability of Mortality"#,
      #title = "Effect of Covariates on Probability of Mortality"
    ) +
    species_fill + species_color + species_linetype +
    theme_bw(base_size = 14)+theme(panel.grid = element_blank())
  ggsave(height = 6, width = 12, 
         paste0(output.folder, "SPCD_stanoutput_joint_v3/predicted_mort/summary/model_", model.no, "_species_level_model_marginal_growth_site_noBF.png"))
  
  ##############################################################################################
  # make one giant figure for the marginal species effects:
  # generate figures for all of these separately 
  growth.marginal <- ggplot(marginal_response_df %>% filter( predictor.class %in% "Growth & Size"), aes(x = Value, y = 1-mean, color = Species, linetype = Species)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
    facet_grid(cols = vars(predictor.class), rows = vars(Predictor), scales = "free") +
    labs(
      x = "Predictor Value",
      y = "Annual Probability of Mortality",
      #title = "Effect of Predictors on Probability of Mortality"
    ) +
    species_fill + species_color + named_species_linetype +
    theme_bw(base_size = 14)+theme(panel.grid = element_blank(), 
                                   legend.key.size = unit(1, "cm")  
    )
  
  Competition.marginal <- ggplot(marginal_response_df %>% filter( predictor.class %in% "Competition"), aes(x = Value, y = 1-mean, color = Species, linetype = `Shade Tolerance`)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
    facet_grid(cols = vars(predictor.class), rows = vars(Predictor), scales = "free") +
    labs(
      x = "Predictor Value",
      y = "Annual Probability of Mortality",
      #title = "Effect of Predictors on Probability of Mortality"
    ) +
    species_fill + species_color + species_linetype +
    theme_bw(base_size = 14)+theme(panel.grid = element_blank())
  
  site.cond.marginal <- ggplot(marginal_response_df %>% filter( predictor.class %in% "Site Conditions"), aes(x = Value, y = 1-mean, color = Species, linetype = `Shade Tolerance`)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
    facet_grid(cols = vars(predictor.class), rows = vars(Predictor), scales = "free") +
    labs(
      x = "Predictor Value",
      y = "Annual Probability of Mortality",
      #title = "Effect of Predictors on Probability of Mortality"
    ) +
    species_fill + species_color + species_linetype +
    theme_bw(base_size = 14)+theme(panel.grid = element_blank())
  
  Climate.marginal <- ggplot(marginal_response_df %>% filter( predictor.class %in% "Climate"), aes(x = Value, y = 1-mean, color = Species, linetype = `Shade Tolerance`)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
    facet_grid(cols = vars(predictor.class), rows = vars(Predictor), scales = "free") +
    labs(
      x = "Predictor Value",
      y = "Annual Probability of Mortality",
      #title = "Effect of Predictors on Probability of Mortality"
    ) +
    species_fill + species_color + species_linetype +
    theme_bw(base_size = 14)+theme(panel.grid = element_blank())
  
  
  
  species.legend <- cowplot::get_legend(growth.marginal)
  
  cowplot::plot_grid(cowplot::plot_grid(growth.marginal + theme(legend.position = "none"), 
                                        Competition.marginal + theme(legend.position = "none"), 
                                        site.cond.marginal + theme(legend.position = "none"), 
                                        Climate.marginal + theme(legend.position = "none"), align = "hv"), species.legend, rel_widths = c(1,0.2))
  
  ggsave(height = 12, width = 12, 
         paste0(output.folder, "SPCD_stanoutput_joint_v3/predicted_mort/summary/model_",model.no,"_species_level_model_marginal_main_effects_summary.png"), 
         dpi = 550)
  
  # make one for the marginal effects, but without balsam fir:
  ##############################################################################################
  # make one giant figure for the marginal species effects:
  # generate figures for all of these separately 
  marginal_response_nobf <- marginal_response_df %>% filter(!Species %in% "balsam fir")
  growth.marginal <- ggplot(marginal_response_nobf %>% filter( predictor.class %in% "Growth & Size"), aes(x = Value, y = 1-mean, color = Species, linetype = `Shade Tolerance`)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
    facet_grid(cols = vars(predictor.class), rows = vars(Predictor), scales = "free") +
    labs(
      x = "Predictor Value",
      y = "Annual Probability of Mortality",
      #title = "Effect of Predictors on Probability of Mortality"
    ) +
    species_fill + species_color + species_linetype +
    theme_bw(base_size = 12)+theme(panel.grid = element_blank())
  
  Competition.marginal <- ggplot(marginal_response_nobf %>% filter( predictor.class %in% "Competition"), aes(x = Value, y = 1-mean, color = Species, linetype = `Shade Tolerance`)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
    facet_grid(cols = vars(predictor.class), rows = vars(Predictor), scales = "free") +
    labs(
      x = "Predictor Value",
      y = "Annual Probability of Mortality",
      #title = "Effect of Predictors on Probability of Mortality"
    ) +
    species_fill + species_color + species_linetype +
    theme_bw(base_size = 12)+theme(panel.grid = element_blank())
  
  site.cond.marginal <- ggplot(marginal_response_nobf %>% filter( predictor.class %in% "Site Conditions"), aes(x = Value, y = 1-mean, color = Species, linetype = `Shade Tolerance`)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
    facet_grid(cols = vars(predictor.class), rows = vars(Predictor), scales = "free") +
    labs(
      x = "Predictor Value",
      y = "Annual Probability of Mortality",
      #title = "Effect of Predictors on Probability of Mortality"
    ) +
    species_fill + species_color + species_linetype +
    theme_bw(base_size = 12)+theme(panel.grid = element_blank())
  
  Climate.marginal <- ggplot(marginal_response_nobf %>% filter( predictor.class %in% "Climate"), aes(x = Value, y = 1-mean, color = Species, linetype = `Shade Tolerance`)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
    facet_grid(cols = vars(predictor.class), rows = vars(Predictor), scales = "free") +
    labs(
      x = "Predictor Value",
      y = "Annual Probability of Mortality",
      #title = "Effect of Predictors on Probability of Mortality"
    ) +
    species_fill + species_color + species_linetype +
    theme_bw(base_size = 12)+theme(panel.grid = element_blank())
  
  
  
  species.legend <- cowplot::get_legend(growth.marginal)
  
  cowplot::plot_grid(cowplot::plot_grid(growth.marginal + theme(legend.position = "none"), 
                                        Competition.marginal + theme(legend.position = "none"), 
                                        site.cond.marginal + theme(legend.position = "none"), 
                                        Climate.marginal + theme(legend.position = "none"), align = "hv"), species.legend, rel_widths = c(1,0.2))
  
  ggsave(height = 12, width = 12, 
         paste0(output.folder, "SPCD_stanoutput_joint_v3/predicted_mort/summary/model_",model.no,"_species_level_model_marginal_main_effects_summary_noBF.png"), 
         dpi = 350)
  
  
  ################################################################
  # finally, only plot effects significantly different from zero:
  ################################################################
  marginal_response_df
  
  significant.effects <- all.joint.betas %>% filter(significance %in% "significant") %>% 
    select(Covariate, Species) #%>% na.omit()
  
  
  sig_marginal_effects <- left_join(significant.effects, marginal_response_df)
  sig_marginal_effects_noBF <- sig_marginal_effects %>% filter(!Species %in% "balsam fir")
  
  growth.marginal <- ggplot(sig_marginal_effects %>% filter( predictor.class %in% "Growth & Size"), aes(x = Value, y = 1-mean, color = Species, linetype = `Shade Tolerance`)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
    facet_grid(cols = vars(predictor.class), rows = vars(Predictor), scales = "free") +
    labs(
      x = "Predictor Value",
      y = "Annual Probability of Mortality",
      #title = "Effect of Predictors on Probability of Mortality"
    ) +
    species_fill + species_color + species_linetype +
    theme_bw(base_size = 12)+theme(panel.grid = element_blank())
  
  Competition.marginal <- ggplot(sig_marginal_effects %>% filter( predictor.class %in% "Competition"), aes(x = Value, y = 1-mean, color = Species, linetype = `Shade Tolerance`)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
    facet_grid(cols = vars(predictor.class), rows = vars(Predictor), scales = "free") +
    labs(
      x = "Predictor Value",
      y = "Annual Probability of Mortality",
      #title = "Effect of Predictors on Probability of Mortality"
    ) +
    species_fill + species_color + species_linetype +
    theme_bw(base_size = 12)+theme(panel.grid = element_blank())
  
  site.cond.marginal <- ggplot(sig_marginal_effects %>% filter( predictor.class %in% "Site Conditions"), aes(x = Value, y = 1-mean, color = Species, linetype = `Shade Tolerance`)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
    facet_grid(cols = vars(predictor.class), rows = vars(Predictor), scales = "free") +
    labs(
      x = "Predictor Value",
      y = "Annual Probability of Mortality",
      #title = "Effect of Predictors on Probability of Mortality"
    ) +
    species_fill + species_color + species_linetype +
    theme_bw(base_size = 12)+theme(panel.grid = element_blank())
  
  Climate.marginal <- ggplot(sig_marginal_effects %>% filter( predictor.class %in% "Climate"), aes(x = Value, y = 1-mean, color = Species, linetype = `Shade Tolerance`)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
    facet_grid(cols = vars(predictor.class), rows = vars(Predictor), scales = "free") +
    labs(
      x = "Predictor Value",
      y = "Annual Probability of Mortality",
      #title = "Effect of Predictors on Probability of Mortality"
    ) +
    species_fill + species_color + species_linetype +
    theme_bw(base_size = 12)+theme(panel.grid = element_blank())
  
  
  
  species.legend <- cowplot::get_legend(growth.marginal)
  
  cowplot::plot_grid(cowplot::plot_grid(growth.marginal + theme(legend.position = "none"), 
                                        Competition.marginal + theme(legend.position = "none"), 
                                        site.cond.marginal + theme(legend.position = "none"), 
                                        Climate.marginal + theme(legend.position = "none"), align = "hv"), species.legend, rel_widths = c(1,0.2))
  
  ggsave(height = 12, width = 12, 
         paste0(output.folder, "SPCD_stanoutput_joint_v3/predicted_mort/summary/model_",model.no,"_joint_marginal_main_effects_summary_sig_only.png"), 
         dpi = 350)
  
  
  # do the same but without balsam fir:
  growth.marginal <- ggplot(sig_marginal_effects_noBF %>% filter( predictor.class %in% "Growth & Size"), aes(x = Value, y = 1-mean, color = Species, linetype = `Shade Tolerance`)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
    facet_grid(cols = vars(predictor.class), rows = vars(Predictor), scales = "free") +
    labs(
      x = "Predictor Value",
      y = "Annual Probability of Mortality",
      #title = "Effect of Predictors on Probability of Mortality"
    ) +
    species_fill + species_color + species_linetype +
    theme_bw(base_size = 12)+theme(panel.grid = element_blank())
  
  Competition.marginal <- ggplot(sig_marginal_effects_noBF %>% filter( predictor.class %in% "Competition"), aes(x = Value, y = 1-mean, color = Species, linetype = `Shade Tolerance`)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
    facet_grid(cols = vars(predictor.class), rows = vars(Predictor), scales = "free") +
    labs(
      x = "Predictor Value",
      y = "Annual Probability of Mortality",
      #title = "Effect of Predictors on Probability of Mortality"
    ) +
    species_fill + species_color + species_linetype +
    theme_bw(base_size = 12)+theme(panel.grid = element_blank())
  
  site.cond.marginal <- ggplot(sig_marginal_effects_noBF %>% filter( predictor.class %in% "Site Conditions"), aes(x = Value, y = 1-mean, color = Species, linetype = `Shade Tolerance`)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
    facet_grid(cols = vars(predictor.class), rows = vars(Predictor), scales = "free") +
    labs(
      x = "Predictor Value",
      y = "Annual Probability of Mortality",
      #title = "Effect of Predictors on Probability of Mortality"
    ) +
    species_fill + species_color + species_linetype +
    theme_bw(base_size = 12)+theme(panel.grid = element_blank())
  
  Climate.marginal <- ggplot(sig_marginal_effects_noBF %>% filter( predictor.class %in% "Climate"), aes(x = Value, y = 1-mean, color = Species, linetype = `Shade Tolerance`)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
    facet_grid(cols = vars(predictor.class), rows = vars(Predictor), scales = "free") +
    labs(
      x = "Predictor Value",
      y = "Annual Probability of Mortality",
      #title = "Effect of Predictors on Probability of Mortality"
    ) +
    species_fill + species_color + species_linetype +
    theme_bw(base_size = 12)+theme(panel.grid = element_blank())
  
  
  
  species.legend <- cowplot::get_legend(growth.marginal)
  
  cowplot::plot_grid(cowplot::plot_grid(growth.marginal + theme(legend.position = "none"), 
                                        Competition.marginal + theme(legend.position = "none"), 
                                        site.cond.marginal + theme(legend.position = "none"), 
                                        Climate.marginal + theme(legend.position = "none"), align = "hv"), species.legend, rel_widths = c(1,0.2))
  
  ggsave(height = 12, width = 12, 
         paste0(output.folder, "SPCD_stanoutput_joint_v3/predicted_mort/summary/model_",model.no,"_joint_marginal_main_effects_summary_sig_only_noBF.png"), 
         dpi = 350)
  
  
  #########################################################################
  # Marginal effects, but plotting effects for species in affected by different pests
  
  marginal_response_df$predictor.class2
  interaction_response_df<- left_join(interaction_response_df, Covariate.types.df)
  
  spruce.fir <- c("balsam fir", "red spruce", "paper birch", "northern white-cedar")
  mixed <- c("American beech", "eastern hemlock", "yellow birch")
  oak.hickory <- c("chestnut oak", "white oak", "hickory", "northern red oak")
  spongy.immune <- c("white ash", "yellow-poplar", "black cherry")
  spongy.resist <- c("sugar maple", "red maple", "eastern white pine")
  
  
  ggplot(marginal_response_df %>% filter( Species %in% "northern white-cedar" &
                                            predictor.class2 %in% c("N deposition", 
                                                                    "Ndep x G & S")), aes(x = Value, y = 1-mean, color = Species, linetype = Species)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
    facet_wrap(~Predictor) +
    labs(
      x = "Predictor Value",
      y = "Annual Probability of Mortality",
      #title = "Effect of Predictors on Probability of Mortality"
    ) +
    species_fill + species_color + named_species_linetype +
    theme_bw(base_size = 14)+theme(panel.grid = element_blank(), 
                                   legend.key.size = unit(1, "cm")  
    )
  
  
  # important effects for american beech:
  
  ggplot(marginal_response_df %>% filter( Species %in% "American beech" &
                                            pred.type %in% "Main Effects"), aes(x = Value, y = 1-mean, color = Species, linetype = Species)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
    facet_wrap(~Predictor) +
    labs(
      x = "Predictor Value",
      y = "Annual Probability of Mortality",
      #title = "Effect of Predictors on Probability of Mortality"
    ) +
    species_fill + species_color + named_species_linetype +
    theme_bw(base_size = 14)+theme(panel.grid = element_blank(), 
                                   legend.key.size = unit(1, "cm")  
    )
  
  # beech interaction plots
  # important main effects to show:
  # Tmax
  # Dia
  # Dia diff
  # bal.scaled
  beech.main <- c("DIA_DIFF_scaled",
                  "DIA_scaled", 
                  "MATmax.scaled", 
                  "BAL.scaled")
  # important interactions to show:
  # Tmax x dia
  # Tmax x dia diff
  beech.interactions <- c("MATmax.scaled_growth.int", 
                          "MATmax.scaled_DIA.int",
                          "DIA_scaled_growth.int", 
                          "ba.scaled_growth.int")
  
  species <- "American beech"
  covar <- "MATmax.scaled_growth.int"
  color.pred.class.2 <- c(
    
    "Size" = "darkgreen",
    "Change in Size" = "#66a61e",
    "Growth x Size" = "#a1d99b" ,
    
    "Climate" = "darkblue",
    "Climate x G & S" = "#1d91c0",
    
    "N deposition" = "#67000d" , 
    "Ndep x G & S"= "#bd0026",
    
    "% Damage" = "#762a83", 
    "Damage x G & S" = "#9970ab",
    
    
    "Competition" = "#8c510a",
    "Competition x G & S" = "#d95f02" ,
    
    "Site Conditions" = "#e6ab02",
    "Site x G & S" = "yellow3")
  
  var.part.fill <- scale_fill_manual(values = color.pred.class.2, name = "",
                                     labels = c("Change in Size" = expression(Delta ~ "D" ), 
                                                "Growth x Size" = expression(Delta ~ "D x D"), 
                                                "Size" = "Diameter (D)", 
                                                "Climate x G & S" = expression("Climate x " ~Delta ~ "D or D"), 
                                                "Ndep x G & S" = expression("N dep. x " ~Delta ~ "D or D"),
                                                "Damage x G & S" = expression("Damage x " ~Delta ~ "D or D"), 
                                                "Competition x G & S" = expression("Compeition x " ~Delta ~ "D or D"),
                                                "Site x G & S" = expression("Site x " ~Delta ~ "D or D")
                                     ))
  
  var.part.color <- scale_color_manual(values = color.pred.class.2, name = "",
                                       labels = c("Change in Size" = expression(Delta ~ "D" ), 
                                                  "Growth x Size" = expression(Delta ~ "D x D"), 
                                                  "Size" = "Diameter (D)", 
                                                  "Climate x G & S" = expression("Climate x " ~Delta ~ "D or D"), 
                                                  "Ndep x G & S" = expression("N dep. x " ~Delta ~ "D or D"),
                                                  "Damage x G & S" = expression("Damage x " ~Delta ~ "D or D"), 
                                                  "Competition x G & S" = expression("Compeition x " ~Delta ~ "D or D"),
                                                  "Site x G & S" = expression("Site x " ~Delta ~ "D or D")
                                       ))
  
  ggplot(marginal_response_df %>% filter( Species %in% "American beech" &
                                            pred.type %in% "Main Effects"), aes(x = Value, y = 1-mean, color = Species, linetype = Species)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
    facet_wrap(~Predictor) +
    labs(
      x = "Predictor Value",
      y = "Annual Probability of Mortality",
      #title = "Effect of Predictors on Probability of Mortality"
    ) +
    species_fill + species_color + named_species_linetype +
    theme_bw(base_size = 14)+theme(panel.grid = element_blank(), 
                                   legend.key.size = unit(1, "cm")  
    )
  
  
  # function to plot out the case study species:
  plot.main.interaction.effect <- function(species, covar){
    marge <- marginal_response_df %>% filter(Species %in% species)
    inter <- interaction_response_df %>% filter(Species %in% species)
    ymax.spp <- max(max(1-inter$ci.hi.90), 1- marge$ci.hi)
    
    if(!covar %in% interaction_response_df$covariate){
      df.species <- marginal_response_df %>% filter(Species %in% species)%>%
        filter(Covariate %in% covar)
      unique(df.species$predictor.class2)
      strip.fill <- as.character(color.pred.class.2[unique(df.species$predictor.class2)])
      
      
      p1 <-  ggplot(df.species , aes(x = Value, y = 1-mean, color = predictor.class2)) +
        geom_line(size = 1) +
        geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = predictor.class2), alpha = 0.2, color = NA) +
        facet_wrap(~Predictor) +
        labs(
          x = paste("Scaled", df.species$Predictor),
          y = "Annual Probability of Mortality",
          #title = "Effect of Predictors on Probability of Mortality"
        ) +
        var.part.fill + var.part.color+
        theme_bw(base_size = 14)+theme(panel.grid = element_blank(), 
                                       legend.key.size = unit(1, "cm"), 
                                       legend.position = "none", 
                                       strip.background = element_rect(fill = strip.fill)
        )+ylim(0, ymax.spp)
    }else{
      
      
      df.species <- interaction_response_df %>% filter(Species %in% species)%>%
        filter(covariate %in% covar)
      df.species$p2.rank <- factor(df.species$p2.rank, levels = c("high", "median", "low"))
      strip.fill <- as.character(color.pred.class.2[unique(df.species$predictor.class2)])
      pred.name <-  marginal_response_df %>% filter(Covariate %in% unique(df.species$Pred.1))%>% select(Predictor)%>% distinct()
      
      if(unique(df.species$Pred.2) %in% "DIA_DIFF_scaled"){
        
        p1 <-  ggplot(df.species, aes(x = p1.value, y = 1-mean, color = p2.rank)) +
          geom_line(size = 1) +
          geom_ribbon(aes(ymin = 1-ci.lo.10, ymax = 1-ci.hi.90, fill = p2.rank), alpha = 0.2, color = NA) +
          #facet_grid(cols = vars(predictor.class2), rows = vars(Predictor)) +
          facet_wrap(~Predictor)+
          labs(
            x = paste("Scaled", pred.name$Predictor),
            y = "Annual Probability of Mortality",
            #title = "Effect of Predictors on Probability of Mortality"
          ) + scale_color_manual(values = c("low"="#a1dab4" ,
                                            "median" = "#41b6c4",
                                            "high"= "#225ea8" ), 
                                 name = expression(Delta ~ "Diameter"))+
          scale_fill_manual(values = c("low"="#a1dab4" ,
                                       "median" = "#41b6c4",
                                       "high"= "#225ea8" ), 
                            name = expression(Delta ~ "Diameter"))+
          
          #species_fill + species_color + 
          #named_species_linetype +
          theme_bw(base_size = 14)+theme(panel.grid = element_blank(), 
                                         legend.key.size = unit(1, "cm"), 
                                         strip.background = element_rect(fill = strip.fill)
          )+ylim(0, ymax.spp)
      }else{
        
        p1 <-  ggplot(df.species, aes(x = p1.value, y = 1-mean, color = p2.rank)) +
          geom_line(size = 1) +
          geom_ribbon(aes(ymin = 1-ci.lo.10, ymax = 1-ci.hi.90, fill = p2.rank), alpha = 0.2, color = NA) +
          #facet_grid(cols = vars(predictor.class2), rows = vars(Predictor)) +
          facet_wrap(~Predictor)+
          labs(
            x = paste("Scaled", pred.name$Predictor),
            y = "Annual Probability of Mortality",
            #title = "Effect of Predictors on Probability of Mortality"
          ) + scale_color_manual(values = c("low"="#fbb4b9" ,
                                            "median" = "#f768a1",
                                            "high"= "#ae017e" ), 
                                 name =  "Diameter", 
                                 labels = c("low" = "small", 
                                            "median" = "medium", 
                                            "high" = "large"))+
          
          
          scale_fill_manual(values = c("low"="#fbb4b9" ,
                                       "median" = "#f768a1",
                                       "high"= "#ae017e" ), 
                            name = "Diameter", 
                            labels = c("low" = "small", 
                                       "median" = "medium", 
                                       "high" = "large"))+
          
          #species_fill + species_color + 
          #named_species_linetype +
          theme_bw(base_size = 14)+theme(panel.grid = element_blank(), 
                                         legend.key.size = unit(1, "cm"),
                                         strip.background = element_rect(fill = strip.fill)
          )+ylim(0, ymax.spp)
        
        
      }
    }
    
    return(p1)
  }
  
  a <- plot.main.interaction.effect(species = "American beech", 
                                    covar = "DIA_scaled")
  b<- plot.main.interaction.effect(species = "American beech", 
                                   covar = "DIA_DIFF_scaled")+
    ylim(0, 0.035)
  c<-plot.main.interaction.effect(species = "American beech", 
                                  covar = "MATmax.scaled")+
    ylim(0, 0.035)
  
  d<-plot.main.interaction.effect(species = "American beech", 
                                  covar = "MATmax.scaled_growth.int")+
    ylim(0, 0.03)
  
  e<- plot.main.interaction.effect(species = "American beech", 
                                   covar = "MATmax.scaled_DIA.int")+
    ylim(0, 0.035)
  
  f<-plot.main.interaction.effect(species = "American beech", 
                                  covar = "DIA_scaled_growth.int")+
    ylim(0, 0.035)
  
  
  beech.top <- plot_grid(a,b,f, c,d,e, ncol = 3, align = "hv")
  save_plot(paste0(output.folder,"images/beech_top_conditionals.svg"), beech.top, base_width = 12, base_height = 6) 
  save_plot(paste0(output.folder,"images/beech_top_conditionals.png"), beech.top, base_width = 12, base_height = 8) 
  
  # make plots of all the species inteacation terms:
  unique(marginal_response_df$Covariate)
  spec <- "balsam fir"
  plt.all.spp.conditional <- function(spec){
    spp.plt.list <- list()
    for(h in 1:12){
      spp.plt.list[[h]] <- plot.main.interaction.effect(species = spec, 
                                                        covar = unique(marginal_response_df$Covariate)[h])
    }
    
    
    main.effects <- plot_grid(plotlist = spp.plt.list, align = "hv")
    save_plot(paste0(output.folder,"images/",spec,"_main_effects_conditionals.png"), 
              main.effects, base_width = 12, base_height = 12) 
    save_plot(paste0(output.folder,"images/",spec,"_main_effects_conditionals.svg"), 
              main.effects, base_width = 12, base_height = 12) 
    
    
    spp.plt.int.list <- list()
    for(h in 13:33){
      spp.plt.int.list[[h-12]] <- plot.main.interaction.effect(species = spec, 
                                                               covar = unique(marginal_response_df$Covariate)[h])
    }
    
    
    int.effects <- plot_grid(plotlist =spp.plt.int.list , align = "hv")
    save_plot(paste0(output.folder,"images/",spec,"_interaction_conditionals.png"), 
              int.effects, base_width = 24, base_height = 24) 
    
    save_plot(paste0(output.folder,"images/",spec,"_interaction_conditionals.svg"), 
              int.effects, base_width = 24, base_height = 24) 
  }
  
  lapply(unique(marginal_response_df$Species), plt.all.spp.conditional)  
  
  # now get the top responses for a few species:
  
  
  
}
#}

# read in all of the AUC outputs for each model and make one big plot and one big summary file:
AUC.oos.list <- AUC.is.list <- list()
for(model.no in 1:9){
  
  AUC.oos.df <-  readRDS(paste0(input.folder, "SPCD_stanoutput_joint_v3_model_", model.no, "/AUC_oos_with_uncertainty.rds"))
  AUC.oos.df$Model <- paste0("Model ", model.no)
  AUC.oos.df$AUC.type <- "Out-of-sample"
  
  AUC.is.df <-  readRDS(paste0(input.folder, "SPCD_stanoutput_joint_v3_model_", model.no, "/AUC_is_with_uncertainty.rds"))
  AUC.is.df$Model <- paste0("Model ", model.no)
  AUC.is.df$AUC.type <- "In-sample"
  
  AUC.oos.list[[model.no]]<- AUC.oos.df
  AUC.is.list[[model.no]]<-  AUC.is.df
  
}

summary.spp.table <- rbind(main.table %>% select(-Species), spp.table)
AUC.oos.df <- do.call(rbind, AUC.oos.list) %>% left_join(., summary.spp.table)
AUC.is.df <- do.call(rbind, AUC.is.list) %>% left_join(., summary.spp.table)
AUC.is.df$COMMON <- factor(AUC.is.df$COMMON, unique(summary.spp.table$COMMON))
AUC.oos.df$COMMON <- factor(AUC.is.df$COMMON, unique(summary.spp.table$COMMON))


is.auc <- ggplot(AUC.is.df, aes(x = Model, y = median, shape = Model %in% "Model 6"))+geom_point()+
  geom_errorbar(data = AUC.is.df, aes(x = Model, ymin = auc.ci.lo, ymax = auc.ci.hi))+
  facet_wrap(~COMMON, scales =  "free_y")+
  scale_color_manual( values = c( "black" , 
                                  "red" ))+
  theme_bw(base_size = 12)+
  theme(axis.text.x = element_text(angle = 45,  hjust=1), legend.position = "none", 
        panel.grid = element_blank()) +
  xlab("")+ylab("In Sample AUC")



ggsave(paste0(output.folder, "/model_summary_full/All_joint_model_all9models_compare-auc-insample_joint_species.png"), 
       is.auc,
       width = 10, height = 6, 
       dpi = 350)


oos.auc <- ggplot(AUC.is.df, aes(x = Model, y = median, shape = Model %in% "Model 6"))+geom_point()+
  geom_errorbar(data = AUC.is.df, aes(x = Model, ymin = auc.ci.lo, ymax = auc.ci.hi))+
  facet_wrap(~COMMON, scales =  "free_y")+
  scale_color_manual( values = c( "black" , 
                                  "red" ))+
  theme_bw(base_size = 12)+
  theme(axis.text.x = element_text(angle = 45,  hjust=1), legend.position = "none", 
        panel.grid = element_blank()) +
  xlab("")+ylab("Out of Sample AUC")

ggsave(paste0(output.folder, "/model_summary_full/All_joint_model_all9models_compare-auc-outofsample_joint_species.png"), 
       oos.auc,
       width = 10, height = 6)
##############
# same plots but without all the free axes
is.auc <- ggplot(AUC.is.df, aes(x = Model, y = median, shape = Model %in% "Model 6"))+geom_point()+
  geom_errorbar(data = AUC.is.df, aes(x = Model, ymin = auc.ci.lo, ymax = auc.ci.hi))+
  facet_wrap(~COMMON)+
  scale_color_manual( values = c( "black" , 
                                  "red" ))+
  theme_bw(base_size = 12)+
  theme(axis.text.x = element_text(angle = 45,  hjust=1), legend.position = "none", 
        panel.grid = element_blank()) +
  xlab("")+ylab("In Sample AUC")



ggsave(paste0(output.folder, "/model_summary_full/All_joint_model_all9models_compare-auc-insample_joint_species_same_scale.png"), 
       is.auc,
       width = 10, height = 6, 
       dpi = 350)


oos.auc <- ggplot(AUC.is.df, aes(x = Model, y = median, shape = Model %in% "Model 6"))+geom_point()+
  geom_errorbar(data = AUC.is.df, aes(x = Model, ymin = auc.ci.lo, ymax = auc.ci.hi))+
  facet_wrap(~COMMON)+
  scale_color_manual( values = c( "black" , 
                                  "red" ))+
  theme_bw(base_size = 12)+
  theme(axis.text.x = element_text(angle = 45,  hjust=1), legend.position = "none", 
        panel.grid = element_blank()) +
  xlab("")+ylab("Out of Sample AUC")

ggsave(paste0(output.folder, "/model_summary_full/All_joint_model_all9models_compare-auc-outofsample_joint_species_same_scale.png"), 
       oos.auc,
       width = 10, height = 6)

saveRDS(AUC.is.df, paste0(output.folder, "/model_summary_full/AUC_is_all_models_joint.rds"))
saveRDS(AUC.oos.df, paste0(output.folder, "/model_summary_full/AUC_is_all_models_joint.rds"))



########################################################################
# also plot the computational resources
########################################################################
# read in all of the core.hours outputs for each model and make one big plot and one big summary file:
core.hours.list <- list()
for(model.no in 1:9){
  
  core.hours.oos.df <-  read.csv(paste0(input.folder, "SPCD_stanoutput_joint_v3_model_", model.no, "/joint_model_time_diag_SPCD_joint_model_",model.no, "_remper_0.5.csv"))
  core.hours.oos.df$Model <- paste0("Model ", model.no)
  #core.hours.oos.df$core.hours.type <- "Out-of-sample"
  
  
  
  core.hours.list[[model.no]]<- core.hours.oos.df
  
  
}

summary.spp.table <- rbind(main.table %>% select(-Species), spp.table)
core.hours.df <- do.call(rbind, core.hours.list) %>% left_join(., summary.spp.table)
core.hours.df$COMMON <- FIESTA::ref_species[match(core.hours.df$SPCD.id, FIESTA::ref_species$SPCD),]$COMMON
core.hours.df$COMMON<- factor(core.hours.df$COMMON, unique(summary.spp.table$COMMON))
core.hours.df$Species <- FIESTA::ref_species[match(core.hours.df$SPCD.id, FIESTA::ref_species$SPCD),]$Species

#colnames(core.hours.df)[13] <- "Species"
# plot up the core hours vs in sample and out of sample AUC scores:
AUC.oos.time <- left_join(core.hours.df %>% 
                            select(n..:cores, cores:Model)%>% 
                            rename("SPCD" = "SPCD.id"), AUC.oos.df)


AUC.is.time <- left_join(core.hours.df %>% 
                           select(n..:cores, cores:Model)%>% 
                           rename("SPCD" = "SPCD.id"), AUC.is.df)

ggplot()+geom_text(data = AUC.oos.time, aes( x = core.hours, y =  median.oos, label = model ))+
  facet_wrap(~COMMON, scales = "free")+
  theme_bw()+theme(panel.grid = element_blank())+ylab("Out-of-sample AUC")+xlab("Core Hours")
ggsave(height = 7, width = 10, dpi = 350, 
       paste0(output.folder, "/model_summary_full/core_hours_oosAUC_all_models_joint.png"))

ggplot()+geom_text(data = AUC.oos.time, aes( x = core.hours, y =  median.oos, label = model ))+
  facet_wrap(~COMMON, scales = "free_x")+
  theme_bw()+theme(panel.grid = element_blank())+ylab("Out-of-sample AUC")+xlab("Core Hours")

ggsave(height = 7, width = 10, dpi = 350, 
       paste0(output.folder, "/model_summary_full/core_hours_oosAUC_all_models_joint_common_y.png"))


ggplot()+geom_text(data = AUC.is.time, aes( x = core.hours, y =  median, label = model ))+
  facet_wrap(~COMMON, scales = "free")+
  theme_bw()+theme(panel.grid = element_blank())+ylab("In-sample AUC")+xlab("Core Hours")
ggsave(height = 7, width = 10, dpi = 350, 
       paste0(output.folder, "/model_summary_full/core_hours_isAUC_all_models_joint.png"))

ggplot()+geom_text(data = AUC.is.time, aes( x = core.hours, y =  median, label = model ))+
  facet_wrap(~COMMON, scales = "free_x")+
  theme_bw()+theme(panel.grid = element_blank())+ylab("In-sample AUC")+xlab("Core Hours")

ggsave(height = 7, width = 10, dpi = 350, 
       paste0(output.folder, "/model_summary_full/core_hours_isAUC_all_models_joint_common_y.png"))


saveRDS(AUC.is.time, paste0(output.folder, "/model_summary_full/core_hours_isAUC_all_models_joint.RDS"))
saveRDS(AUC.oos.time, paste0(output.folder, "/model_summary_full/core_hours_oosAUC_all_models_joint.RDS"))

# save all outputs to the cyverse output directory now:
file.rename(
  "SPCD_stanoutput_joint_v3",
  "joint_model_summary_v3"
)


# copy to the data-store output
system(paste(
  "cp -r",
  "joint_model_summary_v3",
  "data-store/data/output/"
))


system(paste(
  "cp -r",
  "model_summary_full",
  "data-store/data/output/"
))


# file.rename(
#   "joint_model_summary_v3",
#   "SPCD_stanoutput_joint_v3"
#   
# )
