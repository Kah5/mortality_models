##############################################################################
# Script to take the posterior beta estimates and generate a figure of effects
##############################################################################
library(rstan)
library(tidyverse)
library(posterior)
library(FIESTA)
################################################################################
source("R/modelAssessment/connect2cyverse.R") # has folder location for cyverse data & output folder information

# set up species dataframe
nspp <- data.frame(SPCD = c(316, 318, 833, 832, 261, 531, 802, 129, 762,  12, 541,  97, 621, 400, 371, 241, 375))
nspp$Species <- paste(FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD),]$GENUS, FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD),]$SPECIES)
nspp$COMMON <- paste(FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD),]$GENUS, FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD),]$SPECIES)

# model 6 marginal effects:
model.no <- 6

# balsam fir model 6 samples are missing--it showsup on cyverse, but when I download it has 0 KB size?

# set up species.numbers to loop over
spp.num <- c(17:11, 9:1)

for(i in spp.num){
  cat(paste0("generating marginal predictions for all covaraites for model ", model.no," SPCD ",i, " ", nspp[i,]$Species, "\n"))
  
SPCD.id <- nspp[i,]$SPCD
#output.folder = "C:/Users/KellyHeilman/Box/01. kelly.heilman Workspace/mortality/Eastern-Mortality/mortality_models/"

# Extract posterior samples for the model and species
fit <- readRDS(paste0(output.folder,"/SPCD_stanoutput_full_standardized_v3/samples/model_",model.no,"_SPCD_",SPCD.id, "_remper_correction_0.5.RDS"))
#fit <- readRDS(url(paste0(cyverse.folder,"/SPCD_stanoutput_full_standardized_v3/samples/model_",model.no,"_SPCD_",SPCD.id, "_remper_correction_0.5.RDS")))
#saveRDS(fit, paste0(output.folder,"/SPCD_stanoutput_full_standardized_v3/samples/model_",model.no,"_SPCD_",SPCD.id, "_remper_correction_0.5.RDS"))


fit_ssm_df <- as_draws_df(fit) # takes awhile to convert to df

# get all the covariates using posterior package
betas.quant <- subset_draws(fit_ssm_df, variable = "u_beta") %>% summarise_draws(median, ~quantile(., probs = c(0.025, 0.975))) %>%
  rename(`ci.lo` = "2.5%", `ci.hi` = "97.5%") %>% 
  mutate(remper.cor = 0.5, 
         SPCD = SPCD.id)

# Example parameters from the Stan model
# Replace 'beta' with your actual parameter names
# Assume 'beta_0' is the intercept and 'beta' are the slopes for covariates
beta_0 <- data.frame(subset_draws(fit_ssm_df, variable = "alpha_SPP")) %>% select(-.chain, -.iteration, -.draw)  # Intercept
beta <- data.frame(subset_draws(fit_ssm_df, variable = "u_beta")) %>% select(-.chain, -.iteration, -.draw)     # Covariates (matrix: iterations x covariates)

saveRDS(beta, paste0(output.folder,"/SPCD_stanoutput_full_standardized_v3/samples/model_",model.no, "_SPCD_", SPCD.id, "betas.RDS"))
saveRDS(beta_0, paste0(output.folder,"/SPCD_stanoutput_full_standardized_v3/samples/model_",model.no, "_SPCD_", SPCD.id, "alpha.RDS"))

# read in the species data and covariates to get the min and max and get ranges
load(paste0("SPCD_standata_general_full_standardized_v3/SPCD_",SPCD.id, "remper_correction_0.5model_",model.no, ".Rdata")) # load the species code data
mod.data$K <- ncol(mod.data$xM)

covariate_names <- c(colnames(mod.data$xM))  # Replace with your covariate names

var.mins <- as.vector(apply(data.frame(mod.data$xM),2 , function(x)quantile(x, 0.025)))
var.maxes <- as.vector(apply(data.frame(mod.data$xM),2 , function(x) quantile(x, 0.975)))

covariate_ranges_df <- data.frame(
  covariate  = covariate_names, 
  mins = var.mins, 
  maxes = var.maxes
)


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
#cov_name <- "aspect.scaled"


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
    linear_predictor <- beta_0 + rowSums(as.matrix(beta)%*% covariate_matrix[j,] )
  probabilities[[j]] <-  as.vector(inv_logit_fxn(linear_predictor)$alpha_SPP)
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
  
ggsave(height = 10, width = 14, 
       paste0(output.folder, "SPCD_stanoutput_full_standardized_v3/predicted_mort/model_",model.no,"_all_marginal_SPCD_", SPCD.id, ".png"))

# save for this species:
plot_data$SPCD <- SPCD.id
saveRDS(plot_data, paste0(output.folder, "SPCD_stanoutput_full_standardized_v3/predicted_mort/model_",model.no,"_all_marginal_SPCD_", SPCD.id, ".RDS") )

}


################################################################################
# combine all of the species conditional responses together
################################################################################
marginal_response <- list()
for (i in spp.num){
  SPCD.id <- nspp[i,]$SPCD
  marginal_response[[i]] <- readRDS( paste0(output.folder, "SPCD_stanoutput_full_standardized_v3/predicted_mort/model_",model.no,"_all_marginal_SPCD_", SPCD.id, ".RDS") )
  
}
marginal_response_df <- do.call(rbind, marginal_response)

# get species common names
marginal_response_df$Species <- FIESTA::ref_species[match(marginal_response_df$SPCD, FIESTA::ref_species$SPCD),]$COMMON

Covariate.types.df <- data.frame(Covariate = unique(marginal_response_df$Covariate), 
          Predictor = c("Diameter Difference \n (DIA_DIFF)", 
                        "Diameter", 
                        "Plot Basal Area (BA)", 
                        "Basal Area Larger \n (BAL)", 
                        "% Damage", 
                        "MATmax", 
                        "MAP", 
                        "Precip. Anomaly", 
                        "Max. Temp. Anomaly", 
                        "Slope", 
                        "Aspect", 
                        "N Deposition", 
                        
                        # diameter difference interactions
                        "Diameter x DIA_DIFF", 
                        "BA x DIA_DIFF", 
                        "BAL x DIA_DIFF", 
                        "% Damage x DIA_DIFF", 
                        "MATmax x DIA_DIFF", 
                        "MAP x DIA_DIFF", 
                        "Precip. Anomaly \n x DIA_DIFF", 
                        "Max. Temp. Anomaly \n x DIA_DIFF", 
                        "Slope x DIA_DIFF", 
                        "Aspect x DIA_DIFF", 
                        "N Dep. x DIA_DIFF",
                        
                        # interactions with Diameter
                        
                        "BA x Diameter", 
                        "BAL x Diameter", 
                        "% Damage x Diameter", 
                        "MATmax x Diameter", 
                        "MAP x Diameter", 
                        "Precip. Anomaly x \n Diameter", 
                        "Max. Temp. Anomaly x \n Diameter", 
                        "Slope x Diameter", 
                        "Aspect x Diameter", 
                        "N Dep. x Diameter" ), 
          # whether it is a main effect of interaction
          pred.type = c(rep("Main Effects", 12), rep("DIA_DIFF interactions", 11), rep("Diameter interactions", 10)), 
          
          # what type of variable it relates to
          predictor.class = c("Growth & Size", 
                              "Growth & Size", 
                              "Competition", 
                              "Competition", 
                              "Competition", 
                              "Climate", 
                              "Climate", 
                              "Climate", 
                              "Climate", 
                              "Site Conditions", 
                              "Site Conditions", 
                              "Site Conditions", 
                              
                              
                              "Growth & Size",
                              "Competiton x G & S", 
                              "Competiton x G & S",
                              "Competiton x G & S", 
                              "Climate x G & S",
                              "Climate x G & S",
                              "Climate x G & S",
                              "Climate x G & S",
                              "Site x G & S",
                              "Site x G & S",
                              "Site x G & S",
                              
                              
                              
                              "Competiton x G & S", 
                              "Competiton x G & S",
                              "Competiton x G & S",
                              "Climate x G & S",
                              "Climate x G & S",
                              "Climate x G & S",
                              "Climate x G & S",
                              "Site x G & S",
                              "Site x G & S",
                              "Site x G & S"
                              ))

# join to new variable names and make sure they are in order
marginal_response_df <- marginal_response_df %>% left_join(., Covariate.types.df)
marginal_response_df$Predictor <- factor(marginal_response_df$Predictor, levels = unique(marginal_response_df$Predictor))


ggplot(marginal_response_df, aes(x = Value, y = 1-mean, color = Species)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
  facet_wrap(~ Predictor, scales = "free") +
  labs(
    x = "Scaled Covariate Value",
    y = "Annual Probability of Mortality",
    title = "Effect of Covariates on Probability of Mortality"
  ) +
  theme_bw(base_size = 12)+theme(panel.grid = element_blank()) 
ggsave(height = 12, width = 14, paste0(output.folder, "SPCD_stanoutput_full_standardized_v3/images/predicted_mortality/model_",model.no,"_species_level_models_all_species_marginal_effects.png"))



# generate figures for all of these separately 
ggplot(marginal_response_df %>% 
         filter( predictor.class %in% c("Growth & Size", "Competition", "Climate", "Site Conditions") & 
                   !Predictor %in% "Diameter x DIA_DIFF"), aes(x = Value, y = 1-mean, color = Species)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
  facet_wrap(~ Predictor, scales = "free") +
  labs(
    x = "Scaled Covariate Value",
    y = "Annual Probability of Mortality",
    title = "Effect of Covariates on Probability of Mortality"
  ) +
  theme_bw(base_size = 16)+theme(panel.grid = element_blank())
ggsave(height = 8.5, width = 12.5, 
       paste0(output.folder,  "SPCD_stanoutput_full_standardized_v3/images/predicted_mortality/model_",model.no,"_species_level_model_marginal_main_effects.png"))

# balsam fir really domninates the graph so take it out to view
ggplot(marginal_response_df %>% 
         filter( predictor.class %in% c("Growth & Size", "Competition", "Climate", "Site Conditions") & 
                   !Predictor %in% "Diameter x DIA_DIFF" & !Species %in% "balsam fir"), aes(x = Value, y = 1-mean, color = Species)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
  facet_wrap(~ Predictor, scales = "free") +
  labs(
    x = "Scaled Covariate Value",
    y = "Annual Probability of Mortality",
    title = "Effect of Covariates on Probability of Mortality"
  ) +
  theme_bw(base_size = 12)+theme(panel.grid = element_blank())
ggsave(height = 8.5, width = 12.5, 
       paste0(output.folder, "SPCD_stanoutput_full_standardized_v3/images/predicted_mortality/model_6_species_level_model_marginal_main_effects_no_balsam_fir.png"))


ggplot(marginal_response_df %>% filter( predictor.class %in% "Growth & Size" ), aes(x = Value, y = 1-mean, color = Species)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
  facet_wrap(~ Predictor, scales = "free") +
  labs(
    x = "Scaled Covariate Value",
    y = "Annual Probability of Mortality",
    title = "Growth & Size effects on Probability of Mortality"
  ) +
  theme_bw(base_size = 16)+theme(panel.grid = element_blank())
ggsave(height = 5, width = 12, 
       paste0(output.folder, "SPCD_stanoutput_full_standardized_v3/images/predicted_mortality/model_6_species_level_model_marginal_growth_diam.png"))

ggplot(marginal_response_df %>% filter( predictor.class %in% "Competition" ), aes(x = Value, y = 1-mean, color = Species)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
  facet_wrap(~ Predictor, scales = "free", ncol = 5) +
  labs(
    x = "Scaled Covariate Value",
    y = "Annual Probability of Mortality",
    title = "Competition effects on Probability of Mortality"
  ) +
  theme_bw(base_size = 16)+theme(panel.grid = element_blank())
ggsave(height = 5, width = 12, 
       paste0(output.folder, "SPCD_stanoutput_full_standardized_v3/images/predicted_mortality/model_6_species_level_model_marginal_competition.png"))


ggplot(marginal_response_df %>% filter(  predictor.class %in% "Climate" ), aes(x = Value, y = 1-mean, color = Species)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
  facet_wrap(~ Predictor, scales = "free", ncol = 2) +
  labs(
    x = "Scaled Covariate Value",
    y = "Annual Probability of Mortality",
    title = "Climate effects on Probability of Mortality"
  ) +
  theme_bw(base_size = 16)+theme(panel.grid = element_blank())
ggsave(height = 8, width = 12, 
       paste0(output.folder, "SPCD_stanoutput_full_standardized_v3/images/predicted_mortality/model_6_species_level_model_marginal_climate_vars.png"))


ggplot(marginal_response_df %>% filter( predictor.class %in% "Site Conditions" ), aes(x = Value, y = 1-mean, color = Species)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
  facet_wrap(~ Predictor, scales = "free", ncol = 5) +
  labs(
    x = "Scaled Covariate Value",
    y = "Annual Probability of Mortality",
    title = "Site Condition effects on Probability of Mortality"
  ) +
  theme_bw(base_size = 16)+theme(panel.grid = element_blank())
ggsave(height = 5, width = 12, 
       paste0(output.folder, "SPCD_stanoutput_full_standardized_v3/images/predicted_mortality/model_6_species_level_model_marginal_site_vars.png"))

### interaction terms
ggplot(marginal_response_df %>% filter( predictor.class %in% "Competiton x G & S" & !Species %in% "balsam fir"), aes(x = Value, y = 1-mean, color = Species)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
  facet_wrap(~ Predictor, scales = "free", ncol = 3) +
  labs(
    x = "Covariate Value",
     y = "Annual Probability of Mortality" #,
    # title = "Competition x G & S  on Probability of Mortality"
  ) +
  theme_bw(base_size = 16)+theme(panel.grid = element_blank())
ggsave(height = 5, width = 12, 
       paste0(output.folder, "SPCD_stanoutput_full_standardized_v3/images/predicted_mortality/model_6_species_level_model_marginal_growth_competition_noBF.png"))


ggplot(marginal_response_df %>% filter( predictor.class %in% "Climate x G & S" & !Species %in% "balsam fir" ), aes(x = Value, y = 1-mean, color = Species)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
  facet_wrap(~ Predictor, scales = "free", ncol = 4) +
  labs(
    x = "Covariate Value",
    y = "Annual Probability of Mortality"#,
    #title = "Effect of Covariates on Probability of Mortality"
  ) +
  theme_bw(base_size = 16)+theme(panel.grid = element_blank())
ggsave(height = 8, width = 14, 
       paste0(output.folder, "SPCD_stanoutput_full_standardized_v3/images/predicted_mortality/model_6_species_level_model_marginal_growth_cliamte_noBF.png"))

ggplot(marginal_response_df %>% filter(  predictor.class %in% "Site x G & S" & !Species %in% "balsam fir" ), aes(x = Value, y = 1-mean, color = Species)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
  facet_wrap(~ Predictor, scales = "free", ncol = 3) +
  labs(
    x = "Covariate Value",
    y = "Annual Probability of Mortality"#,
    #title = "Effect of Covariates on Probability of Mortality"
  ) +
  theme_bw(base_size = 16)+theme(panel.grid = element_blank())
ggsave(height = 5, width = 12, 
       paste0(output.folder, "SPCD_stanoutput_full_standardized_v3/images/predicted_mortality/model_6_species_level_model_marginal_growth_site_noBF.png"))

##############################################################################################
# make one giant figure for the marginal species effects:
# generate figures for all of these separately 
growth.marginal <- ggplot(marginal_response_df %>% filter( predictor.class %in% "Growth & Size"), aes(x = Value, y = 1-mean, color = Species)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
  facet_grid(cols = vars(predictor.class), rows = vars(Predictor), scales = "free") +
  labs(
    x = "Predictor Value",
    y = "Annual Probability of Mortality",
    #title = "Effect of Predictors on Probability of Mortality"
  ) +
  theme_bw(base_size = 14)+theme(panel.grid = element_blank())

Competition.marginal <- ggplot(marginal_response_df %>% filter( predictor.class %in% "Competition"), aes(x = Value, y = 1-mean, color = Species)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
  facet_grid(cols = vars(predictor.class), rows = vars(Predictor), scales = "free") +
  labs(
    x = "Predictor Value",
    y = "Annual Probability of Mortality",
    #title = "Effect of Predictors on Probability of Mortality"
  ) +
  theme_bw(base_size = 14)+theme(panel.grid = element_blank())

site.cond.marginal <- ggplot(marginal_response_df %>% filter( predictor.class %in% "Site Conditions"), aes(x = Value, y = 1-mean, color = Species)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
  facet_grid(cols = vars(predictor.class), rows = vars(Predictor), scales = "free") +
  labs(
    x = "Predictor Value",
    y = "Annual Probability of Mortality",
    #title = "Effect of Predictors on Probability of Mortality"
  ) +
  theme_bw(base_size = 14)+theme(panel.grid = element_blank())

Climate.marginal <- ggplot(marginal_response_df %>% filter( predictor.class %in% "Climate"), aes(x = Value, y = 1-mean, color = Species)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
  facet_grid(cols = vars(predictor.class), rows = vars(Predictor), scales = "free") +
  labs(
    x = "Predictor Value",
    y = "Annual Probability of Mortality",
    #title = "Effect of Predictors on Probability of Mortality"
  ) +
  theme_bw(base_size = 14)+theme(panel.grid = element_blank())



species.legend <- cowplot::get_legend(growth.marginal)

cowplot::plot_grid(cowplot::plot_grid(growth.marginal + theme(legend.position = "none"), 
                                      Competition.marginal + theme(legend.position = "none"), 
                                      site.cond.marginal + theme(legend.position = "none"), 
                                      Climate.marginal + theme(legend.position = "none"), align = "hv"), species.legend, rel_widths = c(1,0.2))

ggsave(height = 12, width = 12, 
       paste0(output.folder, "SPCD_stanoutput_full_standardized_v3/images/predicted_mortality/model_6_species_level_model_marginal_main_effects_summary.png"))
