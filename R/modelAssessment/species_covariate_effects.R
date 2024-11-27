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

for(i in 1:17){
  
SPCD.id <- nspp[i,]$SPCD
output.folder = "C:/Users/KellyHeilman/Box/01. kelly.heilman Workspace/mortality/Eastern-Mortality/mortality_models/"

# Extract posterior samples for the model and species
fit <- readRDS(paste0(output.folder,"/SPCD_stanoutput_full_standardized/samples/model_",model.no,"_SPCD_",SPCD.id, "_remper_correction_0.5.RDS"))
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
for(i in 1:length (covariate_ranges_df[,2])){
  cov.list[[i]] <- seq(covariate_ranges_df[i,2], covariate_ranges_df[i,3], length.out = 25)
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
  for(i in 1:length(covariate_matrix[,1])){
    linear_predictor <- beta_0 + rowSums( as.matrix(beta)%*% covariate_matrix[i,] )
  probabilities[[i]] <-  as.vector(inv_logit_fxn(linear_predictor))
  }
  

  
  # Summarize probabilities (mean and 90% credible interval)
  prob_summary <- lapply(probabilities, function(p) {
    c(mean = median(p$alpha_SPP), ci.lo = quantile(p$alpha_SPP, 0.1), ci.hi = quantile(p$alpha_SPP, 0.9))
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
ggsave(height = 10, width = 12, paste0(output.folder, "images/predicted_mortality/model_6_all_marginal_SPCD_", SPCD.id, ".png"))

# save for this species:
plot_data$SPCD <- SPCD.id
saveRDS(plot_data, paste0(output.folder, "SPCD_stanoutput_full_standardized/predicted_mort/model_6_all_marginal_SPCD_", SPCD.id, ".RDS") )
}


################################################################################
# combine all of the species conditional responses together
################################################################################
marginal_response <- list()
for (i in 1:17){
  SPCD.id <- nspp[i,]$SPCD
  marginal_response[[i]] <- readRDS( paste0(output.folder, "SPCD_stanoutput_full_standardized/predicted_mort/model_6_all_marginal_SPCD_", SPCD.id, ".RDS") )
  
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
  theme_minimal() 
ggsave(height = 12, width = 12, paste0(output.folder, "images/predicted_mortality/model_6_species_level_models_all_species_marginal_effects.png"))


## split up by variable type:
## main effects:
growth.diam <- c(unique(marginal_response_df$Covariate)[1:2], unique(marginal_response_df$Covariate)[19])
competition <- c(unique(marginal_response_df$Covariate)[3:7])
climate <- c(unique(marginal_response_df$Covariate)[8:13])
site.vars  <- c(unique(marginal_response_df$Covariate)[14:18])

# interaction terms
# growth 
growth.competition.int <- c(unique(marginal_response_df$Covariate)[20:24])
growth.climate.int <- c(unique(marginal_response_df$Covariate)[25:30])
growth.site.int <- c(unique(marginal_response_df$Covariate)[31:35])

# diameter
diameter.competition.int <- c(unique(marginal_response_df$Covariate)[36:40])
diameter.climate.int <- c(unique(marginal_response_df$Covariate)[41:46])
diameter.site.int <- c(unique(marginal_response_df$Covariate)[47:51])

# generate figures for all of these separately 
ggplot(marginal_response_df %>% filter( Covariate %in% growth.diam ), aes(x = Value, y = 1-mean, color = Species)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
  facet_wrap(~ Covariate, scales = "free") +
  labs(
    x = "Covariate Value",
    y = "Annual Probability of Mortality",
    title = "Effect of Covariates on Probability of Mortality"
  ) +
  theme_bw(base_size = 12)+theme(panel.grid = element_blank())
ggsave(height = 5, width = 12, 
       paste0(output.folder, "images/predicted_mortality/model_6_species_level_model_marginal_growth_diam.png"))

ggplot(marginal_response_df %>% filter( Covariate %in% competition ), aes(x = Value, y = 1-mean, color = Species)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
  facet_wrap(~ Covariate, scales = "free", ncol = 5) +
  labs(
    x = "Covariate Value",
    y = "Annual Probability of Mortality",
    title = "Effect of Covariates on Probability of Mortality"
  ) +
  theme_bw(base_size = 12)+theme(panel.grid = element_blank())
ggsave(height = 5, width = 12, 
       paste0(output.folder, "images/predicted_mortality/model_6_species_level_model_marginal_competition.png"))


ggplot(marginal_response_df %>% filter( Covariate %in% climate ), aes(x = Value, y = 1-mean, color = Species)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
  facet_wrap(~ Covariate, scales = "free", ncol = 3) +
  labs(
    x = "Covariate Value",
    y = "Annual Probability of Mortality",
    title = "Effect of Covariates on Probability of Mortality"
  ) +
  theme_bw(base_size = 12)+theme(panel.grid = element_blank())
ggsave(height = 8, width = 12, 
       paste0(output.folder, "images/predicted_mortality/model_6_species_level_model_marginal_climate_vars.png"))


ggplot(marginal_response_df %>% filter( Covariate %in% site.vars ), aes(x = Value, y = 1-mean, color = Species)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
  facet_wrap(~ Covariate, scales = "free", ncol = 5) +
  labs(
    x = "Covariate Value",
    y = "Annual Probability of Mortality",
    title = "Effect of Covariates on Probability of Mortality"
  ) +
  theme_bw(base_size = 12)+theme(panel.grid = element_blank())
ggsave(height = 5, width = 12, 
       paste0(output.folder, "images/predicted_mortality/model_6_species_level_model_marginal_site_vars.png"))

### interaction terms
ggplot(marginal_response_df %>% filter( Covariate %in% growth.competition.int ), aes(x = Value, y = 1-mean, color = Species)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
  facet_wrap(~ Covariate, scales = "free", ncol = 5) +
  labs(
    x = "Covariate Value",
    y = "Annual Probability of Mortality",
    title = "Effect of Covariates on Probability of Mortality"
  ) +
  theme_bw(base_size = 12)+theme(panel.grid = element_blank())
ggsave(height = 5, width = 12, 
       paste0(output.folder, "images/predicted_mortality/model_6_species_level_model_marginal_growth_competition.png"))


ggplot(marginal_response_df %>% filter( Covariate %in% growth.climate.int ), aes(x = Value, y = 1-mean, color = Species)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
  facet_wrap(~ Covariate, scales = "free", ncol = 3) +
  labs(
    x = "Covariate Value",
    y = "Annual Probability of Mortality",
    title = "Effect of Covariates on Probability of Mortality"
  ) +
  theme_bw(base_size = 12)+theme(panel.grid = element_blank())
ggsave(height = 8, width = 12, 
       paste0(output.folder, "images/predicted_mortality/model_6_species_level_model_marginal_growth_cliamte.png"))

ggplot(marginal_response_df %>% filter( Covariate %in% growth.site.int ), aes(x = Value, y = 1-mean, color = Species)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
  facet_wrap(~ Covariate, scales = "free", ncol = 5) +
  labs(
    x = "Covariate Value",
    y = "Annual Probability of Mortality",
    title = "Effect of Covariates on Probability of Mortality"
  ) +
  theme_bw(base_size = 12)+theme(panel.grid = element_blank())
ggsave(height = 5, width = 12, 
       paste0(output.folder, "images/predicted_mortality/model_6_species_level_model_marginal_growth_site.png"))


# diameter interaction terms
ggplot(marginal_response_df %>% filter( Covariate %in% diameter.competition.int ), aes(x = Value, y = 1-mean, color = Species)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
  facet_wrap(~ Covariate, scales = "free", ncol = 5) +
  labs(
    x = "Covariate Value",
    y = "Annual Probability of Mortality",
    title = "Effect of Covariates on Probability of Mortality"
  ) +
  theme_bw(base_size = 12)+theme(panel.grid = element_blank())
ggsave(height = 5, width = 12, 
       paste0(output.folder, "images/predicted_mortality/model_6_species_level_model_marginal_diameter_competition.png"))


ggplot(marginal_response_df %>% filter( Covariate %in% diameter.climate.int ), aes(x = Value, y = 1-mean, color = Species)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA) +
  facet_wrap(~ Covariate, scales = "free", ncol = 3) +
  labs(
    x = "Covariate Value",
    y = "Annual Probability of Mortality",
    title = "Effect of Covariates on Probability of Mortality"
  ) +
  theme_bw(base_size = 12)+theme(panel.grid = element_blank())
ggsave(height = 8, width = 12, 
       paste0(output.folder, "images/predicted_mortality/model_6_species_level_model_marginal_diameter_cliamte.png"))

ggplot(marginal_response_df %>% filter( Covariate %in% diameter.site.int ), aes(x = Value, y = 1-mean, color = Species)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = 1-ci.lo, ymax = 1-ci.hi, fill = Species), alpha = 0.2, color = NA, color = NA) +
  facet_wrap(~ Covariate, scales = "free", ncol = 5) +
  labs(
    x = "Covariate Value",
    y = "Annual Probability of Mortality",
    title = "Effect of Covariates on Probability of Mortality"
  ) +
  theme_bw(base_size = 12)+theme(panel.grid = element_blank())
ggsave(height = 5, width = 12, 
       paste0(output.folder, "images/predicted_mortality/model_6_species_level_model_marginal_diameter_site.png"))


