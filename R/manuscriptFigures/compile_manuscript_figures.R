library(tidyverse)
library(ggplot2)
library(cowplot)
library(posterior)
library(FIESTA)
#library(spdep)
library(sfdep)
library(sf)
library(spdep)
library(gt)
# Setting up species, colorschemes, etc: ----
output.dir <- output.folder <- "C:/Users/KellyHeilman/Box/01. kelly.heilman Workspace/mortality/Eastern-Mortality/mortality_models/"

csv.subdir <- "SPCD_stanoutput_cmdstan/fittedmodels"   # <- parent folder holding one subfolder per model
json.dir   <- file.path("SPCD_standata_json") # <- persistent location of the model input jsons

# # get the complete species list
nspp <- data.frame(SPCD = c(316, 318, 833, 832, 261, 531, 802, 129, 762, 12, 541, 97, 621, 400, 371, 241, 375))
nspp$Species <- paste(FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD), ]$GENUS,
                      FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD), ]$SPECIES)
nspp$COMMON <- FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD), ]$COMMON_NAME

spp.table <- data.frame(SPCD.id = nspp[1:17, ]$SPCD, spp = 1:17, COMMON = nspp[1:17, ]$COMMON)

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

species_fill <- scale_fill_manual(values = sppColors, name = "Species")
species_color <- scale_color_manual(values = sppColors, name = "Species")

# common species ordering scheme:
disturb.species.order <- c(
  # spruce fir forests
  "balsam fir", 
  "red spruce", 
  "northern white-cedar", 
  # HWA & beech
  "eastern hemlock", 
  "American beech", 
  
  # spongy moth susceptible
 # "black oak", 
  "chestnut oak", 
  "northern red oak",
  "white oak", 
  "yellow birch", 
  "paper birch",

  
  # spongy moth resistant
  "hickory spp.", 
  "eastern white pine", 
  "red maple", 
  "sugar maple", 
  
  # spongy moth immune
  "black cherry", 
  "white ash", 
  "yellow-poplar"
)

# update table 1 in the manuscript: ----

df <- tibble(
  Model = c("Model 1", 
                  "Model 2", 
                  "Model 3", 
                  "Model 4", 
                  "Model 5", 
                  "Model 6", 
                  "Model 7", 
                  "Model 8",
                  "Model 9"
  ),
  Covariates = c("Diameter Difference", 
                 "Diameter Difference & Diameter)", 
                 "Diameter Difference, Diameter & Neighborhood attributes", 
                 "Diameter Difference, Diameter, Neighborhood attributes & Climate Variables", 
                 "Diameter Difference, Diameter, Neighborhood attributes, Climate Variables & Site Conditions", 
                 "Model 5 + all Diameter difference and Diameter interactions", 
                 "Model 6 + all competition interactions",  
                  "Model 7 + all climate interactions",
                  "Model 8 + all site interactions (all fixed efffects and two way interactions)"
  ),
  Formula = c("$\\beta_{DBH}DBH$", 
              
              "$\\bm{\\sum\\limits_{s = 1}^{2}\\beta_{Size}Size_{s}} = \\beta_{DBH}DBH + \\beta_{DBH}DBHdiff$", 
              
              "$ \\bm{\\beta_{Size}Size_{s}} + \\bm{\\sum\\limits_{h = 1}^{3}\\beta_{Neighbor}*Neighbor_{h}}$", 
              
              "$\\bm{\\beta_{Size}Size_{s}} + \\bm{\\sum\\limits_{h = 1}^{3}\\beta_{Neighbor}*Neighbor_{h}} +
              \\bm{\\sum\\limits_{c = 1}^{3}\\beta_{Climate}*Climate_{c}}$", 
              
              "$\\bm{\\sum\\limits_{j = 1}^{12}\\beta_{j}X_{j}} = \\bm{\\beta_{Size}*Size{s}} + \\bm{\\sum\\limits_{h = 1}^{3}\\beta_{Neighbor}*Neighbor_{h}} +
              \\bm{\\sum\\limits_{c = 1}^{3}\\beta_{Climate}*Climate_{c}} + \\bm{\\sum\\limits_{d = 1}^{4}\\beta_{Site}*Site_{d}}$",
              
              "$\\bm{\\sum\\limits_{j = 1}^{12}\\beta_{j}X_{j}} + \\sum\\limits_{j > s} \\beta_{j,Size}(X_j Size_{s})$",
              
              "$\\bm{\\sum\\limits_{j = 1}^{12}\\beta_{j}X_{j}} + \\sum\\limits_{j > s} \\beta_{j,Size}(X_j Size_{s}) + \\sum\\limits_{j > h} \\beta_{j,Neighbor}(X_j Neighbor_{h})$",
              
              "$\\bm{\\sum\\limits_{j = 1}^{12}\\beta_{j}X_{j}} + \\sum\\limits_{j > s} \\beta_{j,Size}(X_j Size_{s}) + \\sum\\limits_{j > h} \\beta_{j,Neighbor}(X_j Neighbor_{h}) +
              \\sum\\limits_{j > c} \\beta_{j, Climate}(X_j Climate_{c})$",
              
              "$\\bm{\\sum\\limits_{j = 1}^{12}\\beta_{j}X_{j}} + \\sum\\limits_{j > s} \\beta_{j,Size}(X_j Size_{s}) + \\sum\\limits_{j > h} \\beta_{j,Neighbor}(X_j Neighbor_{h}) +
              \\sum\\limits_{j > c} \\beta_{j, Climate}(X_j Climate_{c}) +
              \\sum\\limits_{j > d} \\beta_{j, Site}(X_j Site_{d})$"
              
             
  ), 

  Number.betas = c(1, 
                   2, 
                   5,
                   9, 
                   12, 
                   33, 
                   57, 
                   75, 
                   78)
)

table1 <- df |> gt()|> 
  fmt_markdown(columns = Formula) |>
  cols_label(Covariates = "Description", 
             Formula = md("$\\bm{\\beta_{K}X_{K}} = $"),
             Number.betas = md("Number of $\\beta$'s (K)"))|>
  cols_width(Model ~ px(80), 
             Covariates ~ px(300),
             Formula ~ px(700), 
             Number.betas ~ px(80))|>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()
  )|>
  tab_style(
    style = cell_text(font = "Times New Roman"),
    locations = cells_body()
  ) |>
  # Apply to column labels
  tab_style(
    style = cell_text(font = "Times New Roman"),
    locations = cells_column_labels()
  )|>
  tab_options(
    # Add horizontal lines between rows
    table.border.bottom.color = "black"
  )#|>

# save png and tex format
# pretty version:
gtsave(data = table1, filename = paste0(output.dir, "images/Table_1_model_summary.png"), vwidth = 1500)
# tables have to be inserted into MS word format: (ugh!)
gtsave(data = table1, filename = paste0(output.dir, "images/Table_1_model_summary.docx"))
#--------------------------------------------------------------------------------------
# Observed Mortality (number of trees, rates, etc)--------
  
# Counts of trees included and not included:-----

# unfiltered tree remeasurement data
TREE.remeas <- readRDS( "data/unfiltered_TREE.remeas.rds")
TREE.remeas$SCIENTIFIC_NAME <- ref_species[match(TREE.remeas$SPCD, ref_species$SPCD),]$SCIENTIFIC_NAME


# Table S1: Count of trees above and below a diameter threshold
TREE.remeas %>% filter(!is.na(dbhold))%>%
  filter(SPCD %in% spp.table[1:17,]$SPCD.id) %>%
  filter(remper > 0 & 
  !is.na(remper)) %>% 
  
  group_by(SCIENTIFIC_NAME, Species, SPCD, dbhold >= 5& dbhcur >=5) %>%
  summarise(n()) %>% 
  group_by(SCIENTIFIC_NAME, Species, SPCD) %>%
  spread(`dbhold >= 5 & dbhcur >= 5`, `n()`) %>%
  ungroup()%>%
  rename("DBH >= 12.7 cm" = "FALSE", 
         "DBH < 12.7 cm" = "TRUE", 
         "Scientific Name" = "SCIENTIFIC_NAME")|>
  gt() %>%
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(
      columns = `Scientific Name`
    )
  ) %>% grand_summary_rows(
    columns = c(`DBH >= 12.7 cm`, `DBH < 12.7 cm`),
    fns = list(
     Total ~ sum(.)
    ),
    fmt = ~ fmt_number(., use_seps = FALSE,decimals = 0)
  )|>gtsave(filename = paste0(output.dir, "images/tables/Ntrees_by_diameter_threshold.png"))

# Table S2: Count of trees > 12.7 cm by Status code and 
# diameter change between measurements (DIAMETER_diff).
TREE.remeas %>% filter(!is.na(dbhold))%>%
  filter(SPCD %in% spp.table[1:17,]$SPCD.id) %>%
  filter(remper > 0 & 
           !is.na(remper)) %>% 
  filter(dbhold >=5 & dbhcur >=5) %>%
  filter(!is.na(DIA_DIFF))%>%
  mutate(status = ifelse(status == 6, 5, status))%>%
  
  group_by(status, DIA_DIFF <= 0) %>%
  summarise(n()) %>% 
  group_by(`DIA_DIFF <= 0`) %>%
  spread(status, `n()`) %>%
  ungroup()%>%
  mutate(`DIA_DIFF <= 0` = ifelse(`DIA_DIFF <= 0` == FALSE, "Greater than 0", "less than or equal to 0"))%>%
  mutate(Totals = `1`+`2`+`3`+`4`+`5`)%>%
  mutate(`Total (non-cut)` = `1`+`2`+`4`+`5`)%>%
  rename("Live" = `1`, 
         "Salvagable Dead" = `2`, 
         "Cut" = `3`,
         "Non-Salvagable Dead" = `4`, 
         "Snags" = `5`,
         "Diameter Difference" = `DIA_DIFF <= 0`) %>%
  #select(`Scientific Name`, Species, SPCD, Live,`Salvagable Dead`,`Non-Salvagable Dead`, `Snags`, Cut)|>
  gt() |>gtsave(filename = paste0(output.dir, "images/tables/Ntrees_DIA_DIFF_by_status.png"))




# Table S3: Total number of trees used in this analysis
TREE.remeas %>% filter(!is.na(dbhold))%>%
  filter(SPCD %in% spp.table[1:17,]$SPCD.id) %>%
  filter(remper > 0 & 
           !is.na(remper)) %>% 
  filter(dbhold >=5 & dbhcur >=5) %>%
  filter(DIA_DIFF > 0) %>%
  mutate(status = ifelse(status == 6, 5, status))%>%
  
  group_by(SCIENTIFIC_NAME, Species, SPCD, status) %>%
  summarise(n()) %>% 
  group_by(SCIENTIFIC_NAME, Species, SPCD) %>%
  spread(status, `n()`) %>%
  ungroup()%>%
  rename("Live" = `1`, 
         "Salvagable Dead" = `2`, 
         "Cut" = `3`,
         "Non-Salvagable Dead" = `4`, 
         "Snags" = `5`,
         "Scientific Name" = "SCIENTIFIC_NAME") %>%
  mutate(`Non-cut Dead` = `Salvagable Dead`+`Non-Salvagable Dead`+`Snags`) %>%
  select(`Scientific Name`, Species, SPCD, Cut, `Salvagable Dead`,`Non-Salvagable Dead`, `Snags`, `Non-cut Dead`, Live)|>
  gt() %>%
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(
      columns = `Scientific Name`
    )
  ) %>% grand_summary_rows(
    columns = c("Live","Non-cut Dead", "Salvagable Dead", "Non-Salvagable Dead", "Snags","Cut"),
    fns = list(
      Total ~ sum(.)
    ),
    fmt = ~ fmt_number(., use_seps = FALSE,decimals = 0)
  )|>
  gtsave(filename = paste0(output.dir, "images/tables/Ntrees_in_analysis_5in_by_status_.png"),vwidth = 1500, vheight = 1000)

# Table S4: covariate predictor variables explored in GLM models

tribble(
  ~Class, ~Covariate, ~Definition, ~Source, ~Scaling,  ~Included,
  "Size-Related", "DIA_DIFF", "Diameter difference betwen remeasurement", "Periodic NFI",  "by species","included",
  "Size-Related", "Diameter", "Inital tree diameter at breast height", "Periodic NFI", "by species","included",
 
  "Site Conditions", "slope", "slope of plot", "Periodic NFI", "by species","included",
  "Site Conditions","aspect", "aspect of plot", "Periodic NFI", "by species","included",
  "Site Conditions","elev", "Elevation",  "PRISM", "all plots","included",
  "Site Conditions","SI",        "Site Index",  "Periodic NFI", "by species","excluded",
  "Site Conditions","PHYSIO", "Physiographic class of plot", "Periodic NFI", "none","excluded",
  
  "Climate Normals","MAP", "Annual Precipitation Normal",  "PRISM", "all plots","included",
  "Climate Normals", "MATmax", "Annual Maximum Temperature Normal", "PRISM", "all plots","included",
  "Climate Normals","MATmin", "Annual Minimum Temperature Normal",  "PRISM", "all plots","excluded",
  "Climate Anomalies","MAPanom", "Annual Precipitation anomaly over the remeasurement period",  "PRISM", "by species","included",
  "Climate Anomalies","MATmaxanom", "Maximum temperature anomaly over the remeasurement period",  "PRISM", "all plots","included",
  "Climate Anomalies","MATminanom", "Minimum temperature anomaly over the remeasurement period",  "PRISM", "all plots","excluded",
  "Pollution","Ndep", "Nitrogen Deposition change over the remeasurement period",  "all plots","TREND-Nitrogen data","included",
  
  "Neighborhood Characteristics","damage", "Percent of trees with damage or mortality agent listed on the plot", "Periodic NFI - derived", "by species","included",
  "Neighborhood Characteristics","BAL", "Basal Area Larger Than (tree-level)", "Periodic NFI - derived", "by species","included",
  "Neighborhood Characteristics","BA", "Basal Area (Plot)", "Periodic NFI - derived", "by species","included",
  "Neighborhood Characteristics","RD", "Relative Density of plot", "Periodic NFI - derived", "by species","excluded",
  
  
  "Neighborhood Characteristics","SPCD.BA", "Basal area of focal species on plot", "Periodic NFI - derived", "by species","excluded",
  "Neighborhood Characteristics","non_SPCD.BA", "Basal area of non-focal species on plot", "Periodic NFI - derived", "by species","excluded",
  "Neighborhood Characteristics","prop.focal.BA", "Porportion basal area on plot made up by focal species", "Periodic NFI - derived", "by species","excluded"
)%>%
  group_by(Class)|>
  gt()%>%
  tab_style(
    style = cell_text(weight = "bold", color = "black"),
    locations = cells_row_groups()
  ) %>%tab_style(
    style = cell_fill(color = "lightgrey"),
    locations = cells_row_groups()
   ) %>%
  # 
  # # 5. Additional Styling
  opt_stylize(style = 6, color = 'gray') %>%
   tab_options(
     row_group.padding = px(5)
   )|>gtsave(filename = paste0(output.dir, "images/tables/Covariates_considered_in_analysis.png"),vwidth = 1500)

# Table S5: GLM models fit:
read.csv("SPCD_glm_output/GLM_reduced_table.csv") %>%
  rename("Model Number" = "Model.Number", 
         "Model Type" = "Model.Type")|>
  gt()%>%
  tab_style(
    style = cell_fill(color = "#00BFC4"),
    locations = cells_body(rows = c(1:14)) # Colors 2nd, 4th, and 6th rows
  )%>%
  tab_style(
    style = cell_fill(color = "#7CAE00"),
    locations = cells_body(rows = c(15:27)) # Colors 2nd, 4th, and 6th rows
  )%>%
  tab_style(
    style = cell_fill(color = "#F8766D"),
    locations = cells_body(rows = c(28:39)) # Colors 2nd, 4th, and 6th rows
  ) %>%
  tab_options(
    table_body.vlines.color  = "white",
    column_labels.border.top.color = "black",
    column_labels.border.bottom.color = "black",
    table_body.border.bottom.color = "black"
  )|>gtsave(filename = paste0(output.dir, "images/tables/GLM_model_building.png"),vheight = 1500)

#####################################################################################
# re-make AUC score figure with weighted model summary ----------------------

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
  
  weighted.file <- list.files(path = paste0(output.dir,"SPCD_stanoutput_cmdstan/AUC/"), 
                                   pattern =paste0(SPCD.id,"_stacking_wts.qs"), full.names = T)
  
  
  all.AUC.files <- c(spp.AUC.files, hierarchical.AUC.files, weighted.file)
  
  AUC_results_all <- lapply(all.AUC.files, qs_read)
  
  lapply(AUC_results_all, colnames)
  
  AUC_confusion_summary <- do.call(rbind, lapply(AUC_results_all, function(x){
    x|> group_by(model.number, type, model.type, SPCD)%>%
      summarise(AUC_median = median(AUC),
                AUC_ci.lo = quantile(AUC, 0.025),
                AUC_ci.hi = quantile(AUC, 0.975),
                
                True_surv_rate = median(`True survival rate`),
                True_surv_rate_ci.lo = quantile(`True survival rate`, 0.025),
                True_surv_rate_ci.hi = quantile(`True survival rate`, 0.975),
                
                True_mort_rate = median(`True mortality rate`),
                True_mort_rate_ci.lo = quantile(`True mortality rate`, 0.025),
                True_mort_rate_ci.hi = quantile(`True mortality rate`, 0.975), .groups = "drop_last")%>%
      mutate(model.number = as.character(model.number))
  }))
  return(AUC_confusion_summary)
}


AUC.df <- do.call(rbind, lapply(nspp$SPCD, AUC_summarise_SPCD))
AUC.df$COMMON_NAME <- FIESTA::ref_species[match(AUC.df$SPCD, FIESTA::ref_species$SPCD),]$COMMON_NAME
AUC.df$COMMON_NAME <- factor(AUC.df$COMMON_NAME, levels = unique(nspp$COMMON))
AUC.df <- AUC.df %>% 
  mutate(`Model Type` = ifelse(model.type %in% c("stacking_wts", "weighting"), "Stacked Posteriors", model.type))%>%
  mutate(model.number = ifelse(model.number %in% "BMA", "MA", model.number))

# make new figures with AUC statistics for BMA AUC scores---
AUC.df %>% filter(type == "in-sample")|> 
  ggplot()+
  geom_pointrange(aes(x = as.character(model.number), y = AUC_median, ymin = AUC_ci.lo, ymax = AUC_ci.hi, shape = `Model Type`, color =`Model Type`), position = position_dodge(width = 1))+
  facet_wrap(~COMMON_NAME)+theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  ylab("In-sample AUC score")+
  xlab("Model Number")
  

ggsave(filename = paste0(output.dir, "SPCD_stanoutput_cmdstan/summary/AUC_is_all_models_stacking.png"), 
       height = 7, width = 10)

AUC.df %>% filter(type == "out-of-sample")|> 
  ggplot()+geom_pointrange(aes(x = as.character(model.number), y = AUC_median, ymin = AUC_ci.lo, ymax = AUC_ci.hi, shape = `Model Type`, color = `Model Type`), position = position_dodge(width = 1))+
  facet_wrap(~COMMON_NAME)+theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  ylab("Held-out AUC score")+
  xlab("Model Number")

ggsave(filename = paste0(output.dir, "SPCD_stanoutput_cmdstan/summary/AUC_oos_all_models_stacking.png"), 
       height = 7, width = 10)

# #####################################################################################
# # Brier scores-- tree level predictions
# 
# # calculate Brier score for each draw:------
# brier_calc <- function(prob, y_obs){
#   mean((prob - y_obs)^2)
# }
#   
# # get posterior predictive probabilities:
# 
# # 2. Get the actual observed response vector
# y_obs <- as.numeric(fit$data$your_outcome_variable) - 1 # ensure 0 and 1
# 
# # 3. Calculate Brier score for each posterior draw
# # This yields a distribution of Brier scores capturing parameter uncertainty
# brier_dist <- apply(p_draws, 2, function(p) mean((p - y_obs)^2))
# 
# # 4. Summarize the posterior Brier score (mean and 95% credible interval)
# mean(brier_dist)
# quantile(brier_dist, probs = c(0.025, 0.975))
# 
# AUC.is.samples.df  <- apply(pSurv_hat_samps, 1, function(p) as.numeric(pROC::auc(actuals, p, quiet = TRUE)))
# AUC.oos.samples.df <- apply(pSurv_rep_samps, 1, function(p) as.numeric(pROC::auc(actuals.oos, p, quiet = TRUE)))
# #


#####################################################################################
# Evaluated predicted vs observed mortality rates estimated for species at the county, state, and overall scale -----
#scale = c("Regional", "State", "ST_CTY") # naming convention for the files

scale <- "State"
SPCD.id <- nspp$SPCD[17]

# get the predicted mortality rates across the region, weighted by the number of trees each tree represented:
# list all files for the scale
summarise_mort_scale <- function(SPCD.id, scale){
  
  cat(paste0("reading and summarising ", scale," estimates for SPCD ", SPCD.id, "\n"))
        all.scale.files <- list.files(paste0(output.folder,"SPCD_stanoutput_cmdstan/predicted_mort/"), 
                                      pattern = scale, full.names = T)
        # get the speices files
        species.scale.files <- grep( paste0("_", SPCD.id, "(_|\\.)"),  
                                                  all.scale.files, 
                                                      value = TRUE, 
                                                      perl = TRUE)
        
        
        scale_summary_df <- do.call(rbind, lapply(species.scale.files, FUN = function(file_scale){
          # add special summary grouping by scale:
          
          # state-level scaling 
          if(scale %in% "State"){
            
          mortality_rate_summary <- qs2::qs_read(file_scale) %>% 
            group_by(SPCD, model.number, model.type, state)%>%
            
            summarise(obs_M_median = median(Obs_mort_rate, na.rm =TRUE), 
                      n_obs = median(total_obs, na.rm =TRUE),
                      pred_M_median = median(Pred_mort_rate, na.rm =TRUE), 
                      pred_M_2.5.ci.lo = quantile(Pred_mort_rate, 0.025, na.rm =TRUE), 
                      pred_M_97.5.ci.hi = quantile(Pred_mort_rate, 0.975, na.rm =TRUE),
                      pred_M_5.ci.lo = quantile(Pred_mort_rate, 0.05, na.rm =TRUE), 
                      pred_M_95.ci.hi = quantile(Pred_mort_rate, 0.95, na.rm =TRUE),
                      pred_M_25.ci.lo = quantile(Pred_mort_rate, 0.25, na.rm =TRUE), 
                      pred_M_75.ci.hi = quantile(Pred_mort_rate, 0.75, na.rm =TRUE), 
                      
                      obs_M_expn_median = median(Obs_mort_rate_expn, na.rm =TRUE), 
                      
                      pred_M_expn_median = median(Pred_mort_rate_expn, na.rm =TRUE), 
                      pred_M_expn_2.5.ci.lo = quantile(Pred_mort_rate_expn, 0.025, na.rm =TRUE), 
                      pred_M_expn_97.5.ci.hi = quantile(Pred_mort_rate_expn, 0.975, na.rm =TRUE),
                      pred_M_expn_5.ci.lo = quantile(Pred_mort_rate_expn, 0.05, na.rm =TRUE), 
                      pred_M_expn_95.ci.hi = quantile(Pred_mort_rate_expn, 0.95, na.rm =TRUE),
                      pred_M_expn_25.ci.lo = quantile(Pred_mort_rate_expn, 0.25, na.rm =TRUE), 
                      pred_M_expn_75.ci.hi = quantile(Pred_mort_rate_expn, 0.75, na.rm =TRUE), 
                      .groups = "drop")%>% 
            mutate(scale = "State")
          return(mortality_rate_summary)
          }
          
          # county-level scaling 
          if(scale %in% "ST_CTY"){
            mortality_rate_summary <- qs2::qs_read(file_scale) %>% 
              group_by(SPCD, model.number, model.type, ST_CTY)%>%
              
              summarise(obs_M_median = median(Obs_mort_rate, na.rm =TRUE), 
                        n_obs = median(total_obs, na.rm =TRUE),
                        pred_M_median = median(Pred_mort_rate, na.rm =TRUE), 
                        pred_M_2.5.ci.lo = quantile(Pred_mort_rate, 0.025, na.rm =TRUE), 
                        pred_M_97.5.ci.hi = quantile(Pred_mort_rate, 0.975, na.rm =TRUE),
                        pred_M_5.ci.lo = quantile(Pred_mort_rate, 0.05, na.rm =TRUE), 
                        pred_M_95.ci.hi = quantile(Pred_mort_rate, 0.95, na.rm =TRUE),
                        pred_M_25.ci.lo = quantile(Pred_mort_rate, 0.25, na.rm =TRUE), 
                        pred_M_75.ci.hi = quantile(Pred_mort_rate, 0.75, na.rm =TRUE), 
                        
                        obs_M_expn_median = median(Obs_mort_rate_expn, na.rm =TRUE), 
                        
                        pred_M_expn_median = median(Pred_mort_rate_expn, na.rm =TRUE), 
                        pred_M_expn_2.5.ci.lo = quantile(Pred_mort_rate_expn, 0.025, na.rm =TRUE), 
                        pred_M_expn_97.5.ci.hi = quantile(Pred_mort_rate_expn, 0.975, na.rm =TRUE),
                        pred_M_expn_5.ci.lo = quantile(Pred_mort_rate_expn, 0.05, na.rm =TRUE), 
                        pred_M_expn_95.ci.hi = quantile(Pred_mort_rate_expn, 0.95, na.rm =TRUE),
                        pred_M_expn_25.ci.lo = quantile(Pred_mort_rate_expn, 0.25, na.rm =TRUE), 
                        pred_M_expn_75.ci.hi = quantile(Pred_mort_rate_expn, 0.75, na.rm =TRUE), 
                        .groups = "drop") %>% 
              mutate(scale = "County")
            return(mortality_rate_summary)
          }
          
          # regional overall scaling
          if(scale %in% "Regional"){
            mortality_rate_summary <- qs2::qs_read(file_scale) %>% 
              group_by(SPCD, model.number, model.type)%>%
              
              summarise(obs_M_median = median(Obs_mort_rate, na.rm =TRUE), 
                        n_obs = median(total_obs, na.rm =TRUE),
                        pred_M_median = median(Pred_mort_rate, na.rm =TRUE), 
                        pred_M_2.5.ci.lo = quantile(Pred_mort_rate, 0.025, na.rm =TRUE), 
                        pred_M_97.5.ci.hi = quantile(Pred_mort_rate, 0.975, na.rm =TRUE),
                        pred_M_5.ci.lo = quantile(Pred_mort_rate, 0.05, na.rm =TRUE), 
                        pred_M_95.ci.hi = quantile(Pred_mort_rate, 0.95, na.rm =TRUE),
                        pred_M_25.ci.lo = quantile(Pred_mort_rate, 0.25, na.rm =TRUE), 
                        pred_M_75.ci.hi = quantile(Pred_mort_rate, 0.75, na.rm =TRUE), 
                        
                        obs_M_expn_median = median(Obs_mort_rate_expn, na.rm =TRUE), 
                        
                        pred_M_expn_median = median(Pred_mort_rate_expn, na.rm =TRUE), 
                        pred_M_expn_2.5.ci.lo = quantile(Pred_mort_rate_expn, 0.025, na.rm =TRUE), 
                        pred_M_expn_97.5.ci.hi = quantile(Pred_mort_rate_expn, 0.975, na.rm =TRUE),
                        pred_M_expn_5.ci.lo = quantile(Pred_mort_rate_expn, 0.05, na.rm =TRUE), 
                        pred_M_expn_95.ci.hi = quantile(Pred_mort_rate_expn, 0.95, na.rm =TRUE),
                        pred_M_expn_25.ci.lo = quantile(Pred_mort_rate_expn, 0.25, na.rm =TRUE), 
                        pred_M_expn_75.ci.hi = quantile(Pred_mort_rate_expn, 0.75, na.rm =TRUE), 
                        .groups = "drop") %>% 
              mutate(scale = "Regional")
            return(mortality_rate_summary)
          }
          
        })
        )
        
        
        
        scale_summary_df$COMMON_NAME <- ref_species[match(scale_summary_df$SPCD, ref_species$SPCD),]$COMMON_NAME
        return(scale_summary_df)

}

# example usage:
#summarise_mort_scale(SPCD.id = 97, scale = "ST_CTY")


# get rmse for all species models at regional scale
#rmse <- sqrt(mean((actual - predicted)^2))
regional_summary_df <- do.call(rbind, lapply(spp.table$SPCD.id, FUN = function(sppcode){summarise_mort_scale(SPCD.id = sppcode, scale = "Regional")}))
regional_summary_df$COMMON_NAME <- ref_species[match(regional_summary_df$SPCD, ref_species$SPCD),]$COMMON_NAME

regional_summary_df %>% 
  group_by(COMMON_NAME, model.type, model.number)%>%
  summarise(RMSE_pct_mort = sqrt(mean((obs_M_expn_median*100 - pred_M_expn_median*100)^2)))|>
  ggplot()+geom_bar(aes(x = model.number, y = RMSE_pct_mort, group = model.type, fill = model.type), position = "dodge", stat = "identity")+
  facet_wrap(~COMMON_NAME, scales = "free")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  ylab("Root Mean Squared Error in Species-scale mortality rate (% tree mortality/year)")


regional_summary_df %>%
  mutate(residual = obs_M_expn_median*100 - pred_M_expn_median*100) %>% View()

# save summaries:
qs2::qs_save(regional_summary_df, paste0(output.dir, "SPCD_stanoutput_cmdstan/summary/Regional_summary_p_o_mort.qs"))
regional_summary_df <- qs2::qs_read( paste0(output.dir, "SPCD_stanoutput_cmdstan/summary/Regional_summary_p_o_mort.qs"))


# get state-level summaries
State_summary_df <- do.call(rbind, lapply(spp.table$SPCD.id, FUN = function(sppcode){summarise_mort_scale(SPCD.id = sppcode, scale = "State")}))
State_summary_df$COMMON_NAME <- ref_species[match(State_summary_df$SPCD, ref_species$SPCD),]$COMMON_NAME
# save summaries:
qs2::qs_save(State_summary_df, paste0(output.dir, "SPCD_stanoutput_cmdstan/summary/State_summary_p_o_mort.qs"))
State_summary_df <- qs2::qs_read( paste0(output.dir, "SPCD_stanoutput_cmdstan/summary/State_summary_p_o_mort.qs"))

State_summary_df %>% 
  group_by(COMMON_NAME, model.type, model.number)%>% 
  mutate(n_state = n_distinct(state))%>%
  summarise(MSE = mean((obs_M_expn_median*100 - pred_M_expn_median*100)^2, na.rm = TRUE))%>%
  group_by(COMMON_NAME, model.type, model.number)%>% 
  summarise(RMSE_pct_mort = sqrt(MSE))|>
  
  ggplot()+geom_bar(aes(x = model.number, y = RMSE_pct_mort, group = model.type, fill = model.type), position = "dodge", stat = "identity")+
  facet_wrap(~COMMON_NAME, scales = "free")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  ylab("Root Mean Squared Error in State-scale mortality rate (% tree mortality/year)")




# get ST_CTY-level summaries
ST_CTY_summary_df <- do.call(rbind, lapply(spp.table$SPCD.id, FUN = function(sppcode){summarise_mort_scale(SPCD.id = sppcode, scale = "ST_CTY")}))
ST_CTY_summary_df$COMMON_NAME <- ref_species[match(ST_CTY_summary_df$SPCD, ref_species$SPCD),]$COMMON_NAME
# save summaries:
qs2::qs_save(ST_CTY_summary_df, paste0(output.dir, "SPCD_stanoutput_cmdstan/summary/ST_CTY_summary_p_o_mort.qs"))
ST_CTY_summary_df <- qs2::qs_read(paste0(output.dir, "SPCD_stanoutput_cmdstan/summary/ST_CTY_summary_p_o_mort.qs"))



ST_CTY_summary_df %>%
  #filter(COMMON_NAME %in% "sugar maple")%>%
  group_by(COMMON_NAME, model.type, model.number)%>%
  summarise(MSE_pct = mean((obs_M_expn_median*100 - pred_M_expn_median*100)^2, na.rm =TRUE))%>%
  mutate(RMSE_pct_mort = sqrt(MSE_pct))|>
  #summarise(RMSE_pct_mort = sqrt(mean((obs_M_expn_median*100 - pred_M_expn_median*100)^2, na.rm =TRUE), na.rm =TRUE))|>
  ggplot()+geom_bar(aes(x = model.number, y = RMSE_pct_mort, group = model.type, fill = model.type), position = "dodge", stat = "identity")+
  facet_wrap(~COMMON_NAME, scales = "free")

ST_CTY_summary_df %>%
  filter(n_obs >= 25)%>%
  group_by(COMMON_NAME, model.type, model.number)%>%
  summarise(MSE_pct = mean((obs_M_expn_median*100 - pred_M_expn_median*100)^2, na.rm =TRUE))%>%
  mutate(RMSE_pct_mort = sqrt(MSE_pct))|>
  #summarise(RMSE_pct_mort = sqrt(mean((obs_M_expn_median*100 - pred_M_expn_median*100)^2, na.rm =TRUE), na.rm =TRUE))|>
  ggplot()+geom_bar(aes(x = model.number, y = RMSE_pct_mort, group = model.type, fill = model.type), position = "dodge", stat = "identity")+
  facet_wrap(~COMMON_NAME, scales = "free")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  ylab("Root Mean Squared Error in County-scale mortality rate (% tree mortality/year)")



ST_CTY_summary_df %>% filter(n_obs > 25)%>%
  #filter(model.number == "Stacked")|>
  ggplot()+geom_point(aes(x = obs_M_expn_median, y = pred_M_expn_median, color = model.type))+
  geom_errorbar(aes(x = obs_M_expn_median, ymin = pred_M_expn_5.ci.lo, ymax = pred_M_expn_95.ci.hi, color = model.type))+
  geom_errorbar(aes(x = obs_M_expn_median, ymin = pred_M_expn_25.ci.lo, ymax = pred_M_expn_75.ci.hi, color = model.type))+
  theme_bw()+
  geom_abline(aes(intercept = 0, slope = 1))+
  facet_wrap(~model.number)+ylab("Species predicted mortality rate (county-level)")+
  xlab("Species observed Mortality rate (county-level)")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggsave(paste0(output.dir, "SPCD_stanoutput_cmdstan/summary/ST_CTY_Species_pred_obs_25obs.png"), 
       height = 6, width = 9)


ST_CTY_summary_df %>% filter(n_obs > 50)%>%
  #filter(model.number == "Stacked")|>
  ggplot()+geom_point(aes(x = obs_M_expn_median, y = pred_M_expn_median, color = model.type))+
  geom_errorbar(aes(x = obs_M_expn_median, ymin = pred_M_expn_5.ci.lo, ymax = pred_M_expn_95.ci.hi, color = model.type))+
  geom_errorbar(aes(x = obs_M_expn_median, ymin = pred_M_expn_25.ci.lo, ymax = pred_M_expn_75.ci.hi, color = model.type))+
  theme_bw()+
  geom_abline(aes(intercept = 0, slope = 1))+
  facet_wrap(~model.number)+ylab("Species predicted mortality rate (county-level)")+
  xlab("Species observed Mortality rate (county-level)")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggsave(paste0(output.dir, "SPCD_stanoutput_cmdstan/summary/ST_CTY_Species_pred_obs_50obs.png"), 
       height = 6, width = 9)


ST_CTY_summary_df %>% 
  #filter(model.number == "Stacked")|>
  ggplot()+geom_point(aes(x = obs_M_expn_median, y = pred_M_expn_median, color = model.type))+
  geom_errorbar(aes(x = obs_M_expn_median, ymin = pred_M_expn_5.ci.lo, ymax = pred_M_expn_95.ci.hi, color = model.type))+
  geom_errorbar(aes(x = obs_M_expn_median, ymin = pred_M_expn_25.ci.lo, ymax = pred_M_expn_75.ci.hi, color = model.type))+
  theme_bw()+
  geom_abline(aes(intercept = 0, slope = 1))+
  facet_wrap(~model.number)+ylab("Species predicted mortality rate (county-level)")+
  xlab("Species observed Mortality rate (county-level)")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggsave(paste0(output.dir, "SPCD_stanoutput_cmdstan/summary/ST_CTY_Species_pred_obs.png"), 
       height = 6, width = 9)

State_summary_df %>% filter(n_obs > 50)%>%
  #filter(model.number == "Stacked")|>
  ggplot()+geom_point(aes(x = obs_M_expn_median, y = pred_M_expn_median, color = model.type))+
  geom_errorbar(aes(x = obs_M_expn_median, ymin = pred_M_expn_5.ci.lo, ymax = pred_M_expn_95.ci.hi, color = model.type))+
  geom_errorbar(aes(x = obs_M_expn_median, ymin = pred_M_expn_25.ci.lo, ymax = pred_M_expn_75.ci.hi, color = model.type))+
  theme_bw()+
  geom_abline(aes(intercept = 0, slope = 1))+
  facet_wrap(~model.number)+ylab("Species predicted mortality rate (state-level)")+
  xlab("Species observed Mortality rate (state-level)")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggsave(paste0(output.dir, "SPCD_stanoutput_cmdstan/summary/State_Species_pred_obs.png"), 
       height = 6, width = 9)



regional_summary_df %>% filter(n_obs > 50)%>%
  #filter(model.number == "Stacked")|>
  ggplot()+geom_point(aes(x = obs_M_expn_median, y = pred_M_expn_median, color = model.type))+
  geom_errorbar(aes(x = obs_M_expn_median, ymin = pred_M_expn_5.ci.lo, ymax = pred_M_expn_95.ci.hi, color = model.type))+
  geom_errorbar(aes(x = obs_M_expn_median, ymin = pred_M_expn_25.ci.lo, ymax = pred_M_expn_75.ci.hi, color = model.type))+
  theme_bw()+
  geom_abline(aes(intercept = 0, slope = 1))+
  facet_wrap(~model.number)+ylab("Species overall predicted mortality rate")+
  xlab("Species regional Mortality rate")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggsave(paste0(output.dir, "SPCD_stanoutput_cmdstan/summary/Regional_Species_pred_obs.png"), 
       height = 6, width = 9)


# save the stacked pred vs observed by species for all regional values
all_scale_stacked <- regional_summary_df %>%
  filter(model.number == "Stacked") %>%
  mutate(group = "region") %>%
  rbind(., 
        State_summary_df %>% 
          filter(model.number == "Stacked")%>%
          rename("group" = "state")%>%
          select(everything(), group)) %>%
  rbind(., ST_CTY_summary_df %>% filter(n_obs > 50)%>%
          filter(model.number == "Stacked")%>%
          rename("group" = "ST_CTY")%>%
          select(everything(), group))

all_scale_stacked$scale <- factor(all_scale_stacked$scale, levels = c("County", "State", "Regional"))

p.o_stacked_expanded <- all_scale_stacked %>% filter(n_obs >25)|>
  ggplot()+geom_point(aes(x = obs_M_expn_median*100, y = pred_M_expn_median*100, color = COMMON_NAME), size = 2)+
  geom_errorbar(aes(x = obs_M_expn_median*100, ymin = pred_M_expn_5.ci.lo*100, ymax = pred_M_expn_95.ci.hi*100, color = COMMON_NAME), size = 0.5)+
  geom_errorbar(aes(x = obs_M_expn_median*100, ymin = pred_M_expn_25.ci.lo*100, ymax = pred_M_expn_75.ci.hi*100, color = COMMON_NAME), size = 1)+
  theme_bw(base_size = 16)+
  facet_wrap(~scale)+
  geom_abline(aes(intercept = 0, slope = 1))+
  ylab("Predicted mortality rate (% trees/year)")+
  xlab("Observed Mortality rate (% trees/year)")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  species_color+
  theme(panel.grid = element_blank())

ggsave(paste0(output.dir, "SPCD_stanoutput_cmdstan/summary/Stacked_model_pred_obs_EXPN.png"), 
       p.o_stacked_expanded,
       height = 6, width = 11)

p.o_stacked_per_acre <- all_scale_stacked %>% filter(n_obs >25)|>
  ggplot()+geom_point(aes(x = obs_M_median*100, y = pred_M_median*100, color = COMMON_NAME), size = 2)+
  geom_errorbar(aes(x = obs_M_median*100, ymin = pred_M_5.ci.lo*100, ymax = pred_M_95.ci.hi*100, color = COMMON_NAME), size = 0.5, alpha = 0.5)+
  geom_errorbar(aes(x = obs_M_median*100, ymin = pred_M_25.ci.lo*100, ymax = pred_M_75.ci.hi*100, color = COMMON_NAME), size = 1, alpha = 0.75)+
  theme_bw(base_size = 16)+
  facet_wrap(~scale)+
  geom_abline(aes(intercept = 0, slope = 1))+
  ylab("Predicted mortality rate (% trees/acre/year)")+
  xlab("Observed Mortality rate (% trees/acre/year)")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  species_color+
  coord_cartesian(ylim = c(0, 2.5), xlim = c(0, 2.5))+
  theme(panel.grid = element_blank())

ggsave(paste0(output.dir, "SPCD_stanoutput_cmdstan/summary/Stacked_model_pred_obs_per_acre.png"), 
       p.o_stacked_per_acre,
       height = 6, width = 11)


################################################################################
# calculate RMSE & BIAS for each model, species, scale ---
# regional scale RMSE and BIAS
Regional_RMSE_BIAS <- regional_summary_df %>% filter(n_obs > 50)%>%
  mutate(residual_diff = obs_M_expn_median*100 - pred_M_expn_median*100)%>%
  group_by(model.number, model.type, COMMON_NAME)%>%
  summarise(BIAS = sum(residual_diff, na.rm = TRUE), 
            RMSE = sqrt(mean(residual_diff^2, na.rm = TRUE)))

ggplot(data = Regional_RMSE_BIAS)+
  geom_point(aes(model.number, BIAS, color = model.type, shape = model.type))+
  facet_wrap(~COMMON_NAME)+
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey")+
  theme_bw(base_size = 16)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        panel.grid = element_blank())+
  ylab("Bias in Regional-scale mortality rate (% trees/year)")+
  xlab(" Model Number ")
ggsave(paste0(output.dir, "SPCD_stanoutput_cmdstan/summary/Regional_Species_BIAS.png"), 
       height = 8, width = 11)


ggplot(data = Regional_RMSE_BIAS)+
  geom_point(aes(model.number, RMSE, color = model.type, shape = model.type))+
  facet_wrap(~COMMON_NAME)+
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey")+
  theme_bw(base_size = 16)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        panel.grid = element_blank())+
  ylab("RMSE in Regional-scale mortality rate (% trees/year)")+
  xlab(" Model Number ")

ggsave(paste0(output.dir, "SPCD_stanoutput_cmdstan/summary/Regional_Species_RMSE.png"), 
       height = 6, width = 9)


ggplot(data = Regional_RMSE_BIAS)+
  #geom_point(aes( RMSE, BIAS, color = model.number, shape = model.type))+
  geom_text(aes( RMSE, BIAS, label = model.number, color = model.type))+
  facet_wrap(~COMMON_NAME, scales = "free")+
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey")+
  theme_bw(base_size = 16)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        panel.grid = element_blank())+
  ylab("BIAS in Regional-scale mortality rate (% trees/year)")+
  xlab("RMSE in Regional-scale mortality rate (% trees/year)")

ggsave(paste0(output.dir, "SPCD_stanoutput_cmdstan/summary/State_Species_RMSE_BIAS.png"), 
       height = 6, width = 9)


# state scale RMSE and BIAS---
State_RMSE_BIAS <- State_summary_df %>% filter(n_obs > 25)%>%
  mutate(residual_diff = obs_M_expn_median*100 - pred_M_expn_median*100)%>%
  group_by(model.number, model.type, COMMON_NAME)%>%
  summarise(BIAS = sum(residual_diff, na.rm = TRUE), 
            RMSE = sqrt(mean(residual_diff^2, na.rm = TRUE)))

ggplot(data = State_RMSE_BIAS)+
  geom_point(aes(RMSE, BIAS, shape = model.number, color = model.type))+
  facet_wrap(~COMMON_NAME, scales = "free")

ggplot(data = State_RMSE_BIAS)+
  geom_point(aes(model.number, BIAS, color = model.type, shape = model.type))+
  facet_wrap(~COMMON_NAME)+
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey")+
  theme_bw(base_size = 16)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        panel.grid = element_blank())+
  ylab("Bias in State-scale mortality rate (% trees/year)")+
  xlab(" Model Number ")

ggsave(paste0(output.dir, "SPCD_stanoutput_cmdstan/summary/State_Species_BIAS.png"), 
       height = 8, width = 11)

ggplot(data = State_RMSE_BIAS)+
  geom_point(aes(model.number, RMSE, color = model.type, shape = model.type))+
  facet_wrap(~COMMON_NAME)+
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey")+
  theme_bw(base_size = 16)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        panel.grid = element_blank())+
  ylab("RMSE in State-scale mortality rate (% trees/year)")+
  xlab(" Model Number ")

ggsave(paste0(output.dir, "SPCD_stanoutput_cmdstan/summary/State_Species_RMSE.png"), 
       height = 6, width = 9)


ggplot(data = State_RMSE_BIAS)+
  #geom_point(aes( RMSE, BIAS, color = model.number, shape = model.type))+
  geom_text(aes( RMSE, BIAS, label = model.number, color = model.type))+
  facet_wrap(~COMMON_NAME, scales = "free")+
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey")+
  theme_bw(base_size = 16)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        panel.grid = element_blank())+
  ylab("BIAS in State-scale mortality rate (% trees/year)")+
  xlab("RMSE in State-scale mortality rate (% trees/year)")

ggsave(paste0(output.dir, "SPCD_stanoutput_cmdstan/summary/State_Species_RMSE_BIAS.png"), 
       height = 6, width = 9)


#  County scale RMSE and BIAS---
# calculate RMSE & BIAS for each model, species, scale
ST_CTY_RMSE_BIAS <- ST_CTY_summary_df %>% filter(n_obs > 25)%>%
  mutate(residual_diff = obs_M_expn_median*100 - pred_M_expn_median*100)%>%
  group_by(model.number, model.type, COMMON_NAME)%>%
  summarise(BIAS = sum(residual_diff, na.rm = TRUE), 
            RMSE = sqrt(mean(residual_diff^2, na.rm = TRUE)))



ggplot(data = ST_CTY_RMSE_BIAS)+
  geom_point(aes(model.number, BIAS, color = model.type, shape = model.type))+
  facet_wrap(~COMMON_NAME)+
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey")+
  theme_bw(base_size = 16)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        panel.grid = element_blank())+
  ylab("Bias in County-scale mortality rate (% trees/year)")+
  xlab(" Model Number ")

ggsave(paste0(output.dir, "SPCD_stanoutput_cmdstan/summary/ST_CTY_Species_BIAS.png"), 
       height = 6, width = 9)

ggplot(data = ST_CTY_RMSE_BIAS)+
  geom_point(aes(model.number, RMSE, color = model.type, shape = model.type))+
  facet_wrap(~COMMON_NAME)+
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey")+
  theme_bw(base_size = 16)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        panel.grid = element_blank())+
  ylab("RMSE in County-scale mortality rate (% trees/year)")+
  xlab(" Model Number ")

ggsave(paste0(output.dir, "SPCD_stanoutput_cmdstan/summary/ST_CTY_Species_RMSE.png"), 
       height = 6, width = 9)


ggplot(data = ST_CTY_RMSE_BIAS)+
  #geom_point(aes( RMSE, BIAS, color = model.number, shape = model.type))+
  geom_text(aes( RMSE, BIAS, label = model.number, color = model.type))+
  facet_wrap(~COMMON_NAME, scales = "free")+
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey")+
  theme_bw(base_size = 16)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        panel.grid = element_blank())+
  ylab("BIAS in ST_CTY-scale mortality rate (tree mortality/year)")+
  xlab("RMSE in ST_CTY-scale mortality rate (tree mortality/year)")

ggsave(paste0(output.dir, "SPCD_stanoutput_cmdstan/summary/ST_CTY_Species_RMSE_BIAS.png"), 
       height = 6, width = 9)


# output a table for stacking

#--------------------------------------------------------------------------------
# Map of predicted mortality rates with residual from calculated ---
#--------------------------------------------------------------------------------
# join up with state and county -level data
ST_CTY_FP <- ST_CTY_summary_df %>% 
  separate(ST_CTY, into = c("state", "county"), remove = F)%>%
mutate(COUNTYFP = str_pad(county, side = "left",3, pad = 0), 
       STATEFP = str_pad(state, side = "left", 2, pad = 0))

State_FP <- State_summary_df %>% 
  mutate(STATEFP = str_pad(state, side = "left", 2, pad = 0))


# join up with state and county shapefiles for plotting:
library(tigris)
# Get the most recent county data for all the states
counties <- counties(state = c("DE","RI","ME", "MA","MD", "NH", "NJ", "NY", "OH", "PA", "VT", "WV"), cb = TRUE) %>% 
  data.frame()%>% 
  rename(
    "CONAME" = "NAME")  %>%
  dplyr::select("STATEFP",
                "COUNTYFP","COUNTYNS" , "GEOID" ,"CONAME",    
                "NAMELSAD","STUSPS","STATE_NAME", "LSAD","ALAND","AWATER",    
                "geometry" )%>% as.data.frame()
# separate handling for CT, which changed their county names around 2017
CT.counties <- counties(state = "CT", cb = FALSE, year = 2015) %>% 
  rename("CONAME" = "NAME") %>%
  mutate(STATE_NAME = "Connecticut", 
         STUSPS = "CT") %>%
  dplyr::select("STATEFP",
                "COUNTYFP","COUNTYNS" , "GEOID" ,"CONAME",    
                "NAMELSAD","STUSPS","STATE_NAME", "LSAD","ALAND","AWATER",    
                "geometry" )%>% as.data.frame()

# join CT and join up mortality data
counties <- counties %>% rbind(., CT.counties)
id.var.counties <- c("STATEFP",
                     "COUNTYFP","COUNTYNS" , "GEOID" ,"CONAME",    
                     "NAMELSAD","STUSPS","STATE_NAME", "LSAD","ALAND","AWATER",    
                     "geometry" )

# get state-level geometry form tigris

# Get the most recent county data for all the states
states_sf <- states(cb = TRUE) 

states_sf


###################################################################################################
# plot out predicted maps by species, with 95% CI, observed rates, and residuals ---
#

SPCD.id <- 531

n_obs_threshold <- 5

map_stacking_pred_obs <- function(SPCD.id, n_obs_threshold = 5, use_EXPN = TRUE){
  
  Common_name_spp <- spp.table[spp.table$SPCD.id %in% SPCD.id, ]$COMMON
      # get species-level stacked estimates for all counties with more than the n_obs_threshold trees 
      ST_CTY_FP_SPECIES <- ST_CTY_FP %>% 
        filter(SPCD %in% SPCD.id) %>% 
        filter(model.number %in% "Stacked")%>%
        filter(n_obs >= n_obs_threshold)
      
      counties_SPP <- left_join(counties, ST_CTY_FP_SPECIES) %>% st_as_sf()
      
      # if use_EXPN == TRUE, then use the predicted and observed values on the trees/year scale
    if(use_EXPN == TRUE){
      
      counties_SPP_long <- counties_SPP%>%
        select(STATEFP:county, geometry, obs_M_expn_median, pred_M_expn_median, pred_M_expn_2.5.ci.lo, pred_M_expn_97.5.ci.hi) %>%
        pivot_longer(
          cols = c(obs_M_expn_median, pred_M_expn_median, pred_M_expn_2.5.ci.lo, pred_M_expn_97.5.ci.hi),
          names_to = "Pred_obs_type",
          values_to = "Mort_rate"
        ) %>% 
        mutate(`Pct_mort_year` = Mort_rate*100) %>%
        mutate(`Pred_obs_names` = ifelse(Pred_obs_type %in% "obs_M_expn_median", "Observed", 
                                         ifelse(Pred_obs_type %in% "pred_M_expn_median", "Predicted (median)", 
                                                ifelse(Pred_obs_type %in% "pred_M_expn_2.5.ci.lo", "Predicted (2.5% C.I.)", 
                                                       ifelse(Pred_obs_type %in% "pred_M_expn_97.5.ci.hi", "Predicted (97.5% C.I.)", NA)))))
      
      
       #labels for maps
      units_label <- "(% trees/year)"
      fn_add_on <- "_pct_Trees_yr"
      
      
      
      # get residuals:
      
      counties_residual <- counties_SPP %>%
        select(STATEFP:county, n_obs, geometry, obs_M_expn_median, pred_M_expn_median, pred_M_expn_2.5.ci.lo, pred_M_expn_97.5.ci.hi) %>%
        mutate(residual_M = obs_M_expn_median*100 - pred_M_expn_median*100, 
               residual_M_2.5.ci.lo = obs_M_expn_median*100 - pred_M_expn_2.5.ci.lo*100, 
               residual_M_97.5.ci.hi = obs_M_expn_median*100 - pred_M_expn_97.5.ci.hi*100)%>%
        mutate(n_obs_group = ifelse(n_obs > 25, ">= 25", " < 25"))%>%
        mutate(obs_pct = obs_M_expn_median*100, 
               pred_pct = pred_M_expn_median*100, 
               pred_ci.lo = pred_M_expn_2.5.ci.lo*100, 
               pred_ci.hi = pred_M_expn_97.5.ci.hi*100)
      
      
      
      }else{
        # if use_EXPN == FALSE, then use the predicted and observed values on the per acre scale
        # trees/acre/year scale
        counties_SPP_long <- counties_SPP%>%
          select(STATEFP:county, geometry, obs_M_median, 
                 pred_M_median, pred_M_2.5.ci.lo, pred_M_97.5.ci.hi) %>%
          pivot_longer(
            cols = c(obs_M_median, pred_M_median, pred_M_2.5.ci.lo, pred_M_97.5.ci.hi),
            names_to = "Pred_obs_type",
            values_to = "Mort_rate"
          ) %>% 
          mutate(`Pct_mort_year` = Mort_rate*100) %>%
          mutate(`Pred_obs_names` = ifelse(Pred_obs_type %in% "obs_M_median", "Observed", 
                                           ifelse(Pred_obs_type %in% "pred_M_median", "Predicted (median)", 
                                                  ifelse(Pred_obs_type %in% "pred_M_2.5.ci.lo", "Predicted (2.5% C.I.)", 
                                                         ifelse(Pred_obs_type %in% "pred_M_97.5.ci.hi", "Predicted (97.5% C.I.)", NA)))))
        
        # labels for maps
        units_label <- "(% trees/acre/year)"
        fn_add_on <- "_pct_TPA_yr"
        
        
        
        # get residuals:
        
        counties_residual <- counties_SPP %>%
          select(STATEFP:county, n_obs, geometry, obs_M_median, pred_M_median, pred_M_2.5.ci.lo, pred_M_97.5.ci.hi) %>%
          mutate(residual_M = obs_M_median*100 - pred_M_median*100, 
                 residual_M_2.5.ci.lo = obs_M_median*100 - pred_M_2.5.ci.lo*100, 
                 residual_M_97.5.ci.hi = obs_M_median*100 - pred_M_97.5.ci.hi*100)%>%
          mutate(n_obs_group = ifelse(n_obs > 25, ">= 25", " < 25"))%>%
          mutate(obs_pct = obs_M_median*100, 
                 pred_pct = pred_M_median*100, 
                 pred_ci.lo = pred_M_2.5.ci.lo*100, 
                 pred_ci.hi = pred_M_97.5.ci.hi*100)
        
      }
      # reorder the facet labels
      counties_SPP_long$Pred_obs_names <- factor(counties_SPP_long$Pred_obs_names, levels = unique(counties_SPP_long$Pred_obs_names))
      
      p_o_maps <- ggplot(data = counties_SPP_long) +
                        geom_sf(aes(fill = Pct_mort_year)) +
                        facet_wrap(~ Pred_obs_names, ncol = 4) +
                        scale_fill_viridis_b(option = "inferno", 
                          breaks = c(0, 0.1, 0.25, 0.5, 0.75, 1, 1.5, 2),
                          name = paste0("Mortality Rate\n",units_label)) + 
                        theme_minimal(base_size = 12)+
                        theme(legend.position = "bottom", 
                              legend.direction = "horizontal", 
                              legend.text = element_text(angle = 90))+
                        ggtitle(paste(Common_name_spp))
      
      ggsave(filename = paste0(output.dir, "SPCD_stanoutput_cmdstan/summary/pred_obs_maps/Stacking_p.o.maps_SPCD_", SPCD.id, fn_add_on,".png"), 
             plot = p_o_maps ,
             height = 4, width = 11)
      
    
      res_map <- ggplot(data = counties_residual) +
        geom_sf(aes(fill = residual_M)) +
        scale_fill_fermenter(palette = "RdBu",
          breaks = c( -1.5, -1, -0.5, -0.25, 0, 0.25, 0.5, 1, 1.5), 
          name = paste0("Residuals " ,units_label,"\n(Observed - Predicted)"))+ 
        theme_minimal(base_size = 12)+
        theme(legend.position = "bottom", 
              legend.direction = "horizontal", 
              legend.text = element_text(angle = 90, vjust = 0.5))+
        ggtitle(paste(Common_name_spp))
      
      
      
      
      RMSE <- counties_residual %>% filter(!is.na(residual_M)) %>% st_drop_geometry()%>%
        mutate(Squared_residuals = residual_M^2) %>% 
        group_by(n_obs_group)%>%
        summarise(#MAE = mean(abs(residual_expn), na.rm = TRUE), 
                  #MSE = mean(Squared_residuals, na.rm = TRUE), 
                  RMSE = sqrt(mean(Squared_residuals, na.rm = TRUE)), 
                  n_counties = n())%>%
        bind_rows(
          counties_residual %>% filter(!is.na(residual_M)) %>% st_drop_geometry()%>%
            mutate(Squared_residuals = residual_M^2, 
                   n_obs_group = "Overall") %>% 
            group_by(n_obs_group)%>%
            summarise(#MAE = mean(abs(residual_expn), na.rm = TRUE), 
              #MSE = mean(Squared_residuals, na.rm = TRUE), 
              RMSE = sqrt(mean(Squared_residuals, na.rm = TRUE)), 
              n_counties = n())
        ) %>%
        mutate(plt_label = paste0("RMSE (", n_obs_group, ") = ", round(RMSE, digits = 2)))%>%
        arrange(desc(n_obs_group))
      
      
      
      # counties_residual %>% st_drop_geometry()%>%
      #   filter(!is.na(residual_M))|>
      #   ggplot()+geom_dots(aes(x = residual_M, y = STATE_NAME),layout = "weave", size =2)
      # 
      # 
      # counties_residual %>% st_drop_geometry()%>%
      #   filter(!is.na(residual_M))|>
      #   ggplot()+geom_dots(aes(x = residual_M),layout = "weave", size = 2)
      
      p.o.plot <- ggplot(data = counties_residual %>% filter(!is.na(pred_pct))) +
        geom_pointrange(aes(x = obs_pct, 
                            y = pred_pct,
                            ymin = pred_ci.lo, 
                            ymax = pred_ci.hi, 
                            color = n_obs_group), alpha = 0.75) +
        scale_color_manual(values = c(">= 25" = "darkred", 
                                      " < 25" = "darkgrey"), name = "# trees included")+
        geom_abline(aes(intercept = 0, slope = 1), linetype = "dashed")+
        
        annotate("text", x = -Inf, y = Inf, label = RMSE[RMSE$n_obs_group %in% "Overall",]$plt_label, color = "black", 
                 hjust = -0.1, vjust = 1.5, size = 4)+
        annotate("text", x = -Inf, y = Inf, label = RMSE[RMSE$n_obs_group %in% " < 25",]$plt_label, color = "darkgrey", 
                 hjust = -0.1, vjust = 3, size = 4)+
        annotate("text", x = -Inf, y = Inf, label = RMSE[RMSE$n_obs_group %in% ">= 25",]$plt_label, color = "darkred", 
                 hjust = -0.1, vjust = 4.5, size = 4)+
        
        #theme_minimal(base_size = 12)+
        theme_bw(base_size = 12)+
        theme(legend.position = "bottom",  
              legend.direction = "horizontal", 
              axis.line.x.top = element_blank() )+
        ylab(paste0("Predicted County Mortality ", units_label))+
        xlab(paste0("Observed County Mortality ", units_label))+
        ggtitle(paste(Common_name_spp))
      
      
      
      ggsave(filename = paste0(output.dir, "SPCD_stanoutput_cmdstan/summary/pred_obs_maps/Stacking_residual.map_SPCD_", SPCD.id, fn_add_on,".png"), 
             plot =  res_map , height = 4, width = 4)
      
      ggsave(filename = paste0(output.dir, "SPCD_stanoutput_cmdstan/summary/pred_obs_maps/Stacking_p.o_scatter_SPCD_", SPCD.id, fn_add_on,".png"), 
             plot =  p.o.plot , height = 4, width = 6)

}

# make maps for trees/yr expn factor
lapply(spp.table$SPCD.id, map_stacking_pred_obs)

# make maps for trees/acre/yr expn factor
lapply(spp.table$SPCD.id, FUN = function(x)map_stacking_pred_obs(x, n_obs_threshold = 5, use_EXPN = FALSE))


# create paired maps and barplots for state-level estimates
map_stacking_STATE_pred_obs <- function(SPCD.id, n_obs_threshold = 5, use_EXPN = TRUE){
  
  Common_name_spp <- spp.table[spp.table$SPCD.id %in% SPCD.id, ]$COMMON
  # get species-level stacked estimates for all states with more than the n_obs_threshold trees 
  State_FP_SPECIES <- State_FP %>% 
    filter(SPCD %in% SPCD.id) %>% 
    filter(model.number %in% "Stacked")%>%
    filter(n_obs >= n_obs_threshold)
  
  states_SPP <- left_join(states_sf, State_FP_SPECIES) %>% st_as_sf()
  
  # if use_EXPN == TRUE, then use the predicted and observed values on the trees/year scale
  if(use_EXPN == TRUE){
    
    states_SPP_long <- states_SPP%>%
      select(STATEFP:state, geometry, obs_M_expn_median, pred_M_expn_median, pred_M_expn_2.5.ci.lo, pred_M_expn_97.5.ci.hi) %>%
      pivot_longer(
        cols = c(obs_M_expn_median, pred_M_expn_median, pred_M_expn_2.5.ci.lo, pred_M_expn_97.5.ci.hi),
        names_to = "Pred_obs_type",
        values_to = "Mort_rate"
      ) %>% 
      mutate(`Pct_mort_year` = Mort_rate*100) %>%
      mutate(`Pred_obs_names` = ifelse(Pred_obs_type %in% "obs_M_expn_median", "Observed", 
                                       ifelse(Pred_obs_type %in% "pred_M_expn_median", "Predicted (median)", 
                                              ifelse(Pred_obs_type %in% "pred_M_expn_2.5.ci.lo", "Predicted (2.5% C.I.)", 
                                                     ifelse(Pred_obs_type %in% "pred_M_expn_97.5.ci.hi", "Predicted (97.5% C.I.)", NA)))))
    
    
    #labels for maps
    units_label <- "(% trees/year)"
    fn_add_on <- "_pct_Trees_yr"
    
    
    
    # get residuals:
    
    states_residual <- states_SPP %>%
      #select(STATEFP:state, n_obs, geometry, obs_M_expn_median, pred_M_expn_median, pred_M_expn_2.5.ci.lo, pred_M_expn_97.5.ci.hi) %>%
      mutate(residual_M = obs_M_expn_median*100 - pred_M_expn_median*100, 
             residual_M_2.5.ci.lo = obs_M_expn_median*100 - pred_M_expn_2.5.ci.lo*100, 
             residual_M_97.5.ci.hi = obs_M_expn_median*100 - pred_M_expn_97.5.ci.hi*100)%>%
      mutate(n_obs_group = ifelse(n_obs > 25, ">= 25", " < 25"))%>%
      mutate(obs_pct = obs_M_expn_median*100, 
             pred_pct = pred_M_expn_median*100, 
             pred_ci.lo = pred_M_expn_2.5.ci.lo*100, 
             pred_ci.hi = pred_M_expn_97.5.ci.hi*100, 
             pred_25ci.lo = pred_M_expn_25.ci.lo*100, 
             pred_75ci.hi = pred_M_expn_75.ci.hi*100)
    
    
    
  }else{
    # if use_EXPN == FALSE, then use the predicted and observed values on the per acre scale
    # trees/acre/year scale
    states_SPP_long <- states_SPP %>%
      select(STATEFP:state, geometry, obs_M_median, 
             pred_M_median, pred_M_2.5.ci.lo, pred_M_97.5.ci.hi) %>%
      pivot_longer(
        cols = c(obs_M_median, pred_M_median, pred_M_2.5.ci.lo, pred_M_97.5.ci.hi),
        names_to = "Pred_obs_type",
        values_to = "Mort_rate"
      ) %>% 
      mutate(`Pct_mort_year` = Mort_rate*100) %>%
      mutate(`Pred_obs_names` = ifelse(Pred_obs_type %in% "obs_M_median", "Observed", 
                                       ifelse(Pred_obs_type %in% "pred_M_median", "Predicted (median)", 
                                              ifelse(Pred_obs_type %in% "pred_M_2.5.ci.lo", "Predicted (2.5% C.I.)", 
                                                     ifelse(Pred_obs_type %in% "pred_M_97.5.ci.hi", "Predicted (97.5% C.I.)", NA)))))
    
    # labels for maps
    units_label <- "(% trees/acre/year)"
    fn_add_on <- "_pct_TPA_yr"
    
    
    
    # get residuals:
    
    states_residual <- states_SPP %>%
     # select(STATEFP:state, n_obs, geometry, obs_M_median, pred_M_median, pred_M_2.5.ci.lo, pred_M_97.5.ci.hi) %>%
      mutate(residual_M = obs_M_median*100 - pred_M_median*100, 
             residual_M_2.5.ci.lo = obs_M_median*100 - pred_M_2.5.ci.lo*100, 
             residual_M_97.5.ci.hi = obs_M_median*100 - pred_M_97.5.ci.hi*100)%>%
      mutate(n_obs_group = ifelse(n_obs > 25, ">= 25", " < 25"))%>%
      mutate(obs_pct = obs_M_median*100, 
             pred_pct = pred_M_median*100, 
             pred_ci.lo = pred_M_2.5.ci.lo*100, 
             pred_ci.hi = pred_M_97.5.ci.hi*100, 
             pred_25ci.lo = pred_M_25.ci.lo*100, 
             pred_75ci.hi = pred_M_75.ci.hi*100)
    
  }
  # reorder the facet labels
  states_SPP_long$Pred_obs_names <- factor(states_SPP_long$Pred_obs_names, levels = unique(states_SPP_long$Pred_obs_names))
  
  map_bbox <- st_bbox(counties_SPP)
  
  
  p_o_maps <- ggplot(data = states_SPP_long) +
    geom_sf(aes(fill = Pct_mort_year)) +
    facet_wrap(~ Pred_obs_names, ncol = 4) +
    scale_fill_viridis_b(option = "inferno", 
                         breaks = c(0, 0.1, 0.25, 0.5, 0.75, 1, 1.5, 2),
                         name = paste0("Mortality Rate\n",units_label)) + 
    theme_minimal(base_size = 12)+
    theme(legend.position = "bottom", 
          legend.direction = "horizontal", 
          legend.text = element_text(angle = 90))+
    ggtitle(paste(Common_name_spp))+
    coord_sf(xlim = c(map_bbox["xmin"],map_bbox["xmax"]),
             ylim = c(map_bbox["ymin"], map_bbox["ymax"]))
  
  ggsave(filename = paste0(output.dir, "SPCD_stanoutput_cmdstan/summary/pred_obs_maps/Stacking_p.o.STATE_maps_SPCD_", SPCD.id, fn_add_on,".png"), 
         plot = p_o_maps ,
         height = 4, width = 11)
  
  
  res_map <- ggplot(data = states_residual) +
    geom_sf(aes(fill = residual_M)) +
    scale_fill_fermenter(palette = "RdBu",
                         breaks = c( -1.5, -1, -0.5, -0.25, 0, 0.25, 0.5, 1, 1.5), 
                         name = paste0("Residuals " ,units_label,"\n(Observed - Predicted)"))+ 
    theme_minimal(base_size = 12)+
    theme(legend.position = "bottom", 
          legend.direction = "horizontal", 
          legend.text = element_text(angle = 90, vjust = 0.5))+
    ggtitle(paste(Common_name_spp))+
    coord_sf(xlim = c(map_bbox["xmin"],map_bbox["xmax"]),
             ylim = c(map_bbox["ymin"], map_bbox["ymax"]))
  
  
  
  
  RMSE <- states_residual %>% filter(!is.na(residual_M)) %>% st_drop_geometry()%>%
    mutate(Squared_residuals = residual_M^2) %>% 
    group_by(n_obs_group)%>%
    summarise(#MAE = mean(abs(residual_expn), na.rm = TRUE), 
      #MSE = mean(Squared_residuals, na.rm = TRUE), 
      RMSE = sqrt(mean(Squared_residuals, na.rm = TRUE)), 
      n_states = n())%>%
    bind_rows(
      states_residual %>% filter(!is.na(residual_M)) %>% st_drop_geometry()%>%
        mutate(Squared_residuals = residual_M^2, 
               n_obs_group = "Overall") %>% 
        group_by(n_obs_group)%>%
        summarise(
          RMSE = sqrt(mean(Squared_residuals, na.rm = TRUE)), 
          n_states = n())
    ) %>%
    mutate(plt_label = paste0("RMSE (", n_obs_group, ") = ", round(RMSE, digits = 2)))%>%
    arrange(desc(n_obs_group))
  
  
  
 
  p.o.plot <- ggplot(data = states_residual %>% filter(!is.na(pred_pct))) +
    geom_pointrange(aes(x = obs_pct, 
                        y = pred_pct,
                        ymin = pred_ci.lo, 
                        ymax = pred_ci.hi, 
                        color = NAME), alpha = 0.75) +
    scale_color_manual(values = c("New York" = "#a6cee3" ,
     "Maine" = "#1f78b4" ,
     "New Hampshire" = "#b2df8a"  ,
      "Vermont" = "#33a02c" ,
      "Pennsylvania" = "#fb9a99" ,
       "Ohio" = "#e31a1c",
       "New Jersey" = "#fdbf6f" ,
      "West Virginia" =  "#ff7f00" ,
       "Connecticut" = "#cab2d6" ,
      "Maryland" = "#6a3d9a" ), name = "")+
    geom_abline(aes(intercept = 0, slope = 1), linetype = "dashed")+
    
    annotate("text", x = -Inf, y = Inf, label = RMSE[RMSE$n_obs_group %in% "Overall",]$plt_label, color = "black", 
             hjust = -0.1, vjust = 1.5, size = 4)+
    theme_bw(base_size = 16)+
    theme(#legend.position = "bottom",  
          #legend.direction = "horizontal", 
          axis.line.x.top = element_blank(), 
          panel.grid.minor = element_blank())+
    ylab(paste0("Predicted State Mortality\n", units_label))+
    xlab(paste0("Observed State Mortality ", units_label))+
    ggtitle(paste(Common_name_spp))
  
  
  
  ggsave(filename = paste0(output.dir, "SPCD_stanoutput_cmdstan/summary/pred_obs_maps/Stacking_residual_STATE.map_SPCD_", SPCD.id, fn_add_on,".png"), 
         plot =  res_map , height = 4, width = 4)
  
  ggsave(filename = paste0(output.dir, "SPCD_stanoutput_cmdstan/summary/pred_obs_maps/Stacking_p.o_STATE_scatter_SPCD_", SPCD.id, fn_add_on,".png"), 
         plot =  p.o.plot , height = 6, width = 8)
  
}

lapply(spp.table$SPCD.id, FUN = function(x)map_stacking_STATE_pred_obs(x, n_obs_threshold = 5, use_EXPN = FALSE))

rm(p_o_maps, p.o.plot, p.o_stacked_expanded, p.o_stacked_per_acre, res.map)
gc()

######################################################################################
# Plot up beta posterior predictions from all models ----------

# read_qs
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
# get covariate names
cov.names <- data.frame(cov.number = unique(u_beta_names$cov.number),
                        Covariate = colnames(mod.data$xM))

spp.table <- read.csv(file = paste0(output.dir, "/data/Hierarchical_obs_model_7.csv"))

u_betas.quant <- u_betas %>% rename("ci.lo"="q5", "ci.hi" = "q95")%>% 
  left_join(., u_beta_names)%>% left_join(., spp.table) %>%
  left_join(., cov.names)
u_betas.quant$Covariate <- factor(u_betas.quant$Covariate, levels = cov.names$Covariate)


# get overlapping zero to color the error bars
u_betas.quant$`significance` <- ifelse(u_betas.quant$ci.lo < 0 & u_betas.quant$ci.hi < 0, "significant", 
                                             ifelse(u_betas.quant$ci.lo > 0 & u_betas.quant$ci.hi > 0, "significant", "overlapping zero"))


# get the species.level u_betas----

species_draws_path <- function(spcd, k){
 
  paste0(output.dir, "SPCD_stanoutput_cmdstan/betas/u_beta_alpha_samps_mort_model_",k,"_SPCD_",spcd,"_remper_correction_0.5_niter_1000_nchain_4.qs")
}


spp_beta_sumary_df <- do.call(rbind, lapply(1:9, function(i){
  
  
  
    ubetas <- do.call(rbind, lapply(spp.table$SPCD, FUN = function(spcd, k = i){
      qs2::qs_read(species_draws_path(spcd = spcd, k)) %>% subset_draws(., "u_beta")%>%  
        summarise_draws() %>% 
        mutate(model.number = i, 
             SPCD = spcd) %>% 
      mutate(param = "u_beta")
    }))
  ubetas

})
)

# species u_beta summaries
u_betas_spp <- spp_beta_sumary_df %>% filter(param %in% "u_beta")

u_beta_spp_names <- data.frame(variable = unique(u_betas_spp$variable), 
                           #SPP = rep(1:17, 78),
                           cov.number = 1:78)
# get covariate names
cov.names_spp <- data.frame(cov.number = unique(u_beta_spp_names$cov.number),
                        Covariate = colnames(mod.data$xM))


u_betas_spp.quant <- u_betas_spp %>% rename("ci.lo"="q5", "ci.hi" = "q95")%>% 
  left_join(., u_beta_spp_names)%>% left_join(., spp.table) %>%
  left_join(., cov.names_spp)
u_betas_spp.quant$Covariate <- factor(u_betas_spp.quant$Covariate, levels = cov.names_spp$Covariate)


# get overlapping zero to color the error bars
u_betas_spp.quant$`significance` <- ifelse(u_betas_spp.quant$ci.lo < 0 & u_betas_spp.quant$ci.hi < 0, "significant", 
                                       ifelse(u_betas_spp.quant$ci.lo > 0 & u_betas_spp.quant$ci.hi > 0, "significant", "overlapping zero"))



# combine the species and hierarchical summaries together:
u_betas_all_models <- rbind(
    u_betas_spp.quant %>%  select(SPCD, COMMON, SPP, model.number, Covariate, median, ci.lo, ci.hi, significance) %>%
        mutate(model.type = "Species"), 


      u_betas.quant %>% 
      select(SPCD, COMMON, SPP, model.number, Covariate, median, ci.lo, ci.hi, significance) %>%
           mutate(model.type = "Hierarchical") 
)

u_betas_signifcance_summary <- u_betas_all_models  %>% 
  group_by(SPCD, COMMON, SPP, Covariate) %>% 
  summarise(total_num = n(), 
            num_sig = sum(significance %in% "significant"), 
            num_pos.sig = sum(significance %in% "significant" & median > 0), 
            num_neg.sig = sum(significance %in% "significant" & median < 0))%>%
  mutate(frac_sig = num_sig/total_num, 
         frac_pos.sig = num_pos.sig/total_num, 
         frac_neg.sig = num_neg.sig/total_num)

u_betas_signifcance_summary |>
  ggplot()+
  geom_bar(aes(x = Covariate, y = num_sig, fill = COMMON), stat = "identity", position = "stack")+
  species_fill+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


u_betas_signifcance_summary |>
  ggplot()+
  geom_bar(aes(x = Covariate, y = frac_sig/17, fill = COMMON), stat = "identity", position = "stack")+
  species_fill+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

u_betas_signifcance_summary %>% filter(Covariate %in% main_effects)|>
  ggplot()+
  geom_bar(aes(x = Covariate, y = num_sig, fill = COMMON), stat = "identity", position = "stack")+
  species_fill+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

u_betas_signifcance_summary %>% filter(Covariate %in% main_effects)|>
  ggplot()+
  geom_bar(aes(x = Covariate, y = frac_sig/17, fill = COMMON), stat = "identity", position = "stack")+
  species_fill+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))



u_betas_signifcance_summary |>
  ggplot()+
  geom_bar(aes(x = Covariate, y = num_pos.sig, fill = COMMON), stat = "identity", position = "stack")+
  geom_bar(aes(x = Covariate, y = -num_neg.sig, fill = COMMON), stat = "identity", position = "stack")+
  species_fill+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

u_betas_all_models  %>% filter(COMMON %in% "sugar maple") |>
  ggplot()+geom_point(aes(x = model.number, y = median, color = model.type, group = COMMON))+
  facet_wrap(~Covariate)
 



u_betas_all_models %>% filter(Covariate %in% c("DIA_DIFF_scaled", "DIA_scaled")) |>
  ggplot()+geom_pointrange(aes(x = as.character(model.number), y = median, ymin = ci.lo, ymax = ci.hi, group = COMMON, color = significance, shape = model.type))+
  facet_grid(vars(Covariate), vars(COMMON))

u_betas_all_models %>% filter(Covariate %in% c("DIA_DIFF_scaled", "DIA_scaled")) |>
  ggplot()+geom_pointrange(aes(x = as.character(model.number), y = median, ymin = ci.lo, ymax = ci.hi, group = COMMON, color = significance))+
  facet_grid(vars(Covariate), vars(COMMON))


main_effects <- c("DIA_DIFF_scaled", "DIA_scaled", "ba.scaled", "BAL.scaled",
                  "damage.scaled", "MATmax.scaled", "MAP.scaled", "ppt.anom",
                  "tmax.anom", "slope.scaled", "aspect.scaled", "Ndep.scaled")

u_betas_all_models %>% filter(Covariate %in% main_effects) |>
  ggplot()+geom_pointrange(aes(x = as.character(model.number), y = median, ymin = ci.lo, ymax = ci.hi, group = COMMON, color = significance))+
  facet_grid(vars(Covariate), vars(COMMON), scales = "free_y")



ggplot(data = na.omit(u_betas_all_models), aes(x = COMMON, y = median, color = significance))+geom_point()+
  geom_errorbar(data = na.omit(u_betas_all_models), aes(x = COMMON , ymin = ci.lo, ymax = ci.hi, color = significance), width = 0.1)+
  geom_abline(aes(slope = 0, intercept = 0), color = "grey", linetype = "dashed")+
  facet_wrap(~Covariate, scales= "free_y")+
  theme_bw(base_size = 14)+
  theme( axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.grid  = element_blank(), legend.position = "none")+
  ylab("Effect on Survival")+xlab("Parameter")+
  scale_color_manual(values = c("not overlapping zero"="darkgrey", "significant"="black"))



ggplot(data = na.omit(u_betas_all_models), aes(x = Covariate, y = median, color = significance))+geom_jitter()+
  #geom_errorbar(data = na.omit(u_betas_all_models), aes(x = COMMON , ymin = ci.lo, ymax = ci.hi, color = significance), width = 0.1)+
  geom_abline(aes(slope = 0, intercept = 0), color = "grey", linetype = "dashed")+
  facet_wrap(~COMMON, scales= "free_y")+
  theme_bw(base_size = 14)+
  theme( axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.grid  = element_blank(), legend.position = "none")+
  ylab("Effect on Survival")+xlab("Parameter")+
  scale_color_manual(values = c("not overlapping zero"="darkgrey", "significant"="black"))




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
#ggsave(height = 10, width = 15,dpi = 350, units = "in",paste0(output.folder,"SPCD_stanoutput_joint_v3/images/Estimated_effects_on_mortality_model_model_",model.no,"_all_species_betas.png"))



######################################################################################
# Marginal stacked posterior predictions ----------

marg_draws_all_list <- lapply(list.files(paste0(output.dir, "SPCD_stanoutput_cmdstan/model_avg_draws/stacking_wts/"), pattern = "Draws_marginal_main", full.names = T), 
       qs2::qs_read)

# expanding into a very big dataframe
marg_draws_all_df <- do.call(rbind, marg_draws_all_list)

# arrange the species--
marg_draws_all_df$COMMON_NAME <- factor(marg_draws_all_df$COMMON_NAME, levels = c(SP.TRAITS$COMMON_NAME))

# make draws summary plots for different effects to identify which models
# pull out species groups:
n_conifers <- c("balsam fir", "red spruce", "northern white-cedar")
Oak_birch_hickory <- c("white oak", "northern red oak", "chestnut oak","hickory spp.",
                      "American beech", "yellow birch","paper birch", "yellow-poplar")

Maples_ash_bcherry <- c("sugar maple", "red maple", "eastern white pine", 
                        "eastern hemlock","white ash", "black cherry")

plot_draws_predictors <- function( species_names, group_name,  y_max = 1){
      
        size.p <- marg_draws_all_df %>% filter(predictor %in% c("DIA_DIFF_scaled", "DIA_scaled")  &
                                     COMMON_NAME %in% species_names)|>
        ggplot()+
        geom_line(aes(x = grid_val, y = decadal_pMort, group = draw, color = draw_source), alpha = 0.05)+
        theme_bw(base_size = 12)+
        xlab("Scaled Predictor")+
        ylab("10-year predicted mortality probability")+
        facet_grid(vars(predictor), vars(COMMON_NAME), scales = "free_y")+
        theme(panel.grid = element_blank(),
              strip.text.y = element_text(angle = 0), 
              strip.text.x = element_text(angle = 90), 
              legend.position = "bottom", 
              legend.direction = "horizontal")+
        guides( color = guide_legend(override.aes = list(linewidth = 2, alpha = 1)))+
        labs(color = "Weighted Draw Source")+
          coord_cartesian(ylim = c(0, y_max))
      
      ggsave(paste0(output.dir, "SPCD_stanoutput_cmdstan/summary/Size_effects_pMort10_draws_",group_name,".png"), 
              size.p,
             height = 6, width = (length(species_names) + 1))
      
      
      heighbor.p <- marg_draws_all_df %>% filter(COMMON_NAME %in% species_names)%>% 
        filter(predictor %in% c("ba.scaled", "BAL.scaled", "damage.scaled"))|>
        ggplot()+
        geom_line(aes(x = grid_val, y = decadal_pMort, group = draw, color = draw_source), alpha = 0.05)+
        theme_bw(base_size = 12)+
        xlab("Scaled Predictor")+
        ylab("10-year predicted mortality probability")+
        facet_grid(vars(predictor), vars(COMMON_NAME), scales = "free_y")+
        theme(panel.grid = element_blank(),
              strip.text.y = element_text(angle = 0), 
              strip.text.x = element_text(angle = 90), 
              legend.position = "bottom", 
              legend.direction = "horizontal")+
        guides( color = guide_legend(override.aes = list(linewidth = 2, alpha = 1)))+
        labs(color = "Weighted Draw Source")+
        coord_cartesian(ylim = c(0, y_max))
      
      ggsave(paste0(output.dir, "SPCD_stanoutput_cmdstan/summary/Neighborhood_effects_pMort10_draws_",group_name,".png"), 
             heighbor.p,
             height = 8, width = (length(species_names) + 1))
      
      
      # for the climate effects
      clim.p <- marg_draws_all_df %>% filter(predictor %in% c("MATmax.scaled", "MAP.scaled", 
                                                    "ppt.anom",
                                                    "tmax.anom"))%>%
        filter(COMMON_NAME %in% species_names)|>
        ggplot()+
        geom_line(aes(x = grid_val, y = decadal_pMort, group = draw, color = draw_source), alpha = 0.05)+
        theme_bw(base_size = 12)+
        xlab("Scaled Predictor")+
        ylab("10-year predicted mortality probability")+
        facet_grid(vars(predictor), vars(COMMON_NAME), scales = "free_y")+
        theme(panel.grid = element_blank(),
              strip.text.y = element_text(angle = 0), 
              strip.text.x = element_text(angle = 90), 
              legend.position = "bottom", 
              legend.direction = "horizontal")+
        guides( color = guide_legend(override.aes = list(linewidth = 2, alpha = 1)))+
        labs(color = "Weighted Draw Source")+
        coord_cartesian(ylim = c(0, y_max))
      
      ggsave(paste0(output.dir, "SPCD_stanoutput_cmdstan/summary/Climate_effects_pMort10_draws_",group_name,".png"), 
             clim.p,
             height = 9, width = (length(species_names) + 1))
      
      # for the site condition effects
      marg_draws_all_df %>% filter(predictor %in% c("Ndep.scaled","slope.scaled", "aspect.scaled"))%>%
        filter(COMMON_NAME %in% species_names)|>
        ggplot()+
        geom_line(aes(x = grid_val, y = decadal_pMort, group = draw, color = draw_source), alpha = 0.05)+
        theme_bw(base_size = 12)+
        xlab("Scaled Predictor")+
        ylab("10-year predicted mortality probability")+
        facet_grid(vars(predictor), vars(COMMON_NAME), scales = "free_y")+
        theme(panel.grid = element_blank(),
              strip.text.y = element_text(angle = 0), 
              strip.text.x = element_text(angle = 90), 
              legend.position = "bottom", 
              legend.direction = "horizontal")+
        guides( color = guide_legend(override.aes = list(linewidth = 2, alpha = 1)))+
        labs(color = "Weighted Draw Source")+
        coord_cartesian(ylim = c(0, y_max))
      
      ggsave(paste0(output.dir, "SPCD_stanoutput_cmdstan/summary/Site_condition_effects_pMort10_draws_",group_name,".png"), 
             height = 8,  width = (length(species_names) + 1))
      

}

plot_draws_predictors(species_names = n_conifers, group_name = "n_conifers")
plot_draws_predictors(species_names = Oak_birch_hickory, group_name = "Hardwoods", y_max = 0.3)
plot_draws_predictors(species_names = Maples_ash_bcherry, group_name = "Maple_Pine_Hemlock", y_max = 0.1)

########################################################################################
# Marginal effects plots with different CI around them using geom_ribbon------

marg_summary_all_df <- do.call(rbind, 
                             lapply(list.files(path = paste0(output.dir, "SPCD_stanoutput_cmdstan/model_avg_draws/stacking_wts/"), 
                                                pattern ="Summary_marginal_main", full.names = T), 
                                    qs2::qs_read))

main_effects_df <- data.frame(predictor = c("DIA_DIFF_scaled", "DIA_scaled", "ba.scaled", "BAL.scaled",
                  "damage.scaled", "MATmax.scaled", "MAP.scaled", "ppt.anom",
                  "tmax.anom", "slope.scaled", "aspect.scaled", "Ndep.scaled"),
                  
                  Pretty_name = c("Diameter difference", "Diameter", "Plot Basal Area", "Tree BA Larger than", 
                                  "% Damage", "MATmax", "MAP", "Precip anomaly", "Tmax anomaly", "slope", "aspect", "N deposition"))


marg_summary_all_df <- marg_summary_all_df %>% left_join(., main_effects_df)
marg_summary_all_df$predictor <- factor(marg_summary_all_df$predictor, levels = main_effects_df$predictor)
marg_summary_all_df$Pretty_name <- factor(marg_summary_all_df$Pretty_name, levels = main_effects_df$Pretty_name)


# convert grid_value back to raw values for each species:
train.data.list <- test.data.list <- list()

all.spcd.data.files <- list.files(paste0("SPCD_standata_general_full_standardized_v3/"), pattern = "_9.Rdata", full.names = TRUE)
for(i in 1:length(all.spcd.data.files)){
  load(all.spcd.data.files[i])
  train.data.list[[i]] <- train.data
  test.data.list[[i]] <- test.data
}


# covariates that are scaled by species:
train.data <- do.call(rbind, train.data.list)

test.data <-do.call(rbind, test.data.list)
plot.medians <- readRDS("data/plot.medians_SPCD_all.rds")


# get the species scaling information:

# n deposition
species.scaling <- 
 
  #DIA_DIFF_scaled
 train.data %>% select(Species, SPCD, logDIA.DIFF.median, logDIA.DIFF.sd) %>% distinct()%>%
          mutate(Covariate = "DIA_DIFF_scaled", 
                 Clean_Name = "Diameter Difference",
                 Transformation = "log(x)",
                 Units = "Inches/year")%>% 
          rename("Val.mean" = "logDIA.DIFF.median",
                 "Val.sd" = "logDIA.DIFF.sd")%>%
  #DIA
  rbind(., train.data %>% select(Species, SPCD, logDIA.median, logDIA.sd) %>% distinct()%>%
          mutate(Covariate = "DIA_scaled", 
                 Clean_Name = "Diameter",
                 Transformation = "log(x)",
                 Units = "Inches")%>% 
          rename("Val.mean" = "logDIA.median",
                 "Val.sd" = "logDIA.sd") 
  )%>%
  
  #plt_ba_sq_ft_old
  rbind(., train.data %>% select(Species, SPCD, logplt_ba_sq_ft_cur.median, logplt_ba_sq_ft_cur.sd) %>% 
          distinct()%>%
          mutate(Covariate = "ba.scaled", 
                 Clean_Name = "Plot Basal Area",
                 Transformation = "log(1+x)",
                 Units = "sq.ft. per acre")%>% 
          rename("Val.mean" = "logplt_ba_sq_ft_cur.median",
                 "Val.sd" = "logplt_ba_sq_ft_cur.sd")
  )%>%
 
  #BAL ratio
  rbind(.,train.data %>% select(Species, SPCD, BAL.ratio.median, BAL.ratio.sd) %>% distinct()%>%
          mutate(Covariate = "BAL.scaled",
                 Clean_Name = "Basal Area Larger Than",
                 Transformation = "None",
                 Units = "Ratio")%>% 
          rename("Val.mean" = "BAL.ratio.median",
                 "Val.sd" = "BAL.ratio.sd")
  )%>%
  
  #damage.scaled
  rbind(.,train.data %>% select(Species, SPCD) %>% distinct()%>%
          mutate(damage.median = plot.medians$damage.median, 
                 damage.sd = plot.medians$damage.sd) %>%
          mutate(Covariate = "damage.scaled",
                 Clean_Name = "% Damage",
                 Transformation = "log(1+x)",
                 Units = "Percent")%>% 
          rename("Val.mean" = "damage.median",
                 "Val.sd" = "damage.sd")
  )%>%
  
  #MATmax.scaled
  rbind(.,train.data %>% select(Species, SPCD) %>% distinct()%>%
          mutate(MATmax.median = plot.medians$MATmax.median, 
                 MATmax.sd = plot.medians$MATmax.sd) %>%
          mutate(Covariate = "MATmax.scaled",
                 Clean_Name = "Mean Annual Tmax",
                 Transformation = "None",
                 Units = "Degrees")%>% 
          rename("Val.mean" = "MATmax.median",
                 "Val.sd" = "MATmax.sd")
  )%>%
  
  #MAP.scaled
  rbind(.,train.data %>% select(Species, SPCD) %>% distinct()%>%
          mutate(MAP.median = plot.medians$MAP.median, 
                 MAP.sd = plot.medians$MAP.sd) %>%
          mutate(Covariate = "MAP.scaled",
                 Clean_Name = "Mean Annual Precipitation",
                 Transformation = "None",
                 Units = "mm")%>% 
          rename("Val.mean" = "MAP.median",
                 "Val.sd" = "MAP.sd")
  )%>%
  
  
  #ppt.anom
  rbind(.,train.data %>% select(Species, SPCD, ppt.anom.median, ppt.anom.sd) %>% distinct()%>%
          mutate(Covariate = "ppt.anom",
                 Clean_Name = "Precipitation anomaly",
                 Transformation = "None",
                 Units = " ")%>% 
          rename("Val.mean" = "ppt.anom.median",
                 "Val.sd" = "ppt.anom.sd")
  )%>%
  
  #tmax.anom
  rbind(.,train.data %>% select(Species, SPCD, tmax.anom.median, tmax.anom.sd) %>% distinct()%>%
          mutate(Covariate = "tmax.anom",
                 Clean_Name = "Max. Temperature anomaly",
                 Transformation = "None",
                 Units = " ")%>% 
          rename("Val.mean" = "tmax.anom.median",
                 "Val.sd" = "tmax.anom.sd")
  )%>%
  
  
  #slope
  rbind(.,train.data %>% select(Species, SPCD, slope.sin.median, slope.sin.sd) %>% distinct()%>%
          mutate(Covariate = "slope.scaled",
                 Clean_Name = "slope",
                 Transformation = "arcsin(x)",
                 Units = "%")%>% 
          rename("Val.mean" = "slope.sin.median",
                 "Val.sd" = "slope.sin.sd")
  )%>%
  
  #aspect
  rbind(.,train.data %>% select(Species, SPCD, aspect.cos.median, aspect.cos.sd) %>% distinct()%>%
          mutate(Covariate = "aspect.scaled",
                 Clean_Name = "cos(aspect)",
                 Transformation = "None",
                 Units = "Northing")%>% 
          rename("Val.mean" = "aspect.cos.median",
                 "Val.sd" = "aspect.cos.sd")
  )%>%
  
  # N deposition
  rbind(., train.data %>% select(Species, SPCD) %>% distinct()%>%
          mutate(Ndep.median = plot.medians$Ndep.median, 
                 Ndep.sd = plot.medians$Ndep.sd) %>%
  mutate(Covariate = "Ndep.scaled", 
         Clean_Name = "N deposition rate",
         Transformation = "None",
         Units = "kg/m^2/year") %>% 
  rename("Val.mean" = "Ndep.median",
         "Val.sd" = "Ndep.sd")
  )%>%
  select(Species, SPCD, Covariate, Val.mean, Val.sd, Clean_Name, Units, Transformation) %>% 
  rename("predictor" = "Covariate")
  
# convert scaled grid values to the raw values for each species:
marg_with_raw_values <- marg_summary_all_df %>% left_join(., species.scaling)%>% 
  mutate(Descaled = grid_val*Val.sd + Val.mean)%>%
  mutate(Raw_grid_val = ifelse(Transformation %in% "log(x)", exp(Descaled), 
                               ifelse(Transformation %in% "log(1+x)", expm1(Descaled), 
                                      ifelse(Transformation %in% "arcsin(x)", asin(Descaled)/(pi/180),Descaled))))

# We just used a blanket -2 to 2 grid for calculating marginal effects, so 
# find actual ranges for each species and filter:

unique(marg_with_raw_values$predictor)


train.test.data <- rbind(train.data %>% select(Species, SPCD, annual.growth, dbhold, plt_ba_sq_ft_cur, BAL.ratio, damage, 
                                                MATmax, MAP, ppt.anom, tmax.anom, slope, slope.sin, aspect, aspect.cos, Ndep.remper.avg),
                         test.data %>% select(Species, SPCD, annual.growth, dbhold, plt_ba_sq_ft_cur, BAL.ratio, damage, 
                                                MATmax, MAP, ppt.anom, tmax.anom, slope, slope.sin, aspect, aspect.cos, Ndep.remper.avg)) %>%
  mutate(slope_degrees = slope)

raw_value_ranges <- train.test.data %>% pivot_longer(
  cols = c(- Species, -SPCD), 
  names_to = c("raw_name"),
  values_to = "raw_value"
) %>% 
  group_by(Species, SPCD, raw_name) %>% 
  summarise(median_val = median(raw_value), 
            min_val = min(raw_value),
            max_val = max(raw_value), 
            q2.5_val = quantile(raw_value, 0.025),
            q97.5_val = quantile(raw_value, 0.975))

unique(raw_value_ranges$raw_name)
unique(marg_with_raw_values$Clean_Name)
unique(marg_with_raw_values$predictor)
# create a table mapping the raw names to the scaled predictor names

raw2pred_name <- tibble::tribble(
  ~predictor,          ~raw_name,
  "DIA_DIFF_scaled",   "annual.growth",
  "DIA_scaled",        "dbhold", 
  "ba.scaled",         "plt_ba_sq_ft_cur", 
  "BAL.scaled",        "BAL.ratio", 
  "damage.scaled",     "damage",
  "MATmax.scaled",     "MATmax" ,
  "MAP.scaled",        "MAP", 
  "ppt.anom",          "ppt.anom" ,
  "tmax.anom",         "tmax.anom", 
  "slope.scaled",      "slope", 
  "aspect.scaled",     "aspect.cos", 
  "Ndep.scaled",       "Ndep.remper.avg" 
)

imperial2metric <- tibble::tribble(
  ~Units,          ~Units_metric_name,      ~multiplier,
  "Inches",            "cm",            2.54,
  "Inches/year",       "cm/year",       2.54,
  "sq.ft per acre",    "m^2/hectare",   (1/4.356), 
  
)

raw_val_rng_df <- raw_value_ranges %>% left_join(., raw2pred_name) %>% filter(!is.na(predictor))
marg_spec_rng_raw_df <- marg_with_raw_values %>% left_join(., raw_val_rng_df)%>%
  mutate(val_in_spp_rng = ifelse(Raw_grid_val >= min_val & Raw_grid_val <= max_val, "in_spp_rng", "out_spp_range")) %>%
  mutate(Units_metric = ifelse(Units %in% "Inches", "cm", 
                               ifelse(Units %in% "Inches/year", "cm/year", 
                                      ifelse(Units %in% "sq.ft. per acre", "m^2/hectare", Units)))) %>% 
  mutate(raw_grid_val_metric = ifelse(Units %in% c("Inches", "Inches/year"), Raw_grid_val*2.54, 
                                    ifelse(Units %in% "sq.ft. per acre", Raw_grid_val*(1/4.356), Raw_grid_val))) 

  
  marg_spec_rng_raw_df %>% filter(is.na(median_val))

marg_spec_rng_raw_df %>% View()



###################################################################################
# function to plot all the 50% CI by predictor name

plot_50CI_summary_clean <- function(predictor_name, 
                                    species_names = unique(marg_spec_rng_raw_df$Species),  
                                    y_max = 0.25, 
                                    facet_species = TRUE){
  
  SPP.marginals_predictor <- marg_spec_rng_raw_df %>% filter(predictor %in% predictor_name  &
                                    COMMON_NAME %in% species_names) %>% 
                                    filter(val_in_spp_rng %in% "in_spp_rng")
  
  #unique(SPP.marginals_predictor$Clean_Name)
  unit_parentheses <- ifelse(unique(SPP.marginals_predictor$Units_metric) %in% " ", " ", paste0(" (",unique(SPP.marginals_predictor$Units_metric), ")" ))
  if(facet_species == TRUE){
    SPP.marginals_predictor |>
        ggplot()+
        #geom_ribbon(aes(x = Raw_grid_val, ymin = decadal_pMort.05.ci.lo, ymax = decadal_pMort.95.ci.hi, fill = COMMON_NAME), alpha = 0.25)+
        geom_ribbon(aes(x = raw_grid_val_metric, ymin = decadal_pMort.25.ci.lo, ymax = decadal_pMort.75.ci.hi, fill = COMMON_NAME), alpha = 0.65)+
        geom_line(aes(x = raw_grid_val_metric, y = decadal_pMort.median, color = COMMON_NAME))+
        species_fill + species_color+
        theme_bw()+
        facet_grid(vars(Pretty_name),  vars(COMMON_NAME), scales = "free_x")+
        theme(panel.grid = element_blank(),
              strip.text.y = element_text(angle = 0), 
              strip.text.x = element_text(angle = 90), 
              legend.position = "bottom", 
              legend.direction = "horizontal")+
        coord_cartesian(ylim = c(0, y_max))+
        xlab(paste0(unique(SPP.marginals_predictor$Clean_Name), unit_parentheses))+
        ylab("10-year predicted mortality probability")
    
  }else{
    
    SPP.marginals_predictor |>
      ggplot()+
      #geom_ribbon(aes(x = Raw_grid_val, ymin = decadal_pMort.05.ci.lo, ymax = decadal_pMort.95.ci.hi, fill = COMMON_NAME), alpha = 0.25)+
      geom_ribbon(aes(x = raw_grid_val_metric, ymin = decadal_pMort.25.ci.lo, ymax = decadal_pMort.75.ci.hi, fill = COMMON_NAME), alpha = 0.65)+
      geom_line(aes(x = raw_grid_val_metric, y = decadal_pMort.median, color = COMMON_NAME))+
      species_fill + species_color+
      theme_bw()+
      facet_wrap(~Pretty_name, scales = "free_x")+
      theme(panel.grid = element_blank(),
            #strip.text.y = element_text(angle = 0), 
            #strip.text.x = element_text(angle = 90), 
            legend.position = "bottom", 
            legend.direction = "horizontal")+
      coord_cartesian(ylim = c(0, y_max))+
      xlab(paste0(unique(SPP.marginals_predictor$Clean_Name), unit_parentheses))+
      ylab("10-year p(mortality)")
  }
 
}


# default plot is faceted by species
plot_50CI_summary_clean(predictor_name = c("DIA_DIFF_scaled"), 
                        species_names = unique(marg_with_raw_values$Species),  
                        y_max = 0.5, facet_species = FALSE)

# can also remove the facet by species
plot_50CI_summary_clean(predictor_name = "DIA_DIFF_scaled", 
                        species_names = unique(marg_with_raw_values$Species),  
                        y_max = 0.15, 
                        facet_species = FALSE)

plot_50CI_summary_clean(predictor_name = "DIA_DIFF_scaled", 
                        species_names = unique(marg_with_raw_values$Species),  
                        y_max = 0.25)

plot_50CI_summary_clean(predictor_name = "DIA_scaled", 
                        species_names = unique(marg_with_raw_values$Species),  
                        y_max = 0.25, 
                        facet_species = FALSE)

plot_50CI_summary_clean(predictor_name = "ba.scaled", 
                        species_names = unique(marg_with_raw_values$Species),  
                        y_max = 0.25)

plot_50CI_summary_clean(predictor_name = "ba.scaled", 
                        species_names = unique(marg_with_raw_values$Species),  
                        y_max = 0.10, 
                        facet_species = FALSE)

plot_50CI_summary_clean(predictor_name = "BAL.scaled", 
                        species_names = unique(marg_with_raw_values$Species),  
                        y_max = 0.15)

plot_50CI_summary_clean(predictor_name = "damage.scaled", 
                        species_names = unique(marg_with_raw_values$Species),  
                        y_max = 0.25)

plot_50CI_summary_clean(predictor_name = "Ndep.scaled", 
                        species_names = unique(marg_with_raw_values$Species),  
                        y_max = 0.25)

plot_50CI_summary_clean(predictor_name = "MAP.scaled", 
                        species_names = unique(marg_with_raw_values$Species),  
                        y_max = 0.25)

plot_50CI_summary_clean(predictor_name = "MATmax.scaled", 
                        species_names = unique(marg_with_raw_values$Species),  
                        y_max = 0.25)


plot_50CI_summary_clean(predictor_name = "ppt.anom", 
                        species_names = unique(marg_with_raw_values$Species),  
                        y_max = 0.25)

plot_50CI_summary_clean(predictor_name = "tmax.anom", 
                        species_names = unique(marg_with_raw_values$Species),  
                        y_max = 0.25)


plot_50CI_summary_clean(predictor_name = "slope.scaled", 
                        species_names = unique(marg_with_raw_values$Species),  
                        y_max = 0.25)

plot_50CI_summary_clean(predictor_name = "aspect.scaled", 
                        species_names = unique(marg_with_raw_values$Species),  
                        y_max = 0.25)







# TODO: 
# classify general shapes of the curves and plot up together?

marg_with_raw_values %>% group_by(SPCD, COMMON_NAME, predictor, Clean_Name) %>%
  summarise(min_predictor = min(Raw_grid_val, na.rm =TRUE), 
            max_predictor = max(Raw_grid_val, na.rm =TRUE)) %>% 
  mutate(linear_diff_pred = (max_predictor - min_predictor)/min_predictor)



plot_ribbons_predictors <- function( species_names, group_name,  y_max = 1){
  
  
 dia.diff.p <-  plot_50CI_summary_clean(predictor_name = "DIA_DIFF_scaled", 
                          species_names = species_names,  
                          y_max = y_max, facet_species = FALSE)
 dia.p <-  plot_50CI_summary_clean(predictor_name = "DIA_scaled", 
                                        species_names = species_names,  
                                        y_max = y_max, facet_species = FALSE)
 
  
  size.p <- cowplot::plot_grid(dia.diff.p, dia.p, ncol = 1, align = "hv")
  
  ggsave(paste0(output.dir, "SPCD_stanoutput_cmdstan/summary/Size_effects_pMort10_ribbons_",group_name,".png"), 
         size.p,
         height = 6, width = 4)
  
  
  ba.p <-  plot_50CI_summary_clean(predictor_name = "ba.scaled", 
                                         species_names = species_names,  
                                         y_max = y_max, facet_species = FALSE)
  BAL.p <-  plot_50CI_summary_clean(predictor_name = "BAL.scaled", 
                                    species_names = species_names,  
                                    y_max = y_max, facet_species = FALSE)
  
  damage.p <-  plot_50CI_summary_clean(predictor_name = "damage.scaled", 
                                    species_names = species_names,  
                                    y_max = y_max, facet_species = FALSE)
  
  heighbor.p <-  cowplot::plot_grid(ba.p, BAL.p, damage.p, ncol = 1, align = "hv")
  ggsave(paste0(output.dir, "SPCD_stanoutput_cmdstan/summary/Neighborhood_effects_pMort10_ribbons_",group_name,".png"), 
         heighbor.p,
         height = 8, width = 4)
  
  
  # for the climate effects
  MAT.p <-  plot_50CI_summary_clean(predictor_name = "MATmax.scaled", 
                                   species_names = species_names,  
                                   y_max = y_max, facet_species = FALSE)
  MAP.p <-  plot_50CI_summary_clean(predictor_name = "MAP.scaled", 
                                    species_names = species_names,  
                                    y_max = y_max, facet_species = FALSE)
  
  ppt_anom.p <-  plot_50CI_summary_clean(predictor_name = "ppt.anom", 
                                       species_names = species_names,  
                                       y_max = y_max, facet_species = FALSE)
  
  tmax_anom.p <-  plot_50CI_summary_clean(predictor_name = "tmax.anom", 
                                       species_names = species_names,  
                                       y_max = y_max, facet_species = FALSE)
  
  clim.p <- cowplot::plot_grid(MAT.p, MAP.p, tmax_anom.p,ppt_anom.p,  ncol = 1, align = "hv")
  
  ggsave(paste0(output.dir, "SPCD_stanoutput_cmdstan/summary/Climate_effects_pMort10_ribbons_",group_name,".png"), 
         clim.p,
         height = 9, width = 4)
  
  # for the site condition effects
  
  
  Ndep.p <-  plot_50CI_summary_clean(predictor_name = "Ndep.scaled", 
                                   species_names = species_names,  
                                   y_max = y_max, facet_species = FALSE)
  slope.p <-  plot_50CI_summary_clean(predictor_name = "slope.scaled", 
                                    species_names = species_names,  
                                    y_max = y_max, facet_species = FALSE)
  
  aspect.p <-  plot_50CI_summary_clean(predictor_name = "aspect.scaled", 
                                       species_names = species_names,  
                                       y_max = y_max, facet_species = FALSE)
  
  
  
  
 site.p <-  cowplot::plot_grid(Ndep.p, slope.p, aspect.p,  ncol = 1, align = "hv")
          
  ggsave(paste0(output.dir, "SPCD_stanoutput_cmdstan/summary/Site_condition_effects_pMort10_ribbons_",group_name,".png"), 
         site.p,
         height = 8,  width = 4)
  
  
}

plot_ribbons_predictors(species_names = n_conifers, group_name = "n_conifers", y_max = 0.5)
plot_ribbons_predictors(species_names = Oak_birch_hickory, group_name = "Hardwoods", y_max = 0.3)
plot_ribbons_predictors(species_names = Maples_ash_bcherry, group_name = "Maple_Pine_Hemlock", y_max = 0.1)


# plot up the 50% CI only across the values:
ggplot(data = marg_summary_all_df)+
  #geom_ribbon(aes(x = grid_val, ymin = decadal_pMort.05.ci.lo, ymax = decadal_pMort.95.ci.hi, fill = COMMON_NAME), alpha = 0.25)+
  geom_ribbon(aes(x = grid_val, ymin = decadal_pMort.25.ci.lo, ymax = decadal_pMort.75.ci.hi, fill = COMMON_NAME), alpha = 0.5)+
  geom_line(aes(x = grid_val, y = decadal_pMort.median, color = COMMON_NAME))+
  facet_wrap(~Pretty_name, scales = "free_y")+
  species_fill + species_color+
  theme_bw()

# omit balsam fir
ggplot(data = marg_summary_all_df %>% filter(! COMMON_NAME %in% c("balsam fir", "northern white-cedar", "red spruce")))+
  geom_ribbon(aes(x = grid_val, ymin = decadal_pMort.05.ci.lo, ymax = decadal_pMort.95.ci.hi, fill = COMMON_NAME), alpha = 0.25)+
  geom_ribbon(aes(x = grid_val, ymin = decadal_pMort.25.ci.lo, ymax = decadal_pMort.75.ci.hi, fill = COMMON_NAME), alpha = 0.75)+
  geom_line(aes(x = grid_val, y = decadal_pMort.median, color = COMMON_NAME))+
  facet_wrap(~Pretty_name, scales = "free_y")+
  species_fill + species_color+
  theme_bw()

######################################################################################
# Plot predicted county mortality case studies ----
# Beech bark disease 
# Hemlock wooley adelgid
# Spongy moth
# Spruce Budworm



########################################################################################
# Old code------


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


joint.betas <- readRDS(paste0(output.folder, "SPCD_stanoutput_joint_v3/samples/u_betas_model_6_5000samples.RDS"))
joint.betas.summary <- joint.betas %>% select(-.iteration, -.chain, -.draw)%>% reshape2::melt() %>%
  group_by(variable)%>% summarise(median = quantile(value, 0.5), 
                                  ci.lo = quantile(value, 0.025), 
                                  ci.hi = quantile(value, 0.975), 
                                  log.odds.median = quantile(exp(value), 0.5), 
                                  log.odds.ci.lo = quantile(exp(value), 0.025), 
                                  log.odds.ci.hi = quantile(exp(value), 0.975))


joint.alphas <- readRDS(paste0(output.folder, "SPCD_stanoutput_joint_v3/samples/alpha.spp_model_6_5000samples.RDS"))
joint.alphas.summary <- joint.alphas %>% select(-.iteration, -.chain, -.draw)%>% reshape2::melt() %>%
  group_by(variable)%>% summarise(median = quantile(value, 0.5), 
                                  ci.lo = quantile(value, 0.025), 
                                  ci.hi = quantile(value, 0.975), 
                                  log.odds.median = quantile(exp(value), 0.5), 
                                  log.odds.ci.lo = quantile(exp(value), 0.025), 
                                  log.odds.ci.hi = quantile(exp(value), 0.975)) 

# read in the data used to fit
mod.data.full <- readRDS ( paste0(output.folder,"SPCD_stanoutput_joint_v3/all_SPCD_model_6.RDS"))
param.id = rep(1:33, each = 17)
spp.id = rep(1:17, 33)
beta.names.df <- data.frame(variable = unique(joint.betas.summary$variable), 
                            spp = spp.id, 
                            Param.id = param.id, 
                            Parameter.name = rep(colnames(mod.data.full$xM), each = 17))

alpha.names.df <- data.frame(variable = unique(joint.alphas.summary$variable), 
                             spp = 1:17, 
                             Param.id = 34, 
                             Parameter.name = rep("alpha", 17))

joint.betas.summary <- joint.betas.summary %>% left_join(.,beta.names.df) %>% left_join(., spp.table)
joint.betas.summary$significance <- ifelse(joint.betas.summary$ci.lo > 0 & joint.betas.summary$ci.hi > 0, "significant", 
                                           ifelse(joint.betas.summary$ci.lo < 0 & joint.betas.summary$ci.hi < 0, "significant", "not significant"))
joint.betas.summary$Parameter.name <- factor(joint.betas.summary$Parameter.name, levels = unique(beta.names.df$Parameter.name))


joint.alphas.summary <- joint.alphas.summary %>% left_join(.,alpha.names.df) %>% left_join(., spp.table)
joint.alphas.summary$significance <- ifelse(joint.alphas.summary$ci.lo > 0 & joint.alphas.summary$ci.hi > 0, "significant", 
                                            ifelse(joint.alphas.summary$ci.lo < 0 & joint.alphas.summary$ci.hi < 0, "significant", "not significant"))
joint.alphas.summary$Parameter.name <- factor(joint.alphas.summary$Parameter.name, levels = unique(alpha.names.df$Parameter.name))

joint.betas.alpha.summary <- rbind(joint.betas.summary, joint.alphas.summary)

model6.betas <- ggplot()+geom_point(data = joint.betas.alpha.summary, aes(x = COMMON, y = median, color = significance))+
  geom_errorbar(data = joint.betas.alpha.summary, aes(x = COMMON, ymin = ci.lo, ymax = ci.hi, color = significance))+
  facet_wrap(~Parameter.name, scales = "free_y")+
  scale_color_manual(values = c("significant" = "black", 
                                "not significant" = "darkgrey"))+
  theme_bw()+
  theme(axis.text.x = element_text(vjust = 0.5, angle = 90, hjust = 1), 
        panel.grid = element_blank())+
  ylab("Effect on Survival")+
  xlab("Species")
ggsave(filename = paste0(output.folder,"images/model_6_regression_betas.png"), 
       model6.betas,
       height = 10, width =12, units = "in")


ggplot()+geom_point(data = joint.betas.summary  %>% filter(significance %in% "significant"), aes(x = Parameter.name, y = log.odds.median, color = log.odds.median >=1))+
  geom_errorbar(data = joint.betas.summary  %>% filter(significance %in% "significant"), aes(x = Parameter.name, ymin = log.odds.ci.lo, ymax = log.odds.ci.hi, color = log.odds.median >=1))+
  geom_hline(yintercept = 1, color = "red", linetype = "dashed")+facet_wrap(~COMMON, scales = "free")+theme(axis.text.x = element_text(hjust = 1, angle = 45), legend.position = "none")
ggsave(filename = paste0(output.folder,"images/joint_model_log_odds_betas_significant.png"), 
       height = 10, width =12, units = "in")


rhat_ESS_model6_hiearachical <- rbind(
  # betas
  summarise_draws(joint.betas) %>% select(variable,  rhat, ess_bulk, ess_tail) %>% 
    left_join(., beta.names.df) %>% left_join(., spp.table),
  # alphas
  summarise_draws(joint.alphas) %>% select(variable,  rhat, ess_bulk, ess_tail) %>% 
    left_join(., alpha.names.df) %>% left_join(., spp.table))

bulk.ess.h <- rhat_ESS_model6_hiearachical |>
  ggplot()+geom_histogram(aes(y = ess_bulk))+
  geom_hline(aes(yintercept = 1000), linetype = "dashed")+
  ylab("Bulk ESS")+theme_bw()

tail.ess.h <- rhat_ESS_model6_hiearachical |>
  ggplot()+geom_histogram(aes(y = ess_tail))+
  geom_hline(aes(yintercept = 1000), linetype = "dashed")+
  ylab("Tail ESS")+theme_bw()

rhat.ess.h <- rhat_ESS_model6_hiearachical |>
  ggplot()+geom_histogram(aes(y = rhat))+
  geom_hline(aes(yintercept = 1.01), linetype = "dashed")+
  ylab("R-hat")+theme_bw()

figureS9<- cowplot::plot_grid(rhat.ess.h,bulk.ess.h , tail.ess.h, 
                              align = "hv", ncol = 3,
                              labels = "AUTO")

ggsave(filename = paste0(output.folder,"images/ESS_Rhat_statistics_model6_heirarchical.png"), 
       figureS9,
       height = 4, width =11, units = "in")

# total maple syrup production value in US
300396*1000/10e6

# in region:
((899 + 27690 + 2808+7972 +28933 + 3965 + 6989 + 95416 +570))*1000
#/10e6


species.name <- "yellow birch"
hotspot.variable <- "n_county_pmort"

get.spatial.hotspots <- function(species.name, hotspot.variable ){
        # do a hotspot analysis for eastern hemlock:
  
  
  species.data.all <- county.pmort.sf %>% filter(Species %in% species.name )
  
  species.data.all$hotspot.var <- species.data.all[[hotspot.variable]]
    
  species.data.all <- species.data.all %>% mutate(
    cut.mort.rate = cut(hotspot.var*100, breaks = c(0, 0.25, 0.5, 0.75, 1, 1.5, 2, 4, 100))
  )
  cut.mort.values <- data.frame(
    cut.mort.rate = c("(0,0.25]", "(0.25,0.5]", "(0.5,0.75]", "(0.75,1]", "(1,1.5]", "(1.5,2]", "(2,4]", "(4,100]"), 
    mortality.rate = c("< 0.25", "0.25 - 0.5", "0.5 - 0.75", "0.75 - 1", "1 - 1.5", "1.5 - 2", "2 - 4", "> 4"), 
    red.hex.colors = c(
      "#ffffcc",
      "#ffeda0",
      "#fed976",
      "#feb24c",
      "#fd8d3c",
      "#fc4e2a",
      "#e31a1c",
      "#b10026"),
    hex.colors = c(
                  "#fff7f3",
                   "#fde0dd",
                   "#fcc5c0",
                   "#fa9fb5",
                   "#f768a1",
                   "#dd3497",
                   "#ae017e",
                   "#7a0177"))
  
  hex.vector <- as.vector(cut.mort.values$red.hex.colors)
  names( hex.vector) <- cut.mort.values$mortality.rate
  fill_mort_rate_red <- scale_fill_manual(values = hex.vector, name = "Mortality\nRate\n(%/year)", drop = FALSE)
  
  
  hex.vector <- as.vector(cut.mort.values$hex.colors)
  names( hex.vector) <- cut.mort.values$mortality.rate
  fill_mort_rate_pink <- scale_fill_manual(values = hex.vector, name = "Mortality\nRate\n(%/year)", drop = FALSE)
  
  species.data.all <-  species.data.all %>% left_join(., cut.mort.values)
  species.data.all$mortality.rate <- factor( species.data.all$mortality.rate, levels = c("< 0.25", "0.25 - 0.5", "0.5 - 0.75", "0.75 - 1", "1 - 1.5", "1.5 - 2", "2 - 4", "> 4"))
  
  # how widespread is mortality probability across the range?---
  
  # if high mortality is >1% per year, sum up the # of counties with greater than that
  pct.counties.high <- species.data.all%>% as.data.frame()%>% filter(!is.na(hotspot.var))%>%
    mutate(over.threshold = ifelse(hotspot.var*100 >= 0.75, 1, 0))%>%
    summarise(ncounties = n(), 
              n.over.threshold = sum(over.threshold, na.rm =TRUE))%>%
    mutate(pct.counties.high = round(sum(n.over.threshold/ncounties)*100, 1))
    
  
  
  
  
  # Visualize the variable on the map
  #species.data.all |>
  spp.map.pmort <- ggplot() +
    geom_sf(data = species.data.all, aes(fill = mortality.rate),color = "black", lwd = 0.1, show.legend = TRUE) +
    
    fill_mort_rate_pink+
    theme_void() +
    labs(
      fill = paste("mortality probability\n(%/year)"),
      title = paste(species.name, "-", hotspot.variable), 
      subtitle = paste0(pct.counties.high$pct.counties.high, "% of counties > 0.75%/year" ))
  

      county.pmort.sf <- county.pmort.sf[!st_is_empty(county.pmort.sf), ]
  
        hemlock.counties <- county.pmort.sf[!is.na(county.pmort.sf[,hotspot.variable]),] %>% 
          filter(Species %in% species.name ) %>% filter(!is.na(geometry))
        
        
        # Create a neighbor list based on queen contiguity
        list_nb <- poly2nb(hemlock.counties, queen = TRUE)
        
        empty_nb <- which(card(list_nb) == 0)
        empty_nb       
        
        if(length(empty_nb)>0){
        hemlock_subset <- hemlock.counties[-empty_nb, ]
        }else{
          hemlock_subset <- hemlock.counties
        }
        # Now that we removed empty neighbor sets (tes_subset)
        # Identify neighbors with queen contiguity (edge/vertex touching)
        hemlock_nb <- poly2nb(hemlock_subset, queen = TRUE)
        
        # Binary weighting assigns a weight of 1 to all neighboring features 
        # and a weight of 0 to all other features
        hemlock_w_binary <- nb2listw(hemlock_nb, style="B")
        
        # Calculate spatial lag of county probability of mortality
        # rename the variable of interest 
        hemlock_subset$hotspot.var <- hemlock_subset[[hotspot.variable]]
        
        hemlock_lag <- lag.listw(hemlock_w_binary, hemlock_subset$hotspot.var)
        
        #  global G statistic of county probability of mortality
       global.G.results <-  globalG.test( hemlock_subset$hotspot.var, hemlock_w_binary)
        
        # local gi test:
        # Identify neighbors, create weights, calculate spatial lag
        hemlock_nbs <- hemlock_subset |> 
          mutate(
            nb = st_contiguity(geometry),        # neighbors share border/vertex
            wt = st_weights(nb),                 # row-standardized weights
            tes_lag = st_lag(hotspot.var, nb, wt)    # calculate spatial lag of TreEqty
          ) 
        
        # Calculate the Gi using local_g_perm
        hemlock_hot_spots <- hemlock_nbs |> 
          mutate(
            Gi = local_g_perm(hotspot.var, nb, wt, nsim = 999)
            # nsim = number of Monte Carlo simulations (999 is default)
          ) |> 
          # The new 'Gi' column itself contains a dataframe 
          # We can't work with that, so we need to 'unnest' it
          unnest(Gi) 
        
        # plot Gi star
        hemlock_hot_spots |> 
          ggplot((aes(fill = gi))) +
          geom_sf(color = "black", lwd = 0.15) +
          scale_fill_gradient2() 
        
        # Create a new data frame with pretty hotspots
     local.hotspot.sf <-   hemlock_hot_spots |> 
          # with the columns 'gi' and 'p_folded_sim"
          # 'p_folded_sim' is the p-value of a folded permutation test
          dplyr::select(gi, p_folded_sim) |> 
          mutate(
            # Add a new column called "classification"
            classification = case_when(
              # Classify based on the following criteria:
              gi > 0 & p_folded_sim <= 0.01 ~ "Very hot",
              gi > 0 & p_folded_sim <= 0.05 ~ "Hot",
              gi > 0 & p_folded_sim <= 0.1 ~ "Somewhat hot",
              gi < 0 & p_folded_sim <= 0.01 ~ "Very cold",
              gi < 0 & p_folded_sim <= 0.05 ~ "Cold",
              gi < 0 & p_folded_sim <= 0.1 ~ "Somewhat cold",
              TRUE ~ "Insignificant"
            ),
            # Convert 'classification' into a factor for easier plotting
            classification = factor(
              classification,
              levels = c("Very hot", "Hot", "Somewhat hot",
                         "Insignificant",
                         "Somewhat cold", "Cold", "Very cold")
            )
          ) %>%
       mutate(globalG_st_dev = global.G.results$statistic[1,1], 
              globalG_pval = global.G.results$p.value[1,1], 
              globalG_estimate = as.numeric(global.G.results$estimate[1]), 
              globalG_expectation = as.numeric(global.G.results$estimate[2]), 
              globalG_Variance = as.numeric(global.G.results$estimate[2]),
              
              Species = species.name)
     
     local.hotspot.map <-  
          # Visualize the results with ggplot2
          ggplot() +
          geom_sf(data = species.data.all, color = "black", fill = "grey")+
          geom_sf(data = local.hotspot.sf, aes(fill = classification),color = "black", lwd = 0.1) +
          scale_fill_brewer(type = "div", palette = 5) +
          theme_void() +
          labs(
            fill = "Hot Spot Classification",
            title = paste(species.name, "mortality hotspots -", hotspot.variable),
            subtitle = paste("global G =", round(global.G.results$statistic[1,1], 2), "\np-val =", round(global.G.results$p.value[1,1], 3))
          )
     
     
     
local.hotspot.variable <- cowplot::plot_grid(local.hotspot.map, spp.map.pmort, align = "hv")
# save the hotspot map
ggsave(filename = paste0(output.dir, "/images/hotspots/hotspot_", species.name, "_", hotspot.variable, ".png"), 
       plot = local.hotspot.variable, 
       height = 4, width = 10, units = "in")
        
# return the local Gi test data
return(local.hotspot.sf)
        

}

# make maps of all spatial hotspots for pmort weighted by # trees in the county
hotspots.n_pmort <- lapply(unique(county.pmort.sf$Species), 
       function(x){get.spatial.hotspots(species.name = x, 
                                        hotspot.variable = "n_county_pmort")})
# do the same with the volfac estimates:
hotspots.volfac_pmort <- lapply(unique(county.pmort.sf$Species), 
                           function(x){get.spatial.hotspots(species.name = x, 
                                                            hotspot.variable = "pmort_weighted_1")})

hotspots.volfac_pmort.df <- do.call(rbind, hotspots.volfac_pmort)
# get summaries of gedis ord gi* and % of counties with > 1% probability mortality per year:


n_counties_over_1pct <-  county.pmort.sf %>% as.data.frame()%>% 
  filter(!is.na(pmort_weighted_1))%>%
  group_by(Species, COUNTYFP, STATEFP, STATE_NAME, geometry)%>%
  
  mutate(over.0.75.pct = ifelse(pmort_weighted_1*100 >= 0.75,  1, 0), 
         over.1.pct = ifelse(pmort_weighted_1*100 >= 1,  1, 0), 
         nplot.over.5 = ifelse(n_plots >= 2, 1, 0))%>%
  ungroup()%>%
  group_by(Species)%>%
  
  summarise(n_high_1 = sum(over.1.pct, na.rm = TRUE), 
            n_high_0.75 = sum(over.0.75.pct, na.rm = TRUE),
            n_presence_plt = sum(nplot.over.5, na.rm =TRUE),
            n_presence = n())%>%
  mutate(pct_high_1 = (n_high_1/n_presence_plt)*100, 
         pct_high_0.75 = (n_high_0.75/n_presence_plt)*100)%>%
  ungroup()%>% left_join(.,hotspots.volfac_pmort.df %>% as.data.frame()%>%
                           ungroup()%>%
                           dplyr::select(Species, globalG_pval, globalG_st_dev, globalG_estimate, globalG_expectation)%>%
                           distinct()) %>% 
  mutate(Global.G.sig = ifelse(globalG_pval <= 0.05, "significant", "n.s."), 
         Global.G.clustering = ifelse(globalG_estimate - globalG_expectation > 0 , "high", "low"))%>% 
  arrange(desc(pct_high_1))

species.hotspot.summaries <- hotspots.volfac_pmort.df %>% as.data.frame()%>%
  ungroup()%>%
  dplyr::select(Species, globalG_pval, globalG_st_dev, globalG_estimate, globalG_expectation, globalG_Variance)%>%
  distinct() %>% left_join(.,

hotspots.volfac_pmort.df %>% as.data.frame()%>%
  ungroup()%>%
  mutate(hotspot_type_sig = ifelse(p_folded_sim <=0.05 & gi >=0, "Hot", 
                                   ifelse(p_folded_sim <= 0.05 & gi <= 0, "Cold", "N.S.")))%>%
  
  group_by(Species)%>%
  mutate(n_counties = n())%>%
  group_by(Species, hotspot_type_sig)%>%
  summarise(n_sig = n(), 
            n_counties = mean(n_counties),
            pct_sig = (n()/ mean(n_counties))*100)%>%
  ungroup() %>% 
  filter(hotspot_type_sig %in% "Hot")
)

uniformity.summary <- left_join(species.hotspot.summaries, n_counties_over_1pct)%>%
  arrange(desc(pct_high_1))

uniformity.summary$Species <- factor(uniformity.summary$Species, levels = disturb.species.order)

  ggplot(uniformity.summary , aes(x = globalG_st_dev, y = pct_sig, color = Species))+
    geom_point()+
    species_color
  
  ggplot(uniformity.summary , aes(x = globalG_st_dev , y = pct_sig, 
                                  color = Global.G.sig))+
    geom_point()
  
  ggplot(uniformity.summary , aes(x = pct_sig, y = globalG_st_dev, color = Species))+
    geom_point()+
    species_color +
    theme_bw()+
    xlab("% of counties included in hotspots")+
    xlab("Global G statistic")
  
  
 uniformity.bar.plot <-  ggplot(uniformity.summary, aes(x = Species, y = pct_high_1, fill = Global.G.sig, color = Global.G.sig))+
    geom_bar(stat = "identity", alpha = 0.5)+
    geom_point()+
    scale_fill_manual(values = c("n.s." = "grey", "significant" = "red"), name = "Global G significance")+
    scale_color_manual(values = c("n.s." = "grey", "significant" = "red"), name = "Global G significance")+
    #species_fill+
    theme_bw(base_size = 14)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
          panel.grid = element_blank(), 
          legend.position = c(0.8,0.8))+ylab("% of counties with high posterior\nmortality probabilities (> 1%/year)")
  
  ggsave(filename = paste0(output.dir, "/images/hotspots/high_mort_uniformity.png"), 
         plot = uniformity.bar.plot, 
         height = 6, width = 8, units = "in")
  
  
  ggplot(uniformity.summary, aes(x = Species, y = pct_sig, fill = Species))+
    geom_bar(stat = "identity")+
    species_fill+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  

# do the hotspots for species overlap?-----
 gi.matrix <-  hotspots.volfac_pmort.df %>% as.data.frame()%>%
    ungroup()%>% dplyr::select(gi, Species, geometry)%>%
    group_by(geometry)%>%
    spread(Species, gi)
 
 library(Hmisc)
 library(reshape2)
 county.corr.gi <- rcorr(as.matrix( gi.matrix[,2:ncol( gi.matrix)]), type = "pearson" ) # or "spearman"
 colnames( county.corr.gi$r)
 # Get upper triangle of the correlation matrix
 get_upper_tri <- function(cormat){
   cormat[lower.tri(cormat)]<- NA
   return(cormat)
 }
 
 get_lower_tri<-function(cormat){
   cormat[upper.tri(cormat)] <- NA
   return(cormat)
 }
 melted.correlation <- get_lower_tri( county.corr.gi$r[disturb.species.order,disturb.species.order]) %>% reshape2::melt()%>%
   rename("r"="value") %>% left_join(., 
                                     get_lower_tri( county.corr.gi$n[disturb.species.order,disturb.species.order]) %>% reshape2::melt())%>%
   rename("n" = "value")%>%
   left_join(., 
             get_lower_tri( county.corr.gi$P[disturb.species.order,disturb.species.order]) %>% reshape2::melt())%>%
   rename("pval" = "value")%>%
   mutate(R_revised = ifelse(n>=25, round(r, digits = 1), NA), 
          P_revised = ifelse(n >=25, pval, NA))%>%
   mutate(R_revised = ifelse(R_revised == 1, NA, R_revised))%>%
   mutate(sig.label = ifelse(P_revised <= 0.1, R_revised, NA))%>%
   filter(!is.na(R_revised))# if the species are colocated in less than 10 plots omit the correlation
 
 
 gi_correlation_county <-  ggplot(data = melted.correlation, aes(Var2, Var1, fill = R_revised))+
   geom_tile(color = "white")+
   scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                        midpoint = 0, 
                        limit = c(-1,1), 
                        space = "Lab", 
                        name="Pearson\nCorrelation") +
   theme_bw(base_size = 16)+ 
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                    hjust = 0))+
   coord_fixed()+
   scale_x_discrete(position = "top",
                    limits = levels(melted.correlation$Var2))+
   scale_y_discrete(position = "left",
                    limits = rev(levels(melted.correlation$Var1)))+
   
   geom_text(aes(Var2, Var1, label = sig.label), color = "black", size = 4) +
   theme(
     axis.title.x = element_blank(),
     axis.title.y = element_blank(),
     panel.grid.major = element_blank(),
     #panel.border = element_blank(),
     panel.background = element_blank(),
     axis.ticks = element_blank(),
     legend.justification = c(1, 0),
     legend.position = c(0.9, 0.65),
     legend.direction = "horizontal")+
   guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                                title.position = "top", title.hjust = 0.5))#+coord_flip()
 
 
 ggsave(filename = paste0(output.dir, "images/posterior_county_Gi_star_species_correlations.png"), 
        plot = gi_correlation_county, 
        height = 6, width = 6, units = "in", 
        dpi = 350)
 
 
 # is mortality correlated with one another----------------------
 pmort.matrix <-  county.pmort.sf %>% as.data.frame()%>%
   ungroup()%>% dplyr::select(pmort_weighted_1, Species, geometry)%>%
   group_by(geometry)%>%
   spread(Species, pmort_weighted_1)
 
 county.corr.pmort <- rcorr(as.matrix( pmort.matrix[,2:ncol( pmort.matrix)]), type = "pearson" ) # or "spearman"
 colnames( county.corr.pmort$r)
 # Get upper triangle of the correlation matrix
 get_upper_tri <- function(cormat){
   cormat[lower.tri(cormat)]<- NA
   return(cormat)
 }
 
 get_lower_tri<-function(cormat){
   cormat[upper.tri(cormat)] <- NA
   return(cormat)
 }
 melted.correlation <- get_lower_tri( county.corr.pmort$r[disturb.species.order,disturb.species.order]) %>% reshape2::melt()%>%
   rename("r"="value") %>% left_join(., 
                                     get_lower_tri( county.corr.pmort$n[disturb.species.order,disturb.species.order]) %>% reshape2::melt())%>%
   rename("n" = "value")%>%
   left_join(., 
             get_lower_tri( county.corr.pmort$P[disturb.species.order,disturb.species.order]) %>% reshape2::melt())%>%
   rename("pval" = "value")%>%
   mutate(R_revised = ifelse(n>=25, round(r, digits = 1), NA), 
          P_revised = ifelse(n >=25, pval, NA))%>%
   mutate(R_revised = ifelse(R_revised == 1, NA, R_revised))%>%
   mutate(sig.label = ifelse(P_revised <= 0.1, R_revised, NA))%>%
   filter(!is.na(R_revised))# if the species are colocated in less than 10 plots omit the correlation
 
 
 pmort_correlation_county <-  ggplot(data = melted.correlation, aes(Var2, Var1, fill = R_revised))+
   geom_tile(color = "white")+
   scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                        midpoint = 0, 
                        limit = c(-1,1), 
                        space = "Lab", 
                        name="Pearson\nCorrelation") +
   theme_bw(base_size = 16)+ 
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                    hjust = 0))+
   coord_fixed()+
   scale_x_discrete(position = "top",
                    limits = levels(melted.correlation$Var2))+
   scale_y_discrete(position = "left",
                    limits = rev(levels(melted.correlation$Var1)))+
   
   geom_text(aes(Var2, Var1, label = sig.label), color = "black", size = 4) +
   theme(
     axis.title.x = element_blank(),
     axis.title.y = element_blank(),
     panel.grid.major = element_blank(),
     #panel.border = element_blank(),
     panel.background = element_blank(),
     axis.ticks = element_blank(),
     legend.justification = c(1, 0),
     legend.position = c(0.9, 0.65),
     legend.direction = "horizontal")+
   guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                                title.position = "top", title.hjust = 0.5))#+coord_flip()
 
 
# Are hotspots related to temporal sampling?---
 
# get remper information by county/state--
 county.meas.years <- TREE.remeas %>% filter(dbhold >= 5 & dbhcur >=5)%>%
   dplyr::select(state, county, PLOT.ID, date, remper)%>%
   distinct()%>%
   mutate(T1 = date - remper, 
          T2 = date)%>%
   group_by(county, state)%>%
   summarise(T1.avg.year = round(mean(T1)), 
             T2.avg.year = mean(T2))%>%
   mutate(COUNTYFP = str_pad(county, side = "left",3, pad = 0), 
          STATEFP = str_pad(state, side = "left", 2, pad = 0))
   
 
 # join to county-species information on measurements:
county.mort.gi.remper <- county.pmort.sf %>% as.data.frame()%>%
   ungroup()%>% left_join(., county.meas.years) %>%
  left_join(., hotspots.volfac_pmort.df) %>%
  mutate(midpoint.remper = T2.avg.year -((T2.avg.year - T1.avg.year)/2))
county.mort.gi.remper$Species <- factor(county.mort.gi.remper$Species, levels = disturb.species.order)

ggplot(data = county.mort.gi.remper, aes(x = midpoint.remper, y = pmort_weighted_1*100, color = Species)) +
  geom_pointrange(aes(ymin = pmort_weighted_1.ci.lo*100,ymax = pmort_weighted_1.ci.hi*100), 
                  position=position_jitter(width=0.5), 
                  linetype='dotted', alpha = 0.5) +
  stat_smooth(aes(x = midpoint.remper, y = pmort_weighted_1*100), method = "lm", color = "black")+
  theme_bw()+facet_wrap(~Species, scales = "free_y")

temporal_pmort <- ggplot(data = county.mort.gi.remper, aes(x = midpoint.remper, y = pmort_weighted_1*100, color = Species)) +
  geom_pointrange(aes(ymin = pmort_weighted_1.ci.lo*100,ymax = pmort_weighted_1.ci.hi*100), 
                  position=position_jitter(width=0.5), 
                  linetype='dotted', alpha = 0.5) +
  stat_smooth(aes(x = midpoint.remper, y = pmort_weighted_1*100), method = "lm", color = "black")+
  facet_wrap(~Species, scales = "free_y")+
  species_color+
  scale_x_continuous(breaks = seq(min(county.mort.gi.remper$midpoint.remper, na.rm =TRUE), 
                                  max(county.mort.gi.remper$midpoint.remper, na.rm =TRUE), by = 5))+
  theme_bw(base_size = 12)+
  theme(legend.position = "none", 
        panel.grid = element_blank())+
  ylab("County average mortality probability (%/year)")+
  xlab("Midpoint of Remeasurement Year")
   
ggsave(filename = paste0(output.dir, "images/posterior_county_pmort_over_time.png"), 
       plot = temporal_pmort, 
       height = 6, width = 8, units = "in", 
       dpi = 350)

# Do hotspots correspond to temporal trends in sampling
# plot gistar values for hotspots against midpoint of remeausrement year:
#-Yes, for beech, sugar maple, eastern white pine, possibly black cherry and White oak

temporal_gistar <- ggplot(data = county.mort.gi.remper %>% filter(!is.na(gi)))+
  geom_jitter(aes(x = midpoint.remper, y = gi, color = classification), alpha = 0.75)+
  stat_smooth(aes(x = midpoint.remper, y = gi), method = "lm", color = "black")+
  facet_wrap(~Species, scales = "free_y")+
  scale_color_brewer(type = "div", palette = 5) +
  #species_color+
  scale_x_continuous(breaks = seq(min(county.mort.gi.remper$midpoint.remper, na.rm =TRUE), 
                                  max(county.mort.gi.remper$midpoint.remper, na.rm =TRUE), by = 5))+
  theme_bw(base_size = 12)+
  theme( 
        panel.grid = element_blank())+
  ylab("County Gi* (z-score)")+
  xlab("Midpoint of Remeasurement Year")

ggsave(filename = paste0(output.dir, "images/posterior_county_Gi_star_over_time.png"), 
       plot = temporal_gistar, 
       height = 6, width = 10, units = "in", 
       dpi = 350)


# plot mortality rates over time and space
ggplot(data = county.mort.gi.remper %>% filter(!is.na(gi)))+
  geom_segment(aes(x = T1.avg.year, xend = T2.avg.year, y = pmort_weighted_1*100, color = Species), alpha = 0.75)+
  species_color+
  facet_wrap(~STATE_NAME, scales = "free_y")+
  theme_bw()

ggplot(data = county.mort.gi.remper %>% filter(!is.na(gi)))+
  geom_segment(aes(x = T1.avg.year, xend = T2.avg.year, y = pmort_weighted_1*100, color = STATE_NAME), alpha = 0.75)+
  scale_color_viridis_d()+
  facet_wrap(~Species, scales = "free_y")+
  theme_bw()

# spongy moth susceptible species:

ggplot(data = county.mort.gi.remper %>% filter(!is.na(gi)) %>% filter(Species %in% c("chestnut oak", "white oak","northern red oak", "yellow birch", "paper birch")))+
  geom_segment(aes(x = T1.avg.year, xend = T2.avg.year, y = pmort_weighted_1*100, color = Species), alpha = 0.75)+
  #scale_color_viridis_d()+
  species_color+
  facet_wrap(~STATE_NAME, scales = "free_y")+
  theme_bw()

# spruce budworm susceptible species:

ggplot(data = county.mort.gi.remper %>% filter(!is.na(gi)) %>% filter(Species %in% c("balsam fir",  "red spruce")))+
  geom_segment(aes(x = T1.avg.year, xend = T2.avg.year, y = pmort_weighted_1*100, color = Species), alpha = 0.75)+
  #scale_color_viridis_d()+
  species_color+
  facet_wrap(~STATE_NAME, scales = "free_y")+
  theme_bw()

# beech, hemlock, eastern white pine:
ggplot(data = county.mort.gi.remper %>% filter(!is.na(gi)) %>% 
         filter(Species %in% c("American beech",  "eastern hemlock", "eastern white pine")))+
  geom_segment(aes(x = T1.avg.year, xend = T2.avg.year, y = pmort_weighted_1*100, color = Species), alpha = 0.75)+
  #scale_color_viridis_d()+
  species_color+
  facet_wrap(~STATE_NAME)+
  theme_bw()

ggplot(data = county.mort.gi.remper %>% filter(!is.na(gi)) %>% 
         filter(Species %in% c("eastern hemlock")))+
  geom_segment(aes(x = T1.avg.year, xend = T2.avg.year, y = pmort_weighted_1*100, color = Species), alpha = 0.75)+
  #scale_color_viridis_d()+
  species_color+
  facet_wrap(~STATE_NAME)+
  theme_bw()

# look at the # of years with HWA present
HWA <- read.csv("data/HWA_invested_counties_2023_3_6_24.csv")  # just gets the remper states

# since this is just the time it was first observed, we want county area or something:
# Get county data for all the states
HWA.cos <- counties %>% mutate(CNTYNAME = toupper(CONAME)) %>% 
  mutate(ST = STUSPS)%>% left_join(.,HWA %>% distinct() ) # some duplicates in the data

hemlock.co.mort.gi <- county.mort.gi.remper %>% 
  filter(Species %in% "eastern hemlock")%>% 
  left_join(., HWA.cos)%>%
  ungroup()%>%
  mutate(no.years.present.T1 = T1.avg.year - YEARFIRST, 
         no.years.present.T2 = T2.avg.year - YEARFIRST)%>%
  mutate(present = ifelse(no.years.present.T2 >=-1, "present", 
                          ifelse(is.na(no.years.present.T2) & !is.na(pmort_weighted_1), "not present", "not present")))
  

hemlock.co.mort.gi


mort.hwa.detect <- ggplot(data = hemlock.co.mort.gi, aes(x = no.years.present.T2, y = pmort_weighted_1*100 ))+
  geom_point()+
  geom_errorbar(aes(x = no.years.present.T2, ymin = pmort_weighted_1.ci.lo*100, ymax = pmort_weighted_1.ci.hi*100), width = 0.25, alpha = 0.5)+theme_bw()+ylab("County average Hemlock \n mortality probability (%/year)")+
  xlab("# of years of detected HWA")

ggsave(filename = paste0(output.dir, "images/HWA_county_pmort_Hemlock.png"), 
       plot = mort.hwa.detect, 
       height = 4, width = 6, units = "in", 
       dpi = 350)

hemlock.co.mort.gi %>% mutate(present.in.periodic = ifelse( is.na(no.years.present.T2) | no.years.present.T2 < 0, "not present", "present"))%>%
  group_by(present.in.periodic) %>%
  summarise(median_annual_mort = median(pmort_weighted_1*100, na.rm =TRUE), 
            sd_annual_mort = sd(pmort_weighted_1*100, na.rm =TRUE))

hemlock.co.mort.gi %>% mutate(present.in.periodic = ifelse( is.na(no.years.present.T2) | no.years.present.T2 < 0, "not present", "present"))%>%
  group_by(present.in.periodic)|>
  ggplot()+geom_boxplot(aes(x = present.in.periodic, y = pmort_weighted_1*100))

ggplot(hemlock.co.mort.gi, aes(x = no.years.present.T2, y = gi, color = classification ))+
  geom_point()+
  scale_color_brewer(palette = "Spectral")

ggplot(hemlock.co.mort.gi, aes(x = gi, y = pmort_weighted_1,color = classification  ))+
  geom_point()+
  scale_color_brewer(palette = "Spectral")


ggplot(data = st_as_sf(hemlock.co.mort.gi))+
  geom_sf(aes(fill = YEARFIRSTC))

volfac.map.hemlock <- ggplot(data = st_as_sf(hemlock.co.mort.gi))+
  geom_sf(aes(fill = pmort_weighted_1*100))+
  scale_fill_distiller(palette = "Spectral")

years.present.HWA <- ggplot(data = st_as_sf(hemlock.co.mort.gi))+
  geom_sf(aes(fill = no.years.present.T2))+
  scale_fill_distiller(palette = "Spectral")


hemlock.co.mort.gi[is.na(hemlock.co.mort.gi$present) & !is.na(hemlock.co.mort.gi$pmort_weighted_1),]$present <- "not present"

HWA.presence.mort.barplot <- ggplot(data = hemlock.co.mort.gi %>% filter(!is.na(present)))+
  geom_jitter(aes(x = present, y = pmort_weighted_1*100, color = Species))+
  geom_boxplot(aes(x = present, y = pmort_weighted_1*100, color = Species), fill = NA, outliers = FALSE)+
  species_color+
  ylab("Eastern Hemlock\nmortality probability (%/year)")+
  xlab("Hemlock wooley adelgid presence in county")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = "none")


present.hwa <- hemlock.co.mort.gi %>% filter(present %in% "present")
notpresent.hwa <- hemlock.co.mort.gi %>% filter(present %in% "not present")

t.test(notpresent.hwa$pmort_weighted_1, present.hwa$pmort_weighted_1)

ggsave(filename = paste0(output.dir, "images/box_plot_HWA_county_pmort_Hemlock.png"), 
       plot = HWA.presence.mort.barplot  , 
       height = 6, width = 6, units = "in", 
       dpi = 350)


hemlock.co.mort.gi %>% filter(pmort_weighted_1*100 >= 0.5) %>% 
  group_by(present)%>% summarise(n())

hemlock.co.mort.gi %>% filter(!is.na(pmort_weighted_1*100)) %>% 
  group_by(present)%>% summarise(n())



ggplot()+
  geom_sf(data = st_as_sf(hemlock.co.mort.gi), aes(fill = present))+
  geom_sf(data = st_as_sf(hemlock.co.mort.gi) %>% filter(pmort_weighted_1*100 >= 0.45), aes(fill = present), color =  "red")+
  scale_fill_brewer(palette = "Spectral")+theme_minimal()


hemlock.co.mort.gi %>% filter(is.na(present)& !is.na(pmort_weighted_1))

maps.hwa.hemlock <- cowplot::plot_grid(volfac.map.hemlock, years.present.HWA, align = "hv")
ggsave(filename = paste0(output.dir, "images/map_HWA_county_pmort_Hemlock.png"), 
       plot = maps.hwa.hemlock , 
       height = 6, width = 12, units = "in", 
       dpi = 350)

# compare beech mortality to BBD data:
beech.bark <- read.csv("data/beech_bark_spread/Cale_Morin-BeechScaleDatesCanadaUS.csv") %>%
  filter(COUNTRY %in% "USA") %>%
  rename("STATE_NAME" = "PRVSTTNAME") 

beech.bark.cos <- counties %>% left_join(.,beech.bark %>% dplyr::select(-REFERENCE)) #%>% dplyr::select(-STUSPS) 
beech.co.mort.gi <- county.mort.gi.remper %>% 
  filter(Species %in% "American beech")%>% 
  left_join(., beech.bark.cos)%>%
  ungroup()%>%
  mutate(no.years.present.T1 = T1.avg.year - SCALEYR, 
         no.years.present.T2 = T2.avg.year - SCALEYR)




mort.bbd.detect <- ggplot(beech.co.mort.gi, aes(x = no.years.present.T2, y = pmort_weighted_1*100 ))+
  geom_point()+
  geom_point()+
  geom_errorbar(aes(x = no.years.present.T2, ymin = pmort_weighted_1.ci.lo*100, ymax = pmort_weighted_1.ci.hi*100), width = 0.25, alpha = 0.5)+
  theme_bw()+ylab("County average Beech \n mortality probability (%/year)")+
  xlab("# of years of detected Beech Scale")

ggsave(filename = paste0(output.dir, "images/BBD_county_pmort_Beech.png"), 
       plot = mort.bbd.detect, 
       height = 4, width = 6, units = "in", 
       dpi = 350)



ggplot(beech.co.mort.gi, aes(x = no.years.present.T2, y = gi, color = classification ))+
  geom_point()+
  scale_color_brewer(palette = "Spectral")

ggplot(beech.co.mort.gi, aes(x = gi, y = pmort_weighted_1,color = classification  ))+
  geom_point()+
  scale_color_brewer(palette = "Spectral")

ggplot(beech.co.mort.gi, aes(x = no.years.present.T2, y = pmort_weighted_1,color = STATE_NAME  ))+
  geom_point()+
  scale_color_brewer(palette = "Spectral")



volfac.map.beech <- ggplot(data = st_as_sf(beech.co.mort.gi))+
  geom_sf(aes(fill = pmort_weighted_1*100))+
  scale_fill_distiller(palette = "Spectral")

years.present.BBD <- ggplot(data = st_as_sf(beech.co.mort.gi))+
  geom_sf(aes(fill = no.years.present.T2))+
  scale_fill_distiller(palette = "Spectral")

first.scaleYR.present.BBD <- ggplot(data = st_as_sf(beech.co.mort.gi))+
  geom_sf(aes(fill = SCALEYR))+
  scale_fill_distiller(palette = "Spectral")

maps.bbd.beech <- cowplot::plot_grid(volfac.map.beech, years.present.BBD, align = "hv")
ggsave(filename = paste0(output.dir, "images/map_BBD_county_pmort_beech.png"), 
       plot = maps.bbd.beech , 
       height = 6, width = 12, units = "in", 
       dpi = 350)

# looking at spatial variation in county level mortality for the maples----
sugar.bush.county <- county.mort.gi.remper %>% 
  filter(Species %in% c("sugar maple", "red maple"))%>%
  filter(!is.na(pmort_weighted_10))%>%
  mutate(high_sugarbush = ifelse(STUSPS %in% c("VT", "NH", "ME", "NY"), "sugarbush", "not high sugarbush"))


sugar.bush.county.high.prod <- sugar.bush.county %>% filter(high_sugarbush %in% "sugarbush")
sugar.bush.county.nothigh.prod <- sugar.bush.county %>% filter(high_sugarbush %in% "not high sugarbush")

t.test(sugar.bush.county.high.prod$pmort_weighted_10*100, sugar.bush.county.nothigh.prod$pmort_weighted_10*100)

sugar.bush.county|>
  ggplot()+geom_pointrange(aes(x = STATE_NAME, y = pmort_weighted_10*100, ymin = pmort_weighted_10.ci.lo*100, ymax = pmort_weighted_10.ci.hi*100, color = Species), 
                           position = "jitter", 
                           alpha = 0.5)+species_color+theme_bw()+
  ylab("Predicted 10-year mortality risks")+
  facet_wrap(~high_sugarbush, scales = "free_x")

county.mort.gi.remper %>% filter(STUSPS %in% "OH") %>%
  filter(Species %in% c("sugar maple", "red maple"))%>%
  filter(!is.na(pmort_weighted_10))%>% arrange(desc(pmort_weighted_10))

# compare spongy moth records to mortality:--
state.area <- counties %>% mutate(ALAND.ha = ALAND/10000) %>%
  mutate(State = STATE_NAME)%>%
  group_by(State, STUSPS)%>% summarise(total.land.ha = sum(ALAND.ha))

# join county area up to spongy moth and convert spongy moth defoliation to hectares too

spongy <- read.csv( "data/NE_spongy_moth_outbreaks.csv") %>%
  mutate(HA.Defoliated = Acres.Defoliated/2.47105381)%>%
  left_join(., state.area) %>% 
  mutate(fraction.Defoliated = HA.Defoliated/total.land.ha)


spongy %>% group_by(State)%>%
  summarise(state.max.fraction.defolated = max(fraction.Defoliated*100))%>%
  arrange(desc(state.max.fraction.defolated))

library(pracma)
st.name.id <- "New York"
state.defoliation.peaks.list <- list()
for(i in 1:length(unique(spongy$State))){
  
  st.name.id <- unique(spongy$State)[i]
  
  
  
  
  
  # make and save some plots of the spongy moth defoliation records and the mortality by species:
  state.species.spongy <- county.mort.gi.remper %>% 
    filter(!is.na(gi) & STATE_NAME %in% st.name.id) %>% 
    filter(Species %in% c("chestnut oak", "white oak","northern red oak", 
                          "yellow birch", "paper birch"))
  
  if(nrow(state.species.spongy) >0){
  
  min(state.species.spongy$T1.avg.year)
  max(state.species.spongy$T2.avg.year)
  
  state.spongy <- ggplot(spongy %>% filter(State %in% st.name.id), 
                         aes(x = year, y = HA.Defoliated/1000000))+
    geom_line()+
    geom_vline(aes(xintercept =  min(state.species.spongy$T1.avg.year)), color = "red", linetype = "dashed")+
    geom_vline(aes(xintercept = max(state.species.spongy$T2.avg.year)), color = "red", linetype = "dashed")+
    facet_wrap(~State, scales = "free_y")+
    theme_bw()+
    xlim(1925, 2015)+
    ylim(0,2)+
    ylab("State Defoliation\n(millions of Hectares")
  
  
  
  state.mort <- ggplot(data =state.species.spongy )+
    geom_segment(aes(x = T1.avg.year, xend = T2.avg.year, y = pmort_weighted_1*100, color = Species), alpha = 0.75)+
    #scale_color_viridis_d()+
    species_color+
    facet_wrap(~STATE_NAME, scales = "free_y")+
    theme_bw()+
    xlim(1925, 2015)+
    ylim(0, 3.5)+
    ylab("mortality probability\n(%/year)")+
    xlab("year")
  
  ggsave(filename = paste0(output.dir, "images/hotspots/hotspots_spongy_species_", st.name.id, ".png"), 
          plot =   cowplot::plot_grid(state.mort, state.spongy, align = "hv", ncol = 1), 
         height = 4, width = 5, units = "in")

  }
  
 
    NY.data <- spongy %>% filter(State %in% st.name.id)%>%
      arrange(year)
    peaks_info <- findpeaks(NY.data$fraction.Defoliated)
    
    state.defoliation.peaks <- data.frame(
      STATE_NAME = st.name.id, 
      # height peak defoliation
      height = peaks_info[,1],
      
      # year of peak defoliation
      year.peak = NY.data[peaks_info[,2],]$year,
      
      # start and end of the peak defoliation
      start.year = NY.data[peaks_info[,3],]$year,
      end.year = NY.data[peaks_info[,4],]$year
    )
    
    state.defoliation.peaks.list[[i]] <- state.defoliation.peaks
}
spongy.defoliation.state <- do.call(rbind, state.defoliation.peaks.list)

# compare # of defoliation peaks in the remper to pmort
co.remper.years <- county.mort.gi.remper %>% 
  dplyr::select(STATE_NAME, CONAME, ALAND, T1.avg.year, T2.avg.year)%>% distinct()%>%
  filter(!is.na(T1.avg.year))%>%
  mutate(number.spongy.peaks = NA, 
         max.defoliation = NA)

for(i in 1:length(co.remper.years$STATE_NAME)){
  peaks.in.remper <- spongy.defoliation.state %>% filter(STATE_NAME %in% co.remper.years[i,]$STATE_NAME) %>%
    filter( year.peak >= co.remper.years[i,]$T1.avg.year & 
           year.peak <= co.remper.years[i,]$T2.avg.year)
  npeaks = nrow(peaks.in.remper)
  
  if(nrow(peaks.in.remper)> 0){
    co.remper.years[i,]$number.spongy.peaks = nrow(peaks.in.remper)
    co.remper.years[i,]$max.defoliation = max(peaks.in.remper$height, na.rm =TRUE)
  }else{
    co.remper.years[i,]$number.spongy.peaks = 0
  }
 
}

co.remper.years <- co.remper.years %>% group_by(STATE_NAME)%>%
  mutate(Total_state_ALAND = sum(ALAND, na.rm =TRUE)) %>% ungroup()

county.mort.gi.remper %>% left_join(.,co.remper.years) %>% 
  filter(Species %in% c("chestnut oak", "northern red oak", "white oak", "yellow birch", "paper birch"))%>%
  ungroup()%>%
  mutate(CO.LAND.FRACTION = ALAND/Total_state_ALAND)%>%
  filter(!is.na(pmort_weighted_1) & pmort_weighted_1*100 >= 1)%>%
  group_by(STATE_NAME, Species)%>%
  summarise(total.land.fraction = sum(CO.LAND.FRACTION, na.rm =TRUE), 
            maxstate.defoliation = mean(max.defoliation))|>
  ggplot()+geom_point(aes(x = total.land.fraction, y = maxstate.defoliation, color = STATE_NAME))+
  facet_wrap(~Species)+
  ylab("Maximum Defoliation Fraction")+
  xlab("Total fraction of state with > 1%/year mortality risk")


county.mort.gi.remper %>% left_join(.,co.remper.years) %>% 
  filter(Species %in% c("chestnut oak", "northern red oak", "white oak", 
                        "yellow birch", "paper birch", 
                        "eastern hemlock", "eastern white pine", 
                        "red maple", "sugar maple", 
                        "red spruce"))%>%
  ungroup()%>%
  mutate(CO.LAND.FRACTION = ALAND/Total_state_ALAND)%>%
  filter(!is.na(pmort_weighted_1))|>
  ggplot()+geom_boxplot(aes(x = number.spongy.peaks, y = pmort_weighted_1*100, group = number.spongy.peaks))+
  facet_wrap(~Species, scales = "free_y")

# from table 2 here:https://www.fs.usda.gov/foresthealth/docs/fidls/FIDL-162-SpongyMoth.pdf
spongy.susceptible.df <- data.frame(Species = unique(county.mort.gi.remper$Species), 
                                    Spongy.susceptibility = c(
                                      "Resistant", 
                                      "Resistant", 
                                      "Immune", 
                                      "Susceptible", 
                                      "Resistant", 
                                      "Resistant", 
                                      "Resistant", 
                                      "Susceptible", 
                                      "Immune", 
                                      "Susceptible", 
                                      "Resistant", 
                                      "Resistant", 
                                      "Resistant", 
                                      "Immune", 
                                      "Susceptible", 
                                      "Immune", 
                                      "Susceptible"
                                  
                                    ))

co.with.spongy <- county.mort.gi.remper %>% left_join(.,co.remper.years) %>% 
  left_join(., spongy.susceptible.df)%>%
  
  ungroup()%>%
  mutate(CO.LAND.FRACTION = ALAND/Total_state_ALAND)%>%
  filter(!is.na(pmort_weighted_1)) %>%
  mutate(CO.withpeak = ifelse(number.spongy.peaks == 0, "no spongy defoliation peaks", "1+ defoliation peaks"))

co.with.spongy %>% filter(Spongy.susceptibility %in% "Susceptible")|> 
  ggplot()+
  geom_boxplot(aes(x = CO.withpeak, y = pmort_weighted_1*100, group = CO.withpeak, fill = Species))+
  facet_wrap(~ Species, scales = "free_y")

co.with.spongy %>% filter(Spongy.susceptibility %in% "Resistant")|> 
  ggplot()+
  geom_boxplot(aes(x = CO.withpeak, y = pmort_weighted_1*100, group = CO.withpeak, fill = Species))+
  facet_wrap(~ Species, scales = "free_y")

co.with.spongy %>% filter(Spongy.susceptibility %in% "Immune")|> 
  ggplot()+
  geom_boxplot(aes(x = CO.withpeak, y = pmort_weighted_1*100, group = CO.withpeak, fill = Species))+
  facet_wrap(~ Species, scales = "free_y")


co.with.spongy %>% group_by(Species, CO.withpeak)%>%
  summarise(median.volfac = median(pmort_weighted_1*100))




defoliation.peaks <- co.with.spongy %>% filter(CO.withpeak %in% "1+ defoliation peaks")
no.defoliation.peaks <- co.with.spongy %>% filter(!CO.withpeak %in% "1+ defoliation peaks")

# differences in mortality:
t.test(defoliation.peaks$pmort_weighted_1*100, no.defoliation.peaks$pmort_weighted_1*100)

# filter for counties with at least one primary host detected:
spongy.host.present <- co.with.spongy %>% filter(Spongy.susceptibility %in% "Susceptible") %>% 
  select(STATEFP:STATE_NAME)%>% distinct()%>%
  mutate(spongy.host.presence = "Yes")

co.with.spongy <- co.with.spongy %>% left_join(., spongy.host.present)

# differences in mortality for chestnut oak:
t.test.spongy <- function(Species.name){
    defoliation.peaks.chestnut <- co.with.spongy %>% 
      filter(CO.withpeak %in% "1+ defoliation peaks" &
               spongy.host.presence %in% "Yes")%>%
      filter(Species %in% Species.name)
    no.defoliation.peaks.chestnut <- co.with.spongy %>% 
      filter(!CO.withpeak %in% "1+ defoliation peaks" &
               spongy.host.presence %in% "Yes")%>%
      filter(Species %in% Species.name)
    
    t.test(defoliation.peaks.chestnut$pmort_weighted_1*100, no.defoliation.peaks.chestnut$pmort_weighted_1*100)
}

# for the susceptible species:
t.test.spongy("chestnut oak")# significant
t.test.spongy("northern red oak") # NS
t.test.spongy("white oak")# signfficant
t.test.spongy("yellow birch")# NS
#t.test.spongy("paper birch")# not enough counties witn no defoliation observations

# for the resistant species:
#t.test.spongy("balsam fir") # no counties with no defoliatino
t.test.spongy("red spruce") # NS
t.test.spongy("eastern hemlock") # signifcantly higher
t.test.spongy("American beech")# significantly higher
t.test.spongy("hickory spp.") # N.S
t.test.spongy("eastern white pine") # significantly higher
t.test.spongy("red maple") # signifcantly higher
t.test.spongy("sugar maple") # signficantly higher

# for the immune species:
#t.test.spongy("northern white-cedar") # not enough data
t.test.spongy("black cherry") # significantly hihger
t.test.spongy("white ash")
t.test.spongy("yellow-poplar")


# differences in mortality for white oak:
defoliation.peaks.white <- co.with.spongy %>% 
  filter(CO.withpeak %in% "1+ defoliation peaks")%>%
  filter(Species %in% "white oak")
no.defoliation.peaks.white <- co.with.spongy %>% 
  filter(!CO.withpeak %in% "1+ defoliation peaks")%>%
  filter(Species %in% "white oak")

t.test(defoliation.peaks.white$pmort_weighted_1, 
       no.defoliation.peaks.white$pmort_weighted_1)
# white oaks significantly different

# differences in mortality for yellow birch:
defoliation.peaks.yb <- co.with.spongy %>% 
  filter(CO.withpeak %in% "1+ defoliation peaks")%>%
  filter(Species %in% "yellow birch")
no.defoliation.peaks.yb <- co.with.spongy %>% 
  filter(!CO.withpeak %in% "1+ defoliation peaks")%>%
  filter(Species %in% "yellow birch")

t.test(defoliation.peaks.yb$pmort_weighted_1, 
       no.defoliation.peaks.yb$pmort_weighted_1)
# yellow birch not significantly different


spongy.mort.num.peaks.plt <- co.with.spongy %>% mutate(spongy.peaks = ifelse(CO.withpeak %in% "1+ defoliation peaks", "1+", 
                                                ifelse(CO.withpeak %in% "no spongy defoliation peaks", "none", NA)))|> ggplot()+
  geom_jitter(aes(x =spongy.peaks, y = pmort_weighted_1*100, group = spongy.peaks, color = Species))+
  geom_boxplot(aes(x =spongy.peaks, y = pmort_weighted_1*100, group = spongy.peaks), fill = NA, outliers = FALSE)+
  facet_wrap(~Species, scales = "free_y", ncol = 5)+
  species_color+
  ylab("Mortality probability (%/year)")+
  xlab("Number of spongy moth defoliation peaks")+
  theme_bw()+
  theme(panel.grid = element_blank())

ggsave(filename = paste0(output.dir, "images/county_pMort_spongy_species_num_peaks.png"), 
       plot = spongy.mort.num.peaks.plt , 
       height = 4, width = 9, units = "in", 
       dpi = 350)

spongysusceptible.mort.num.peaks.plt <- co.with.spongy %>% filter(!is.na(CO.withpeak)) %>% filter(Species %in% c("chestnut oak", "white oak", "northern red oak", "yellow birch", "paper birch"))%>%
  mutate(spongy.peaks = ifelse(CO.withpeak %in% "1+ defoliation peaks", "1+", 
                                                                             ifelse(CO.withpeak %in% "no spongy defoliation peaks", "none", NA)))|> ggplot()+
  #geom_jitter(aes(x =spongy.peaks, y = pmort_weighted_1*100, group = spongy.peaks, color = Species), alpha = 0.65)+
  geom_pointrange(aes(x =spongy.peaks, y = pmort_weighted_1*100, group = spongy.peaks, color = Species,
                      ymin = pmort_weighted_1.ci.lo*100, ymax = pmort_weighted_1.ci.hi*100), 
                  position=position_jitter(width=0.25), 
                  linetype='dotted', alpha = 0.5) +
  geom_boxplot(aes(x = spongy.peaks, y = pmort_weighted_1*100, group = spongy.peaks), fill = NA, outliers = FALSE)+
  facet_wrap(~Species, scales = "free_y", ncol = 5)+
  species_color+
  ylab("Mortality probability (%/year)")+
  xlab("Number of spongy moth defoliation peaks")+
  theme_bw()+
  theme(panel.grid = element_blank())
spongysusceptible.mort.num.peaks.plt

ggsave(filename = paste0(output.dir, "images/county_pMort_spongy_species_num_peaks.png"), 
       plot = spongysusceptible.mort.num.peaks.plt , 
       height = 3.5, width = 9, units = "in", 
       dpi = 350)

county.mort.gi.remper %>% left_join(.,co.remper.years) %>% 
  filter(Species %in% c("chestnut oak", "northern red oak", "white oak", "yellow birch", "paper birch"))%>%
  ungroup()%>%
  mutate(CO.LAND.FRACTION = ALAND/Total_state_ALAND)%>%
  filter(!is.na(pmort_weighted_1) & pmort_weighted_1*100 >= 1)%>%
  group_by(STATE_NAME, Species)%>%
  summarise(total.land.fraction = sum(CO.LAND.FRACTION, na.rm =TRUE), 
            maxstate.peaks = max(number.spongy.peaks))|>
  ggplot()+geom_point(aes(x = total.land.fraction, y = maxstate.peaks, color = STATE_NAME))+
  facet_wrap(~Species)+
  ylab("Max. # of peaks during remper")+
  xlab("Total fraction of state with > 1%/year mortality risk")

county.mort.gi.remper %>% left_join(.,co.remper.years)%>%
  filter(Species %in% c("chestnut oak"))%>%
  st_as_sf()|>
  ggplot()+geom_sf(aes(fill = number.spongy.peaks))+
  scale_fill_distiller(palette = "Spectral")

county.mort.gi.remper %>% left_join(.,co.remper.years)%>%
  filter(Species %in% c("paper birch"))%>%
  st_as_sf()|>
  ggplot()+geom_sf(aes(fill = max.defoliation))+
  scale_fill_distiller(palette = "Spectral")

county.mort.gi.remper %>% left_join(.,co.remper.years)%>%
  filter(Species %in% c("chestnut oak"))%>%
  st_as_sf()|>
  ggplot()+geom_sf(aes(fill =pmort_weighted_1*100))+
  scale_fill_distiller(palette = "Spectral")

county.mort.gi.remper %>% left_join(.,co.remper.years)|>
  ggplot()+geom_point(aes(y = pmort_weighted_1, x= max.defoliation))+
  facet_wrap(~Species, scales = "free_y")


county.mort.gi.remper %>% left_join(.,co.remper.years)|>
  ggplot()+geom_point(aes(y = pmort_weighted_1, x= max.defoliation))+
  facet_wrap(~Species, scales = "free_y")

# get total land and convert to hectares
# according to tigris census.gov all units are in square meters
state.area <- counties %>% mutate(ALAND.ha = ALAND/10000) %>%
  group_by(State, STUSPS)%>% summarise(total.land.ha = sum(ALAND.ha))
# join county area up to spongy moth and convert spongy moth defoliation to hectares too
spongy <- spongy %>% left_join(., state.area) %>% 
  mutate(HA.Defoliated = Acres.Defoliated/2.47105381)%>% #conversion
  mutate(fraction.Defoliated = HA.Defoliated/total.land.ha)



# compare spruce budworm records to mortality:--
state.area <- counties %>% mutate(ALAND.ha = ALAND/10000) %>%
  mutate(State = STATE_NAME)%>%
  group_by(State, STUSPS)%>% summarise(total.land.ha = sum(ALAND.ha))

# join county area up to budworm moth and convert budworm moth defoliation to hectares too

budworm <- read.csv("data/NE_spruce_budworm_outbreaks.csv") %>%
  mutate(HA.Defoliated = Spruce.Budworm.Acres.Defoliated/2.47105381)%>%
  left_join(., state.area) %>% 
  mutate(fraction.Defoliated = HA.Defoliated/total.land.ha)%>%
  filter(!is.na(STUSPS))%>%
  rename("year" = "Year")


budworm %>% group_by(State)%>%
  summarise(max(fraction.Defoliated, na.rm =TRUE), 
            max(HA.Defoliated/1000000, na.rm =TRUE))

st.name.id <- "New York"
state.defoliation.peaks.list <- list()

for(i in 1:length(unique(budworm$State))){
  
  st.name.id <- unique(budworm$State)[i]
  
  

  
  
  # make and save some plots of the budworm moth defoliation records and the mortality by species:
  state.species.budworm <- county.mort.gi.remper %>% 
    filter(!is.na(gi) & STATE_NAME %in% st.name.id) %>% 
    filter(Species %in% c("balsam fir", "red spruce","eastern white pine", "eastern hemlock"))
  
  if(nrow(state.species.budworm) >0){
    
    min(state.species.budworm$T1.avg.year)
    max(state.species.budworm$T2.avg.year)
    
    state.budworm <- ggplot(budworm %>% filter(State %in% st.name.id), 
                           aes(x = year, y = HA.Defoliated/1000000))+
      geom_line()+
      geom_vline(aes(xintercept =  min(state.species.budworm$T1.avg.year)), color = "red", linetype = "dashed")+
      geom_vline(aes(xintercept = max(state.species.budworm$T2.avg.year)), color = "red", linetype = "dashed")+
      facet_wrap(~State, scales = "free_y")+
      theme_bw()+
      xlim(1950, 2015)+
      #ylim(0,4)+
      ylab("State Defoliation\n(millions of Hectares")
    
    
    state.fraction.budworm <- ggplot(budworm %>% filter(State %in% st.name.id), 
                            aes(x = year, y = fraction.Defoliated))+
      geom_line()+
      geom_vline(aes(xintercept =  min(state.species.budworm$T1.avg.year)), color = "red", linetype = "dashed")+
      geom_vline(aes(xintercept = max(state.species.budworm$T2.avg.year)), color = "red", linetype = "dashed")+
      facet_wrap(~State, scales = "free_y")+
      theme_bw()+
      xlim(1950, 2015)+
      #ylim(0,4)+
      ylab("State Defoliation\n(millions of Hectares")
    
    
    
    state.mort <- ggplot(data =state.species.budworm )+
      geom_segment(aes(x = T1.avg.year, xend = T2.avg.year, y = pmort_weighted_1*100, color = Species), alpha = 0.75)+
      #scale_color_viridis_d()+
      species_color+
      facet_wrap(~STATE_NAME, scales = "free_y")+
      theme_bw()+
      xlim(1925, 2015)+
      ylim(0, 5.5)+
      ylab("mortality probability\n(%/year)")+
      xlab("year")
    
    ggsave(filename = paste0(output.dir, "images/hotspots/hotspots_budworm_species_", st.name.id, ".png"), 
           plot =   cowplot::plot_grid(state.mort, state.budworm, align = "hv", ncol = 1), 
           height = 4, width = 5, units = "in")
    
  }
  
  
  NY.data <- budworm %>% filter(State %in% st.name.id)%>%
    arrange(year) %>% mutate(fraction.Defoliated = ifelse(is.na(fraction.Defoliated), 0, fraction.Defoliated))
  peaks_info <- findpeaks(NY.data$fraction)
  
  state.defoliation.peaks <- data.frame(
    STATE_NAME = st.name.id, 
    # height peak defoliation
    height = peaks_info[,1],
    
    # year of peak defoliation
    year.peak = NY.data[peaks_info[,2],]$year,
    
    # start and end of the peak defoliation
    start.year = NY.data[peaks_info[,3],]$year,
    end.year = NY.data[peaks_info[,4],]$year
  )
  
  state.defoliation.peaks.list[[i]] <- state.defoliation.peaks
}
budworm.defoliation.state <- do.call(rbind, state.defoliation.peaks.list)

# compare # of defoliation peaks in the remper to pmort
co.remper.years <- co.remper.years %>% 
  mutate(number.budworm.peaks = NA, 
         max.budworm.defoliation = NA, 
         number.budworm.peaks.10.prior = NA, 
         max.budworm.defoliation.10.prior = NA)

for(i in 1:length(co.remper.years$STATE_NAME)){
  peaks.in.remper <- budworm.defoliation.state %>% 
    filter(STATE_NAME %in% co.remper.years[i,]$STATE_NAME) %>%
    filter( year.peak >= co.remper.years[i,]$T1.avg.year & 
              year.peak <= co.remper.years[i,]$T2.avg.year)
  npeaks = nrow(peaks.in.remper)
  
  
  # get the peaks 10 years prior to the remper:
  peaks.10.prior.remper <- budworm.defoliation.state %>% 
    filter(STATE_NAME %in% co.remper.years[i,]$STATE_NAME) %>%
    filter( year.peak >= co.remper.years[i,]$T1.avg.year-10 & 
              year.peak < co.remper.years[i,]$T1.avg.year)
  npeaks.10.prior = nrow(peaks.10.prior.remper)
  
 
  
  
  if(nrow(peaks.in.remper)> 0){
    co.remper.years[i,]$number.budworm.peaks = nrow(peaks.in.remper)
    co.remper.years[i,]$max.budworm.defoliation = max(peaks.in.remper$height, na.rm =TRUE)
  }else{
    co.remper.years[i,]$number.budworm.peaks = 0
    
  }
  
  if(nrow(peaks.10.prior.remper)> 0){
    co.remper.years[i,]$number.budworm.peaks.10.prior = nrow(peaks.10.prior.remper)
    co.remper.years[i,]$max.budworm.defoliation.10.prior = max(peaks.10.prior.remper$height, na.rm =TRUE)
  }else{
    co.remper.years[i,]$number.budworm.peaks.10.prior = 0
    
  }
  
}

co.remper.years <- co.remper.years %>% group_by(STATE_NAME)%>%
  mutate(Total_state_ALAND = sum(ALAND, na.rm =TRUE)) %>% ungroup()

county.mort.gi.remper %>% left_join(.,co.remper.years) %>% 
  filter(Species %in% c("balsam fir", "red spruce", "eastern white pine", "eastern hemlock"))%>%
  ungroup()%>%
  mutate(CO.LAND.FRACTION = ALAND/Total_state_ALAND)%>%
  filter(!is.na(pmort_weighted_1) & pmort_weighted_1*100 >= 0.7)%>%
  group_by(STATE_NAME, Species)%>%
  summarise(total.land.fraction = sum(CO.LAND.FRACTION, na.rm =TRUE), 
            maxstate.defoliation = mean(max.budworm.defoliation))|>
  ggplot()+geom_point(aes(x = total.land.fraction, y = maxstate.defoliation, color = STATE_NAME))+
  facet_wrap(~Species)+
  ylab("Maximum Defoliation Fraction")+
  xlab("Total fraction of state with > 0.7%/year mortality risk")


county.mort.gi.remper %>% left_join(.,co.remper.years) %>% 
  filter(Species %in% c("balsam fir", "red spruce", "eastern white pine", "eastern hemlock"))%>%
  ungroup()%>%
  mutate(CO.LAND.FRACTION = ALAND/Total_state_ALAND)%>%
  filter(!is.na(pmort_weighted_1))|>
  ggplot()+geom_boxplot(aes(x = number.budworm.peaks + number.budworm.peaks.10.prior, y = pmort_weighted_1, group = number.spongy.peaks))+
  facet_wrap(~Species, scales = "free_y")

co.with.budworm <- county.mort.gi.remper %>% left_join(.,co.remper.years) %>% 
  filter(Species %in% c("balsam fir", "red spruce", "eastern white pine", "eastern hemlock"))%>%
  ungroup()%>%
  mutate(CO.LAND.FRACTION = ALAND/Total_state_ALAND)%>%
  filter(!is.na(pmort_weighted_1)) %>%
  mutate(CO.withpeak = ifelse(number.budworm.peaks == 0, "no budworm defoliation peaks", "1+ defoliation peaks"))

co.with.budworm |> ggplot()+geom_boxplot(aes(x = CO.withpeak, y = pmort_weighted_1*100, group = CO.withpeak, fill = Species))+
  facet_wrap(~Species, scales = "free_y")


co.with.budworm %>% group_by(Species, CO.withpeak)%>%
  summarise(median.volfac = median(pmort_weighted_1*100))




defoliation.peaks <- co.with.budworm %>% filter(CO.withpeak %in% "1+ defoliation peaks")
no.defoliation.peaks <- co.with.budworm %>% filter(!CO.withpeak %in% "1+ defoliation peaks")

# differences in mortality:
t.test(defoliation.peaks$pmort_weighted_1*100, no.defoliation.peaks$pmort_weighted_1*100)

# differences in mortality for balsam fir:
defoliation.peaks.balsam <- co.with.budworm %>% 
  filter(CO.withpeak %in% "1+ defoliation peaks")%>%
  filter(Species %in% "balsam fir")
no.defoliation.peaks.balsam <- co.with.budworm %>% 
  filter(!CO.withpeak %in% "1+ defoliation peaks")%>%
  filter(Species %in% "balsam fir")

t.test(defoliation.peaks.balsam$pmort_weighted_1*100, 
       no.defoliation.peaks.balsam$pmort_weighted_1*100)


# differences in mortality for red spruce:
defoliation.peaks.redspruce <- co.with.budworm %>% 
  filter(CO.withpeak %in% "1+ defoliation peaks")%>%
  filter(Species %in% "red spruce")
no.defoliation.peaks.redspruce <- co.with.budworm %>% 
  filter(!CO.withpeak %in% "1+ defoliation peaks")%>%
  filter(Species %in% "red spruce")

t.test(defoliation.peaks.redspruce$pmort_weighted_1, 
       no.defoliation.peaks.redspruce$pmort_weighted_1)


# differences in mortality for eastern white pine:
defoliation.peaks.wp <- co.with.budworm %>% 
  filter(CO.withpeak %in% "1+ defoliation peaks")%>%
  filter(Species %in% "eastern white pine")
no.defoliation.peaks.wp <- co.with.budworm %>% 
  filter(!CO.withpeak %in% "1+ defoliation peaks")%>%
  filter(Species %in% "eastern white pine")

t.test(defoliation.peaks.wp$pmort_weighted_1, 
       no.defoliation.peaks.wp$pmort_weighted_1)

# differences in mortality for eastern white pine:
defoliation.peaks.hemlock <- co.with.budworm %>% 
  filter(CO.withpeak %in% "1+ defoliation peaks")%>%
  filter(Species %in% "eastern white pine")
no.defoliation.peaks.hemlock <- co.with.budworm %>% 
  filter(!CO.withpeak %in% "1+ defoliation peaks")%>%
  filter(Species %in% "eastern white pine")

t.test(defoliation.peaks.hemlock$pmort_weighted_1*100, 
       no.defoliation.peaks.hemlock$pmort_weighted_1*100)

co.with.budworm %>% filter(Species %in% "balsam fir")%>%st_as_sf()|>
  ggplot()+geom_sf(aes(fill =  number.budworm.peaks + number.budworm.peaks.10.prior))

state.num.peaks <- co.with.budworm %>% group_by(STATE_NAME)%>%
  summarise(remper.peaks = median(number.budworm.peaks), 
            prior.10.peaks = median(number.budworm.peaks.10.prior))%>%
  mutate(total = remper.peaks + prior.10.peaks)%>%
  arrange(total)

co.with.budworm$STATE_NAME <- factor(co.with.budworm$STATE_NAME, levels = state.num.peaks$STATE_NAME)
remper.prior.peaks.volfac <- co.with.budworm %>%
  mutate(total.peak.budworm.defoliation = as.character(as.integer(number.budworm.peaks + number.budworm.peaks.10.prior)))|>
  ggplot()+
  geom_jitter(aes(x = STATE_NAME, y = pmort_weighted_1*100, color = total.peak.budworm.defoliation))+
  geom_boxplot(aes(x = STATE_NAME, y = pmort_weighted_1*100, color =total.peak.budworm.defoliation), fill = NA)+
  scale_color_manual(values = c("4" = "#d7191c",
                                
                                "3" ="#fc8d59",
                                "2" ="#fee090",
                
                                "1" ="#abd9e9",
                                "0" ="#2c7bb6"),
                     name = "budworm defoliation peaks\n(10 years prior to T1 - T2)")+
  facet_wrap(~Species, scales = "free_y", ncol = 4)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        panel.grid = element_blank())+
  ylab("Mortality probability (%/year)")+
  xlab("")


remper.peaks.volfac <- co.with.budworm %>%
  mutate(total.peak.budworm.defoliation = as.character(as.integer(number.budworm.peaks + number.budworm.peaks.10.prior)), 
         peak.budworm.defolation.remper = as.character(as.integer(number.budworm.peaks)), 
         peak.budworm.defoliation.prior = as.character(as.integer(number.budworm.peaks.10.prior)))|>
  ggplot()+
  geom_jitter(aes(x = STATE_NAME, y = pmort_weighted_1*100, color = peak.budworm.defolation.remper))+
  geom_boxplot(aes(x = STATE_NAME, y = pmort_weighted_1*100, color =peak.budworm.defolation.remper), fill = NA)+
  scale_color_manual(values = c("4" = "#d7191c",
                                
                                "3" ="#fc8d59",
                                "2" ="#fee090",
                                
                                "1" ="#abd9e9",
                                "0" ="#2c7bb6"),
                     name = "budworm defoliation peaks\n(T1 - T2)")+
  facet_wrap(~Species, scales = "free_y", ncol = 4)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        panel.grid = element_blank())+
  ylab("Mortality probability (%/year)")+
  xlab("")


prior.peaks.volfac <- co.with.budworm %>%
  mutate(total.peak.budworm.defoliation = as.character(as.integer(number.budworm.peaks + number.budworm.peaks.10.prior)), 
         peak.budworm.defolation.remper = as.character(as.integer(number.budworm.peaks)), 
         peak.budworm.defoliation.prior = as.character(as.integer(number.budworm.peaks.10.prior)))|>
  ggplot()+
  geom_jitter(aes(x = STATE_NAME, y = pmort_weighted_1*100, color = peak.budworm.defoliation.prior), size = 0.1)+
  geom_boxplot(aes(x = STATE_NAME, y = pmort_weighted_1*100, color =peak.budworm.defoliation.prior), fill = NA)+
  scale_color_manual(values = c("4" = "#d7191c",
                                
                                "3" ="#fc8d59",
                                "2" ="#fee090",
                                
                                "1" ="#abd9e9",
                                "0" ="#2c7bb6"),
                     name = "budworm defoliation peaks\n(T1 - T2)")+
  facet_wrap(~Species, scales = "free_y", ncol = 4)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        panel.grid = element_blank())+
  ylab("Mortality probability (%/year)")+
  xlab("")


ggsave(filename = paste0(output.dir, "images/county_pMort_budworm_species_num_peaks_remper_prior.png"), 
       plot = remper.prior.peaks.volfac , 
       height = 4, width = 12, units = "in", 
       dpi = 350)

budworm.mort.num.peaks.plt <- co.with.budworm %>% 
  mutate(budworm.peaks = ifelse(CO.withpeak %in% "1+ defoliation peaks", "1+", 
                                ifelse(CO.withpeak %in% "no budworm defoliation peaks", "none", NA)))|> ggplot()+
  #geom_jitter(aes(x =budworm.peaks, y = pmort_weighted_1*100, group = budworm.peaks, color = Species), alpha = 0.65)+
  
  geom_pointrange(aes(x = budworm.peaks, y = pmort_weighted_1*100, group = budworm.peaks, color = Species,
                      ymin = pmort_weighted_1.ci.lo*100, ymax = pmort_weighted_1.ci.hi*100), 
                  position=position_jitter(width=0.25), 
                  linetype='dotted', alpha = 0.5) +
  
  geom_boxplot(aes(x =budworm.peaks, y = pmort_weighted_1*100, group = budworm.peaks), fill = NA, outliers = FALSE)+
  facet_wrap(~Species, scales = "free_y", ncol = 5)+
  species_color+
  ylab("Mortality probability (%/year)")+
  xlab("Number of budworm moth defoliation peaks")+
  theme_bw()+
  theme(panel.grid = element_blank())

ggsave(filename = paste0(output.dir, "images/county_pMort_budworm_species_num_peaks.png"), 
       plot = budworm.mort.num.peaks.plt , 
       height = 3.5, width = 9, units = "in", 
       dpi = 350)


# figure 6:
figure6 <- plot_grid(budworm.mort.num.peaks.plt, mort.bbd.detect,
          spongysusceptible.mort.num.peaks.plt, mort.hwa.detect, ncol =2, 
          rel_widths = c(1, 0.5, 1, 0.5), 
          labels = c("A", "C", "B", "C"))

ggsave(filename= paste0(output.dir, "images/Figure_6_mortality_estimates_vs_data.png"), 
       figure6, 
       height = 6, width = 12)

co.with.budworm %>% group_by(STATE_NAME)%>%
  summarise(max(number.budworm.peaks), 
            max(max.budworm.defoliation, na.rm =TRUE))


co.with.budworm %>% 
  mutate(budworm.peaks = ifelse(CO.withpeak %in% "1+ defoliation peaks", "1+", 
                                ifelse(CO.withpeak %in% "no budworm defoliation peaks", "none", NA)))|> 
  ggplot()+
  geom_jitter(aes(x =STATE_NAME, y = pmort_weighted_1*100, group = Species, color = Species))+
  geom_boxplot(aes(x =STATE_NAME, y = pmort_weighted_1*100, group = STATE_NAME), fill = NA, outliers = FALSE)+
  facet_wrap(~Species, scales = "free_y", ncol = 5)+
  species_color+
  ylab("Mortality probability (%/year)")+
  xlab("Number of budworm moth defoliation peaks")+
  theme_bw()+
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1))


# Quantifying Specificity:--------
# Is there spatial overlap across species on county pmort?--- 

species.data.all <- county.pmort.sf %>% mutate(
  cut.mort.rate = cut(pmort_weighted_1*100, breaks = c(0, 0.25, 0.5, 0.75, 1, 1.5, 2, 4, 100))
)
cut.mort.values <- data.frame(
  cut.mort.rate = c("(0,0.25]", "(0.25,0.5]", "(0.5,0.75]", "(0.75,1]", "(1,1.5]", "(1.5,2]", "(2,4]", "(4,100]"), 
  mortality.rate = c("< 0.25", "0.25 - 0.5", "0.5 - 0.75", "0.75 - 1", "1 - 1.5", "1.5 - 2", "2 - 4", "> 4"), 
  red.hex.colors = c(
    "#ffffcc",
    "#ffeda0",
    "#fed976",
    "#feb24c",
    "#fd8d3c",
    "#fc4e2a",
    "#e31a1c",
    "#b10026"),
  hex.colors = c(
    "#fff7f3",
    "#fde0dd",
    "#fcc5c0",
    "#fa9fb5",
    "#f768a1",
    "#dd3497",
    "#ae017e",
    "#7a0177"))

hex.vector <- as.vector(cut.mort.values$red.hex.colors)
names( hex.vector) <- cut.mort.values$mortality.rate
fill_mort_rate_red <- scale_fill_manual(values = hex.vector, name = "Mortality\nRate\n(%/year)", drop = FALSE)


hex.vector <- as.vector(cut.mort.values$hex.colors)
names( hex.vector) <- cut.mort.values$mortality.rate
fill_mort_rate_pink <- scale_fill_manual(values = hex.vector, name = "Mortality\nRate\n(%/year)", drop = FALSE)

species.data.all <-  species.data.all %>% left_join(., cut.mort.values)
species.data.all$mortality.rate <- factor( species.data.all$mortality.rate, levels = c("< 0.25", "0.25 - 0.5", "0.5 - 0.75", "0.75 - 1", "1 - 1.5", "1.5 - 2", "2 - 4", "> 4"))

# how widespread is high mortality probability across the range?

# if high mortality is >1% per year, sum up the # of counties with greater than that
county.high.pmort.species <- species.data.all %>% as.data.frame()%>% 
  filter(!is.na(pmort_weighted_1))%>%
  ungroup()%>%
  # calculate species mean pmort 
  group_by(Species)%>%
  mutate(mean.pmort = mean(pmort_weighted_1, na.rm =TRUE), 
         median.pmort = median(pmort_weighted_1, na.rm =TRUE), 
         pmort.75 = quantile(pmort_weighted_1,0.75, na.rm =TRUE))%>%
  ungroup()%>%
  mutate(over.threshold = ifelse(pmort_weighted_1*100 >= pmort.75*100, 1, 0))%>%
  ungroup()

county.high.pmort.species  %>% 
  group_by(STATEFP, STATE_NAME, COUNTYFP, CONAME, geometry)%>%
  summarise(nspecies_high = sum(over.threshold, na.rm =TRUE))%>%
  st_as_sf()|>
  ggplot() +
  geom_sf(aes(fill = nspecies_high),color = "black", lwd = 0.1, show.legend = TRUE) +
  scale_fill_viridis_c()

high.species <- county.high.pmort.species  %>% as.data.frame()%>%
 
  group_by(STATEFP, STATE_NAME, COUNTYFP, CONAME, geometry, Species)%>%
  summarise(high_species= sum(over.threshold, na.rm =TRUE))%>%
  ungroup()%>%
  dplyr::select(STATEFP, STATE_NAME, COUNTYFP, CONAME, geometry, Species, high_species)%>%
  group_by(STATEFP, STATE_NAME, COUNTYFP, CONAME, geometry)%>%
  spread(Species, high_species, fill = 0)%>%
  ungroup()%>%
  mutate(site = 1:length(STATEFP))%>%
  as.data.frame()


all.species.presence <- county.high.pmort.species  %>% as.data.frame()%>%
  
  group_by(STATEFP, STATE_NAME, COUNTYFP, CONAME, geometry, Species)%>%
  summarise(spp_presence = n())%>%
  ungroup()%>%
  dplyr::select(STATEFP, STATE_NAME, COUNTYFP, CONAME, geometry, Species, spp_presence)%>%
  group_by(STATEFP, STATE_NAME, COUNTYFP, CONAME, geometry)%>%
  spread(Species, spp_presence, fill = 0)%>%
  ungroup()%>%
  mutate(site = 1:length(STATEFP))%>%
  as.data.frame()


count_high_mort_matrix <- t(as.matrix(high.species[,6:22])) %*% as.matrix(high.species[,6:22])

count_total_cooccurance_matrix <- t(as.matrix(all.species.presence[,6:22])) %*% as.matrix(all.species.presence[,6:22])

# Convert to a data frame and make row names a column
count_high_df <- as.data.frame(as.matrix(count_high_mort_matrix[disturb.species.order, disturb.species.order]))
count_high_df$Species1 <- rownames(count_high_df)

count_cooccurance_df <- as.data.frame(as.matrix(count_total_cooccurance_matrix[disturb.species.order, disturb.species.order]))
count_cooccurance_df$Species1 <- rownames(count_cooccurance_df)


get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

count_high_long <- count_high_df %>% get_upper_tri()%>% reshape2::melt(., id.vars = "Species1")%>%
  rename("Species2" = "variable", 
         "Counties_high" = "value")%>%
  #filter(!Species1 == Species2)%>%
  filter(!is.na(Counties_high))
count_high_long$Species1 <- factor(count_high_long$Species1, levels = disturb.species.order)
count_high_long$Species2 <- factor(count_high_long$Species2, levels = disturb.species.order)

count_cooccurance_long <- count_cooccurance_df %>% get_upper_tri()%>% reshape2::melt(., id.vars = "Species1")%>%
  rename("Species2" = "variable", 
         "co.occurance" = "value")
count_cooccurance_long$Species1 <- factor(count_cooccurance_long$Species1, levels = disturb.species.order)
count_cooccurance_long$Species2 <- factor(count_cooccurance_long$Species2, levels = disturb.species.order)


co.occurance.mort.species <- count_high_long %>% 
  left_join(.,count_cooccurance_long)%>% ungroup()%>%
  mutate(prop.cooccurance = Counties_high/co.occurance)%>%
  mutate(prop.cooccurance.v2 = ifelse(co.occurance < 10, NA, prop.cooccurance))%>%
  filter(!Species1 == Species2)

co.occurance.high_mort_rates <- co.occurance.mort.species %>% 
  ggplot()+geom_tile(aes(y = Species2, x = Species1, fill = prop.cooccurance.v2), color = "black")+
  scale_fill_distiller(palette = "OrRd", direction = 1, na.value = "grey")+
  theme_minimal()+
  coord_equal()+
  scale_x_discrete(position = "top",
                   limits = levels(count_high_long$Species2))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                   hjust = 0), 
        axis.title = element_blank())+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    #panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.85, 0.15),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))#+coord_flip()


  
  
ggsave(filename = paste0(output.dir, "images/cooccurance_high_species_county_pMort.png"), 
       plot = co.occurance.high_mort_rates, 
       height = 6, width = 6, units = "in", 
       dpi = 350)

  
  


county.high.pmort.species %>% group_by(Species)%>%
  summarise(ncounties = n(), 
            n.over.threshold = sum(over.threshold, na.rm =TRUE))%>%
  mutate(pct.counties.high = round(n.over.threshold/ncounties*100, 1))




species.ranked.percentiles <- species.data.all %>% as.data.frame()%>% 
  filter(!is.na(pmort_weighted_1))%>%
  ungroup()%>%
  # calculate species mean pmort 
  group_by(Species) %>%
  mutate(mean.mort = mean(pmort_weighted_1*100, na.rm =TRUE), 
         sd.mort = sd(pmort_weighted_1*100, na.rm =TRUE))%>%
  #mutate(mort.percentile = percent_rank(pmort_weighted_1, na.rm =TRUE))%>%
  ungroup()%>%
  group_by(STATEFP, STATE_NAME, COUNTYFP, CONAME, geometry, Species)%>%
  mutate(standardized.pmort = ((pmort_weighted_1*100)-mean.mort)/sd.mort, 
         log.pmort = log(pmort_weighted_1*100))%>%
  #View()
  ungroup()%>%
  dplyr::select(STATEFP, STATE_NAME, COUNTYFP, CONAME, geometry, Species, standardized.pmort)%>%
  ungroup()%>%
  group_by(STATEFP, STATE_NAME, COUNTYFP, CONAME, geometry)%>%
  spread(Species, standardized.pmort, fill = NA)%>%
  ungroup()#%>% data.frame()
  

# Correlation of county standardized species mortality probabilities----
# get correlations by species:
county.corr.log.pmort <- rcorr(as.matrix(species.ranked.percentiles[,6:ncol(species.ranked.percentiles)]), type = "spearman" ) # or "spearman"
colnames(county.corr.log.pmort$r)


# assign any NA values to 0:

county.corr.log.pmort$r[is.na(county.corr.log.pmort$r)] <- 0
# also assign non-significant values to zero
county.corr.log.pmort$r[county.corr.log.pmort$P > 0.05] <- 0

distance_matrix <- as.dist(1 - county.corr.log.pmort$r)

hc_result <- hclust(distance_matrix, method = "ward.D2") 
hc_result_cormat <- hclust(as.dist(county.corr.log.pmort$r), 
                           method = "ward.D2")
# Plot the dendrogram
plot(hc_result, main = "Hierarchical Clustering based on Correlation")
plot(hc_result_cormat, main = "Hierarchical Clustering based on Correlation")

# ordering by shared hosts?
library(vegan)

pest.metadata <- read.csv(paste0(output.dir, "data/pest_metadata_new_england.csv"))
host.pests.general <- read.csv(paste0(output.dir, "data/host_pests_general.csv"), check.names = F)
host.pests.mat <- host.pests.general%>% column_to_rownames("Species") %>% as.matrix()


# cluster the species (rows)
host.pests.general.dist <- vegdist(host.pests.mat, method ="euclidian")
species_clust <- hclust(host.pests.general.dist, method = "ward.D2") 
plot(species_clust, main = "Hierarchical Clustering based on forest pests")
species_order <- species_clust$labels[species_clust$order]

species_order <-  c("northern red oak",
                    "white oak",
                    "chestnut oak" ,
                    "yellow birch", 
                    "paper birch" ,
                    "sugar maple" , 
                    "red maple" , 
                    "black cherry",
                    "American beech",
                    "white ash",
                    "hickory spp.",
                    "yellow-poplar",
                    "eastern white pine",
                    "eastern hemlock",
                    "northern white-cedar", 
                    "red spruce", 
                    "balsam fir") 

pests.dist <- vegdist(t(host.pests.mat), method ="euclidian")
pest_clust <- hclust(pests.dist , method = "ward.D2") 
plot(pest_clust, main = "Hierarchical Clustering based on forest pests")
pest_order <- pest_clust$labels[pest_clust$order]




county.corr.log.pmort <- rcorr(as.matrix(species.ranked.percentiles[,6:ncol(species.ranked.percentiles)]), type = "spearman" ) # or "spearman"


# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
species_order
disturb.species.order




#id which pests may be at play:

host.pests.m <- host.pests.general %>% melt(., id.vars = "Species") %>% 
  mutate(susceptible = ifelse(value == 0, NA, ifelse(value == 1,"occasional host", ifelse(value == 2, "preferred host", NA)))) %>%
  mutate(pest = tolower(variable)) 

correlations.m <- correlations %>% melt(.) %>% rename("mort.cor"="value") 

host.pests.m$Species <- factor(host.pests.m$Species, levels =species_order)
host.pests.m$variable <- factor(host.pests.m$variable, levels =pest_order)


forest.host.bar <- ggplot(data = host.pests.m)+geom_col(aes(x = variable, y = value, fill = Species))+
  species_fill+
  theme_bw()+
  ylab("Impact on species")+xlab("")+
  theme(axis.text.x = element_text(angle = 90, hjust =1, vjust = 0.5))

ggsave(filename = paste0(output.dir, "images/forest_pests_by_host.png"), 
       plot = forest.host.bar, 
       height = 6, width = 6, units = "in", 
       dpi = 350)



forest.host.bar.preferred <- ggplot(data = host.pests.m %>% filter(value == 2))+geom_col(aes(x = variable, y = value, fill = Species))+
  species_fill+
  theme_bw()+
  ylab("Impact on species")+xlab("")+
  theme(axis.text.x = element_text(angle = 90, hjust =1, vjust = 0.5))

ggsave(filename = paste0(output.dir, "images/forest_pests_by_host_preferred.png"), 
       plot = forest.host.bar.preferred, 
       height = 6, width = 6, units = "in", 
       dpi = 350)





host.species.tile <- ggplot(data =host.pests.m )+geom_tile(aes(x = Species, y = variable, fill = susceptible))+
  
  theme_bw()+
  ylab("")+xlab("")+
  theme(axis.text.x = element_text(angle = 90, hjust =1, vjust = 0.5))

ggsave(filename = paste0(output.dir, "images/host_pest_tile.png"), 
       plot = host.species.tile, 
       height = 6, width = 6, units = "in")



melted.correlation <- get_lower_tri(county.corr.log.pmort$r[species_order,species_order]) %>% reshape2::melt()%>%
  rename("r"="value") %>% left_join(., 
                                    get_lower_tri(county.corr.log.pmort$n[species_order,species_order]) %>% reshape2::melt())%>%
  rename("n" = "value")%>%
  left_join(., 
            get_lower_tri(county.corr.log.pmort$P[species_order,species_order]) %>% reshape2::melt())%>%
  rename("pval" = "value")%>%
  mutate(R_revised = ifelse(n>=50, round(r, digits = 1), NA), 
         P_revised = ifelse(n >=50, pval, NA))%>%
  mutate(R_revised = ifelse(R_revised == 1, NA, R_revised))%>%
  mutate(sig.label = ifelse(P_revised <= 0.05, R_revised, NA))%>%
  filter(!is.na(R_revised))# if the species are colocated in less than 10 plots omit the correlation


pmort_correlation_county <-  ggplot(data = melted.correlation, aes(Var2, Var1, fill = R_revised))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-0.6,0.6), space = "Lab", 
                       name="Spearman\nCorrelation") +
  theme_bw(base_size = 16)+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                   hjust = 1))+
  coord_fixed()+
  scale_x_discrete(position = "bottom",
                   limits = levels(melted.correlation$Var2))+
  scale_y_discrete(position = "left",
                   limits = rev(levels(melted.correlation$Var1)))+
  
  geom_text(aes(Var2, Var1, label = sig.label), color = "black", size = 3) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    #panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.9, 0.65),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))#+coord_flip()
pmort_correlation_county

ggsave(filename = paste0(output.dir, "images/posterior_county_pMort_species_correlations.png"), 
       plot = pmort_correlation_county, 
       height = 6, width = 6, units = "in", 
       dpi = 350)


# save figure 5:----

figure5 <- plot_grid(plot_grid(uniformity.bar.plot+theme(axis.title.x = element_blank(), legend.position = c(0.5, 0.6),
                                                         legend.key.size = unit(0.5, "cm"),
                                                         legend.text = element_text(size = 10),
                                                         legend.title = element_text(size = 12)), 
                               pmort_correlation_county + theme(axis.title = element_blank(), legend.position = c(0.88, 0.6), 
                                                                legend.direction = "horizontal",legend.key.size = unit(0.25, "cm"),
                                                                legend.text = element_text(size = 10),
                                                                legend.title = element_text(size = 12), 
                                                                legend.background = element_rect(fill = "transparent", color = NA)), 
                                ncol = 2, 
                               #align = "h",
                               rel_widths = c(1, 1.5),
                               rel_heights = c(0.5, 1.5),
                               labels = c("A", "B")),
temporal_pmort, labels = c(NA, "C"), 
rel_heights = c(1,1.1),
rel_widths = c(1,0.8),
ncol = 1)

figure5

ggsave(filename = paste0(output.dir, "images/figure_5_hotspots_specificity_synchronicity.png"), 
       plot = figure5, 
       height = 11, width = 8.5, units = "in", bg = "white", 
       dpi = 350)

uniformity.bar.plot/pmort_correlation_county

# Visualize correlations for highly correlated species on the map:

# identify significantly positive correlated species:
cor.pairs <- melted.correlation %>% filter(R_revised >=0.2 & P_revised <=0.05)

for(i in 1:nrow(cor.pairs)){
  species1 <- cor.pairs[i,]$Var1
  species2 <- cor.pairs[i,]$Var2
  
  overlap <- species.ranked.percentiles[,species1] > 1 & species.ranked.percentiles[,species2] >1
  cat(as.character(species1), "-", as.character(species2), ":", sum(overlap, na.rm =TRUE), "overlappping high counties \n")
  
  }

# get index for all high values for 
hot <- species.ranked.percentiles[,species_order] >= 1
storage.mode(hot) <- "integer"

pair_overlap_long <- purrr::map_dfr(seq_len(nrow(cor.pairs)), function(i){
  Species1 <- cor.pairs[i,]$Var1
  Species2 <- cor.pairs[i,]$Var2
  
  tibble(species.ranked.percentiles %>% select(STATEFP, COUNTYFP, STATE_NAME, CONAME, geometry), 
         Spp1 = Species1, 
         Spp2 = Species2, 
         integer.overlap = as.integer(hot[, Species1] & hot[,Species2]),
         continous.overlap = (species.ranked.percentiles[,Species1]*species.ranked.percentiles[,Species2])[,1])
})

# map the burden for each county:

burden_county.sf <- pair_overlap_long %>% group_by(STATEFP, COUNTYFP, STATE_NAME, CONAME, geometry) %>%
  summarise(n_pairs_overlapping = sum(integer.overlap, na.rm =TRUE), 
            burden_overlapping = sum(continous.overlap, na.rm =TRUE)) %>% st_as_sf()

any.overlaps <- burden_county.sf %>% arrange(desc(n_pairs_overlapping)) %>% 
  filter(n_pairs_overlapping >=1)

ggplot(burden_county.sf)+geom_sf(aes(fill = n_pairs_overlapping))+
  scale_fill_viridis_c()

ggplot(burden_county.sf)+geom_sf(aes(fill = burden_overlapping))+
  scale_fill_viridis_c()

# get the overlaps between sugar maple and other species:

plot.overlapping <- function(species){

  overlapping.pair <- pair_overlap_long %>% 
    filter(Spp1 %in% species | Spp2 %in% species) %>%
    mutate(pair = paste0(Spp1, " x ", Spp2)) %>% 
    st_as_sf() 
  
  overlapping.pair|>
    ggplot()+geom_sf(aes(fill = integer.overlap))+
    facet_wrap(~pair, drop = TRUE)+scale_fill_viridis_c()
  
  # overlapping.pair %>% group_by(STATEFP, COUNTYFP, STATE_NAME, CONAME, geometry) %>%
  #   summarise(n_pairs_overlapping = sum(integer.overlap, na.rm =TRUE)) %>% st_as_sf()|>
  #   ggplot()+geom_sf(aes(fill = n_pairs_overlapping))+
  #  scale_fill_viridis_c()
  
  # overlapping.pair %>% as.data.frame()%>% filter(integer.overlap == 1) %>% 
  #   group_by(STATEFP, COUNTYFP, STATE_NAME, CONAME, geometry) %>% 
  #   summarise(n_high_species = sum(integer.overlap, na.rm =TRUE), 
  #             high_species = str_c(unique(unlist(across(c(Spp1, Spp2)))), collapse = ", ")) %>% 
  #   ungroup()%>%
  #   arrange(desc(n_high_species)) %>% select(STATE_NAME, CONAME, n_high_species, high_species)|>
  #   gt()

}

# oaks: 
plot.overlapping(species = "white oak")
plot.overlapping(species = "chestnut oak")
plot.overlapping(species = "northern red oak")

# maples:
plot.overlapping(species = "sugar maple")
plot.overlapping(species = "red maple")

# eastern hemlock
plot.overlapping(species = "eastern hemlock")
plot.overlapping(species = "black cherry") # over laps with red maple & beech: half wing geometer (r maple, )

plot.overlapping(species = "paper birch") 

plot.overlapping(species = "white ash") 

# alternative approach: id the common species for each pest and see which counties these have high mortality in


pair_overlap_long %>% filter(Spp2 %in% "red maple") %>%
  st_as_sf() |>
  ggplot()+geom_sf(aes(fill = integer.overlap))+
  facet_wrap(~Spp1)


pair_overlap_long %>% filter(Spp1 %in% "chestnut oak") %>%
  st_as_sf() |>
  ggplot()+geom_sf(aes(fill = integer.overlap))+
  facet_wrap(~Spp2)

pair_overlap_long %>% filter(Spp1 %in% "chestnut oak") %>%
  st_as_sf() |>
  ggplot()+geom_sf(aes(fill = integer.overlap))+
  facet_wrap(~Spp2)

#species.data.all |>
spp.map.pmort <- ggplot() +
  geom_sf(data = species.data.all, aes(fill = mortality.rate),color = "black", lwd = 0.1, show.legend = TRUE) +
  
  fill_mort_rate_pink+
  theme_void() +
  labs(
    fill = paste("mortality probability\n(%/year)"),
    title = paste(species.name, "-", hotspot.variable), 
    subtitle = paste0(pct.counties.high$pct.counties.high, "% of counties > 0.7%/year" ))



ggplot(data = county.pmort.sf %>% filter(Species %in% "white oak")%>% 
         mutate(ncount_pmort_plot_3 = ifelse(n_plots > 3, n_county_pmort, NA)))+
  geom_sf(aes(fill = ncount_pmort_plot_3), color = NA)+scale_fill_viridis_c()+
  facet_wrap(~Species)

ggplot(data = county.pmort.sf %>% filter(Species %in% "eastern hemlock")%>% 
         mutate(ncount_pmort_plot_3 = ifelse(n_plots > 3, n_county_pmort, NA)))+
  geom_sf(aes(fill = ncount_pmort_plot_3), color = NA)+scale_fill_viridis_c()+
  facet_wrap(~Species)

ggplot(data = county.pmort.sf %>% filter(Species %in% "eastern hemlock")%>% 
         mutate(ncount_pmort_plot_3 = ifelse(n_plots > 5, n_county_pmort, NA)))+
  geom_sf(aes(fill = ncount_pmort_plot_3), color = NA)+scale_fill_viridis_c()+
  facet_wrap(~Species)

ggplot(data = county.pmort.sf %>% filter(Species %in% "eastern hemlock")
       %>% mutate(ncount_pmort_plot_3 = ifelse(n_plots > 5, n_county_pmort, 0)))+
  geom_sf(aes(fill = pmort_weighted_1), color = NA)+scale_fill_viridis_c()+
  facet_wrap(~Species)






avg.plt.mort.annual <- spp.tree.mort.probs %>% group_by(state, PLOT.ID, county, Species, SPCD, LONG_FIADB, LAT_FIADB)%>%
  summarise(pMort_annual_plot = median(p1year.mortality, na.rm =TRUE), 
            pMort_annual_plot.sd = sd(p1year.mortality, na.rm =TRUE), 
            nTrees = n())

ggplot(data = avg.plt.mort.annual)+
  geom_density(aes(y = pMort_annual_plot*100, fill = Species))+
  facet_wrap(~Species, scales = "free")+species_fill

ggplot(data = avg.plt.mort.annual)+
  geom_boxplot(aes(x = Species, y = pMort_annual_plot*100, fill = Species))+
  species_fill

ggplot(data = avg.plt.mort.annual)+
  geom_density(aes(y = pMort_annual_plot*100, fill = Species))+
  facet_wrap(~Species, scales = "free")+species_fill

ggplot(data = avg.plt.mort.annual)+
  geom_point(aes(x = LONG_FIADB, y = LAT_FIADB, color = pMort_annual_plot*100), size = 0.5)+
  facet_wrap(~Species)+
  scale_color_viridis_c()
hist(avg.plt.mort.annual$pMort_annual_plot*100)

# Create discrete categories from continuous_var
avg.plt.mort.annual$pMort_annual_plot_discrete <- cut(avg.plt.mort.annual$pMort_annual_plot*100,
                                breaks = c(0, 0.5, 1, 1.5, 2, 5, 7.5, 10, 20  ),
                                labels = c("< 0.5", "0.5 - 1", "1 - 1.5", "1.5 - 2", "2 - 5", "5 - 7.5", "7.5 - 10", "> 10"))
nspp$COMMON

disturb.species.order <- c("balsam fir" ,"red spruce","northern white-cedar","eastern hemlock",     
                          "American beech" ,"chestnut oak","northern red oak",    
                           "white oak", "yellow birch" ,"paper birch" ,"hickory spp.",        
                            "eastern white pine","red maple","sugar maple" ,"black cherry",      
                            "white ash" ,"yellow-poplar"   )

species.plot.mortality.list <- list()
for(i in 1:length(disturb.species.order)){
  
species.plot.mortality.list[[i]] <- ggplot() + 
  geom_polygon(data = canada, 
               aes(x=long, y=lat, group = group), 
               color = "black", fill = "white") +
  geom_polygon(data = lakes_df , 
               aes(x = long, y = lat, group = group), 
               color = "black", fill = "lightblue") +
  geom_polygon(data = state_sub, 
               aes(x=long, y=lat, group = group), 
               color = "black", fill = "white") +
 
  geom_point(data = avg.plt.mort.annual %>% filter(Species %in% disturb.species.order[i]) , aes(x = LONG_FIADB, y = LAT_FIADB, color = pMort_annual_plot_discrete), size = 0.1)+
  scale_color_manual(values =  c(
    "#74add1",
    "#abd9e9",
   
    "#fcbba1",
    "#fc9272",
    "#fb6a4a",
    "#ef3b2c",
    "#cb181d",
    "#99000d"), name = "Mortality Probability\n(%/year)")+
  

  coord_sf(xlim = c(-85, -67.5), ylim = c(37, 47.5))+theme(panel.grid = element_blank(), #panel.background = element_rect(fill = 'lightblue'), 
                                                          legend.position = "none",
                                                           axis.title  = element_blank(),
                                                           #legend.title = element_blank(),
                                                           legend.background = element_rect(fill = "white", color = "black") 
  ) + ggtitle(paste0(disturb.species.order[i]))
}

mort.pred.legend <- get_legend(ggplot() + 
                                 geom_point(data = avg.plt.mort.annual  , aes(x = LONG_FIADB, y = LAT_FIADB, color = pMort_annual_plot_discrete), size = 3)+
                                 scale_color_manual(values =  c(
                                   "#74add1",
                                   "#abd9e9",
                                   
                                   "#fcbba1",
                                   "#fc9272",
                                   "#fb6a4a",
                                   "#ef3b2c",
                                   "#cb181d",
                                   "#99000d"), name = "Mortality Probability\n(%/year)")+
                                 
                                 
                                 coord_sf(xlim = c(-85, -67.5), ylim = c(37, 47.5))+theme(panel.grid = element_blank(), #panel.background = element_rect(fill = 'lightblue'), 
                                                                                          #legend.position = "none",
                                                                                          axis.title  = element_blank(),
                                                                                          #legend.title = element_blank(),
                                                                                          legend.background = element_rect(fill = "white", color = "black") 
                                 ) )
mort.pred.maps.species <- plot_grid(
  plot_grid(plotlist = species.plot.mortality.list, ncol = 3, align = "hv"), 
  mort.pred.legend, rel_widths = c(1, 0.2), ncol = 2
)

ggsave(paste0(output.dir, "images/map_annual_mortality_prob_plot.png"), plot = mort.pred.maps.species, 
       width = 12, height = 13, units = "in")

ggsave(paste0(output.dir, "images/map_annual_mortality_prob_plot.svg"), plot = mort.pred.maps.species, 
       width = 12, height = 13, units = "in")


# ideas for improving this figure:
# summarise by county (need to scale by volfac)
# alternative colors?

# Looking at four attibutes of tree mortality-----
# after Davis et al. hemlock pollen decline:
# Quaternary history and the stability of forest communities. Pages 132–153 in D. C. West, H. H. Shugart, and D. B. Botkin, eds. Forest succession: concepts and applications. Springer, New York.

# Specificity: Is it only one species or a group?--
# Synchrony: Is it happening at the same time across the range?--
# Rapidity: ??--
  # dia_diff as metric? lower diameter differences mean less growth prior to mortality?
  # may not be able to answer
# Uniformity: Is it happening across the range at similar rates?--

# Specificity: Are Species proability of mortality correlated with one another?
nmort.df <- spp.tree.mort.probs %>% group_by(state, PLOT.ID, Species, SPCD, LONG_FIADB, LAT_FIADB)%>%
  summarise(pMort_annual_plot = median(p1year.mortality, na.rm =TRUE), 
            pMort_annual_plot.sd = sd(p1year.mortality, na.rm =TRUE), 
            nTrees = n(), 
            Nmort = sum(mortality.draw))


# set up a dataframe with plots as rows and columns the pmort for different species
pmort.spread <- avg.plt.mort.annual %>% ungroup()%>% 
  mutate(pMort.pct  = pMort_annual_plot*100)%>%
  select(LAT_FIADB, LONG_FIADB, PLOT.ID, county,Species,pMort.pct)%>%
  group_by(LAT_FIADB, LONG_FIADB, county,PLOT.ID) %>% 
  spread(Species,pMort.pct,  fill = NA)
#p_correlations <- cor(pmort.spread[,4:ncol(pmort.spread)], use = "pairwise.complete")
# library(ggcorrplot)
# library(Hmisc)
# library(GGally)
# 
# 
# ggpairs(pmort.spread, columns = 4:ncol(pmort.spread))

library(ggcorrplot)
library(Hmisc)
p_correlations <- rcorr(as.matrix(pmort.spread[,4:ncol(pmort.spread)]), type = "spearman" ) # or "spearman"

# lets look at correlations for species that co-occur in at least 10 plots
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
melted.correlation <- get_lower_tri(p_correlations$r[disturb.species.order,disturb.species.order]) %>% reshape2::melt()%>%
  rename("r"="value") %>% left_join(., 
                                    get_lower_tri(p_correlations$n[disturb.species.order,disturb.species.order]) %>% reshape2::melt())%>%
  rename("n" = "value")%>%
  left_join(., 
            get_lower_tri(p_correlations$P[disturb.species.order,disturb.species.order]) %>% reshape2::melt())%>%
  rename("pval" = "value")%>%
  mutate(R_revised = ifelse(n>=50, round(r, digits = 1), NA), 
         P_revised = ifelse(n >=50, pval, NA))%>%
  mutate(R_revised = ifelse(R_revised == 1, NA, R_revised))%>%
  mutate(sig.label = ifelse(P_revised <= 0.05, R_revised, NA))%>%
  filter(!is.na(R_revised))# if the species are colocated in less than 10 plots omit the correlation


pmort_correlation <-  ggplot(data = melted.correlation, aes(Var2, Var1, fill = R_revised))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-0.6,0.6), space = "Lab", 
                       name="Spearman\nCorrelation") +
  theme_bw(base_size = 16)+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                    hjust = 0))+
  coord_fixed()+
  scale_x_discrete(position = "top",
                   limits = levels(melted.correlation$Var2))+
  scale_y_discrete(position = "left",
                   limits = rev(levels(melted.correlation$Var1)))+
 
  geom_text(aes(Var2, Var1, label = sig.label), color = "black", size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    #panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.9, 0.65),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))#+coord_flip()
  

ggsave(filename = paste0(output.dir, "images/posterior_plot_pMort_species_correlations.png"), 
       plot = pmort_correlation, 
       height = 6, width = 6, units = "in", 
       dpi = 350)

# fourteen species pairs have R_revised >= 0.3
melted.correlation %>% filter(R_revised >=0.3)
# northern region:
# balsam fir - Red spruce
# balsam fir - eastern hemlock

# eastern hemlock - yellow birch
# American beech - red maple

# American beech - red maple

# oaks:
# chestnut oak - N red oak
# chestnut oak - white oak

# chestnut oak - E. white pine
# chestnut oak - red maple

# N red oak - paper birch
# hickor spp. - northern red oak

# white oak - yellow birch 
# white oak - red maple

# red maple - paper birch
# sugar maple - red maple

# four species pairs have correlations >= 0.4
melted.correlation %>% filter(R_revised >=0.4)


# Synchrony and uniformity across space:
# moran's global 
# moran's local 
# Ripley's K-function:



data_built <- ggplot_build(p)

# The data for the filled contours is in the second layer
data_2d <- data_built$data[[2]]

# The head() function shows the structure
head(data_2d)

library(spdep)
polygon.contour.list <- polygon.buffer.list <- buffer.map.list <- moran.local.list <- polygon_hull.list <- list()
for(i in 1:17){
  
plt_mort_sf <- st_as_sf(avg.plt.mort.annual %>% filter(Species %in% nspp[i,]$COMMON), coords = c("LONG_FIADB", "LAT_FIADB"))
# get the spatial weights matrix 

plot_neighbors <- knn2nb(knearneigh(plt_mort_sf, k = 10)) # 10 nearest neighbors

# get spatial weights 
plot_weights <- nb2listw(plot_neighbors, style = "W")

# global Moran's I
moran_global <- moran.test(plt_mort_sf$pMort_annual_plot, plot_weights)
print(moran_global)

# Moran plot 
moran.plot(plt_mort_sf$pMort_annual_plot, plot_weights, #labels = as.character(tree_data$id),
           main = "Moran Plot for Tree Mortality Risk")


#Calculate Local Moran's I
local_moran_result <- localmoran(plt_mort_sf$pMort_annual_plot, plot_weights)

# Add the local Moran statistics to the sf object
plt_mort_sf$local_moran_I <- local_moran_result[, 1]
plt_mort_sf$p_value <- local_moran_result[, 5]

# get the type of local cluster
# standardization and  filter by significance
plt_mort_sf$cluster_type <- "Not Significant"
plt_mort_sf$standardized_risk <- scale(plt_mort_sf$pMort_annual_plot)
plt_mort_sf$quadrant <- local_moran_result[, 3] # Quadrant for the moran plot

# assign labels
significant_clusters <- which(plt_mort_sf$p_value <= 0.05)
plt_mort_sf$cluster_type[significant_clusters] <- ifelse(
  plt_mort_sf$standardized_risk[significant_clusters] > 0 & plt_mort_sf$quadrant[significant_clusters] > 0, "High-High",
  ifelse(plt_mort_sf$standardized_risk[significant_clusters] < 0 & plt_mort_sf$quadrant[significant_clusters] < 0, "Low-Low",
         ifelse(plt_mort_sf$standardized_risk[significant_clusters] > 0 & plt_mort_sf$quadrant[significant_clusters] < 0, "High-Low",
                ifelse(plt_mort_sf$standardized_risk[significant_clusters] < 0 & plt_mort_sf$quadrant[significant_clusters] > 0, "Low-High", "Not Significant"))))

# plot the  clusters
# ggplot(plt_mort_sf) +
#   geom_sf(aes(color = cluster_type), size = 3) +
#   scale_color_manual(values = c("High-High" = "red", "Low-Low" = "blue", "High-Low" = "pink", "Low-High" = "lightblue", "Not Significant" = "grey")) +
#   theme_minimal()
# 

# identify hotspots of mortality risk and create a polygon
high_high_clusters <- plt_mort_sf %>%
  filter(p_value < 0.05 & local_moran_I > 0)# & pMort_annual_plot >= 0.01)



# Convex hull around all HH cluster points
if (nrow(high_high_clusters) > 0) {
  hh_polygon_chull <- st_convex_hull(st_union(high_high_clusters))
} else {
  message("No HH clusters found to create a polygon.")
  hh_polygon_chull <- NULL
}
if (!is.null(hh_polygon_chull)) {
  ggplot() +
    geom_sf(data = plt_mort_sf, aes(color = pMort_annual_plot), size = 2) +
    geom_sf(data = high_high_clusters, color = "red", size = 3) +
    geom_sf(data = hh_polygon_chull, fill = "red", alpha = 0.3, color = NA) +
    labs(title = "HH Clusters and Convex Hull Polygon") +
    theme_minimal()
}

# Create a buffer around each High-High cluster point, then union
# You might want to adjust the buffer distance (dist)
if (nrow(high_high_clusters) > 0) {
  hh_polygon_buffer <- st_union(st_buffer(high_high_clusters, dist = 0.5))
} else {
  message("No HH clusters found to create a polygon.")
  hh_polygon_buffer <- NULL
}

gg.buffer.map <- ggplot() + 
  geom_polygon(data = canada, 
               aes(x=long, y=lat, group = group), 
               color = "black", fill = "white") +
  geom_polygon(data = lakes_df , 
               aes(x = long, y = lat, group = group), 
               color = "black", fill = "lightblue") +
  geom_polygon(data = state_sub, 
               aes(x=long, y=lat, group = group), 
               color = "black", fill = "white")+
    geom_sf(data = plt_mort_sf, aes(color = pMort_annual_plot_discrete), size = 2) +
    #geom_sf(data = high_high_clusters, color = "red", size = 3) +
    geom_sf(data = hh_polygon_buffer, fill = "red", alpha = 0.3, color = NA) +
    #labs(title = "HH Clusters and Convex Hull Polygon") +
    theme_minimal()+
  scale_color_manual(values =  c(
    "#74add1",
    "#abd9e9",
    
    "#fcbba1",
    "#fc9272",
    "#fb6a4a",
    "#ef3b2c",
    "#cb181d",
    "#99000d"), name = "Mortality Probability\n(%/year)")+
  
  
  coord_sf(xlim = c(-85, -67.5), ylim = c(37, 47.5))+theme(panel.grid = element_blank(), #panel.background = element_rect(fill = 'lightblue'), 
                                                           #legend.position = "none",
                                                           axis.title  = element_blank(),
                                                           #legend.title = element_blank(),
                                                           legend.background = element_rect(fill = "white", color = "black") 
  ) +ggtitle(paste0("mortality risk hotspots for ", nspp[i,]$COMMON))


# # use spatstat package to get KDE
library(spatstat)
spp.mort.rate <- avg.plt.mort.annual %>% filter(Species %in% nspp[i,]$COMMON)
# spp_mort_high_sf <- st_as_sf(spp.mort.rate %>% filter(pMort_annual_plot >= mean(spp.mort.rate$pMort_annual_plot)), coords = c("LONG_FIADB", "LAT_FIADB"), crs = 4269)%>%
#   st_transform(., crs = 5070)
# 
# # use state spatial polygons as the windown for the density object using as.owin 
# us_polygon_df <-map("usa", fill = TRUE, plot = FALSE) %>% st_as_sf() %>%
#   st_transform(., crs = 5070)
# us_state_owin <- as.owin(us_polygon_df)
# 
# ppp_obj_us <- as.ppp(st_coordinates(spp_mort_high_sf), W = us_state_owin)
# kde_result_us <- density(ppp_obj_us, sigma = 12)


# Convert the im object to a data frame for ggplot2
# kde_df <- as.data.frame(kde_result_us)

# Plot with filled contours
# ggplot(kde_df %>% filter(value >=0), aes(x = x, y = y, z = value)) +
#   geom_contour_filled() +
#   scale_fill_viridis_d() # Or other color scales
# 
# plot(kde_result_us)

# Define contour levels (e.g., at 25%, 50%, 75% of the maximum density)
# levels <- quantile(kde_result_us$v, probs = c(0.25, 0.50, 0.75), na.rm =TRUE)
# 
# # Extract contours
# contours_list <- contour(kde_result_us, levels = levels, draw = FALSE)


#mean(spp.mort.rate$pMort_annual_plot)
# if there are any points with > 1% mortality per year:
avg.plt.mort.annual %>% filter(pMort_annual_plot >= 0.005) %>% ungroup()%>% dplyr::select(Species) %>% distinct()

# get 2D KDE of above average mortality rates form stat_density2d
p <- ggplot(data = spp.mort.rate %>% filter(pMort_annual_plot >= quantile(spp.mort.rate$pMort_annual_plot, 0.75)), aes(x = LONG_FIADB, y = LAT_FIADB)) +
  geom_point(alpha = 0.5) + # Optional: show the original points
  coord_sf(xlim = c(-85, -60), ylim = c(30, 50))+
  
  stat_density2d(aes(fill = ..level..,), alpha = 0.5, geom = "polygon", contour_var = "count",#, contour_var = "count",
                 breaks =  c(0,0.5, 1, 2, 5, 10, 15, 20, 25, 30, 60)) +
  scale_fill_viridis_c() + # A color scale for the density levels
  #coord_equal() + # Ensures correct aspect ratio for spatial data
  labs(title = "Density of higher than average mortality",
       x = "Longitude",
       y = "Latitude") +
  theme_minimal()
p


p2 <- ggplot(data = spp.mort.rate %>% filter(pMort_annual_plot >= quantile(spp.mort.rate$pMort_annual_plot, 0.5)), aes(x = LONG_FIADB, y = LAT_FIADB)) +
  geom_point(alpha = 0.5) + # Optional: show the original points
  coord_sf(xlim = c(-85, -60), ylim = c(30, 50))+
  #geom_hdr(probs = c(0.25, 0.5, 0.75, 0.95))
  
   geom_density2d(aes(color = ..level..,), alpha = 0.5)+ 
  #                breaks = density_levels) +
  scale_fill_viridis_c() + # A color scale for the density levels
  #coord_equal() + # Ensures correct aspect ratio for spatial data
  labs(title = "Density of higher than average mortality",
       x = "Longitude",
       y = "Latitude") +
  theme_minimal()







data_built <- ggplot_build(p2)

# The data for the filled contours is in the second layer
data_2d <- data_built$data[[2]]

# set a polygon IDs for each contour piece
data_2d$pol_id <- paste0(data_2d$group)

# Split and create sf polygons based on the IDs
polygons_list <- lapply(unique(data_2d$pol_id), function(x) {
  # Subset data for the current polygon
  polygon_data <- data_2d[data_2d$pol_id == x, ]
  
  # close polygon
  closed_polygon <- rbind(polygon_data, polygon_data[1, ])
  
  # sf polygon 
  polygon_sf <- st_polygon(list(as.matrix(closed_polygon[, c("x", "y")])))
  
  # density level 
  level_data <- unique(polygon_data[, grepl("level", names(polygon_data))])
  
  # sf object with the level data
  st_as_sf(level_data, geometry = st_sfc(polygon_sf))
})

# combine species contours into a df
final_polygons_sf <- do.call(rbind, polygons_list) %>% 
  mutate(Species = nspp[i,]$COMMON, 
         mean.plt.mort.spp = mean(spp.mort.rate$pMort_annual_plot))

final_polygons_sf %>% filter(nlevel >=0.5) %>% ggplot()+geom_sf(aes(fill = nlevel))+scale_fill_viridis_c()+
  geom_point(data = spp.mort.rate %>% filter(pMort_annual_plot >= 0.01), aes(x = LONG_FIADB, y = LAT_FIADB), size = 0.1)




# save each species buffers:
polygon.contour.list[[i]] <- final_polygons_sf
polygon.buffer.list[[i]] <- hh_polygon_buffer
polygon_hull.list[[i]] <- hh_polygon_chull
buffer.map.list[[i]] <- gg.buffer.map
moran.local.list[[i]] <- plt_mort_sf

rm(hh_polygon_buffer, hh_polygon_chull, gg.buffer.map, final_polygons_sf)
}

names(polygon.contour.list) <- nspp$COMMON
names(polygon.buffer.list) <- nspp$COMMON
names(polygon_hull.list) <- nspp$COMMON
#spatial_buffer_points <- st_sfc(polygon.buffer.list[[1]], polygon.buffer.list[[2]], polygon.buffer.list[[3]])


poly.contours <- do.call(rbind, polygon.contour.list)%>% #data.frame()%>%
  group_by(Species) %>% mutate(max.level = max(level))

poly.contours %>% as.data.frame ()%>% group_by(Species)%>% summarise(max.level = max(level), min (level))%>% 
  arrange(max.level)


polygon.contour.list[[1]] %>% filter(nlevel > 0.5) %>% ggplot()+geom_sf(aes(fill = level))+scale_fill_viridis_c()
polygon.contour.list[[2]] %>% filter(nlevel > 0.5) %>% ggplot()+geom_sf(aes(fill = level))+scale_fill_viridis_c()
polygon.contour.list[[3]]  %>% filter(level >= 15) %>% ggplot()+geom_sf(aes(fill = level))+scale_fill_viridis_c()
polygon.contour.list[[4]]   %>% filter(level >= 15)%>% ggplot()+geom_sf(aes(fill = level))+scale_fill_viridis_c()
polygon.contour.list[[5]]  %>% filter(level >= 15)%>% ggplot()+geom_sf(aes(fill = level))+scale_fill_viridis_c()
polygon.contour.list[[6]]  %>% filter(level >= 15)%>% ggplot()+geom_sf(aes(fill = level))+scale_fill_viridis_c()
polygon.contour.list[[7]]   %>% filter(level >= 15)%>% ggplot()+geom_sf(aes(fill = level))+scale_fill_viridis_c()
polygon.contour.list[[8]]   %>% filter(level >= 15)%>% ggplot()+geom_sf(aes(fill = level))+scale_fill_viridis_c()
polygon.contour.list[[9]]   %>% filter(level >= 15)%>% ggplot()+geom_sf(aes(fill = level))+scale_fill_viridis_c()
polygon.contour.list[[10]] %>% filter(level >= 15)%>% ggplot()+geom_sf(aes(fill = level))+scale_fill_viridis_c()
polygon.contour.list[[11]]  %>% filter(level >= 15)%>% ggplot()+geom_sf(aes(fill = level))+scale_fill_viridis_c()
polygon.contour.list[[12]]   %>% filter(level >= 15)%>% ggplot()+geom_sf(aes(fill = level))+scale_fill_viridis_c()
polygon.contour.list[[13]]   %>% filter(level >= 15)%>% ggplot()+geom_sf(aes(fill = level))+scale_fill_viridis_c()
polygon.contour.list[[14]]  %>% filter(level >= 15) %>% ggplot()+geom_sf(aes(fill = level))+scale_fill_viridis_c()
polygon.contour.list[[15]]   %>% filter(level >= 15)%>% ggplot()+geom_sf(aes(fill = level))+scale_fill_viridis_c()
polygon.contour.list[[16]]   %>% filter(level >= 15)%>% ggplot()+geom_sf(aes(fill = level))+scale_fill_viridis_c()
polygon.contour.list[[17]]   %>% filter(level >= 15)%>% ggplot()+geom_sf(aes(alpha = level), fill = "forestgreen")#+scale_fill_viridis_c()


# plot the point-based buffers all together:
base.map <- ggplot() + 
  geom_polygon(data = canada, 
               aes(x=long, y=lat, group = group), 
               color = "black", fill = "white") +
  geom_polygon(data = lakes_df , 
               aes(x = long, y = lat, group = group), 
               color = "black", fill = "lightblue") +
  geom_polygon(data = state_sub, 
               aes(x=long, y=lat, group = group), 
               color = "black", fill = "white")+
  theme_minimal()

# plot kde based polygons for each species:
all.kde <- base.map +
  geom_sf(data = poly.contours %>% filter(nlevel >= 0.75) , aes(fill = Species, group = Species), alpha = 0.7) +
  species_fill+coord_sf(xlim = c(-85, -67.5), ylim = c(37, 47.5))+theme(panel.grid = element_blank(), #panel.background = element_rect(fill = 'lightblue'), 
                                                                       legend.position = "none",
                                                                       axis.title  = element_blank(),
                                                                       #legend.title = element_blank(),
                                                                       legend.background = element_rect(fill = "white", color = "black") )


TREE.remeas %>% dplyr::select(state, stname, date, remper) %>% distinct() %>%
  filter(!remper == 0) %>%
  group_by(state, stname, date)%>%
  summarise(avg.remper = median(remper))%>%
  mutate(date2 = date - avg.remper)

# most mortality here is 1981-1995 (1997)
northern.kde <- base.map +
  geom_sf(data = poly.contours %>% 
            filter(Species %in% c("balsam fir", "northern white-cedar",  "red spruce"))%>%
                     filter(nlevel > 0.75), aes(fill = Species, group = Species), alpha = 0.75, color = NA) +
  species_fill+coord_sf(xlim = c(-85, -67.5), ylim = c(37, 47.5))+theme(panel.grid = element_blank(), #panel.background = element_rect(fill = 'lightblue'), 
                                                                        legend.position = "none",
                                                                        axis.title  = element_blank(),
                                                                        #legend.title = element_blank(),
                                                                        legend.background = element_rect(fill = "white", color = "black") )



  


# peaks of mortality in PA for Oaks, in NY, NH, and VT, and ME for birches
# most mortality here is 1976-1989 for PA and WV and 1980-1997ish for birches
# Oaks:
# PA: remper 1976-1989 # peak SM defoliation = 1990
# NJ: remper 1972-1987 # peak SM defoliation = 1981
# Birches:
# NY: remper = 1980-1993 # peak SM defoliation = 1980
# VT: remper = 1982-1997 # peak SM defoliation = 1953, but some later on
# NH: remper = 1983 - 1997 # peak SM defoliation = 1981, and large peaks in 1980s


oak.kde <- base.map +
  geom_sf(data = poly.contours %>% 
            filter(Species %in% c( "chestnut oak","northern red oak", "white oak", "paper birch", "yellow birch"))%>%
            filter(nlevel >= 0.70), aes(fill = Species, group = Species), alpha = 0.75, color = NA) +
  species_fill+coord_sf(xlim = c(-85, -67.5), ylim = c(37, 47.5))+theme(panel.grid = element_blank(), #panel.background = element_rect(fill = 'lightblue'), 
                                                                        legend.position = "none",
                                                                        axis.title  = element_blank(),
                                                                        #legend.title = element_blank(),
                                                                        legend.background = element_rect(fill = "white", color = "black") )


  
  

# most mortality is in PA, NY, VT, NH, an Maine for both
# could be HWA related for hemlock in PA and NY
# PA = 1976 - 1989 (HWA first detected 1979)
# NY = 1980 - 1993 (HWA first detected 1984)
# VT = 1982 - 1997 (HWA first detected 2008)-likely Hemlock looper (Forest conditions report)
# NH = 1983 - 1997 (HWA first detected 2004)-Hemlock looper?

# Beech bark disease in vermont increased in n ares
# beech scale present in pennsyalvannia in forest conditions report 1990

beech.hemlock<- base.map +
  geom_sf(data = poly.contours %>% 
            filter(Species %in% c("eastern hemlock", "American beech"))%>%
            filter(nlevel > 0.70), aes(fill = Species, group = Species), alpha = 0.75, color = NA) +
  species_fill+coord_sf(xlim = c(-85, -67.5), ylim = c(37, 47.5))+theme(panel.grid = element_blank(), #panel.background = element_rect(fill = 'lightblue'), 
                                                                        legend.position = "none",
                                                                        axis.title  = element_blank(),
                                                                        #legend.title = element_blank(),
                                                                        legend.background = element_rect(fill = "white", color = "black") )





# hickory mortality is WV, OH = 
# maples and other hardwoods: 
spongy.resistant.kde <- base.map +
  geom_sf(data = poly.contours %>% 
            filter(Species %in% c("sugar maple", "red maple", "eastern white pine", "hickory spp."))%>%
            filter(nlevel >= 0.7), aes(fill = Species, group = Species), alpha = 0.75, color = NA) +
  species_fill+coord_sf(xlim = c(-85, -67.5), ylim = c(37, 47.5))+theme(panel.grid = element_blank(), #panel.background = element_rect(fill = 'lightblue'), 
                                                                        legend.position = "none",
                                                                        axis.title  = element_blank(),
                                                                        #legend.title = element_blank(),
                                                                        legend.background = element_rect(fill = "white", color = "black") )






spongy.immune.kde <- base.map +
  geom_sf(data = poly.contours %>% 
            filter(Species %in% c("black cherry", "yellow-poplar", "white ash"))%>%
            filter(nlevel >=0.7), aes(fill = Species, group = Species), alpha = 0.75, color = NA) +
  species_fill+coord_sf(xlim = c(-85, -67.5), ylim = c(37, 47.5))+theme(panel.grid = element_blank(), #panel.background = element_rect(fill = 'lightblue'), 
                                                                        legend.position = "none",
                                                                        axis.title  = element_blank(),
                                                                        #legend.title = element_blank(),
                                                                        legend.background = element_rect(fill = "white", color = "black") )








# all species point-based buffers together:
base.map +
  geom_sf(data = polygon.buffer.list[["balsam fir"]], aes(fill = "balsam fir"), alpha = 0.5, color = NA) +
  geom_sf(data = polygon.buffer.list[["red spruce"]], aes(fill = "red spruce"), alpha = 0.5, color = NA) +
  geom_sf(data = polygon.buffer.list[["northern white-cedar"]], aes(fill = "northern white-cedar"), alpha = 0.5, color = NA) +
  
  geom_sf(data = polygon.buffer.list[["eastern hemlock"]], aes(fill = "eastern hemlock"), alpha = 0.5, color = NA) +
  geom_sf(data = polygon.buffer.list[["American beech"]], aes(fill = "American beech"), alpha = 0.5,color = NA) +
  
  
  geom_sf(data = polygon.buffer.list[["northern red oak"]], aes(fill = "northern red oak"), alpha = 0.5, color = NA) +
  geom_sf(data = polygon.buffer.list[["chestnut oak"]], aes(fill = "chestnut oak"), alpha = 0.5,color = NA) +
  geom_sf(data = polygon.buffer.list[["white oak"]], aes(fill = "white oak"), alpha = 0.5, color = NA) +
  geom_sf(data = polygon.buffer.list[["yellow birch"]], aes(fill = "yellow birch"), alpha = 0.5,color = NA) +
  geom_sf(data = polygon.buffer.list[["paper birch"]], aes(fill = "paper birch"), alpha = 0.5,color = NA) +
  geom_sf(data = polygon.buffer.list[["hickory spp."]], aes(fill = "hickory spp."), alpha = 0.5,color = NA) +
  geom_sf(data = polygon.buffer.list[["eastern white pine"]], aes(fill = "eastern white pine"), alpha = 0.5,color = NA) +
  geom_sf(data = polygon.buffer.list[["sugar maple"]], aes(fill = "sugar maple"), alpha = 0.5,color = NA) +
  geom_sf(data = polygon.buffer.list[["red maple"]], aes(fill = "red maple"), alpha = 0.5,color = NA) +
  geom_sf(data = polygon.buffer.list[["black cherry"]], aes(fill = "black cherry"), alpha = 0.5,color = NA) +
  geom_sf(data = polygon.buffer.list[["yellow-poplar"]], aes(fill = "yellow-poplar"), alpha = 0.5,color = NA) +
  geom_sf(data = polygon.buffer.list[["white ash"]], aes(fill = "white ash"), alpha = 0.5,color = NA) +
  species_fill+coord_sf(xlim = c(-85, -67.5), ylim = c(37, 47.5))+theme(panel.grid = element_blank(), #panel.background = element_rect(fill = 'lightblue'), 
                                                                        #legend.position = "none",
                                                                        axis.title  = element_blank(),
                                                                        #legend.title = element_blank(),
                                                                        legend.background = element_rect(fill = "white", color = "black") 
  )


spruce.fir.hotspots <- base.map +
  geom_sf(data = polygon.buffer.list[["balsam fir"]], aes(fill = "balsam fir"), alpha = 0.9, color = NA) +
  geom_sf(data = polygon.buffer.list[["red spruce"]], aes(fill = "red spruce"), alpha = 0.9, color = NA) +
  geom_sf(data = polygon.buffer.list[["northern white-cedar"]], aes(fill = "northern white-cedar"), alpha = 0.5, color = NA) +
  
  #geom_sf(data = polygon.buffer.list[["eastern hemlock"]], aes(fill = "eastern hemlock"), alpha = 0.5, color = NA) +
  #geom_sf(data = polygon.buffer.list[["American beech"]], aes(fill = "American beech"), alpha = 0.5,color = NA) +
  
  
 
  species_fill+coord_sf(xlim = c(-85, -67.5), ylim = c(37, 47.5))+theme(panel.grid = element_blank(), #panel.background = element_rect(fill = 'lightblue'), 
                                                                        #legend.position = "none",
                                                                        axis.title  = element_blank(),
                                                                        #legend.title = element_blank(),
                                                                        legend.background = element_rect(fill = "white", color = "black") 
  )

oak.hotspots <- base.map +
  geom_sf(data = polygon.buffer.list[["northern red oak"]], aes(fill = "northern red oak"), alpha = 0.5, color = NA) +
  geom_sf(data = polygon.buffer.list[["chestnut oak"]], aes(fill = "chestnut oak"), alpha = 0.5,color = NA) +
  geom_sf(data = polygon.buffer.list[["white oak"]], aes(fill = "white oak"), alpha = 0.5, color = NA) +
  geom_sf(data = polygon.buffer.list[["yellow birch"]], aes(fill = "yellow birch"), alpha = 0.5,color = NA) +
  geom_sf(data = polygon.buffer.list[["paper birch"]], aes(fill = "paper birch"), alpha = 0.5,color = NA) +
  
  species_fill+coord_sf(xlim = c(-85, -67.5), ylim = c(37, 47.5))+theme(panel.grid = element_blank(), #panel.background = element_rect(fill = 'lightblue'), 
                                                                        #legend.position = "none",
                                                                        axis.title  = element_blank(),
                                                                        #legend.title = element_blank(),
                                                                        legend.background = element_rect(fill = "white", color = "black") 
  )



beech.hemlock.hotspots <- base.map +
  geom_sf(data = polygon.buffer.list[["eastern hemlock"]], aes(fill = "eastern hemlock"), alpha = 0.75, color = "black") +
  geom_sf(data = polygon.buffer.list[["American beech"]], aes(fill = "American beech"), alpha = 0.75,color = "black") +
  geom_sf(data = polygon.buffer.list[["eastern white pine"]], aes(fill = "eastern white pine"), alpha = 0.75,color = "black") +
  species_fill+coord_sf(xlim = c(-85, -67.5), ylim = c(37, 47.5))+theme(panel.grid = element_blank(), #panel.background = element_rect(fill = 'lightblue'), 
                                                                        #legend.position = "none",
                                                                        axis.title  = element_blank(),
                                                                        #legend.title = element_blank(),
                                                                        legend.background = element_rect(fill = "white", color = "black") 
  )

beech.hemlock.hotspots



maple.hotspots <- base.map +
  geom_sf(data = polygon.buffer.list[["sugar maple"]], aes(fill = "sugar maple"), alpha = 0.5, color = "black") +
  geom_sf(data = polygon.buffer.list[["red maple"]], aes(fill = "red maple"), alpha = 0.5,color = "black") +
 # geom_sf(data = polygon.buffer.list[["eastern white pine"]], aes(fill = "eastern white pine"), alpha = 0.5,color = "black") +
  species_fill+coord_sf(xlim = c(-85, -67.5), ylim = c(37, 47.5))+theme(panel.grid = element_blank(), #panel.background = element_rect(fill = 'lightblue'), 
                                                                        #legend.position = "none",
                                                                        axis.title  = element_blank(),
                                                                        #legend.title = element_blank(),
                                                                        legend.background = element_rect(fill = "white", color = "black") 
  )

maple.hotspots

other.hotspots <- base.map +
  geom_sf(data = polygon.buffer.list[["black cherry"]], aes(fill = "black cherry"), alpha = 0.7,color = "black") +
  geom_sf(data = polygon.buffer.list[["yellow-poplar"]], aes(fill = "yellow-poplar"), alpha = 0.7,color = "black") +
  geom_sf(data = polygon.buffer.list[["white ash"]], aes(fill = "white ash"), alpha = 0.7,color = "black") +
  geom_sf(data = polygon.buffer.list[["hickory spp."]], aes(fill = "hickory spp."), alpha = 0.7,color = "black") +
  # geom_sf(data = polygon.buffer.list[["eastern white pine"]], aes(fill = "eastern white pine"), alpha = 0.5,color = "black") +
  species_fill+coord_sf(xlim = c(-85, -67.5), ylim = c(37, 47.5))+theme(panel.grid = element_blank(), #panel.background = element_rect(fill = 'lightblue'), 
                                                                        #legend.position = "none",
                                                                        axis.title  = element_blank(),
                                                                        #legend.title = element_blank(),
                                                                        legend.background = element_rect(fill = "white", color = "black") 
  )

other.hotspots

spruce.fir.hothulls <- base.map +
  geom_sf(data = polygon_hull.list[["balsam fir"]], aes(fill = "balsam fir"), alpha = 0.9, color = NA) +
  geom_sf(data = polygon_hull.list[["red spruce"]], aes(fill = "red spruce"), alpha = 0.9, color = NA) +
  geom_sf(data = polygon_hull.list[["northern white-cedar"]], aes(fill = "northern white-cedar"), alpha = 0.5, color = NA) +
  
  #geom_sf(data = polygon.buffer.list[["eastern hemlock"]], aes(fill = "eastern hemlock"), alpha = 0.5, color = NA) +
  #geom_sf(data = polygon.buffer.list[["American beech"]], aes(fill = "American beech"), alpha = 0.5,color = NA) +
  
  
  
  species_fill+coord_sf(xlim = c(-85, -67.5), ylim = c(37, 47.5))+theme(panel.grid = element_blank(), #panel.background = element_rect(fill = 'lightblue'), 
                                                                        #legend.position = "none",
                                                                        axis.title  = element_blank(),
                                                                        #legend.title = element_blank(),
                                                                        legend.background = element_rect(fill = "white", color = "black") 
  )



# geom_boxplot()# to do this, we need to handle some of the predicted trees have a volfac == 0
all.remeas %>% select(state, volfac)%>% distinct()%>% group_by(state) %>% summarise(nunique_volfac = n())

# states 9, 23, 24, 33, and 50 may have fewer plot designs and volfacs are either zero, sawtimber, poletimber, etc
all.remeas %>% select(state, volfac, dbhcur) %>% filter(state %in% c(9, 23, 24, 33, 50))%>% distinct()%>% 
  group_by(state, volfac) %>% 
  summarise(min_dbh = min(dbhcur))

all.remeas %>% filter(dbhcur>5)%>% select(state, PLOT.ID, point, date, cndtn) %>% 
  distinct()%>% group_by(state, PLOT.ID, date, cndtn)%>% summarise(npoints = n())%>% 
  ungroup()%>% group_by(state, date, cndtn) %>% summarise(max_npoints = max(npoints), min_npoints = min(npoints))

all.remeas %>% select(state, volfac, dbhcur) %>% filter(!state %in% c(9, 23, 24, 33, 50))%>% distinct()%>% 
  group_by(state, volfac) %>% 
  summarise(min_dbh = min(dbhcur), 
            dbhcur/volfac)
# formula for TPA
#TPA = (BAF / 0.005454*DIA^2)/Npoints
#volfac*Npoints = BAF/0.005454*DIA^2
(volfac*Npoints)*(0.005454*DIA^2)

# the trees with volfac == 0 are in mostly fixed radius plot states (9, 23, 33, 50), except for 50
all.remeas %>% filter(volfac == 0 & status %in% c(2,4,5,6))%>% group_by(state, mortfac > 0) %>% summarise(n())

# if status == 1 and volfac == 0, assume it is not counted

BAFestimated <- left_join(all.remeas, all.remeas %>% select(state, PLOT.ID, point) %>% distinct()%>%
  group_by(state, PLOT.ID) %>% summarise(Npoints = n())) %>%
  ungroup()%>%
  #filter(Npoints >1) %>% 
  mutate(BAFest = (volfac*Npoints)*(0.005454*dbhcur^2))                                           

# state 36, new york == 10 point design, 37.5 BAF?
# state 39, ohio == 5 point design, 18.7 BAF?
# state 42, Pennsylvannia == 5 point design, 18.7 BAF?
# state 54, West Virginia == 10 point design, 37.5BAF or 5 point design, 18.7 BAF?

ggplot(BAFestimated, aes( y = BAFest, fill = state, group = state))+geom_histogram()+facet_wrap(~state)+ylim(0,50)
ggplot(BAFestimated, aes( y = Npoints, fill = state, group = state))+geom_histogram()+facet_wrap(~state)

ggplot(BAFestimated %>% filter(state == 54), aes( y = BAFest, fill = state, group = state))+geom_histogram()+facet_wrap(~state)
ggplot(BAFestimated %>% filter(state == 54), aes( y = BAFest, fill = state, group = state))+geom_histogram()+facet_wrap(~Npoints >5)

# state 34, new jersey is a mix of fixed radius and variable radius designs
ggplot(BAFestimated %>% filter(state == 34), aes( y = BAFest, fill = state, group = state))+geom_histogram()+facet_wrap(~state)
ggplot(BAFestimated %>% filter(state == 34), aes( y = BAFest, fill = state, group = state))+geom_histogram()+facet_wrap(~Npoints >5)
ggplot(BAFestimated %>% filter(state == 34), aes( y = volfac, fill = state, group = state))+geom_histogram()+facet_wrap(~Npoints >1)

BAFestimated %>% filter(volfac == 0)%>% group_by(state) %>% summarise(n())
# fixed radius volfacs--
# volfac == 10.0196, for trees < 11.0 in
# volfac == 4.9928, for trees >= 11 inches
# all of these are single point plots:
BAFestimated %>% filter(state == 34) %>% filter(volfac == 10.0196 | volfac == 4.9928) %>% 
  group_by(volfac, point)%>% summarise(maxdbh = max(dbhcur), 
                                mindbh = min(dbhcur))

BAFestimated %>% filter(state == 34) %>% filter(!volfac == 10.0196 & !volfac == 4.9928) #%>% 
  
ggplot(BAFestimated %>% filter(state == 34 & !volfac == 10.0196 & !volfac == 4.9928), aes( y = BAFest, fill = state, group = state))+geom_histogram()+facet_wrap(~Npoints >5)
ggplot(BAFestimated %>% filter(state == 34 & !volfac == 10.0196 & !volfac == 4.9928), aes( y = Npoints, fill = state, group = state))+geom_histogram()+facet_wrap(~Npoints >5)


BAFestimated %>% filter(Npoints > 10)%>% View()

# states 34, 36, 39, 42, and 54 likely have variable radius plot designs
all.remeas %>% group_by(state,volfac == 0) %>% summarise(n = n())

# total trees in each species on the stand and in the plots

  
group_by(SPCD, Species, state, remper, PLOT.ID, LONG_FIADB, LAT_FIADB) %>%
  
  summarise(n_trees_volfac = sum(volfac, na.rm = TRUE), 
            n_count = n())




total.spp.predicted.plt <- spp.tree.pred.mort.volfac %>% rename("Species" = "COMMON")%>% 
  group_by(SPCD, Species, state, PLOT.ID, LAT_FIADB, LONG_FIADB)%>%
  summarise(n_predicted_plot = n(), 
            volfac_predicted_plot = sum(volfac, na.rm =TRUE))

spp.predicted.dead.plot <- spp.tree.pred.mort.volfac %>% rename("Species" = "COMMON")%>% 
  mutate(predicted.status = ifelse(survival.draw == 1, "1", "2"))%>%
  filter(predicted.status == "2")%>%
  group_by(SPCD, Species, state, PLOT.ID, LAT_FIADB, LONG_FIADB)%>%
  
  summarise(predicted_dead_trees_volfac = sum(volfac, na.rm = TRUE), 
            n_predicted_dead = n())%>%
  left_join(., total.spp.predicted.plt )



# get the plot-level predicted mortality rates by each species
predicted.mort.rates.spp.plot <- spp.predicted.dead.plot %>% 
  group_by(SPCD, Species, state, PLOT.ID, LAT_FIADB, LONG_FIADB)%>%
  # for each remper period, get the mortality rate per year 
  summarise(predicted_10yrmort_rate_volfac = ((predicted_dead_trees_volfac/volfac_predicted_plot)), 
            predicted_10yrmort_rate = ((n_predicted_dead/n_predicted_plot)))%>%
  ungroup()%>%
  # 10 year survival rates
  mutate(predicted_10yr_surv_rate_volfac = 1-predicted_10yrmort_rate_volfac, 
         predicted_10yr_surv_rate = 1-predicted_10yrmort_rate)%>%
  # get annualized predicted survival and mortality rates
  mutate(predicted_1year_surv_rate_volfac = predicted_10yr_surv_rate_volfac^(1/10),
         predicted_1year_surv_rate = predicted_10yr_surv_rate^(1/10))%>%
  mutate(predicted_1year_mort_rate_volfac = 1 - predicted_1year_surv_rate_volfac, 
         predicted_1year_mort_rate = 1 - predicted_1year_surv_rate)

hist(predicted.mort.rates.spp.plot$predicted_1year_mort_rate_volfac)

  
# plot up 10 year survival rates by plot, by species:
ggplot(data = predicted.mort.rates.spp.plot %>% filter(Species %in% c("balsam fir", "red spruce", "northern white-cedar")))+
  geom_point(aes(x =LONG_FIADB,  y = LAT_FIADB, color = predicted_10yrmort_rate*100))+
  facet_wrap(~Species)+ scale_color_viridis_c(option = "magma")


ggplot(data = predicted.mort.rates.spp.plot %>% filter(Species %in% c("balsam fir", "red spruce", "northern white-cedar")))+
  geom_density(aes(y = predicted_1year_surv_rate*100, fill = Species))+
  facet_wrap(~Species)#+ylab("predicted_1year_mort_rate_volfac")

ggplot(data = predicted.mort.rates.spp.plot %>% filter(Species %in% c("balsam fir", "red spruce", "northern white-cedar")))+
  geom_density(aes(y = predicted_10yrmort_rate*100, fill = Species))+
  facet_wrap(~Species)+ylab("Predicted 10-year mortality rates")


ggplot(data = predicted.mort.rates.spp.plot)+
  geom_density(aes(y = predicted_10yrmort_rate*100, fill = Species))+
  facet_wrap(~Species)+ylab("Predicted 10-year mortality rates")+
  species_fill




ggplot(data = predicted.mort.rates.spp.plot)+
  geom_density(aes(y = predicted_10yrmort_rate_volfac*100, fill = Species))+
  facet_wrap(~Species)+ylab("Predicted 10-year mortality rates (% of stems/species/plot)")+
  species_fill

predicted.mort.rates.spp.plot %>% filter(Species %in% c("balsam fir", "red spruce", "northern white-cedar"))%>%
  ggplot(aes(x = Species, y = predicted_10yrmort_rate_volfac, fill = Species)) +
  geom_dots(layout = "weave") +
  stat_slabinterval()



ggplot(data = predicted.mort.rates.spp.plot %>% filter(Species %in% c("balsam fir", "red spruce", "northern white-cedar")))+
  geom_point(aes(x =LONG_FIADB,  y = LAT_FIADB, color = predicted_10yrmort_rate*100))+
  facet_wrap(~Species)+ scale_color_viridis_c(option = "magma")

# get state-level mortality rate estimates for the modelled estimates---
# sum up volfac by state
# sum up number of plots by state

st.mort.rate.species <- all.remeas %>% group_by(SPCD, Species, state, remper)%>% filter(status %in% c(2, 4, 5, 6)) %>%
  summarise(nonlog_dead_trees_volfac = sum(volfac, na.rm = TRUE), 
            n_nonlog_dead = n()) %>%
  left_join(., st.total.tree.sums) %>%# join to total trees
  
  # for each remper period, get the mortality rate per year 
  mutate(nonlog_mort_rate_volfac = ((nonlog_dead_trees_volfac/total_trees_volfac)*100)/remper, 
         nonlog_mort_rate = ((n_nonlog_dead/total_n)*100)/remper)%>%
  ungroup()%>%
  
  # average all the remper mortality rates together to get a single species value
  group_by(SPCD, Species, state)%>%
  summarise(species_volfac_mort = mean(nonlog_mort_rate_volfac, na.rm =TRUE), 
            species_n_mort = mean(nonlog_mort_rate, na.rm =TRUE))


# get the total trees predicted for each species: 
# predictions were made on standard 10 year intervals
total.spp.predicted <- spp.tree.pred.mort.volfac %>% rename("Species" = "COMMON")%>% 
  group_by(SPCD, Species, state)%>%
  summarise(n_predicted_total = n(), 
            volfac_predicted_total = sum(volfac, na.rm =TRUE))

spp.predicted.dead <- spp.tree.pred.mort.volfac %>% rename("Species" = "COMMON")%>% 
  mutate(predicted.status = ifelse(survival.draw == 1, "1", "2"))%>%
  filter(predicted.status == "2")%>%
  group_by(SPCD, Species, state)%>%
  
  summarise(predicted_dead_trees_volfac = sum(volfac, na.rm = TRUE), 
            n_predicted_dead = n())%>%
  left_join(., total.spp.predicted )
  

predicted.mort.rates.spp.state <- spp.predicted.dead %>% group_by(SPCD, Species, state)%>%
  # for each remper period, get the mortality rate per year 
  summarise(predicted_mort_rate_volfac = ((predicted_dead_trees_volfac/volfac_predicted_total)*100)/10, 
         predicted_mort_rate = ((n_predicted_dead/n_predicted_total)*100)/10)%>%
  ungroup()


# combine with the observed mortality rates by species:
p.o.mort.rates <- left_join(st.mort.rate.species, predicted.mort.rates.spp.state)

ggplot(data = p.o.mort.rates, aes(x = species_volfac_mort, y = predicted_mort_rate_volfac, color = Species))+
  geom_point()+
  species_color+geom_abline(aes(slope = 1, intercept = 0), linetype = "dashed")+
  xlim(0, 6)+ylim(0,6)

# looking at distibution of mortality rates within the context of disturbances----

PLOT <- read_delim(paste0(output.dir,"data/formatted_older_matching_plts_PLOT.txt"))
colnames(PLOT)
unique(PLOT$typcur)
TREE.typcd <- TREE.remeas %>% left_join(., PLOT %>% select(PLOT.ID, cndtn, typcur, typold, date) %>% distinct())

# northern white cedar mortality is in N white cedar forest type (14), and spruce fir types
# red spruce and balsam fir (13), balsam fir(11), and red spruce (19)

TREE.typcd %>% filter(Species %in% "northern white-cedar" & dbhold > 5)%>%
  filter(status == c(1, 2, 4, 5)) %>% 
  mutate(M = ifelse(status == 1, 0, 1))%>%
  group_by( typcur, M)%>%
  summarise(ntree = n()) %>% 
  spread(M, ntree)%>% arrange(desc(`1`)) 

# get plots with at least 1 dead nwc tree:
Trees.maine <- TREE.typcd %>% filter(dbhold > 5)%>%
  filter(status %in% c(1, 2,4, 5)) %>% 
  mutate(M = ifelse(status == 1, 0, 1))%>%
  group_by( PLOT.ID, stname, Species, typcur, M)%>%
  summarise(ntree = n()) %>% 
  ungroup() %>% filter(stname %in% c("ME", "VT", "NH"))

Dead.NWC.plots <- Trees.maine %>% filter(Species %in% "northern white-cedar" ) %>% 
  group_by(PLOT.ID, M)%>%
  summarise(totals = sum(ntree)) %>% 
  spread(M, totals) %>% filter(!is.na(`1`))

noDead.NWC.plots <- Trees.maine %>% filter(Species %in% "northern white-cedar"  )%>% 
  group_by(PLOT.ID, M)%>%
  summarise(totals = sum(ntree)) %>% 
  spread(M, totals) %>% filter(is.na(`1`) | `1`==0)

# get the mortality fraction for other species for plots with dead N. white cedar
Trees.maine %>% filter(PLOT.ID %in% Dead.NWC.plots$PLOT.ID) %>% 
  mutate(plot.type = "Dead white-cedar") %>% rbind(., 
                                                   Trees.maine %>% filter(PLOT.ID %in% noDead.NWC.plots$PLOT.ID) %>% 
                                                     mutate(plot.type = "No Dead white-cedar")) %>% 
  filter(M == 1)%>%
  ggplot()+geom_bar(aes(x = plot.type, y = ntree, fill = Species), position = "stack", stat = "identity")+facet_wrap(~typcur)

TREE.typcd %>% filter(Species %in% "northern white-cedar" & dbhold > 5 & !status == 3) %>% 
  mutate(M = ifelse(status == 1, "live", "dead"))%>%
  ggplot()+geom_boxplot(aes(x = as.factor(typcur), y = dbhcur,  fill = M), alpha = 0.75, position = "dodge")

TREE.typcd %>% filter(Species %in% "northern white-cedar" & dbhold > 5 & !status == 3) %>% 
  mutate(M = ifelse(status == 1, "live", "dead"))%>%
  group_by(M, damage)%>%
  summarise(n())

TREE.typcd %>% filter(Species %in% "northern white-cedar" & dbhold > 5 & !status == 3)%>%
  #filter(damage == 50)%>%
  group_by(PLOT.ID) %>% 
  ggplot()+geom_jitter(aes(x = LONG_FIADB,y = LAT_FIADB, color = as.factor(damage)))


TREE.typcd %>% filter(Species %in% "balsam fir")%>%
  filter(status %in% c(1, 2, 4, 5)) %>% 
  mutate(M = ifelse(status == 1, 0, 1))%>%
  group_by( typcur, M)%>%
  summarise(ntree = n()) %>% 
  spread(M, ntree)%>% arrange(desc(`1`)) 



TREE.typcd %>%filter(Species %in% "balsam fir") %>%
  filter(status %in% c(1, 2, 4, 5)) %>% 
  mutate(M = ifelse(status == 1, 0, 1))%>%
  group_by( typcur, stname, M)%>%
  summarise(ntree = n()) %>% 
  spread(M, ntree)%>% arrange(desc(`1`)) 


# most of the Northern white cedar tree mortality is in Maine
TREE.typcd %>% filter(Species %in% "northern white-cedar" & dbhold > 5)%>%
  filter(status == c(1, 2,4, 5)) %>% 
  mutate(M = ifelse(status == 1, 0, 1))%>%
  group_by( stname, M)%>%
  summarise(ntree = n()) %>% 
  spread(M, ntree)%>% arrange(desc(`1`)) 


TREE.typcd %>% filter(Species %in% "northern white-cedar"& dbhold > 5)%>%
  filter(status == c(1, 2, 4, 5)) %>% 
  mutate(M = ifelse(status == 1, 0, 1))%>%
  group_by(PLOT.ID, M)%>%
  summarise(ntree = n()) %>% 
  spread(M, ntree) %>%
  mutate(total.NWC = sum(c(`0`,`1`), na.rm =TRUE)) %>% 
  arrange(desc(`1`))




# get summary of damages for later use:
N.DAMAGE <- cleaned.data %>% group_by(SPCD, damage) %>% summarise(n.by.damage = n())
N.DAMAGE$SPECIES <- ref_species[match(N.DAMAGE$SPCD, ref_species$SPCD),]$COMMON
ref_damage<- ref_codes %>% filter(VARIABLE %in% "AGENTCD")
N.DAMAGE$damage_agent <- ref_damage[match(N.DAMAGE$damage, ref_damage$VALUE),]$MEANING
N.DAMAGE$damage_agent <- ifelse(N.DAMAGE$damage == 0, "None", N.DAMAGE$damage_agent)





maine_var <- readRDS(paste0(output.folder, "/SPCD_stanoutput_joint_v3/predicted_mort/","variance_explained_state_23.RDS"))






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

# centralized code for the manuscript figures:
input.folder <- "C:/Users/KellyHeilman/Box/01. kelly.heilman Workspace/mortality/Eastern-Mortality/mortality_models/SPCD_stanoutput_joint_v3/"
output.folder <- "C:/Users/KellyHeilman/Box/01. kelly.heilman Workspace/mortality/Eastern-Mortality/mortality_models/SPCD_stanoutput_joint_v3/"

# re-invisioning the betas effects plots ----
# get all the covariates using posterior package
model.no <- 6
mod.data <- readRDS(paste0(input.folder, "all_SPCD_model_",model.no,".RDS"))
ncovar <- length(colnames(mod.data$xM))



full.model <- data.frame(Covariates = colnames(mod.data$xM), 
                         id = 1:length(colnames(mod.data$xM)))

# read in the betas, marginal responses, and variance parsing summaries ---
betas.df <- readRDS(paste0(input.folder, "samples/u_betas_model_", model.no, "_5000samples.rds"))
marginal_response_df <- readRDS(paste0(input.folder, "all_marginal_responses.RDS"))
interaction_response_df <- readRDS(paste0(input.folder, "interaction_responses.RDS"))
var_summary <- read.csv(paste0(output.folder, "predicted_mort/regional_var/variance_partitioning_summary_by_predictor_all.csv"))
var_summary <- readRDS(paste0(output.folder, "predicted_mort/var_summary_regional.rds"))


betas.quant <- betas.df %>% summarise_draws(median, ~quantile(., probs = c(0.025, 0.975))) %>%
  rename(`ci.lo` = "2.5%", `ci.hi` = "97.5%") %>%
  mutate(remper.cor = 0.5)
# relabel u_betas to meaningful species ids names
betas.quant$spp <- rep(1:17, ncovar)
betas.quant$cov <- rep(1:ncovar, each = 17)


covariate_names <- c(colnames(mod.data$xM))  # Replace with your covariate names
betas.quant$Covariate <- rep(covariate_names, each = 17)
betas.quant$Species <- rep(nspp[1:17,]$COMMON, ncovar)


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
#ggsave(height = 10, width = 15,dpi = 350, units = "in",paste0(output.folder,"SPCD_stanoutput_joint_v3/images/Estimated_effects_on_mortality_model_model_",model.no,"_all_species_betas.png"))



betas.quant$Species <- factor(betas.quant$Species, levels = rev(disturb.species.order))

betas.sig.df <- betas.quant %>% filter(significance %in% "significant") %>%
  left_join(., unique(var_summary[,c("Covariate", "Predictor","predictor.class2")]))

plot_beta_effects <- function(pred.group, betas.sig){
  
  df.covar <- betas.sig %>% filter(predictor.class2 %in%  pred.group)
  strip.fill <- as.character(color.pred.class.2[unique(df.covar$predictor.class2)])
  
  ggplot()+
    geom_errorbar(data = na.omit(df.covar), aes(y = Predictor , xmin = ci.lo, xmax = ci.hi, color = Species, linetype = significance), position = position_dodge(width = 0.9), width = 0)+
    geom_point(data = na.omit(df.covar), aes(y = Predictor, x = median, color = Species, shape = significance), position = position_dodge(width = 0.9), size = 1.5)+
    
    
    facet_wrap(~predictor.class2, scales= "free_y", ncol = 1)+
    geom_vline(xintercept = 0, color = "grey", linetype = "dashed")+
    theme_bw(base_size = 14)+
    theme(axis.text.y = element_text(angle = 90, hjust = 0.5),
          panel.grid  = element_blank(), 
          strip.background = element_rect(fill = c(strip.fill)), 
          strip.text = element_text(color = "white"), 
          axis.title.y = element_blank()
    )+
    xlab("Effect on Survival")+ylab("Parameter")+
    species_color+
    scale_shape_manual(values = c("significant" = 19,  "not overlapping zero" = 1), drop = FALSE, name = "Significance")+
    scale_linetype_manual(values = c("significant" = "solid",  "not overlapping zero" = "dotted"), 
                          drop = FALSE,  name = "Significance")+
    guides(color = guide_legend(reverse = TRUE))
  
  
}


plot_beta_effects(pred.group = "Change in Size", betas.sig = betas.sig.df)
plot_beta_effects(pred.group = "Competition", betas.sig = betas.sig.df)

plot_beta_effects(pred.group = "Competition", betas.sig = betas.sig.df)

beta.sig.plots <- lapply(unique(betas.sig.df$predictor.class2), function(x){
  plot_beta_effects(pred.group = x, betas.sig = betas.sig.df)
})

beta.sig.plots


# plot for all the betas, not just significant ones:
betas.all <- betas.quant %>% 
  left_join(., unique(var_summary[,c("Covariate", "Predictor","predictor.class2")]))

beta.all.plots <- lapply(unique(betas.all$predictor.class2), function(x){
  plot_beta_effects(pred.group = x, betas.sig = betas.all)
})


beta.all.plots[[2]]

# plots for the betas that have at least one significant species effect:
betas.sig.all <- betas.quant %>% filter(Covariate %in% unique(betas.sig.df$Covariate))%>% 
  left_join(., unique(var_summary[,c("Covariate", "Predictor","predictor.class2")]))

# group by variables
# list.vars <- list(c("Change in Size", "Size", "Growth x Size"),
# c("% Damage", "Competition", "Competition x G & S"),
# c("Climate", "Climate x G & S"),
# c("N deposition", "Ndep x G & S", "Site Conditions"))

beta.sig.all.plots <- lapply(unique(betas.sig.all$predictor.class2), function(x){
  plot_beta_effects(pred.group = x, betas.sig = betas.sig.all)
})
library(patchwork)
# set up the layouts using patchwork
growth.plt <- (beta.sig.all.plots[[1]] / beta.sig.all.plots[[2]]/ beta.sig.all.plots[[8]] & coord_cartesian(xlim=c(-1, 3.567)))+ plot_layout(ncol = 1, guides = "collect", axes = "collect_x") 

compete.plt <- ( (beta.sig.all.plots[[3]]) / beta.sig.all.plots[[12]]/beta.sig.all.plots[[4]] /beta.sig.all.plots[[9]]  & coord_cartesian(xlim=c(-0.7, 0.5))) + 
  plot_layout(ncol = 1, guides = "collect", axes = "collect_x", heights = c(2, 1,1,1)) 

climate.plt <- ( beta.sig.all.plots[[5]] / beta.sig.all.plots[[10]] & coord_cartesian(xlim=c(-1, 1))) + plot_layout(ncol = 1, guides = "collect", axes = "collect_x") 

site.ndep.plt <- (  beta.sig.all.plots[[7]]/ beta.sig.all.plots[[11]]/beta.sig.all.plots[[6]]   & coord_cartesian(xlim=c(-0.5, 0.7))) + plot_layout(ncol = 1, guides = "collect", axes = "collect_x") 

# save all the plots
save_plot(paste0(output.folder,"images/significant_betas_growth_dia.png"), 
          growth.plt, base_width = 5, base_height = 8) 

save_plot(paste0(output.folder,"images/significant_betas_growth_dia.svg"),
          growth.plt, base_width = 5, base_height = 8)

save_plot(paste0(output.folder,"images/significant_betas_competition.png"), 
          compete.plt, base_width = 5, base_height = 10) 

save_plot(paste0(output.folder,"images/significant_betas_competition.svg"),
          compete.plt, base_width = 5, base_height = 10)

save_plot(paste0(output.folder,"images/significant_betas_climate.png"), 
          climate.plt, base_width = 5, base_height = 10) 

save_plot(paste0(output.folder,"images/significant_betas_climate.svg"),
          climate.plt, base_width = 5, base_height = 10)


save_plot(paste0(output.folder,"images/significant_betas_siteNdep.png"), 
          site.ndep.plt, base_width = 5, base_height = 8) 

save_plot(paste0(output.folder,"images/significant_betas_siteNdep.svg"),
          site.ndep.plt, base_width = 5, base_height = 8)


#########################################################################
# SVG of variance partitions

# read in regional variance partitioning:
var_summary <- readRDS(paste0(output.dir, "SPCD_stanoutput_joint_v3/predicted_mort/var_summary_regional.rds"))

color.pred.class <- c(
  "Size" = "darkgreen",
  "Change in Size" = "#66a61e",
  "Site x G & S"= "#fdbf6f",
  "Competition x G & S"="#d95f02" ,
  "Climate x G & S"="#7570b3" ,
  "Climate"= "#e7298a" ,
  "Growth x Size" =  "#1b9e77",
  "Site Conditions" = "#e6ab02",
  "Competition" = "#a6761d"
)

color.pred.class.2 <- c(
  
  "Size" = "darkgreen",
  "Change in Size" = "#66a61e",
  "Growth x Size" = "#a1d99b" ,
  
  "Climate" = "#081d58",
  "Climate x G & S" = "#1d91c0",
  
  "N deposition" = "#67000d" , 
  "Ndep x G & S"= "#bd0026",
  
  "% Damage" = "#762a83", 
  "Damage x G & S" = "#9970ab",
  
  
  "Competition" = "#8c510a",
  "Competition x G & S" = "#d95f02" ,
  
  "Site Conditions" = "#e6ab02",
  "Site x G & S" = "yellow3")



var.summary.region <- var_summary
var.summary.region$COMMON <- factor(var.summary.region$COMMON, 
                                    levels = rev(c("balsam fir", "red spruce", "northern white-cedar", 
                                                   "eastern hemlock", "American beech", 
                                                   "black oak", "chestnut oak", "northern red oak", "white oak", "yellow birch", "paper birch", 
                                                   "hickory spp.", "eastern white pine", "red maple", "sugar maple", 
                                                   "black cherry", "white ash", "yellow-poplar")))
var.summary.region$predictor.class <- factor(var.summary.region$predictor.class, 
                                             levels = c(
                                               
                                               "Size" ,
                                               "Change in Size" ,
                                               "Climate" ,
                                               "Site Conditions",
                                               "Competition",
                                               "N deposition", 
                                               "% Damage", 
                                               "Growth x Size",
                                               
                                               "Site x G & S",
                                               "Competition x G & S" ,
                                               
                                               "Ndep x G & S",
                                               "Climate x G & S",
                                               "Damage x G & S"))


var.summary.region$predictor.class2 <- factor(var.summary.region$predictor.class2, 
                                              levels = c(
                                                
                                                "Size",
                                                "Change in Size",
                                                "Growth x Size" ,
                                                
                                                "Climate",
                                                "Climate x G & S",
                                                
                                                "N deposition", 
                                                "Ndep x G & S" ,
                                                
                                                "% Damage", 
                                                "Damage x G & S",
                                                
                                                
                                                "Competition",
                                                "Competition x G & S",
                                                
                                                "Site Conditions",
                                                "Site x G & S"
                                                
                                                
                                              ))




ggplot(data = var.summary.region)+
  geom_bar(aes(x = COMMON, y = mean, fill = predictor.class2), stat = "identity", position = "stack")+
  theme_bw(base_size = 16)+ theme(axis.text.x = element_text(angle = 60, hjust = 1), 
                                  panel.background = element_blank())+
  ylab("Across-tree variance in  \n p(survival) explained")+
  scale_fill_manual(values = color.pred.class.2, name = "",
                    labels = c("Change in Size" = expression(Delta ~ "D" ), 
                               "Growth x Size" = expression(Delta ~ "D x D"), 
                               "Size" = "Diameter (D)", 
                               "Climate x G & S" = expression("Climate x " ~Delta ~ "D or D"), 
                               "Ndep x G & S" = expression("N dep. x " ~Delta ~ "D or D"),
                               "Damage x G & S" = expression("Damage x " ~Delta ~ "D or D"), 
                               "Competition x G & S" = expression("Compeition x " ~Delta ~ "D or D"),
                               "Site x G & S" = expression("Site x " ~Delta ~ "D or D")
                    ),
                    guide = guide_legend(direction = "horizontal",
                                         ncol = 14,nrow = 1,reverse = TRUE,
                                         label.position="top", label.hjust = 0,
                                         label.vjust = 0.5,
                                         label.theme = element_text(angle = 90)))+
  coord_flip()+
  xlab("")+
  theme(panel.grid = element_blank(), 
        legend.position = "top", 
        panel.border = element_blank(), 
        axis.ticks.y = element_blank())


plt.regional.plt.all <- ggplot(data = var.summary.region, aes(x = COMMON, y = mean_logit, fill = predictor.class2))+
  #geom_bar(aes(x = COMMON, y = mean_logit, fill = predictor.class2), stat = "identity", position = "stack", width = 1)+
  stat_summary(geom = "col", fun = sum, width = 0.7,
               position = "stack") +
  theme_bw(base_size = 16)+ theme(axis.text.x = element_text(angle = 60, hjust = 1), 
                                  panel.background = element_blank())+
  ylab("Across-tree variance in  \n p(survival) explained")+
  scale_fill_manual(values = color.pred.class.2, name = "",
                    labels = c("Change in Size" = expression(Delta ~ "D" ), 
                               "Growth x Size" = expression(Delta ~ "D x D"), 
                               "Size" = "Diameter (D)", 
                               "Climate x G & S" = expression("Climate x " ~Delta ~ "D or D"), 
                               "Ndep x G & S" = expression("N dep. x " ~Delta ~ "D or D"),
                               "Damage x G & S" = expression("Damage x " ~Delta ~ "D or D"), 
                               "Competition x G & S" = expression("Compeition x " ~Delta ~ "D or D"),
                               "Site x G & S" = expression("Site x " ~Delta ~ "D or D")
                    ),
                    guide = guide_legend(direction = "horizontal",
                                         ncol = 14,nrow = 1,reverse = TRUE,
                                         label.position="top", label.hjust = 0,
                                         label.vjust = 0.5,
                                         label.theme = element_text(angle = 90)))+
  coord_flip()+
  xlab("")+
  theme(panel.grid = element_blank(), 
        legend.position = "top", 
        panel.border = element_blank(), 
        axis.ticks.y = element_blank())

ggsave(plot = plt.regional.plt.all, height = 8, width = 7.5, units = "in", dpi = 500, 
       paste0(output.dir, "images/regional_species_var_partitioning.png"))

ggsave(plot = plt.regional.plt.all, height = 8, width = 7.5, units = "in", dpi = 500, 
       paste0(output.dir, "images/regional_species_var_partitioning.svg"))


plt.regional.plt.all2 <- ggplot(data = var.summary.region, aes(x = COMMON, y = mean, fill = predictor.class2))+
  #geom_bar(aes(x = COMMON, y = mean_logit, fill = predictor.class2), stat = "identity", position = "stack", width = 1)+
  stat_summary(geom = "col", fun = sum, width = 0.7,
               position = "stack") +
  theme_bw(base_size = 16)+ theme(axis.text.x = element_text(angle = 60, hjust = 1), 
                                  panel.background = element_blank())+
  ylab("Across-tree variance in  \n p(survival) explained")+
  scale_fill_manual(values = color.pred.class.2, name = "",
                    labels = c("Change in Size" = expression(Delta ~ "D" ), 
                               "Growth x Size" = expression(Delta ~ "D x D"), 
                               "Size" = "Diameter (D)", 
                               "Climate x G & S" = expression("Climate x " ~Delta ~ "D or D"), 
                               "Ndep x G & S" = expression("N dep. x " ~Delta ~ "D or D"),
                               "Damage x G & S" = expression("Damage x " ~Delta ~ "D or D"), 
                               "Competition x G & S" = expression("Compeition x " ~Delta ~ "D or D"),
                               "Site x G & S" = expression("Site x " ~Delta ~ "D or D")
                    ),
                    guide = guide_legend(direction = "horizontal",
                                         ncol = 14,nrow = 1,reverse = TRUE,
                                         label.position="top", label.hjust = 0,
                                         label.vjust = 0.5,
                                         label.theme = element_text(angle = 90)))+
  coord_flip()+
  xlab("")+
  theme(panel.grid = element_blank(), 
        legend.position = "top", 
        panel.border = element_blank(), 
        axis.ticks.y = element_blank())

ggsave(plot = plt.regional.plt.all2, height = 8, width = 7.5, units = "in", dpi = 500, 
       paste0(output.dir, "images/regional_species_var_partitioning_meanpsurv.png"))

ggsave(plot = plt.regional.plt.all2, height = 8, width = 7.5, units = "in", dpi = 500, 
       paste0(output.dir, "images/regional_species_var_partitioning_meanpsurv.svg"))


# some summaries for putting numbers in the paper
var.summary.region %>% group_by(COMMON) %>% arrange(COMMON, desc(mean_logit)) %>% 
  filter(predictor.class2 %in% c("Size", "Change in Size", "Growth x Size"))%>%
  group_by(COMMON)%>%
  summarise(sum.logit = sum(mean_logit))

var.summary.region %>% group_by(COMMON) %>% arrange(COMMON, desc(mean_logit)) %>% 
  filter(!predictor.class2 %in% c("Size", "Change in Size", "Growth x Size"))%>%
  group_by(COMMON) %>% 
  group_by(predictor.class2, COMMON)%>%
  summarise(total_cat = sum(mean_logit)*100)%>%
  select(total_cat,  predictor.class2, COMMON) %>%
  
  filter(total_cat > 2.5) %>% 
  arrange( desc(total_cat)) %>%
  group_by(COMMON)%>% View()


var.summary.region %>% group_by(COMMON) %>% arrange(COMMON, desc(mean_logit)) %>% 
  filter(!predictor.class2 %in% c("Size", "Change in Size", "Growth x Size"))%>%
  mutate(logit_percent = round(mean_logit*100, digits = 2))%>%
  group_by(predictor.class2, COMMON)%>%
  mutate(total_cat = sum(mean_logit)*100)%>%
  select(Covariate, logit_percent, total_cat,  predictor.class2, COMMON) %>%
  filter(predictor.class2 %in%  "Competition x G & S")%>%
  #filter(total_cat > 5) %>% 
  arrange(desc(total_cat), desc(logit_percent)) %>% 
  filter(total_cat >1)%>% View()

# tree size:--


var.summary.region %>% group_by(COMMON) %>% arrange(COMMON, desc(mean_logit)) %>% 
  #filter(!predictor.class2 %in% c("Size", "Change in Size", "Growth x Size"))%>%
  mutate(logit_percent = round(mean_logit*100, digits = 2))%>%
  group_by(predictor.class2, COMMON)%>%
  mutate(total_cat = sum(mean_logit)*100)%>%
  select(Covariate, logit_percent, total_cat,  predictor.class2, COMMON) %>%
  filter(predictor.class2 %in%  "Size")%>%
  #filter(total_cat > 5) %>% 
  arrange(desc(total_cat), desc(logit_percent)) %>% 
  filter(logit_percent >1)%>% select(COMMON)

# delta diameter--
var.summary.region %>% group_by(COMMON) %>% arrange(COMMON, desc(mean_logit)) %>% 
  #filter(!predictor.class2 %in% c("Size", "Change in Size", "Growth x Size"))%>%
  mutate(logit_percent = round(mean_logit*100, digits = 2))%>%
  group_by(predictor.class2, COMMON)%>%
  mutate(total_cat = sum(mean_logit)*100)%>%
  select(Covariate, logit_percent, total_cat,  predictor.class2, COMMON) %>%
  filter(predictor.class2 %in%  "Change in Size")%>%
  #filter(total_cat > 5) %>% 
  arrange(desc(total_cat), desc(logit_percent)) %>% 
  filter(logit_percent >1)%>% View()

# diameter x delta diameter--
var.summary.region %>% group_by(COMMON) %>% arrange(COMMON, desc(mean_logit)) %>% 
  #filter(!predictor.class2 %in% c("Size", "Change in Size", "Growth x Size"))%>%
  mutate(logit_percent = round(mean_logit*100, digits = 2))%>%
  group_by(predictor.class2, COMMON)%>%
  mutate(total_cat = sum(mean_logit)*100)%>%
  select(Covariate, logit_percent, total_cat,  predictor.class2, COMMON) %>%
  filter(predictor.class2 %in%  "Growth x Size")%>%
  #filter(total_cat > 5) %>% 
  arrange(desc(total_cat), desc(logit_percent)) %>% 
  filter(logit_percent >1)%>% View()

# climate---
var.summary.region %>% group_by(COMMON) %>% arrange(COMMON, desc(mean_logit)) %>% 
  #filter(!predictor.class2 %in% c("Size", "Change in Size", "Growth x Size"))%>%
  mutate(logit_percent = round(mean_logit*100, digits = 2))%>%
  group_by(predictor.class2, COMMON)%>%
  mutate(total_cat = sum(mean_logit)*100)%>%
  select(Covariate, logit_percent, total_cat,  predictor.class2, COMMON) %>%
  filter(predictor.class2 %in%  "Climate")%>%
  #filter(total_cat > 5) %>% 
  arrange(desc(total_cat), desc(logit_percent)) %>% 
  filter(logit_percent >1)%>% View()

var.summary.region %>% group_by(COMMON) %>% arrange(COMMON, desc(mean_logit)) %>% 
  #filter(!predictor.class2 %in% c("Size", "Change in Size", "Growth x Size"))%>%
  mutate(logit_percent = round(mean_logit*100, digits = 2))%>%
  group_by(predictor.class2, COMMON)%>%
  mutate(total_cat = sum(mean_logit)*100)%>%
  select(Covariate, logit_percent, total_cat,  predictor.class2, COMMON) %>%
  filter(predictor.class2 %in%  "Climate x G & S")%>%
  #filter(total_cat > 5) %>% 
  arrange(desc(total_cat), desc(logit_percent)) %>% 
  filter(logit_percent >1)%>% View()

# damage --
var.summary.region %>% group_by(COMMON) %>% arrange(COMMON, desc(mean_logit)) %>% 
  #filter(!predictor.class2 %in% c("Size", "Change in Size", "Growth x Size"))%>%
  mutate(logit_percent = round(mean_logit*100, digits = 2))%>%
  group_by(predictor.class2, COMMON)%>%
  mutate(total_cat = sum(mean_logit)*100)%>%
  select(Covariate, logit_percent, total_cat,  predictor.class2, COMMON) %>%
  filter(predictor.class2 %in%  "% Damage")%>%
  #filter(total_cat > 5) %>% 
  arrange(desc(total_cat), desc(logit_percent)) %>% 
  filter(logit_percent >1)%>% View()

var.summary.region %>% group_by(COMMON) %>% arrange(COMMON, desc(mean_logit)) %>% 
  #filter(!predictor.class2 %in% c("Size", "Change in Size", "Growth x Size"))%>%
  mutate(logit_percent = round(mean_logit*100, digits = 2))%>%
  group_by(predictor.class2, COMMON)%>%
  mutate(total_cat = sum(mean_logit)*100)%>%
  select(Covariate, logit_percent, total_cat,  predictor.class2, COMMON) %>%
  filter(predictor.class2 %in%  "Damage x G & S")%>%
  #filter(total_cat > 5) %>% 
  arrange(desc(total_cat), desc(logit_percent)) %>% 
  filter(logit_percent >1)%>% View()

# neighborhood effects --
var.summary.region %>% group_by(COMMON) %>% arrange(COMMON, desc(mean_logit)) %>% 
  #filter(!predictor.class2 %in% c("Size", "Change in Size", "Growth x Size"))%>%
  mutate(logit_percent = round(mean_logit*100, digits = 2))%>%
  group_by(predictor.class2, COMMON)%>%
  mutate(total_cat = sum(mean_logit)*100)%>%
  select(Covariate, logit_percent, total_cat,  predictor.class2, COMMON) %>%
  filter(predictor.class2 %in%  "Competition")%>%
  #filter(total_cat > 5) %>% 
  arrange(desc(total_cat), desc(logit_percent)) %>% 
  filter(logit_percent >1)%>% View()

var.summary.region %>% group_by(COMMON) %>% arrange(COMMON, desc(mean_logit)) %>% 
  #filter(!predictor.class2 %in% c("Size", "Change in Size", "Growth x Size"))%>%
  mutate(logit_percent = round(mean_logit*100, digits = 2))%>%
  group_by(predictor.class2, COMMON)%>%
  mutate(total_cat = sum(mean_logit)*100)%>%
  select(Covariate, logit_percent, total_cat,  predictor.class2, COMMON) %>%
  filter(predictor.class2 %in%  "Competition x G & S")%>%
  #filter(total_cat > 5) %>% 
  arrange(desc(total_cat), desc(logit_percent)) %>% 
  filter(logit_percent >1)%>% View()

# N deposition --
var.summary.region %>% group_by(COMMON) %>% arrange(COMMON, desc(mean_logit)) %>% 
  #filter(!predictor.class2 %in% c("Size", "Change in Size", "Growth x Size"))%>%
  mutate(logit_percent = round(mean_logit*100, digits = 2))%>%
  group_by(predictor.class2, COMMON)%>%
  mutate(total_cat = sum(mean_logit)*100)%>%
  select(Covariate, logit_percent, total_cat,  predictor.class2, COMMON) %>%
  filter(predictor.class2 %in%  "N deposition")%>%
  #filter(total_cat > 5) %>% 
  arrange(desc(total_cat), desc(logit_percent)) %>% 
  filter(logit_percent >1)%>% View()

var.summary.region %>% group_by(COMMON) %>% arrange(COMMON, desc(mean_logit)) %>% 
  #filter(!predictor.class2 %in% c("Size", "Change in Size", "Growth x Size"))%>%
  mutate(logit_percent = round(mean_logit*100, digits = 2))%>%
  group_by(predictor.class2, COMMON)%>%
  mutate(total_cat = sum(mean_logit)*100)%>%
  select(Covariate, logit_percent, total_cat,  predictor.class2, COMMON) %>%
  filter(predictor.class2 %in%  "Ndep x G & S")%>%
  #filter(total_cat > 5) %>% 
  arrange(desc(total_cat), desc(logit_percent)) %>% 
  filter(logit_percent >1)%>% View()


# site conditions --
var.summary.region %>% group_by(COMMON) %>% arrange(COMMON, desc(mean_logit)) %>% 
  #filter(!predictor.class2 %in% c("Size", "Change in Size", "Growth x Size"))%>%
  mutate(logit_percent = round(mean_logit*100, digits = 2))%>%
  group_by(predictor.class2, COMMON)%>%
  mutate(total_cat = sum(mean_logit)*100)%>%
  select(Covariate, logit_percent, total_cat,  predictor.class2, COMMON) %>%
  filter(predictor.class2 %in%  "Site Conditions")%>%
  #filter(total_cat > 5) %>% 
  arrange(desc(total_cat), desc(logit_percent)) %>% 
  filter(logit_percent >1)%>% View()

var.summary.region %>% group_by(COMMON) %>% arrange(COMMON, desc(mean_logit)) %>% 
  #filter(!predictor.class2 %in% c("Size", "Change in Size", "Growth x Size"))%>%
  mutate(logit_percent = round(mean_logit*100, digits = 2))%>%
  group_by(predictor.class2, COMMON)%>%
  mutate(total_cat = sum(mean_logit)*100)%>%
  select(Covariate, logit_percent, total_cat,  predictor.class2, COMMON) %>%
  filter(predictor.class2 %in%  "Site x G & S")%>%
  #filter(total_cat > 5) %>% 
  arrange(desc(total_cat), desc(logit_percent)) %>% 
  filter(logit_percent >1)%>% View()
#########################################################################
# Marginal effects, but plotting effects for species in affected by different pests



marginal_response_df$predictor.class2
#interaction_response_df<- left_join(interaction_response_df, Covariate.types.df)

spruce.fir <- c("balsam fir", "red spruce",  "northern white-cedar")
mixed <- c("American beech", "eastern hemlock")
spongy.susceptible<- c("chestnut oak", "white oak", "northern red oak", "paper birch", "yellow birch")
spongy.immune <- c("white ash", "yellow-poplar", "black cherry")
spongy.resist <- c("sugar maple", "red maple", "eastern white pine", "hickory spp.")





# set up groups of species to look at:
spruce.fir <- c("balsam fir", "red spruce",  "northern white-cedar")
mixed <- c("American beech", "eastern hemlock")
spongy.susceptible<- c("chestnut oak", "white oak",  "northern red oak", "paper birch", "yellow birch")
spongy.immune <- c("white ash", "yellow-poplar", "black cherry")
spongy.resist <- c("sugar maple", "red maple", "eastern white pine", "hickory spp.")

interaction.terms <- unique(interaction_response_df$Covariate)
all_marginal_response_df <- marginal_response_df

# set up species colors:
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



marginal_response_df

# covariates that are scaled by species:
train.data <- readRDS (

  paste0(
    output.folder,
    "train_data_all_SPCD_model_",
    model.no,
    ".RDS"
  )
)

test.data <- readRDS (

  paste0(
    output.folder,
    "test_data_all_SPCD_model_",
    model.no,
    ".RDS"
  )
)


#Ndep_Diff_per_yr
species.scaling <- train.data %>% select(Species, SPCD, Ndep_Diff_per_yr.median, Ndep_Diff_per_yr.sd) %>% distinct()%>%
  mutate(Covariate = "Ndep.scaled", 
         Clean_Name = "delta N deposition",
         Units = "kg/m^2/year") %>% 
  rename("Val.mean" = "Ndep_Diff_per_yr.median",
         "Val.sd" = "Ndep_Diff_per_yr.sd")%>%
  #DIA_DIFF_scaled
  rbind(., train.data %>% select(Species, SPCD, DIA.DIFF.median, DIA.DIFF.sd) %>% distinct()%>%
          mutate(Covariate = "DIA_DIFF_scaled", 
                 Clean_Name = "Diameter Difference",
                 Units = "Inches")%>% 
          rename("Val.mean" = "DIA.DIFF.median",
                 "Val.sd" = "DIA.DIFF.sd"))%>%
  #BAL
  rbind(.,train.data %>% select(Species, SPCD, BAL.median, BAL.sd) %>% distinct()%>%
          mutate(Covariate = "BAL.scaled",
                 Clean_Name = "Basal Area Larger Than",
                 Units = "ft^2")%>% 
          rename("Val.mean" = "BAL.median",
                 "Val.sd" = "BAL.sd")
  )%>%
  #DIA
  rbind(., train.data %>% select(Species, SPCD, DIA.median, DIA.sd) %>% distinct()%>%
          mutate(Covariate = "DIA_scaled", 
                 Clean_Name = "Diameter",
                 Units = "Inches")%>% 
          rename("Val.mean" = "DIA.median",
                 "Val.sd" = "DIA.sd") 
  )%>%
  #plt_ba_sq_ft_old
  rbind(., train.data %>% select(Species, SPCD, plt_ba_sq_ft_old.median, plt_ba_sq_ft_old.sd) %>% distinct()%>%
          mutate(Covariate = "ba.scaled", 
                 Clean_Name = "Plot Basal Area",
                 Units = "ft^2")%>% 
          rename("Val.mean" = "plt_ba_sq_ft_old.median",
                 "Val.sd" = "plt_ba_sq_ft_old.sd")
  )%>%
  select(Species, SPCD, Covariate, Val.mean, Val.sd, Clean_Name, Units)


#--------------------------------------------------------------------------------
# generate marginal predictions for specific sizes by species (diameter marking thresholds)
#--------------------------------------------------------------------------------
size.raw.scaling <- train.data %>% select(Species, SPCD, DIA.median, DIA.sd) %>% distinct()%>%
  mutate(Covariate = "DIA_scaled", 
         Clean_Name = "Diameter",
         Units = "Inches")%>% 
  rename("Val.mean" = "DIA.median",
         "Val.sd" = "DIA.sd")

# get the posterior samples of betas and alphas:
alpha.fit <- readRDS(paste0(input.folder, "samples/alpha.spp_model_", model.no, "_5000samples.rds"))
alpha_df <- as_draws_df(alpha.fit) 



# get all the covariates using posterior package
betas.df <- readRDS(paste0(input.folder, "samples/u_betas_model_", model.no, "_5000samples.rds"))
betas.quant <- betas.df %>% summarise_draws(median, ~quantile(., probs = c(0.025, 0.975))) %>%
  rename(`ci.lo` = "2.5%", `ci.hi` = "97.5%") %>%
  mutate(remper.cor = 0.5)
# relabel u_betas to meaningful species ids names
betas.quant$spp <- rep(1:17, ncovar)
betas.quant$cov <- rep(1:ncovar, each = 17)


covariate_names <- c(colnames(mod.data$xM))  # Replace with your covariate names
betas.quant$Covariate <- rep(covariate_names, each = 17)
betas.quant$Species <- rep(nspp[1:17,]$COMMON, ncovar)


nspp <- data.frame(SPCD = c(316, 318, 833, 832, 261, 531, 802, 129, 762,  12, 541,  97, 621, 400, 371, 241, 375))
nspp$Species <- paste(FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD),]$GENUS, FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD),]$SPECIES)
# link up to the species table:
nspp$COMMON <- FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD),]$COMMON

spp.table <- data.frame(SPCD.id = nspp[1:17,]$SPCD, 
                        spp = 1:17, 
                        COMMON = nspp[1:17,]$COMMON)


# this function predicts posterior samples of main effects (only does single main effects):
posterior.predict.single <- function(SPCD.id, Beta.var, raw.value, raw.scaling = size.raw.scaling){
      
      
      i <- spp.table[which(spp.table$SPCD.id %in% SPCD.id),]$spp
     
      scaled.value <- raw.scaling %>% 
        filter(Covariate %in% Beta.var) %>%
        filter(SPCD %in% SPCD.id)%>%
        mutate(scaled.value = (raw.value - Val.mean)/Val.sd)
      
      
      beta.species.names <- betas.quant %>% filter(spp %in% i) %>% select(variable)
      # select only the betas and intercepts for the species of interest:
      beta <- subset_draws(betas.df, variable = beta.species.names$variable) %>% select(-.chain, -.iteration, -.draw) 
      beta_0 <- data.frame(subset_draws(alpha_df, variable = paste0("alpha_SPP[",i,"]"))) %>% select(-.chain, -.iteration, -.draw)  # Intercept
      # read in the species data and covariates to get the min and max and get ranges
      mod.data <-
        readRDS (
          paste0(
            input.folder,
            "all_SPCD_model_",
            model.no,
            ".RDS"
          )
        )
      
      
      
      covariate_names <- c(colnames(mod.data$xM))  
      # set up a matrxi of zeros with the covariate ranges
      covariate_ranges_df <- data.frame(names = covariate_names, 
                                        medians = 0) %>%
                                    mutate(medians = ifelse(names %in% Beta.var, scaled.value$scaled.value, 0))
                                        
      covariate_matrix <- covariate_matrix <- matrix(covariate_ranges_df$medians, 
                                                     nrow = length(raw.value), 
                                                     ncol = length(covariate_names), 
                                                     byrow = TRUE)
      # set up a function calculate probabilities
      inv_logit_fxn <- function(x) {
        1 / (1 + exp(-x))
      }
      # probabilities <- list()
      
        linear_predictor <- beta_0[,1] + rowSums( as.matrix(beta)%*% covariate_matrix[1,] )
        probabilities <-  as.vector(inv_logit_fxn(linear_predictor))
    
      
      # output the simple survival probabilities:
posterior.pred.out <- data.frame(psurv = probabilities, 
                                 iter = 1:length(probabilities),
                                 SPCD = SPCD.id, 
                                 COMMON = spp.table[which(spp.table$SPCD.id %in% SPCD.id),]$COMMON,
                                 Covariate = Beta.var, 
                                 Raw = raw.value, 
                                 scaled = scaled.value$scaled.value)
return(posterior.pred.out)
}


# for each species, thresholds, get the posterior samples:
species.econ.maturity <- data.frame(COMMON = 
                      # those with 16-18 inches
                    c(rep(c("red maple",
                           "American beech",
                       "sugar maple", 
                      
                       "yellow birch",
                       
                       "northern red oak", 
                       
                       "white ash"), each = 2 
), "paper birch", "paper birch", 
"red spruce", "red spruce",
"eastern hemlock", "eastern hemlock",

# those with 18-24 inch maturation for long-lived hardwoods
c( "sugar maple", 
      
      "yellow birch",
      
      "northern red oak", 
      
      
      "white ash"), 

# balsam fir and red spruce maturation ages from silvics manaul:
# balsam fir: 12-18in, red spruce = 12-24 in
# balsam fir economic maturity from 12-14 inches
 "balsam fir", "balsam fir", "balsam fir",

# pulpwood merchantable: 5inch
 "balsam fir",
 "red spruce" ,

# maple tapping ranges:
"sugar maple", "sugar maple", "sugar maple", 
"red maple", "red maple","red maple"
), 
           Diameter.val = c(rep(c(16, 18), times = 6), c(12,14), c(12,16), c(18,24), rep(24, 4), c(12, 14, 18), c(5,5), 
                            rep(c(10, 17, 25), 2)), 
           Diameter.type = c(rep(c("low", "high"), times = 9), rep("high-high", 4), c("low", "high", "high-high"), c("pulp", "pulp"), 
                             rep(c("single tap", "single tap high", "two taps"),2))
) %>% left_join(.,spp.table)


post.samp.dias <- list()
for(k in 1:length(species.econ.maturity$COMMON)){
  post.samp.dias[[k]] <- posterior.predict.single(SPCD.id = species.econ.maturity[k,]$SPCD.id, 
                           Beta.var = "DIA_scaled", 
                           raw.value = species.econ.maturity[k,]$Diameter.val)
}

post.samp.dias.df <- do.call(rbind, post.samp.dias) %>% 
  group_by(SPCD, COMMON, Raw) %>% summarise(pmort_10 = median(1-psurv^10), 
                                    pmort_10.lo = quantile(1-psurv^10, 0.05), 
                                    pmort_10.hi = quantile(1-psurv^10, 0.95))%>% 
  rename("Diameter.val" = "Raw")%>%
  left_join(.,species.econ.maturity)

ggplot(data = post.samp.dias.df, aes(x = Diameter.val, y = pmort_10, fill = COMMON))+
  geom_bar(stat = "identity")+
  geom_errorbar(data = post.samp.dias.df, aes(x = Diameter.val, ymin = pmort_10.lo, ymax = pmort_10.hi))+
  facet_wrap(~COMMON, scales = "free_y")+
  species_fill


# just do it for the any or moderate ranges:
ranges.spp <- species.econ.maturity %>% filter(Diameter.type %in% c("low", "high"))%>%
  mutate(Diameter.val.cm = round(Diameter.val*2.54, digits = 1)) %>%
  select(COMMON, Diameter.type, Diameter.val.cm)%>% 
  group_by(COMMON) %>% 
  spread(Diameter.type, Diameter.val.cm) %>%
  mutate(Range = paste0("(", low,"cm - ", high, " cm)"))%>%
  select(COMMON, Range)

post.samp.dias.df%>% select(COMMON, SPCD, Diameter.type, pmort_10) %>%
  group_by(COMMON, SPCD) %>% spread(Diameter.type, pmort_10) %>%
  mutate(delta.low.high = high - low, 
         pct.delta.low.high = ((high - low)/low)*100, 
         pct.detla.low.high.high = ((`high-high` - low)/low)*100) %>% 
  select(COMMON, SPCD, low, high, pct.delta.low.high) %>%
  left_join(., ranges.spp)%>%
  rename("Species" = "COMMON", 
         "% change (high-low)"= "pct.delta.low.high", 
         "Diameter Range" = "Range") %>% 
  arrange(desc(`% change (high-low)`))

post.samp.dias.df%>% select(COMMON, SPCD, Diameter.type, pmort_10) %>%
  group_by(COMMON, SPCD) %>% filter(COMMON %in% c("balsam fir", "red spruce"))#%>% 
  #filter(Diameter.type %in% "pulp")

post.samp.dias.df%>% select(COMMON, SPCD, Diameter.type, pmort_10) %>%
  group_by(COMMON, SPCD) %>% filter(COMMON %in% c("sugar maple", "red maple"))


post.samp.dias.df.longevity <- do.call(rbind, post.samp.dias) %>% 
  group_by(SPCD, COMMON, Raw) %>% 
  mutate(pmort_annual = 1-psurv) %>% summarise(pmort_1 = median(pmort_annual), 
                                            longevity = median(1/pmort_annual), 
                                            longevity.lo =quantile(1/pmort_annual, 0.05),
                                            longevity.hi = quantile(1/pmort_annual, 0.95),
                                            pmort_10 = median(1-psurv^10), 
                                            pmort_10.lo = quantile(1-psurv^10, 0.05), 
                                            pmort_10.hi = quantile(1-psurv^10, 0.95))%>% 
  rename("Diameter.val" = "Raw")%>%
  left_join(.,species.econ.maturity)

post.samp.dias.df.longevity %>% #select(COMMON, SPCD, Diameter.type, longevity, pmort_10) %>%
  group_by(COMMON, SPCD) %>% filter(COMMON %in% c("sugar maple", "red maple"))


post.samp.dias.df.longevity %>% select(COMMON, SPCD, Diameter.type, longevity, pmort_10) %>%
  group_by(COMMON, SPCD) %>% filter(COMMON %in% c("balsam fir", "red spruce"))



####################################################################################
# marginal response to diameter differences:
size_diff.raw.scaling <- train.data %>% select(Species, SPCD, DIA.DIFF.median, DIA.DIFF.sd) %>% distinct()%>%
  mutate(Covariate = "DIA_DIFF_scaled", 
         Clean_Name = "Diameter_DIFF",
         Units = "Inches")%>% 
  rename("Val.mean" = "DIA.DIFF.median",
         "Val.sd" = "DIA.DIFF.sd")



species.dia.diff.medians <- train.data %>% left_join(., size_diff.raw.scaling) %>% 
  group_by(M, Species, SPCD, remper) %>% 
  mutate(avg.growth.in = DIA_DIFF/remper) %>% 
  ungroup()%>%
  group_by(M, SPCD, Species)%>%
  summarise(avg.growth.cm = mean(avg.growth.in)*2.54,
            mean.dia.diff.cm = mean(DIA_DIFF)*2.54, 
            mean.dia.scaled = mean(DIA_DIFF_scaled),
            #calc.mean.dia.scaled = mean((DIA_DIFF-DIA.DIFF.median)/DIA.DIFF.sd),
            DIA.DIFF.median = mean(DIA.DIFF.median), 
            DIA.DIFF.sd  = mean(DIA.DIFF.sd),
            remper = mean(remper)) %>%
  mutate(mean.dia.in = (mean.dia.scaled*DIA.DIFF.sd) + DIA.DIFF.median) %>%
  mutate(mean.dia.cm = mean.dia.in*2.54) %>%
  rename(
    "SPCD.id" = "SPCD"
  )%>%
  left_join(., spp.table)

spp.table

train.data %>% 
  group_by(M, Species, SPCD, remper) %>% 
  mutate(avg.growth.in = DIA_DIFF/remper) %>% 
  ungroup()%>%
  group_by(M)%>%
  summarise(avg.growth.cm = mean(avg.growth.in)*2.54,
            mean.dia.diff.cm = mean(DIA_DIFF)*2.54, 
            mean(remper))

# specific predictions for diameter differences


post.samp.diadiffs <- list()
for(k in 1:length(species.dia.diff.medians$COMMON)){
  post.samp.diadiffs[[k]] <- posterior.predict.single(SPCD.id = species.dia.diff.medians[k,]$SPCD.id,
                                                  Beta.var = "DIA_DIFF_scaled",
                                                  raw.value = species.dia.diff.medians[k,]$mean.dia.in, 
                                                  raw.scaling = size_diff.raw.scaling)
}

post.samp.diadifs.df <- do.call(rbind, post.samp.diadiffs) %>%
  group_by(SPCD, COMMON, Raw) %>% summarise(pmort_1 = median(1-psurv),
                                            pmort_1.lo = quantile(1-psurv, 0.05),
                                            pmort_1.hi = quantile(1-psurv, 0.95),
                                            
                                            pmort_10 = median(1-psurv^10),
                                            pmort_10.lo = quantile(1-psurv^10, 0.05),
                                            pmort_10.hi = quantile(1-psurv^10, 0.95))%>%
  rename("mean.dia.in" = "Raw")%>%
  left_join(.,species.dia.diff.medians)


# for raw diameter diferences

# 
ggplot(data = post.samp.diadifs.df, aes(x = mean.dia.in, y = pmort_10, fill = COMMON))+
  geom_bar(stat = "identity")+
  geom_errorbar(data = post.samp.diadifs.df, aes(x = mean.dia.in, ymin = pmort_10.lo, ymax = pmort_10.hi))+
  facet_wrap(~COMMON)+
  species_fill

ggplot(data = post.samp.diadifs.df, aes(x = mean.dia.in*2.54, y = pmort_1, color = COMMON))+
  geom_point()+
  geom_errorbar(data = post.samp.diadifs.df, aes(x = mean.dia.in*2.54, ymin = pmort_1.lo, ymax = pmort_1.hi))+
  species_color

ggplot(data = post.samp.diadifs.df, aes(x = mean.dia.in*2.54, y = pmort_10, color = COMMON))+
  geom_point()+
  geom_errorbar(data = post.samp.diadifs.df, aes(x = mean.dia.in*2.54, ymin = pmort_10.lo, ymax = pmort_10.hi))+
  species_color

marginal_response_df_unscaled %>% 
  filter(Clean_Name %in% "Diameter")%>%
  mutate(Raw.value.cm = Raw.value*2.54)%>%
  filter(Species %in% c("sugar maple", "northern red oak", "white ash", "yellow birch"))|>
  ggplot()+geom_line(aes(x = Raw.value.cm, y = 1-(mean)^10, color = Species))+
  geom_ribbon(aes(x = Raw.value.cm, ymin = 1-(ci.lo)^10, ymax = 1-(ci.hi)^10, fill = Species), alpha = 0.1, color = NA)+
  theme_bw()+species_color+species_fill+
  geom_vline(aes(xintercept  = (16*2.54)))+
  geom_vline(aes(xintercept = (18*2.54)))+
  geom_vline(aes(xintercept = (24*2.54)))

# low end of diameter thresholds for these species

# find the closest diameter target value for each species:



marginal_response_df_unscaled %>% 
  filter(Clean_Name %in% "Diameter")%>%
  mutate(Raw.value.cm = Raw.value*2.54)%>%
  filter(Species %in% c("sugar maple", "northern red oak", "white ash", "yellow birch"))%>%
  filter(Raw.value >= 16 & Raw.value < 17)%>%
  mutate(mort_10_yr.med = 1-(mean)^10, 
         mort_10yr.ci.lo = 1-(ci.lo)^10, 
         mort_10yr.ci.hi = 1-(ci.hi)^10)

species.name <- "sugar maple"
target.Diameter.in <- 16
get_est_mort_target <- function(species.name, target.Diameter.in){
  marginal_response_df_unscaled %>% 
    filter(Clean_Name %in% "Diameter")%>%
    mutate(Raw.value.cm = Raw.value*2.54)%>%
    filter(Species %in% species.name)%>%
    group_by(Species)%>%
    filter(abs(Raw.value - target.Diameter.in) == min(abs(Raw.value - target.Diameter.in)))%>%
    mutate(mort_10_yr.med = 1-(mean)^10, 
           mort_10yr.ci.lo = 1-(ci.lo)^10, 
           mort_10yr.ci.hi = 1-(ci.hi)^10)%>%
      mutate(Target.diam = target.Diameter.in)
}

D.16 <- get_est_mort_target(species.name = c( "sugar maple", 
                                      "yellow birch",
                                      "white ash",
                                      "northern red oak",
                                      "red maple",
                                     "American beech", 
                                     "red spruce"), 
                    target.Diameter.in = 16)

D.18 <- get_est_mort_target(species.name = c( "sugar maple", 
                                              "yellow birch",
                                              "white ash",
                                              "northern red oak",
                                              "red maple",
                                              "American beech", 
                                              "eastern hemlock"), 
                            target.Diameter.in = 18)

D.24 <- get_est_mort_target(species.name = c( "sugar maple", 
                                              "yellow birch",
                                              "white ash",
                                              "northern red oak",
                                              "eastern hemlock"), 
                            target.Diameter.in = 24)

D.12 <- get_est_mort_target(species.name = c( "paper birch", "red spruce"), 
                            target.Diameter.in = 12)
D.14 <- get_est_mort_target(species.name = c( "paper birch"), 
                            target.Diameter.in = 14)


# combine all of these together:
D.threshold.pmort = rbind(D.12, D.14, D.16, D.18, D.24)

D.threshold.pmort|>
  ggplot()+geom_bar(aes(x = Target.diam, y = mort_10_yr.med), stat = "identity")+
  geom_errorbar(aes(x = Target.diam, ymin = mort_10yr.ci.lo, ymax = mort_10yr.ci.hi), width = 0.1)+
  facet_wrap(~Species, scales = "free_y")

marginal_response_df_unscaled %>% 
  filter(Clean_Name %in% "Diameter")%>%
  mutate(Raw.value.cm = Raw.value*2.54)%>%
  filter(Species %in% c("sugar maple", "northern red oak", "white ash", "yellow birch"))%>%
  filter(Raw.value >= 24 & Raw.value < 25)%>%
  mutate(mort_10_yr.med = 1-(mean)^10, 
         mort_10yr.ci.lo = 1-(ci.lo)^10, 
         mort_10yr.ci.hi = 1-(ci.hi)^10)

# covariates scaled across space:
plot.medians <- readRDS("data/plot.medians_SPCD_all.rds")
plot.scaled <- plot.medians %>% select(damage.median, damage.sd) %>% 
  mutate(Covariate = "damage.scaled",
         Clean_Name = "Damage",
         Units = "% of stems")%>%
  rename("Val.mean" = "damage.median",
         "Val.sd" = "damage.sd")%>%
  
  #damage.scaled = (damage.total - plot.medians$damage.median)/plot.medians$damage.sd,
  #MAP.scaled = (MAP-plot.medians$MAP.median)/plot.medians$MAP.sd,
  rbind(., plot.medians %>% select(MAP.median, MAP.sd) %>% 
          mutate(Covariate = "MAP.scaled", 
                 Clean_Name = "Mean Annual Precipitation",
                 Units = "mm")%>%
          rename("Val.mean" = "MAP.median",
                 "Val.sd" = "MAP.sd"))%>%
  #MATmax.scaled = (MATmax - plot.medians$MATmax.median)/plot.medians$MATmax.sd)
  rbind(.,plot.medians %>% select(MATmax.median, MATmax.sd) %>% 
          mutate(Covariate = "MATmax.scaled", 
                 Clean_Name = "Average Annual Maximum Temperature",
                 Units = "Degrees C")%>%
          rename("Val.mean" = "MATmax.median",
                 "Val.sd" = "MATmax.sd"))%>%
  
  
  rbind(.,plot.medians %>% select(slope.median, slope.sd) %>% 
          mutate(Covariate = "slope.scaled", 
                 Clean_Name = "Slope",
                 Units = "%")%>%
          rename("Val.mean" = "slope.median",
                 "Val.sd" = "slope.sd"))%>%
  expand_grid(., Species = unique(species.scaling$Species))%>%
  left_join(., unique(species.scaling[,c("Species", "SPCD")]))%>%
  select(Species, SPCD, Covariate, Val.mean, Val.sd, Clean_Name, Units)

# data.frame with the median and sd for each species to rescale the variable by:
combined.scaled.main <- rbind(plot.scaled, species.scaling)

marginal_response_df_unscaled <- marginal_response_df %>% left_join(., combined.scaled.main)%>%
  mutate(Raw.value = ifelse(!is.na(Val.mean), (Value*Val.sd) + Val.mean, Value)) %>%
  mutate(Raw.value.converted = ifelse(Units %in% "Inches", Raw.value*2.54,
                                      ifelse(Units %in% "ft^2", Raw.value*0.0929, Raw.value)), 
         Units.converted = ifelse(Units %in% "Inches", "cm", 
                                  ifelse(Units %in% "ft^2", "m^2", Units)))%>%
  mutate(Predictor.strip = ifelse(Predictor %in% "Diameter Difference \n (DIA_DIFF)", "Diameter Difference", 
                            ifelse(Predictor %in% "Plot Basal Area (BA)", "Plot Basal Area", 
                                   ifelse(Predictor %in% "Basal Area Larger \n (BAL)", "Basal Area Larger", paste0(Predictor)))))

# marginal_response_df_unscaled
#Ndep X dia_diff
#p1.value   p2.value p2.rank             covariate     Pred.1          Pred.2
#1  -1.12091712  0.1224049 

interaction_response_df_unscaled <- interaction_response_df %>% 
  # join the first precdictor with the raw values of the covariate
  left_join(., combined.scaled.main %>% rename("Pred.1" = "Covariate"))%>%
  mutate(Raw.value = ifelse(!is.na(Val.mean), (p1.value*Val.sd) + Val.mean, p1.value)) %>%
  left_join(., combined.scaled.main %>% rename("Pred.2" = "Covariate", 
                                               "Val.mean.2" = "Val.mean", 
                                               "Val.sd.2" = "Val.sd", 
                                               "Clean_Name.2" = "Clean_Name", 
                                               "Units.2" = "Units"))%>%
  mutate(Raw.value.p2 = ifelse(!is.na(Val.mean.2), (p2.value*Val.sd.2) + Val.mean.2, p2.value))%>%
  mutate(Raw.value.converted = ifelse(Units %in% "Inches", Raw.value*2.54,
                                      ifelse(Units %in% "ft^2", Raw.value*0.0929, Raw.value)), 
         Units.converted = ifelse(Units %in% "Inches", "cm", 
                                  ifelse(Units %in% "ft^2", "m^2", Units)))%>%
  mutate(Predictor.strip = ifelse(Predictor %in% "Diameter Difference \n (DIA_DIFF)", "Diameter Difference", 
                                  ifelse(Predictor %in% "Plot Basal Area (BA)", "Plot Basal Area", 
                                         ifelse(Predictor %in% "Basal Area Larger \n (BAL)", "Basal Area Larger", paste0(Predictor)))))


interaction_response_df_unscaled %>% select(Species, p2.rank, Pred.2, Raw.value.p2) %>% distinct()
# need to join the interaction plots by the median and sd values and convert p1 and p2 values
# function to plot up the effects on 10 year mortality probabilities
# using the unscaled "raw" values for covariates




plot_main_region_effects <- function(species.group, covar, ymax.spp){
  
  
  
  if(! covar %in% interaction.terms){  
    
    df.species <- marginal_response_df_unscaled %>% filter(Species %in%  species.group & Covariate %in% covar)
    strip.fill <- as.character(color.pred.class.2[unique(df.species$predictor.class2)])
    
    
    
    p1 <-  ggplot(df.species , aes(x = Raw.value.converted, y = 1-(mean)^10, color = Species)) +
      geom_line(size = 1) +
      geom_ribbon(aes(ymin = 1-(ci.lo)^10, ymax = 1-(ci.hi)^10, fill = Species), alpha = 0.1, color = NA) +
      #geom_label(aes(label = Species, color = Species))+
      facet_wrap(~Predictor.strip) +
      labs(
        x = paste(df.species$predictor.class2, "(", unique(df.species$Units.converted),")"),
        y = "10-year Mortality Probability",
        #title = "Effect of Predictors on Probability of Mortality"
      ) +
      species_color + species_fill+
      #var.part.fill + var.part.color+
      theme_bw(base_size = 14)+theme(panel.grid = element_blank(), 
                                     legend.key.size = unit(1, "cm"), 
                                     legend.position = "none", 
                                     strip.background = element_rect(fill = strip.fill), 
                                     strip.text = element_text(color = "white")
      )+coord_cartesian(ylim = c(0, ymax.spp))
    
    
  }else{
    
    df.species <- interaction_response_df_unscaled %>% filter(Species %in%  species.group & Covariate %in% covar)
    
    
    df.species$p2.rank <- factor(df.species$p2.rank, levels = c("high", "median", "low"))
    strip.fill <- as.character(color.pred.class.2[unique(df.species$predictor.class2)])
    pred.name <-  marginal_response_df %>% filter(Covariate %in% unique(df.species$Pred.1))%>% select(Predictor)%>% distinct()
    
    if(unique(df.species$Pred.2) %in% "DIA_DIFF_scaled"){
      
      p1 <-  ggplot(data = df.species) +
        geom_line(aes(x = Raw.value.converted, y = 1-(mean)^10, color = Species, group = interaction(Species, p2.rank), linetype = p2.rank), size = 1) +
        geom_ribbon(aes(x = Raw.value.converted, ymin = 1-(ci.lo.10)^10, ymax = 1-(ci.hi.90)^10, fill = Species, group = interaction(Species, p2.rank)), color = NA, alpha = 0.1) +
        #facet_grid(cols = vars(predictor.class2), rows = vars(Predictor)) +
        facet_wrap(~Predictor)+
        labs(
          x = paste(pred.name$Predictor, "(", unique(df.species$Units.converted),")"),
          y = "10-year Mortality Probability",
          #title = "Effect of Predictors on Probability of Mortality"
        ) + scale_linetype_manual(values = c("low"= "dotted" ,
                                             "median" = "solid",
                                             "high"= "dashed" ), 
                                  name = expression(Delta ~ "Diameter"))+
        species_color + species_fill+
        # scale_fill_manual(values = c("low"="#a1dab4" ,
        #                              "median" = "#41b6c4",
        #                              "high"= "#225ea8" ), 
        #                   name = expression(Delta ~ "Diameter"))+
        
        #species_fill + species_color + 
        #named_species_linetype +
        theme_bw(base_size = 14)+theme(panel.grid = element_blank(), 
                                       legend.key.size = unit(1, "cm"), 
                                       strip.background = element_rect(fill = strip.fill)
        )+coord_cartesian(ylim = c(0, ymax.spp))
      
      
    }else{
      
      p1 <-  ggplot(data = df.species) +
        geom_line(aes(x = Raw.value.converted, y = 1-(mean)^10, color = Species, group = interaction(Species, p2.rank), linetype = p2.rank, size = p2.rank)) +
        geom_ribbon(aes(x = Raw.value.converted, ymin = 1-(ci.lo.10)^10, ymax = 1-(ci.hi.90)^10, fill = Species, group = interaction(Species, p2.rank)), color = NA, alpha = 0.1) +
        #facet_grid(cols = vars(predictor.class2), rows = vars(Predictor)) +
        facet_wrap(~Predictor)+
        labs(
          x = paste(pred.name$Predictor, "(", unique(df.species$Units.converted),")"),
          y = "10-year Mortality Probability",
          #title = "Effect of Predictors on Probability of Mortality"
        ) + 
        species_color + species_fill+
        scale_size_manual(values = c("low"=0.25 ,
                                     "median" = 1.2,
                                     "high"= 2.2 ), 
                          name = "Diameter", 
                          labels = c("low" = "small", 
                                     "median" = "medium", 
                                     "high" = "large"))+
        scale_linetype_manual(values = c("low"="dashed" ,
                                         "median" = "solid",
                                         "high"= "longdash" ), 
                              name = "Diameter", 
                              labels = c("low" = "small", 
                                         "median" = "medium", 
                                         "high" = "large"))+
        
        
        theme_bw(base_size = 14)+theme(panel.grid = element_blank(), 
                                       legend.key.size = unit(1, "cm"), 
                                       strip.background = element_rect(fill = strip.fill)
        )+coord_cartesian(ylim = c(0, ymax.spp))
      
      
      
      
      
      
    }
  }
  return(p1)
}


unique(marginal_response_df_unscaled$Covariate)


# generate plots for variables with significant contribution to prop variance---


# Diameter
dia.marg.p <- plot_main_region_effects(species.group = 
                           c("hickory spp.", "balsam fir", "chestnut oak", "northern white-cedar", 
                             "northern red oak", "eastern hemlock", "white oak", "eastern white pine", 
                             "red maple", "white ash", "yellow poplar", "sugar maple", "red spruce"), 
                         covar = "DIA_scaled", 
                         ymax.spp = 0.15)

marginal_response_df_unscaled %>%
  filter(Covariate %in% "DIA_scaled") %>%
  filter( Raw.value <= 16 & Raw.value >= 15) %>%
  filter(Species %in% unique(marginal_response_df_unscaled$Species)) %>%
  mutate(ten.mort.mean = 1-(mean)^10) %>%
  select(Species, ten.mort.mean)

# Diameter_diff
dia.diff.marg.p <- plot_main_region_effects(species.group = 
                           unique(marginal_response_df_unscaled$Species), 
                         covar = "DIA_DIFF_scaled", 
                         ymax.spp = 0.25)

# marginal_response_df_unscaled %>% 
#   filter(Covariate %in% "DIA_DIFF_scaled") %>% 
#   filter( Raw.value <= 2 & Raw.value >= 1.8) %>% 
#   filter(Species %in% unique(marginal_response_df_unscaled$Species)) %>%
#   mutate(ten.mort.mean = 1-(mean)^10) %>% 
#   select(Species, ten.mort.mean)

# Diameter_diff x growth interaction
dia_dia.diff.marg.p <- plot_main_region_effects(species.group = 
                           c("northern white-cedar", "black cherry", 
                             "northern red oak", "yellow poplar", "red maple", "sugar maple", 
                             "white ash"), 
                         covar = "DIA_scaled_growth.int", 
                         ymax.spp = 0.075)+theme(legend.position = "none")


dia_dia.diff.balsam.fir.marg.p <- plot_main_region_effects(species.group = 
                                                  c("balsam fir"), 
                                                covar = "DIA_scaled_growth.int", 
                                                ymax.spp = 0.5)+theme(legend.position = "none")

# matmax scaled:
matmax.marg.p <- plot_main_region_effects(species.group = 
                           c("American beech", 
                             "northern red oak", "eastern white pine", 
                             "balsam fir", "red maple"), 
                         covar = "MATmax.scaled", 
                         ymax.spp = 0.15)

# matmax x dia differeces
# more mortality at higher temps
matmax_dia.diff.marg.p <- plot_main_region_effects(species.group = 
                           c( "northern white-cedar", 
                             "yellow-poplar"), 
                         covar = "MATmax.scaled_growth.int", 
                         ymax.spp = 0.075)+theme(legend.position = "none")


# more mortality at cooler temps
matmax_dia.diff.marg.p.cooler <- plot_main_region_effects(species.group = 
                           c("yellow birch", "red spruce", "sugar maple"), 
                         covar = "MATmax.scaled_growth.int", 
                         ymax.spp = 0.1)+theme(legend.position = "none")


tmax_dia.marg.p <- plot_main_region_effects(species.group = 
                           c("hickory spp."), 
                         covar = "tmax.anom_DIA.int", 
                         ymax.spp = 0.025)+theme(legend.position = "none")

# percent damage
pct.damage.marg.p <- plot_main_region_effects(species.group = 
                           c("northern white-cedar", 
                             "white oak","yellow birch", "red spruce", 
                             "balsam fir", "chestnut oak", "northern red oak", "yellow birch", "red maple"), 
                         covar = "damage.scaled", 
                         ymax.spp = 0.075)

marginal_response_df_unscaled %>% 
  filter(Covariate %in% "damage.scaled") %>% 
  filter(Raw.value >= 20 & Raw.value <= 21) %>% 
  filter(Species %in% c("northern white-cedar",  "white oak","yellow birch", "red spruce")) %>%
  mutate(ten.mort.mean = 1-(mean)^10) %>% 
  select(Species, ten.mort.mean)



# Neighborhood 
ba.marg.p <- plot_main_region_effects(species.group = 
                           c("yellow birch"), 
                         covar = "ba.scaled", 
                         ymax.spp = 0.075)

bal.marg.p <- plot_main_region_effects(species.group = 
                           c("yellow birch", "red spruce", "American beech"), 
                         covar = "BAL.scaled", 
                         ymax.spp = 0.05)

ba.dia.marg.p <- plot_main_region_effects(species.group = 
                           c("white oak"), 
                         covar = "ba.scaled_DIA.int", 
                         ymax.spp = 0.015)+theme(legend.position = "none")


# N deposition

ndep.marg.p <-plot_main_region_effects(species.group = 
                           c("chestnut oak"), 
                         covar = "Ndep.scaled", 
                         ymax.spp = 0.05)

ndep.dia.diff.marg.p<- plot_main_region_effects(species.group = 
                           c("eastern white pine"), 
                         covar = "Ndep.scaled_growth.int", 
                         ymax.spp = 0.05)+theme(legend.position = "none")


# site conditions

slop.marg.p <- plot_main_region_effects(species.group = 
                           c("eastern hemlock", "chestnut oak"), 
                         covar = "slope.scaled", 
                         ymax.spp = 0.02)+
  xlab("Slope (%)")

# plot up in a big figure:

# figure 4:-----
cleaned.data$Species <- factor(cleaned.data$Species, disturb.species.order)
barplot.legend <- ggplot()+geom_histogram(data = cleaned.data, aes(fill = Species, y = tree))+species_fill


species.legend <- get_legend(barplot.legend)

figure4 <- plot_grid(plot_grid(dia.marg.p, 
                     dia.diff.marg.p, 
                     matmax.marg.p +xlab("Degrees C"), 
                     pct.damage.marg.p +xlab("% Damage"), 
          ba.marg.p,
          bal.marg.p, 
          slop.marg.p, 
          ndep.marg.p +xlab("N deposition decline \n (kg/m^2/year)"), 
          ncol = 4, align = "hv", labels = "AUTO"), 
          species.legend, ncol = 2, rel_widths = c(1, 0.2))
ggsave(paste0(output.dir,"images/Figure_4_marginal.png"), 
       plot = figure4, width = 12, height = 7, dpi = 300, bg = "white")


dia_plots_marg <- dia.marg.p + dia.diff.marg.p + dia_dia.diff.marg.p  + dia_dia.diff.balsam.fir.marg.p + plot_layout(ncol = 4)
ggsave(paste0(output.dir,"images/diameter_effects_marginal_top_var.png"), 
       plot = dia_plots_marg, width = 10, height = 4, dpi = 300)


#dia_plots_marg <- plot_grid(dia.marg.p, dia.diff.marg.p, dia_dia.diff.marg.p, dia_dia.diff.balsam.fir.marg.p, 
 #                           align = "hv", ncol = 4)



climate_plots_marg <- tmax_dia.marg.p + matmax.marg.p +  matmax_dia.diff.marg.p + matmax_dia.diff.marg.p.cooler  + plot_layout(ncol = 4)

ggsave(paste0(output.dir,"images/climate_effects_marginal_top_var.png"), 
       plot = climate_plots_marg, width = 10, height = 4, dpi = 300)

hood.vars.marg.p <- pct.damage.marg.p + ba.marg.p + bal.marg.p + ba.dia.marg.p + plot_layout(ncol = 4)

ggsave(paste0(output.dir,"images/neighborhood_effects_marginal_top_var.png"), 
       plot = hood.vars.marg.p, width = 10, height = 4, dpi = 300)

site.vars.marg.p <- slop.marg.p + ndep.marg.p + ndep.dia.diff.marg.p+ plot_layout(ncol = 4)

ggsave(paste0(output.dir,"images/site_effects_marginal_top_var.png"), 
       plot = site.vars.marg.p, width = 10, height = 4, dpi = 300)

# just the important main effects

main_plots_marg <- dia.marg.p + dia.diff.marg.p + matmax.marg.p + pct.damage.marg.p + ba.marg.p + bal.marg.p + slop.marg.p + ndep.marg.p + plot_layout(ncol = 4)
ggsave(paste0(output.dir,"images/main_effects_marginal_top_var.png"), 
       plot = main_plots_marg, width = 10, height = 8, dpi = 300)

ggsave(paste0(output.dir,"images/main_effects_marginal_top_var.svg"), 
       plot = main_plots_marg, width = 10, height = 8, dpi = 300)

# just the important interaction effects
interaction_plots_marg <- dia_dia.diff.marg.p + dia_dia.diff.balsam.fir.marg.p + matmax_dia.diff.marg.p + matmax_dia.diff.marg.p.cooler +
  tmax_dia.marg.p + ndep.dia.diff.marg.p + ba.dia.marg.p +plot_layout(ncol = 4)
ggsave(paste0(output.dir,"images/interaction_effects_marginal_top_var.png"), 
       plot = interaction_plots_marg, width = 10, height = 8, dpi = 300)

ggsave(paste0(output.dir,"images/interaction_effects_marginal_top_var.svg"), 
       plot = interaction_plots_marg, width = 10, height = 8, dpi = 300)
# just the important interaction effects:

# Plot up the important marginal effects for different groups of species---
# select main effects and interactions to plot for each species group:
# read in the variance partitioning summary by species
#var_summary <- readRDS(paste0(output.folder, "variance_partitioning_summary_region_species.RDS"))

# plot up all of the 




plot_main_region_effects(species.group = c("northern red oak", "hickory spp."), 
                         covar = c("tmax.anom_DIA.int"),
                         ymax.spp = 0.02)

plot_main_region_effects(species.group = c("chestnut oak"), 
                         covar = c("Ndep.scaled"),
                         ymax.spp = 0.1)

plot_main_region_effects(species.group = c("red maple", "yellow birch"), 
                         covar = c("MATmax.scaled_growth.int"),
                         ymax.spp = 0.25)

plot_main_region_effects(species.group = c("balsam fir"), 
                         covar = c("BAL.scaled_DIA.int"),
                         ymax.spp = 0.25)

plot_main_region_effects(species.group = c("balsam fir", "yellow-poplar", "yellow birch"), 
                         covar = c("ba.scaled_growth.int"),
                         ymax.spp = 0.5)


plot_main_region_effects(species.group = c( "paper birch"), 
                         covar = c("damage.scaled_DIA.int"),
                         ymax.spp = 0.05)

plot_main_region_effects(species.group = c("paper birch"), 
                         covar = c("damage.scaled_growth.int"),
                         ymax.spp = 0.25)

plot_main_region_effects(species.group = c("eastern white pine"), 
                         covar = c("Ndep.scaled_growth.int"),
                         ymax.spp = 0.09)
# spruce fir:----

# get all of the covariates that contribute to 1% or greater of the species
sprucefir.top5 <- var_summary %>% filter(Species %in% spruce.fir) %>% ungroup()%>%
  group_by( Covariate, predictor.class2, Species)%>%
  summarise(pct.region = median(mean_logit)*100)%>% arrange(Species, desc(pct.region)) %>% 
  group_by(Species)%>% #View()
  filter(pct.region >= 1 ) %>% ungroup() #%>% select(Covariate) %>% distinct()

# D X deltaD
# deltaD
# D
# Ndep x deltaD
# damage x deltaD
# damage
# matmax x D
# ba x D

region.plt.list <-list()
for(h in 1:length(unique(sprucefir.top5$Covariate))){
  region.plt.list[[h]] <- plot_main_region_effects(species.group = spruce.fir, 
                                                   covar = unique(sprucefir.top5$Covariate)[h], 
                                                   ymax.spp = 0.5)
}



# make separate legends
dia.diff.legend <- get_legend(region.plt.list[[3]]+scale_color_discrete(guide = "none")+scale_fill_discrete(guide = "none")+theme_bw(base_size = 18)+ theme(legend.key.width = unit(2,"cm")))
#diameter.legend <- get_legend(region.plt.list[[4]]+scale_color_discrete(guide = "none")+scale_fill_discrete(guide = "none")+theme_bw(base_size = 18)+ theme(legend.key.width = unit(2,"cm")))
species.legend <- get_legend(region.plt.list[[1]]+scale_linetype_discrete(guide = "none")+scale_size_discrete(guide = "none")+theme_bw(base_size = 18)+ theme(legend.key.width = unit(2,"cm")))

spruce.fir.effects <-plot_grid(
  plot_grid(region.plt.list[[1]]+theme(legend.position = "none"),
            region.plt.list[[2]]+theme(legend.position = "none"),
            region.plt.list[[5]]+theme(legend.position = "none"),
            
            
            ncol = 3, align = "hv"),
  
  plot_grid(region.plt.list[[3]]+theme(legend.position = "none"),
            region.plt.list[[4]]+theme(legend.position = "none"),
            region.plt.list[[6]]+theme(legend.position = "none"),
            ncol = 3, align = "hv"),
  plot_grid(species.legend, dia.diff.legend, ncol = 2),
  rel_heights = c(1,1,0.3),
  nrow = 4)

save_plot(paste0(output.folder,"images/spruce.fir_regional_marginal_effects_1pct.png"), 
          spruce.fir.effects, base_width = 16, base_height = 15) 

save_plot(paste0(output.folder,"images/spruce.fir_regional_marginal_effects_1pct.svg"),
          spruce.fir.effects, base_width = 16, base_height = 15)




# spongy moth susceptible:----
spongy.susceptible.top5 <- var_summary %>% filter(Species %in% spongy.susceptible) %>% ungroup()%>%
  group_by( Covariate, predictor.class2, Species)%>%
  summarise(pct.region = median(mean_logit)*100)%>% arrange(Species, desc(pct.region)) %>% 
  group_by(Species)%>% #View()
  filter(pct.region >= 1 ) %>% ungroup() #%>% select(Covariate) %>% distinct()

unique(spongy.susceptible.top5$Covariate)

region.plt.list <-list()
for(h in 1:length(unique(spongy.susceptible.top5$Covariate))){
  region.plt.list[[h]] <- plot_main_region_effects(species.group = spongy.susceptible, 
                                                   covar = unique(spongy.susceptible.top5$Covariate)[h], 
                                                   ymax.spp = 0.05)
}



# make separate legends
dia.diff.legend <- get_legend(region.plt.list[[6]]+scale_color_discrete(guide = "none")+scale_fill_discrete(guide = "none")+theme_bw(base_size = 18)+ theme(legend.key.width = unit(2,"cm")))
diameter.legend <- get_legend(region.plt.list[[8]]+scale_color_discrete(guide = "none")+scale_fill_discrete(guide = "none")+theme_bw(base_size = 18)+ theme(legend.key.width = unit(2,"cm")))
species.legend <- get_legend(region.plt.list[[1]]+scale_linetype_discrete(guide = "none")+scale_size_discrete(guide = "none")+theme_bw(base_size = 18)+ theme(legend.key.width = unit(2,"cm")))

spongy.susceptible.effects <-plot_grid(
  plot_grid(region.plt.list[[2]]+theme(legend.position = "none"),
            region.plt.list[[1]]+theme(legend.position = "none"),
            region.plt.list[[6]]+theme(legend.position = "none"),
            region.plt.list[[5]]+theme(legend.position = "none"),
            region.plt.list[[9]]+theme(legend.position = "none"),
            ncol = 6, align = "hv"),
  
  plot_grid(
    region.plt.list[[7]]+theme(legend.position = "none"),
    region.plt.list[[3]]+theme(legend.position = "none"),
    region.plt.list[[4]]+theme(legend.position = "none"),
    region.plt.list[[10]]+theme(legend.position = "none"),
    region.plt.list[[11]]+theme(legend.position = "none"),
    region.plt.list[[8]]+theme(legend.position = "none"),
    
    ncol = 6, align = "hv"),
  plot_grid(species.legend, diameter.legend,  dia.diff.legend, ncol = 2),
  rel_heights = c(1,1,0.3),
  nrow = 4)

save_plot(paste0(output.folder,"images/spongy.susceptible_dominant_marginal_variance_effects.png"), 
          spongy.susceptible.effects, base_width = 16, base_height = 15) 

save_plot(paste0(output.folder,"images/spongy.susceptible_dominant_marginal_variance_effects.svg"),
          spongy.susceptible.effects, base_width = 16, base_height = 15)




# spongy moth resistant:----

spongy.resist.top5 <- var_summary %>% filter(Species %in% spongy.resist) %>% ungroup()%>%
  group_by( Covariate, predictor.class2, Species)%>%
  summarise(pct.region = median(mean_logit)*100)%>% arrange(Species, desc(pct.region)) %>% 
  group_by(Species)%>% #View()
  filter(pct.region >= 1 ) %>% ungroup()

unique(spongy.resist.top5$Covariate)

region.plt.list <-list()
for(h in 1:length(unique(spongy.resist.top5$Covariate))){
  region.plt.list[[h]] <- plot_main_region_effects(species.group = spongy.resist, 
                                                   covar = unique(spongy.resist.top5$Covariate)[h], 
                                                   ymax.spp = 0.05)
}



# make separate legends
dia.diff.legend <- get_legend(region.plt.list[[4]]+scale_color_discrete(guide = "none")+scale_fill_discrete(guide = "none")+theme_bw(base_size = 18)+ theme(legend.key.width = unit(2,"cm")))
diameter.legend <- get_legend(region.plt.list[[5]]+scale_color_discrete(guide = "none")+scale_fill_discrete(guide = "none")+theme_bw(base_size = 18)+ theme(legend.key.width = unit(2,"cm")))
species.legend <- get_legend(region.plt.list[[1]]+scale_linetype_discrete(guide = "none")+scale_size_discrete(guide = "none")+theme_bw(base_size = 18)+ theme(legend.key.width = unit(2,"cm")))

spongy.resist.effects <-plot_grid(
  plot_grid(region.plt.list[[2]]+theme(legend.position = "none"),
            region.plt.list[[1]]+theme(legend.position = "none"),
            region.plt.list[[6]]+theme(legend.position = "none"),
            region.plt.list[[4]]+theme(legend.position = "none"),
            
            
            ncol = 4, align = "hv"),
  
  plot_grid(
    
    region.plt.list[[3]]+theme(legend.position = "none"),
    region.plt.list[[7]]+theme(legend.position = "none"),
    region.plt.list[[5]]+theme(legend.position = "none"),
    
    
    ncol = 4, align = "hv"),
  plot_grid(species.legend, diameter.legend, dia.diff.legend, ncol = 3),
  rel_heights = c(1,1,0.3),
  nrow = 4)

save_plot(paste0(output.folder,"images/spongy.resist_dominant_marginal_variance_effects.png"), 
          spongy.resist.effects, base_width = 16, base_height = 15) 

save_plot(paste0(output.folder,"images/spongy.resist_dominant_marginal_variance_effects.svg"),
          spongy.resist.effects, base_width = 16, base_height = 15)

# spongy moth immune:----

spongy.immune.top5 <- var_summary %>% filter(Species %in% spongy.immune) %>% ungroup()%>%
  group_by( Covariate, predictor.class2, Species)%>%
  summarise(pct.region = median(mean_logit)*100)%>% arrange(Species, desc(pct.region)) %>% 
  group_by(Species)%>% #View()
  filter(pct.region >= 1 ) %>% ungroup()

unique(spongy.immune.top5$Covariate)

region.plt.list <-list()
for(h in 1:length(unique(spongy.immune.top5$Covariate))){
  region.plt.list[[h]] <- plot_main_region_effects(species.group = spongy.immune, 
                                                   covar = unique(spongy.immune.top5$Covariate)[h], 
                                                   ymax.spp = 0.05)
}



# make separate legends
dia.diff.legend <- get_legend(region.plt.list[[2]]+scale_color_discrete(guide = "none")+scale_fill_discrete(guide = "none")+theme_bw(base_size = 18)+ theme(legend.key.width = unit(2,"cm")))
#diameter.legend <- get_legend(region.plt.list[[10]]+scale_color_discrete(guide = "none")+scale_fill_discrete(guide = "none")+theme_bw(base_size = 18)+ theme(legend.key.width = unit(2,"cm")))
species.legend <- get_legend(region.plt.list[[1]]+scale_linetype_discrete(guide = "none")+scale_size_discrete(guide = "none")+theme_bw(base_size = 18)+ theme(legend.key.width = unit(2,"cm")))

spongy.immune.effects <-plot_grid(
  plot_grid(region.plt.list[[3]]+theme(legend.position = "none"),
            region.plt.list[[1]]+theme(legend.position = "none"),
            region.plt.list[[2]]+theme(legend.position = "none"),
            region.plt.list[[4]]+theme(legend.position = "none"),
            
            
            ncol = 4, align = "hv"),
  
  
  
  
  plot_grid(species.legend,  dia.diff.legend, ncol = 3),
  rel_heights = c(1,0.3),
  nrow = 2)

save_plot(paste0(output.folder,"images/spongy.immune_dominant_marginal_variance_effects.png"), 
          spongy.immune.effects, base_width = 16, base_height = 9) 

save_plot(paste0(output.folder,"images/spongy.immune_dominant_marginal_variance_effects.svg"),
          spongy.immune.effects, base_width = 16, base_height = 9)

# hemlock, beech, and :----
mixed.top5 <- var_summary %>% filter(Species %in% mixed) %>% ungroup()%>%
  group_by( Covariate, predictor.class2, Species)%>%
  summarise(pct.region = median(mean_logit)*100)%>% arrange(Species, desc(pct.region)) %>% 
  group_by(Species)%>% #View()
  filter(pct.region >= 1 ) %>% ungroup() #%>% select(Covariate) %>% distinct()


region.plt.list <-list()
for(h in 1:length(unique(mixed.top5$Covariate))){
  region.plt.list[[h]] <- plot_main_region_effects(species.group = mixed, 
                                                   covar = unique(mixed.top5$Covariate)[h], 
                                                   ymax.spp = 0.15)
}


region.plt.list[[3]]

# make separate legends
#dia.diff.legend <- get_legend(region.plt.list[[8]]+scale_color_discrete(guide = "none")+scale_fill_discrete(guide = "none")+theme_bw(base_size = 18)+ theme(legend.key.width = unit(2,"cm")))
#diameter.legend <- get_legend(region.plt.list[[9]]+scale_color_discrete(guide = "none")+scale_fill_discrete(guide = "none")+theme_bw(base_size = 18)+ theme(legend.key.width = unit(2,"cm")))
species.legend <- get_legend(region.plt.list[[1]]+scale_linetype_discrete(guide = "none")+scale_size_discrete(guide = "none")+theme_bw(base_size = 18)+ theme(legend.key.width = unit(2,"cm")))

mixed.effects <- plot_grid(
  plot_grid(region.plt.list[[3]]+theme(legend.position = "none"),
            region.plt.list[[1]]+theme(legend.position = "none"),
            region.plt.list[[2]]+theme(legend.position = "none"),
            region.plt.list[[4]]+theme(legend.position = "none"),
            
            region.plt.list[[5]]+theme(legend.position = "none"),
            
            
            ncol = 5, align = "hv"),
  
  plot_grid(species.legend, ncol = 1),
  rel_heights = c(1,0.3),
  nrow = 2)

save_plot(paste0(output.folder,"images/mixed_dominant_marginal_variance_effects.png"), 
          mixed.effects, base_width = 16, base_height = 9) 

save_plot(paste0(output.folder,"images/mixed_dominant_marginal_variance_effects.svg"),
          mixed.effects, base_width = 16, base_height = 9)





## Select variance partitioning across the region:---------------------
# read in the state-level variance partitioning:
st.var.files <- list.files(paste0(output.folder,"SPCD_stanoutput_joint_v3/predicted_mort/"), pattern = "variance_partitioning_summary_by_predictor_state_", full.names = TRUE)
last_three_str_sub <- str_sub(st.var.files, -3, -1)
Covariate_names <- read.csv(paste0(output.dir, "data/model_covariate_types_v2.csv"))

st.var.all.df <- do.call(rbind, lapply(st.var.files, function(x){
  read.csv(x)%>% mutate(state = str_sub(x, -6, -5)) %>% rename("Covariate" = "predictor")%>% left_join(., Covariate_names)
}))

state.variance.df <- st.var.all.df %>% left_join(.,state.df %>% mutate(state = as.character(state)))
state.variance.df
#var.summary.region <- var_summary
state.variance.df$COMMON <- factor(state.variance.df$Species, 
                                   levels = rev(c("balsam fir", "red spruce", "northern white-cedar", 
                                                  "eastern hemlock", "American beech", 
                                                  "black oak", "chestnut oak", "northern red oak", "white oak", "yellow birch", "paper birch", 
                                                  "hickory spp.", "eastern white pine", "red maple", "sugar maple", 
                                                  "black cherry", "white ash", "yellow-poplar")))
state.variance.df$predictor.class <- factor(state.variance.df$predictor.class, 
                                            levels = c(
                                              
                                              "Size" ,
                                              "Change in Size" ,
                                              "Climate" ,
                                              "Site Conditions",
                                              "Competition",
                                              "N deposition", 
                                              "% Damage", 
                                              "Growth x Size",
                                              
                                              "Site x G & S",
                                              "Competition x G & S" ,
                                              
                                              "Ndep x G & S",
                                              "Climate x G & S",
                                              "Damage x G & S"))


state.variance.df$predictor.class2 <- factor(state.variance.df$predictor.class2, 
                                             levels = c(
                                               
                                               "Size",
                                               "Change in Size",
                                               "Growth x Size" ,
                                               
                                               "Climate",
                                               "Climate x G & S",
                                               
                                               "N deposition", 
                                               "Ndep x G & S" ,
                                               
                                               "% Damage", 
                                               "Damage x G & S",
                                               
                                               
                                               "Competition",
                                               "Competition x G & S",
                                               
                                               "Site Conditions",
                                               "Site x G & S"
                                               
                                               
                                             ))



spruce.fir.state <- ggplot(data = state.variance.df %>% filter(COMMON %in% spruce.fir))+
  geom_bar(aes(x = region, y = mean, fill = predictor.class2), stat = "identity", position = "stack")+
  theme_bw(base_size = 16)+ theme(axis.text.x = element_text(angle = 60, hjust = 1), 
                                  panel.background = element_blank())+
  ylab("Across-tree variance in  \n p(survival) explained")+
  scale_fill_manual(values = color.pred.class.2, name = "",
                    labels = c("Change in Size" = expression(Delta ~ "D" ), 
                               "Growth x Size" = expression(Delta ~ "D x D"), 
                               "Size" = "Diameter (D)", 
                               "Climate x G & S" = expression("Climate x " ~Delta ~ "D or D"), 
                               "Ndep x G & S" = expression("N dep. x " ~Delta ~ "D or D"),
                               "Damage x G & S" = expression("Damage x " ~Delta ~ "D or D"), 
                               "Competition x G & S" = expression("Compeition x " ~Delta ~ "D or D"),
                               "Site x G & S" = expression("Site x " ~Delta ~ "D or D")
                    ),
                    guide = guide_legend(direction = "horizontal",
                                         ncol = 14,nrow = 1,reverse = TRUE,
                                         label.position="top", label.hjust = 0,
                                         label.vjust = 0.5,
                                         label.theme = element_text(angle = 90)))+
  coord_flip()+
  xlab("")+
  theme(panel.grid = element_blank(), 
        legend.position = "top", 
        panel.border = element_blank(), 
        axis.ticks.y = element_blank())+
  facet_wrap(~COMMON, scales = "free_y", ncol = 1)

save_plot(paste0(output.folder,"images/spruce_fir_state_variance_effects.png"), 
          spruce.fir.state, base_width = 16, base_height = 15) 

save_plot(paste0(output.folder,"images/spruce_fir_state_variance_effects.svg"),
          spruce.fir.state, base_width = 16, base_height = 15)


spongy.susceptible.state <- ggplot(data = state.variance.df %>% filter(COMMON %in% spongy.susceptible))+
  geom_bar(aes(x = region, y = mean, fill = predictor.class2), stat = "identity", position = "stack")+
  theme_bw(base_size = 16)+ theme(axis.text.x = element_text(angle = 60, hjust = 1), 
                                  panel.background = element_blank())+
  ylab("Across-tree variance in  \n p(survival) explained")+
  scale_fill_manual(values = color.pred.class.2, name = "",
                    labels = c("Change in Size" = expression(Delta ~ "D" ), 
                               "Growth x Size" = expression(Delta ~ "D x D"), 
                               "Size" = "Diameter (D)", 
                               "Climate x G & S" = expression("Climate x " ~Delta ~ "D or D"), 
                               "Ndep x G & S" = expression("N dep. x " ~Delta ~ "D or D"),
                               "Damage x G & S" = expression("Damage x " ~Delta ~ "D or D"), 
                               "Competition x G & S" = expression("Compeition x " ~Delta ~ "D or D"),
                               "Site x G & S" = expression("Site x " ~Delta ~ "D or D")
                    ),
                    guide = guide_legend(direction = "horizontal",
                                         ncol = 14,nrow = 1,reverse = TRUE,
                                         label.position="top", label.hjust = 0,
                                         label.vjust = 0.5,
                                         label.theme = element_text(angle = 90)))+
  coord_flip()+
  xlab("")+
  theme(panel.grid = element_blank(), 
        legend.position = "right", 
        panel.border = element_blank(), 
        axis.ticks.y = element_blank())+
  facet_wrap(~COMMON, scales = "free_y", ncol = 1)

save_plot(paste0(output.folder,"images/spongy_susceptible_state_variance_effects.png"), 
          spongy.susceptible.state , base_width = 16, base_height = 15) 

save_plot(paste0(output.folder,"images/spongy_susceptible_state_variance_effects.svg"),
          spongy.susceptible.state , base_width = 16, base_height = 15)



beech.hemlock.pine.pop.st <-  ggplot(data = state.variance.df %>% filter(COMMON %in% c("American beech", "eastern hemlock", "yellow-poplar", "eastern white pine")))+
  geom_bar(aes(x = region, y = mean, fill = predictor.class2), stat = "identity", position = "stack")+
  theme_bw(base_size = 16)+ theme(axis.text.x = element_text(angle = 60, hjust = 1), 
                                  panel.background = element_blank())+
  ylab("Across-tree variance in  \n p(survival) explained")+
  scale_fill_manual(values = color.pred.class.2, name = "",
                    labels = c("Change in Size" = expression(Delta ~ "D" ), 
                               "Growth x Size" = expression(Delta ~ "D x D"), 
                               "Size" = "Diameter (D)", 
                               "Climate x G & S" = expression("Climate x " ~Delta ~ "D or D"), 
                               "Ndep x G & S" = expression("N dep. x " ~Delta ~ "D or D"),
                               "Damage x G & S" = expression("Damage x " ~Delta ~ "D or D"), 
                               "Competition x G & S" = expression("Compeition x " ~Delta ~ "D or D"),
                               "Site x G & S" = expression("Site x " ~Delta ~ "D or D")
                    ),
                    guide = guide_legend(direction = "horizontal",
                                         ncol = 14,nrow = 1,reverse = TRUE,
                                         label.position="top", label.hjust = 0,
                                         label.vjust = 0.5,
                                         label.theme = element_text(angle = 90)))+
  coord_flip()+
  xlab("")+
  theme(panel.grid = element_blank(), 
        legend.position = "right", 
        panel.border = element_blank(), 
        axis.ticks.y = element_blank())+
  facet_wrap(~COMMON, scales = "free_y", ncol = 1)

save_plot(paste0(output.folder,"images/mixed_other_state_variance_effects.png"), 
          beech.hemlock.pine.pop.st , base_width = 16, base_height = 15) 

save_plot(paste0(output.folder,"images/mixed_other_state_variance_effects.svg"),
          beech.hemlock.pine.pop.st , base_width = 16, base_height = 15)

state.variance.df %>% filter(COMMON %in% spruce.fir)%>%
  mutate(logit_percent = round(mean_logit*100, digits = 2))%>%
  group_by(predictor.class2, COMMON, region, Covariate)%>%
  summarise(total_cat = sum(mean_logit)*100)%>%
  #select(Covariate, logit_percent, total_cat,  predictor.class2, COMMON, region) %>%
  filter(COMMON %in% "red spruce" )%>%
  filter(total_cat > 1) %>% 
  arrange(region, desc(total_cat))%>% View()

# plot up regional values against the mortality rates
train.data %>% filter(Species %in%  "red spruce") %>% 
  group_by(state) %>% 
  summarise(median(MATmax))

train.data.values <- train.data %>% left_join(., state.df)
ggplot(train.data.values %>% filter(Species %in% "red spruce"), aes(x = region, y = MATmax))+geom_boxplot()
ggplot(train.data.values %>% filter(Species %in% "red spruce"), aes(x = region, y = Difference_per_yr))+geom_boxplot()

ggplot(train.data.values %>% filter(Species %in% spruce.fir), aes(x = region, y = damage.total, fill = Species))+geom_boxplot()
ggplot(train.data.values %>% filter(Species %in% spruce.fir), 
       aes(x = region, y = Difference_per_yr, fill = Species))+geom_boxplot()

ggplot(train.data.values %>% filter(Species %in% spruce.fir), 
       aes(x = region, y = MATmax, fill = Species))+geom_boxplot()
ggplot(train.data.values %>% filter(Species %in% spruce.fir), 
       aes(x = region, y = plt_ba_sq_ft_old, fill = Species))+geom_boxplot()

ggplot(train.data.values %>% filter(Species %in% spongy.susceptible), aes(x = region, y = damage.total, fill = Species))+geom_boxplot()


ggplot(train.data.values %>% filter(Species %in% spruce.fir), aes(x = region, y = damage.total, fill = Species))+geom_boxplot()
ggplot(train.data.values %>% filter(Species %in% spongy.susceptible), aes(x = region, y = damage.total, fill = Species))+geom_boxplot()
