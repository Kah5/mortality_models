### 
# Code to make species distribution maps of the eastern tree species
# ##
# Author: Kelly Heilman
# source of species distribution shapefiles: https://github.com/wpetry/USTreeAtlas
library(sf)
library(ggplot2)
library(tidyverse)
library(here)
library(FIESTA)
library(maps)
library(mapdata)
library(ggrepel)
library(tigris)
library(ggalluvial)


################################################################################
# Read in mortality data for 17 species
################################################################################
cleaned.data <- readRDS( "data/cleaned.data.mortality.TRplots.RDS")

cleaned.data <- cleaned.data %>% filter(!is.na(ba) & !is.na(slope) & ! is.na(physio) & !is.na(aspect))%>% 
  dplyr::select(state, county, pltnum, cndtn, point, tree, PLOT.ID, date,cycle, spp, dbhcur, dbhold, status, damage, Species, SPCD,
                remper, LAT_FIADB, LONG_FIADB, elev, DIA_DIFF, annual.growth, M, relative.growth, si, physio:RD) %>% distinct()
# get summary of damages for later use:
N.DAMAGE <- cleaned.data %>% group_by(SPCD, damage) %>% summarise(n.by.damage = n())
N.DAMAGE$SPECIES <- ref_species[match(N.DAMAGE$SPCD, ref_species$SPCD),]$COMMON
ref_damage<- ref_codes %>% filter(VARIABLE %in% "AGENTCD")
N.DAMAGE$damage_agent <- ref_damage[match(N.DAMAGE$damage, ref_damage$VALUE),]$MEANING
N.DAMAGE$damage_agent <- ifelse(N.DAMAGE$damage == 0, "None", N.DAMAGE$damage_agent)
#saveRDS(N.DAMAGE, "data/N.DAMAGE.table.RDS")


nspp <- cleaned.data %>% group_by(SPCD) %>% summarise(n = n(), 
                                                      pct = n/nrow(cleaned.data)) %>% arrange (desc(`pct`))

nspp$cumulative.pct <- cumsum(nspp$pct)



# link up to the species table:
nspp$COMMON <- FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD),]$COMMON
nspp$Species <- paste(FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD),]$GENUS, FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD),]$SPECIES)

#View(nspp)

nspp[1:17,]$COMMON

library(gt)
nspp[1:17,] %>% mutate(pct = round(pct, 3), 
                       cumulative.pct = round(cumulative.pct, 3)) %>% rename(`# of trees` = "n", 
                                                                             `% of trees` = "pct",
                                                                             `cumulative %` = "cumulative.pct", 
                                                                             `Common name` = "COMMON") %>%
  dplyr::select(Species, `Common name`, SPCD, `# of trees`, `% of trees`, `cumulative %`)|> gt()



cleaned.data$SPGRPCD <- FIESTA::ref_species[match(cleaned.data$SPCD, FIESTA::ref_species$SPCD),]$E_SPGRPCD

SPGRP.df <- FIESTA::ref_codes %>% filter(VARIABLE %in% "SPGRPCD") %>% filter(VALUE %in% unique(cleaned.data$SPGRPCD))
cleaned.data$SPGRPNAME <- SPGRP.df[match(cleaned.data$SPGRPCD, SPGRP.df$VALUE),]$MEANING
cleaned.data$STNAME <- ref_statecd[match(cleaned.data$state, ref_statecd$VALUE),]$MEANING


# set up custom colors for species
# set the species order using the factors:
SP.TRAITS <- read.csv("data/NinemetsSpeciesTraits.csv") %>% filter(COMMON_NAME %in% c(unique(nspp[1:17,]$COMMON), "sweet birch", "Virginia pine"))
# order the trait db by softwood-hardwood, then shade tolerance, then name (this puts all the oaks together b/c hickory and red oak have the same tolerance values)
SP.TRAITS <- SP.TRAITS %>% group_by(SFTWD_HRDWD) %>% arrange(desc(SFTWD_HRDWD), desc(ShadeTol), desc(COMMON_NAME))

SP.TRAITS$Color <- c(# softwoods
  "#b2df8a", # balsam fir
  "#003c30", # hemlock
  "#b2182b", # red spruce
  "#fee090", # n white cedar
  "#33a02c", # eastern white pine
  
  # virginia pine
  "#1a9850",
  
  
  # sugar  maples
  "#a6cee3",
  "#1f78b4",# american beech
  
  # red maple
  "#e31a1c",
  # yellow birch
  "#fdbf6f",
  
  
  # oaks
  "#54278f", # white oak
  "#cab2d6", # chestnut oak
  "#8073ac", # n red oak
  "#6a3d9a",# black oak
  # sweet birch
  "#fccde5",
  
  # hickory
  "#7f3b08",
  # white ash
  "#bababa",
  # black cherry
  "#4d4d4d",
  # yellow poplar
  "#ff7f00"
  
  
  
  
)

SP.TRAITS$`Shade Tolerance`  <- ifelse(SP.TRAITS$ShadeTol >=4, "High", 
                                       ifelse(SP.TRAITS$ShadeTol <=2.5, "Low", "Moderate"))


sppColors <- c( "#f5f5f5", "#d9f0d3","darkgrey", SP.TRAITS$Color )
names(sppColors) <- c("other hardwood", "other softwood","dead", unique(SP.TRAITS$COMMON_NAME))

species_fill <- scale_fill_manual(values = sppColors)
species_color <- scale_color_manual(values = sppColors)




## calculate compositoin from number of trees and tree BA by state
# cleaned.data only has the 17 species of interest, get the full records here:
TREE.remeas <- readRDS("data/unfiltered_TREE.remeas.rds")




################################################################################
# Actual filtering of all the data--
###############################################################################
# caution! this includes logged tree & all species too just for exploration!
TREE.data.all <- TREE.remeas %>% 
  # for each tree caclulate annual growth if dbhold is not NA and the remper is > 0
  group_by(PLOT.ID, point, state, county, pltnum, tree, date) %>% 
  
  # for dead trees (dead, non-salvable dead, and snags), calculate annual growth as diameter diff/ half of remper
  # for live trees, calculate annual growht as diamber diff/remper
  mutate(annual.growth = ifelse(status %in% c(2,  4, 5) & ! is.na(dbhold) & ! remper == 0, 
                                DIA_DIFF/(remper/2), DIA_DIFF/remper),
         M = ifelse(status %in%  c(2, 4, 5), 1, 
                    ifelse(status %in% c(3), 3, 0)), # if a tree is dead, non-salvable deead or a snag, mark as dead
         relative.growth = (annual.growth/dbhold)*100) %>% ungroup() %>% 
  # acutal filtering section
  filter( exprem > 0 & # if exprem == 0, these could be modeled plots?
            dbhold >= 5 & # need an initial dbh greater than 5
            ! remper == 0 & # if remper is listed as zero, filter out
            DIA_DIFF >= 0 ) %>% #& # filter out diameter differences >= 0
            #!status == 3 & # remove cut trees
            #SPCD %in% nspp[1:17,]$SPCD)%>% # filter species in the top 17 of all species
  mutate(Tree.status = ifelse(M == 1, "dead", ifelse(M == 3, "cut", "live")), # add some labels to ID all dead trees
         DIAMETER_diff = ifelse(DIA_DIFF > 0, "positive", 
                                ifelse(DIA_DIFF == 0,"zero", "negative")))
TREE.data.all$STNAME <- ref_statecd[match(TREE.data.all$state, ref_statecd$VALUE),]$MEANING
TREE.data.all$SW_HW <- SP.TRAITS[match(TREE.data.all$SPCD, SP.TRAITS$SPCD),]$SFTWD_HRDWD
TREE.data.all %>% filter(is.na(SW_HW))
# get the top 17 species and make a table
nspp <-  TREE.data.all %>% group_by(SPCD) %>% 
    summarise(n = n(), 
              pct = n/nrow(TREE.data.all)) %>% arrange (desc(`pct`))
nspp$cumulative.pct <- cumsum(nspp$pct)
nspp$COMMON <- FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD),]$COMMON
nspp$Species <- paste(FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD),]$GENUS, FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD),]$SPECIES)

  
nspp %>% mutate(pct = round(pct, 3), 
         cumulative.pct = round(cumulative.pct, 3)) %>% 
  rename(`# of trees` = "n", 
         `% of trees` = "pct",
         `cumulative %` = "cumulative.pct", 
        `Common name` = "COMMON") %>%
  dplyr::select(Species, `Common name`, SPCD, `# of trees`, `% of trees`, `cumulative %`)|> gt()



species.composition <- TREE.data.all %>%
  mutate(Top.Species = ifelse(SPCD %in% nspp[1:17,]$SPCD, Species, 
                              ifelse(SW_HW %in% "S", "other softwood","other hardwood"))) %>% # label for non-top species
  mutate(dbhold.cm = dbhold*2.54) %>%
  mutate(tree_ba_cm = (pi*dbhold.cm^2)/4, 
         tree_ba_m = (pi*(dbhold.cm/100)^2)/4) %>% #select(dbhold, dbhold.cm, tree_ba_cm, tree_ba_m)
  group_by(state, STNAME) %>% mutate(total.trees.st = n(), 
                             total.ba.st_sq_m = sum(tree_ba_m, na.rm =TRUE)) %>%
  group_by(state, STNAME, Top.Species, total.trees.st, total.ba.st_sq_m) %>%
  summarise(Total.trees.sp = n(), 
            Total.ba.sp_sq_m = sum(tree_ba_m, na.rm =TRUE)) %>%
  mutate(`Percent Composition (Density)` = (Total.trees.sp/total.trees.st)*100, 
         `Percent Composition (BA)` = (Total.ba.sp_sq_m/total.ba.st_sq_m)*100)

species.composition.nspp17 <- TREE.data.all %>% filter(SPCD %in% nspp[1:17,]$SPCD)%>%
  mutate(Top.Species = ifelse(SPCD %in% nspp[1:17,]$SPCD, Species, 
                              ifelse(SW_HW %in% "S", "other softwood","other hardwood"))) %>% # label for non-top species
  mutate(dbhold.cm = dbhold*2.54) %>%
  mutate(tree_ba_cm = (pi*dbhold.cm^2)/4, 
         tree_ba_m = (pi*(dbhold.cm/100)^2)/4) %>% #select(dbhold, dbhold.cm, tree_ba_cm, tree_ba_m)
  group_by(state, STNAME) %>% mutate(total.trees.st = n(), 
                                     total.ba.st_sq_m = sum(tree_ba_m, na.rm =TRUE)) %>%
  group_by(state, STNAME, Top.Species, total.trees.st, total.ba.st_sq_m) %>%
  summarise(Total.trees.sp = n(), 
            Total.ba.sp_sq_m = sum(tree_ba_m, na.rm =TRUE)) %>%
  mutate(`Percent Composition (Density)` = (Total.trees.sp/total.trees.st)*100, 
         `Percent Composition (BA)` = (Total.ba.sp_sq_m/total.ba.st_sq_m)*100)

nspp %>% filter(COMMON %in% c("sweet birch", "Virginia pine")) 

species.composition.nsppcored.mort.17 <- TREE.data.all %>% 
  filter(SPCD %in% c(nspp[1:17,]$SPCD, 372,132))%>%
  mutate(Top.Species = ifelse(SPCD %in% c(nspp[1:17,]$SPCD, 372,132), Species, 
                              ifelse(SW_HW %in% "S", "other softwood","other hardwood"))) %>% # label for non-top species
  mutate(dbhold.cm = dbhold*2.54) %>%
  mutate(tree_ba_cm = (pi*dbhold.cm^2)/4, 
         tree_ba_m = (pi*(dbhold.cm/100)^2)/4) %>% #select(dbhold, dbhold.cm, tree_ba_cm, tree_ba_m)
  group_by(state, STNAME) %>% mutate(total.trees.st = n(), 
                                     total.ba.st_sq_m = sum(tree_ba_m, na.rm =TRUE)) %>%
  group_by(state, STNAME, Top.Species, total.trees.st, total.ba.st_sq_m) %>%
  summarise(Total.trees.sp = n(), 
            Total.ba.sp_sq_m = sum(tree_ba_m, na.rm =TRUE)) %>%
  mutate(`Percent Composition (Density)` = (Total.trees.sp/total.trees.st)*100, 
         `Percent Composition (BA)` = (Total.ba.sp_sq_m/total.ba.st_sq_m)*100)


# order the species by shade tolerence
species.composition$Top.Species <- factor(species.composition$Top.Species, levels = c("other softwood", "other hardwood",unique(SP.TRAITS$COMMON_NAME)) )
species.composition$STNAME <- factor(species.composition$STNAME, levels = c("Maine", "New Hampshire", "Vermont", "New York", 
                                                                                 "Connecticut", "New Jersey", "Pennsylvania", "Ohio", "Maryland","West Virginia"))

species.composition.nspp17$Top.Species <- factor(species.composition.nspp17$Top.Species, levels = c(unique(SP.TRAITS$COMMON_NAME)) )
species.composition.nspp17$STNAME <- factor(species.composition.nspp17$STNAME, levels = c("Maine", "New Hampshire", "Vermont", "New York", 
                                                                            "Connecticut", "New Jersey", "Pennsylvania", "Ohio", "Maryland","West Virginia"))

species.composition.nsppcored.mort.17$Top.Species <- factor(species.composition.nsppcored.mort.17$Top.Species, levels = c(unique(SP.TRAITS$COMMON_NAME)) )
species.composition.nsppcored.mort.17$STNAME <- factor(species.composition.nsppcored.mort.17$STNAME, levels = c("Maine", "New Hampshire", "Vermont", "New York", 
                                                                                          "Connecticut", "New Jersey", "Pennsylvania", "Ohio", "Maryland","West Virginia"))

species.composition.nsppcored.mort.17

# plot up compostion charts for all the trees, including other hardwood
ggplot() + 
  geom_bar(data = species.composition, aes(x = STNAME, y = `Percent Composition (Density)`, group = Top.Species, fill = Top.Species), position = "stack", stat = "identity")+
  theme_bw(base_size = 12)+
  xlab("State")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1), 
        panel.grid = element_blank(), legend.title = element_blank())+species_fill
ggsave(height = 5, width = 8, units = "in", dpi = 300, "images/all_state_composition_by_density.png")

st.comp.ba <- ggplot() + 
  geom_bar(data = species.composition %>% filter(!is.na(Top.Species)), aes(x = STNAME, y = `Percent Composition (BA)`, group = Top.Species, fill = Top.Species), position = "stack", stat = "identity")+
  theme_bw(base_size = 14)+
  xlab("State")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1), 
        panel.grid = element_blank(), legend.title = element_blank())+species_fill
ggsave(plot = st.comp.ba, height = 6, width = 8, units = "in", dpi = 300, paste0(output.dir, "images/all_state_composition_by_BA.png"))


# make pie charts:
ggplot() + 
  geom_bar(data = species.composition, aes(x = factor(1), y = `Percent Composition (Density)`, fill = Top.Species), stat = "identity", color = "black", width = 1)+
  coord_polar(theta = "y")+
  facet_wrap(~STNAME)+
  theme_bw(base_size = 14)+
  #theme_void() +
  xlab("State")+
  theme(
         legend.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank())+species_fill+
  labs(x="", y="")
ggsave(height = 5, width = 8, units = "in", dpi = 300, "images/all_state_composition_by_density_pie.png")


# plot up figures of compostion for all the species looked at here
# plot up compostion charts for all the trees, including other hardwood
ggplot() + 
  geom_bar(data = species.composition.nspp17, aes(x = STNAME, y = `Percent Composition (Density)`, group = Top.Species, fill = Top.Species), position = "stack", stat = "identity")+
  theme_bw(base_size = 14)+
  xlab("State")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1), 
        panel.grid = element_blank(), legend.title = element_blank())+species_fill
ggsave(height = 5, width = 8, units = "in", dpi = 300, "images/all_state_composition_by_density_nspp17.png")

ggplot() + 
  geom_bar(data = species.composition.nspp17, aes(x = STNAME, y = `Percent Composition (BA)`, group = Top.Species, fill = Top.Species), position = "stack", stat = "identity")+
  theme_bw(base_size = 14)+
  xlab("State")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1), 
        panel.grid = element_blank(), legend.title = element_blank())+species_fill
ggsave(height = 5, width = 8, units = "in", dpi = 300, "images/all_state_composition_by_BA_nspp17.png")


# make pie charts:
ggplot() + 
  geom_bar(data = species.composition.nspp17, aes(x = factor(1), y = `Percent Composition (Density)`, fill = Top.Species), stat = "identity", color = "black", width = 1)+
  coord_polar(theta = "y")+
  facet_wrap(~STNAME)+
  theme_bw(base_size = 14)+
  #theme_void() +
  xlab("State")+
  theme(
    legend.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid  = element_blank())+species_fill+
  labs(x="", y="")
ggsave(height = 5, width = 8, units = "in", dpi = 300, "images/all_state_composition_by_density_pie_nspp17.png")

######################################################################################
# create a sankey diagram showing tree mortality by species ----
# using ggalluvial


species.mort.summary <- TREE.data.all %>%
  mutate(Top.Species = ifelse(SPCD %in% c(nspp[1:17,]$SPCD, 372,132), Species, 
                              ifelse(SW_HW %in% "S", "other softwood","other hardwood"))) %>% # label for non-top species
  group_by(Top.Species, Tree.status) %>%
  mutate(
    species_t3 = ifelse(Tree.status == "live", Top.Species, NA)  # same species if alive
  )
# Prepare alluvial plot data
tree_alluvial <- species.mort.summary %>%
  mutate(
    T1 = Top.Species,
    T2 = Tree.status,
    T3 = ifelse(Tree.status == "live", species_t3, 
                ifelse(Tree.status == "cut","cut","dead"))
  ) %>%
  
  count(T1, T2, T3, name = "count")
tree_alluvial$T2 <- factor(tree_alluvial$T2, levels = c("live", "cut", "dead"))

tree_alluvial$T1 <- factor(tree_alluvial$T1, levels = c(unique(SP.TRAITS$COMMON_NAME),"other softwood", "other hardwood") )
tree_alluvial$T_pct <- nspp[match(tree_alluvial$T1, nspp$COMMON),]$pct*100
tree_alluvial$T_pct <- ifelse(tree_alluvial$T1 %in% "other hardwood", sum(nspp[18:nrow(nspp),]$pct)*100, tree_alluvial$T_pct) # get the total values for the non-focal species
tree_alluvial$T_pct <- paste0(round(tree_alluvial$T_pct, digits = 2), "%")
tree.shade.order.pct <- tree_alluvial %>% ungroup()%>% select(T_pct, T1) %>% distinct()%>% arrange(by = T1)
tree_alluvial$T_pct <- factor(tree_alluvial$T_pct, levels =c(tree.shade.order.pct$T_pct))
tree.counts <- tree_alluvial %>% ungroup()%>% select(T1, T_pct, count) %>% group_by(T1, T_pct) %>% summarise(total = sum(count)) 
#"linear", "cubic", "quintic", "sine", "arctangent", and "sigmoid". "xspline"


# Plot up the data in a sankey diagram
ggplot(data = tree_alluvial ,
       aes(axis1 = T1, axis2 = T3,  y = count)) +
  geom_alluvium(data = tree_alluvial ,aes(fill = T1), width = 1/12, curve_type = "sigmoid") +
  geom_stratum(data = tree_alluvial ,aes(fill = T1),width = 1/3, color = "black") +
  #geom_label_repel(aes(x = 0.8, label = T_pct, color = T1),  show.legend = FALSE,  position = "stack") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(limits = c("Species Composition at T1", "Species Composition at T2"),
                   expand = c(.05, .05)) +
  labs(title = "Tree Species Composition and Survival Over Time",
       y = "Number of Trees") +
  theme_minimal()+species_fill + species_color#+
  theme(axis.title.y = element_blank(), axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        panel.grid = element_blank(),
        legend.positon = "none")
ggsave(height = 10, width = 7, dpi = 350, "images/Compostion_sankey.png")


ggplot(tree_alluvial %>% filter(T1 %in% nspp[1:17,]$COMMON),
       aes(axis1 = T1, axis2 = T3,  y = count)) +
  geom_alluvium(aes(fill = T1), width = 1/12, curve_type = "sigmoid") +
  geom_stratum(aes(fill = T1),width = 1/3, color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(limits = c("Species Composition at T1", "Species Composition at T2"),
                   expand = c(.05, .05)) +
  labs(title = "Tree Species Composition and Survival Over Time",
       y = "Number of Trees") +
  theme_minimal()+species_fill + species_color+
  theme(axis.title.y = element_blank(), axis.ticks = element_blank(), 
        axis.text.y = element_blank(), 
        panel.grid = element_blank(), 
        legend.positon = "none")
ggsave(height = 10, width = 7, dpi = 350, "images/Compostion_sankey_nspp17.png")


#add percentage to the sankey diagram
# to do this we need to use the lodes formate of the data

tree_alluvial_long <- to_lodes_form(tree_alluvial %>% ungroup()%>% select(T1, T3, T_pct, count), axes = 1:2,
                              key = "variable", value = "value", id = "cohort",
                              diffuse = T1)



alluvial.plt <- ggplot(
  data =tree_alluvial_long %>% filter(!is.na(T1)),
  aes(x = variable, stratum = value, alluvium = cohort, y = count)
) +
  geom_alluvium(aes(fill=T1), width = 1/12, curve_type = "sigmoid") +
  geom_stratum(aes(fill=T1), width = 1/3, color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  theme_minimal() +species_fill+
 
  geom_label(stat = "stratum",
                  aes( fill = T1,
                      label = ifelse(variable == "T1" & !value %in% c(
                        "cut","dead"), as.character(T_pct), NA)),
                  direction = "y", nudge_x = -0.3, min.segment.length = Inf, show.legend = FALSE)+
  # scale_x_discrete(limits = c("Species Composition at T1", "Species Composition at T2"),
  #                  expand = c(.05, .05)) +
  # labs(title = "Tree Species Composition and Survival Over Time",
  #      y = "Number of Trees") +
  theme_minimal()+species_fill + 
  theme(axis.title.y = element_blank(), axis.ticks = element_blank(), 
        axis.text.y = element_blank(), 
        panel.grid = element_blank(), 
        legend.positon = "none")+
  guides(fill=guide_legend(title=NULL))
ggsave(plot = alluvial.plt, height = 10, width = 10, dpi = 350, "images/Compostion_sankey_percentages.png")

############################################################################


cleaned.data.full <- cleaned.data

mortality.summary <- cleaned.data.full %>% filter(SPCD %in% unique(nspp[1:17,]$SPCD)) %>%
  group_by(SPCD, M) %>% summarise(ntrees = n()) %>%
  spread(M, ntrees) %>% rename("live" = `0`, 
                               "dead" = `1`)
mortality.summary$`Common name` <- FIESTA::ref_species[match(mortality.summary$SPCD, FIESTA::ref_species$SPCD),]$COMMON
mortality.summary %>% dplyr::select(`Common name`, SPCD, live, dead) %>% mutate(total = live + dead) %>% ungroup() |> gt()
colnames(cleaned.data.full)

cleaned.data.17 <- cleaned.data %>% dplyr::select(-spp) %>% filter(SPCD %in% unique(nspp[1:17,]$SPCD))

#####################################################################################################state######################################################################################################################
# Plot the species compositition ----
#####################################################################################################state######################################################################################################################



################################################################################
# For each species get the combined species name that matches the little maps----
################################################################################
plt.all <- left_join(cleaned.data.17, ref_species) 
plt.shp <- plt.all %>% mutate(GENUS = tolower(GENUS), 
                              shp.code = paste0(substr(GENUS, 1, 4), substr(SPECIES, 1,4)))

# fix some exceptions to the 4 letter rule:
# caryspp. is separated into multiple carya for the maps (maybe overlay each cayra species on the map later)
# Pinsus strobis is pinustrb instead of pinustrob (pinus strobiformis*)
# sugar maple is acersacr instead of acersacc (silver maple)


plt.shp$shp.code <- ifelse(plt.shp$SCIENTIFIC_NAME %in% "Pinus strobus", "pinustrb", plt.shp$shp.code)
plt.shp$shp.code <- ifelse(plt.shp$SCIENTIFIC_NAME %in% "Acer saccharum", "acersacr", plt.shp$shp.code)


output.dir <- "C:/Users/KellyHeilman/Box/01. kelly.heilman Workspace/mortality/Eastern-Mortality/mortality_models/"


################################################################################
# MAPPING ---figure 1 and supplemental figures
################################################################################
# get the states that we need:
states <- map_data("state")
#9=CT, 25=MA, 33=NH, 23=ME, 50=VT, 44=RI, 42=PA, 39=OH, 54=WV
state_sub <- filter(states, region %in% c("connecticut","maine","new hampshire","vermont","new york", "new jersey",
                                          "rhode island","pennsylvania","ohio","west virginia", "massachusetts", "virginia", "delaware", 
                                          "north carolina", "kentucky", "tennessee", "michigan", "indiana", "district of columbia", "south carolina", 
                                          "georgia", "maryland"))

canada <- map_data("worldHires", "Canada")
spp <- "pinustrb"

# map out all the species on the plot:
ggplot() +
  geom_polygon(data = canada, 
               aes(x=long, y=lat, group = group), 
               color = "black", fill = "white") +
  geom_polygon(data = state_sub, 
               aes(x=long, y=lat, group = group), 
               color = "black", fill = "white") +
  #geom_sf(alpha = 0.75, aes(fill = as.character(Distribution)))+
  #scale_fill_manual(values = c("Species Distribution" = "forestgreen", "Outside Distribution" = "white"))+
  geom_jitter(data = plt.shp, aes(x = LONG_FIADB, y = LAT_FIADB, color = COMMON_NAME), size = 0.5)+theme_bw()+
  coord_sf(xlim = c(-85, -68), ylim = c(35, 48))+theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'lightblue'), 
                                                       #legend.position = c(0.87, 0.25),
                                                       legend.title = element_blank(),
                                                       legend.background = element_rect(fill = "white", color = "black"))+
  scale_color_manual(values = c(
    
    # maples
    "red maple" ="black", 
    "sugar maple" = "#878787", 
    # oaks
    "northern red oak" = "#a6cee3", 
    "white oak" = "#cab2d6", 
    "chestnut oak" = "#1f78b4", 
    "black oak" = "#8dd3c7",
    
    # other decidious
    "hickory spp."= "#6a3d9a",
    "American beech" = "#b15928", 
    "black cherry" = "#c51b7d", 
    "white ash" = "#f1b6da", 
    "yellow-poplar" = "#c2a5cf", 
    "yellow birch" = "#ffff99", 
    
    # conifers
    "red spruce"= "#e31a1c", 
    "balsam fir" = "#ff7f00", 
    "northern white-cedar" = "#fdbf6f", 
    "eastern white pine" = "#66c2a4", 
    "eastern hemlock" = "#00441b"
    
    
  ))

plt.pct.mortality <- plt.shp %>% group_by(state, county, pltnum, LAT_FIADB, LONG_FIADB, M) %>% 
  summarise(ntrees = n()) %>% spread(M, ntrees) %>%
  mutate(mort.ratio = ifelse(is.na(`1`),0,
                             ifelse(is.na(`0`), 1, `1`/(`0`+`1`)))) %>%
  mutate(pct.mortality = ifelse(mort.ratio > 1, 100, mort.ratio*100))

plt.pct.mortality %>% filter(pct.mortality > 100)
hist(plt.pct.mortality$pct.mortality)
ggplot(data = plt.pct.mortality, aes(x = pct.mortality))+geom_histogram()+theme_bw(base_size = 12)+xlab("% mortality on plot")
ggsave(height = 2, width = 3, units = "in", paste0(output.dir, "images/species_comp/all_observed_pct_mortality_hist.png"))

nonzero <- plt.pct.mortality %>% filter(pct.mortality> 0)
quantile(nonzero$pct.mortality)

plt.pct.mortality <- plt.pct.mortality %>% mutate(Mort.quantiles = cut(pct.mortality, 
                                                       breaks = c(0, 1, 2.5, 5, 10, 20, 30, 40, 50, 100), 
                                                       include.lowest=TRUE))

unique(plt.pct.mortality$Mort.quantiles)
large.size.font <- ggplot() +
  geom_polygon(data = canada, 
               aes(x=long, y=lat, group = group), 
               color = "black", fill = "white") +
  geom_polygon(data = state_sub, 
               aes(x=long, y=lat, group = group), 
               color = "black", fill = "white") +
  #geom_sf(alpha = 0.75, aes(fill = as.character(Distribution)))+
  #scale_fill_manual(values = c("Species Distribution" = "forestgreen", "Outside Distribution" = "white"))+
  geom_jitter(data = plt.pct.mortality, aes(x = LONG_FIADB, y = LAT_FIADB, color = Mort.quantiles), size = 5)+theme_bw()+
  coord_sf(xlim = c(-85, -68), ylim = c(35, 48))+theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'lightblue'), 
                                                       #legend.position = c(0.87, 0.25),
                                                       legend.title = element_blank(),
                                                       legend.background = element_rect(fill = "white", color = "black"))+
  scale_color_manual(values = c(
    "[0,1]" ="darkgrey",
    "(1,2.5]" = "#fcc5c0", 
    "(2.5,5]" = "#fa9fb5", 
    "(5,10]" = "#feb24c", 
    "(10,20]" = "#fd8d3c", 
    "(20,30]" = "#fc4e2a", 
    "(30,40]" = "#e31a1c", 
    "(40,50]" =  "#bd0026", 
    "(50,100]" = "#49006a" 
    
    
  ))+ guides(fill=guide_legend(title="% Mortality"))


legend.large <- cowplot::get_legend(large.size.font)

library(rnaturalearth)
library(terra)

nat.earth <- terra::rast(paste0(output.dir,"data/NaturalEarth/NE1_50M_SR_W/NE1_50M_SR_W/NE1_50M_SR_W.tif"))
domain <- ext(-87, -65, 36, 48)

nat.crop <- crop(nat.earth, y=domain)
rast.table <- as.data.frame(nat.crop, xy = TRUE)


raster_data <- as.data.frame(nat.crop, xy = TRUE) 

# Combine relevant layers using rgb (you may need to adjust band indices)
raster_data$rgb <- rgb(raster_data[, 3], raster_data[, 4], raster_data[, 5], maxColorValue = 255) # Assuming RGB bands


lakes = ne_download(scale = 50, type = 'lakes', category = 'physical', returnclass = "sf")

# Convert sf object to a data frame with geometry as coordinates
lakes_df <- lakes %>%
  st_as_sf() %>%
  st_cast("POLYGON") %>% # Ensure geometry is in polygon format
  st_coordinates() %>%
  as.data.frame()

# Add grouping variables for polygons
lakes_df <- lakes_df %>%
  rename(long = X, lat = Y) %>%
  mutate(group = interaction(L1, L2)) # Grouping for polygons

# Plot using geom_polygon
ggplot(nc_df, aes(x = long, y = lat, group = group)) +
  geom_polygon(fill = "lightblue", color = "black") +
  theme_minimal()+coord_sf(xlim = c(-85.5, -67.5), ylim = c(37, 47.5))


# create a basemap from natural earth raster:
NE_basemap <- ggplot(data = raster_data, aes(x = x, y = y)) +
  geom_tile(fill = raster_data$rgb, alpha = 0.75) +
  #geom_sf(data = lakes, fill = "lightblue") + theme_bw()+
  coord_sf(xlim = c(-85.5, -67.5), ylim = c(37, 47.5))+
  theme_void() 


mortality_percent_map <- NE_basemap + 
  geom_polygon(data = canada, 
               aes(x=long, y=lat, group = group), 
               color = "black", fill = NA) +
  geom_polygon(data = lakes_df, 
               aes(x = long, y = lat, group = group), 
               color = "black", fill = "lightblue") +
 
  
  #geom_sf(alpha = 0.75, aes(fill = as.character(Distribution)))+
  #scale_fill_manual(values = c("Species Distribution" = "forestgreen", "Outside Distribution" = "white"))+
  geom_jitter(data = plt.pct.mortality, aes(x = LONG_FIADB, y = LAT_FIADB, color = Mort.quantiles), size = 0.5)+theme_bw()+
  geom_polygon(data = state_sub, 
               aes(x=long, y=lat, group = group), 
               color = "black", fill = NA) +
  coord_sf(xlim = c(-85, -67.5), ylim = c(37, 47.5))+theme(panel.grid = element_blank(), #panel.background = element_rect(fill = 'lightblue'), 
                                                       #legend.position = c(0.87, 0.25),
                                                       axis.title  = element_blank(),
                                                       legend.title = element_blank(),
                                                       legend.background = element_rect(fill = "white", color = "black"), 
                                                       legend.position= "none")+
  scale_color_manual(values = c(
    "[0,1]" ="darkgrey",
    "(1,2.5]" = "#fcc5c0", 
    "(2.5,5]" = "#fa9fb5", 
    "(5,10]" = "#feb24c", 
    "(10,20]" = "#fd8d3c", 
    "(20,30]" = "#fc4e2a", 
    "(30,40]" = "#e31a1c", 
    "(40,50]" =  "#bd0026", 
    "(50,100]" = "#49006a" 
   

  ))

png(height = 5, width = 8.5, res = 350, units = "in", paste0(output.dir, "images/species_comp/all_observed_pct_mortality_map.png") )
cowplot::plot_grid(legend.large, mortality_percent_map,  align = "v", rel_widths = c(0.15, 0.89))
dev.off()





##########################################################################
# function to map out the species presence for the plots -------
plot.cored.spp.distn <- function(spp){
  dir.path = "C:/Users/KellyHeilman/OneDrive - USDA/Documents/USTreeAtlas/shp/" # where the distribution shape files are
  
  
  if(file.exists(paste0(dir.path, spp, "/"))){
    spp.distribution <- st_read(paste0(dir.path, spp, "/"))
    
    plt.spp <- plt.shp %>% filter(shp.code %in% spp)
    
    # plot distribution 
    spp.distribution %>% mutate(Distribution = ifelse(CODE == 1, "Species Distribution", "Outside Distribution"))%>%
      ggplot() +
      geom_polygon(data = canada, 
                   aes(x=long, y=lat, group = group), 
                   color = "black", fill = "white") +
      geom_polygon(data = state_sub, 
                   aes(x=long, y=lat, group = group), 
                   color = "black", fill = "white") +
      geom_sf(alpha = 0.75, aes(fill = as.character(Distribution)))+
      scale_fill_manual(values = c("Species Distribution" = "forestgreen", "Outside Distribution" = "white"))+
      geom_point(data = plt.spp, aes(x = LONG_FIADB, y = LAT_FIADB), size = 0.75)+theme_bw()+ggtitle(unique(plt.spp$COMMON_NAME))+
      coord_sf(xlim = c(-85, -68), ylim = c(35, 48))+theme(panel.grid = element_blank(), legend.position = "none", panel.background = element_rect(fill = 'lightblue'))
    
    ggsave(height = 6, width = 8, units = "in", here("images/species_comp/", paste0(spp, "_distribution_maps.png")))
  }
  else{
    cat(paste0("no range map exists for ", spp))
  }
}







# color by % live and dead for that taxa
plot.mort.spp.distn <- function(spp){
  dir.path = "C:/Users/KellyHeilman/OneDrive - USDA/Documents/USTreeAtlas/shp/" # where the distribution shape files are
  
  
  if(file.exists(paste0(dir.path, spp, "/"))){
    spp.distribution <- st_read(paste0(dir.path, spp, "/"))
    
    plt.spp <- plt.shp %>% filter(shp.code %in% spp)
    
    plt.spp <-  plt.spp %>% group_by(PLOT.ID, LONG_FIADB, LAT_FIADB, COMMON_NAME, M) %>% summarise(ntrees = n()) %>% ungroup()%>%
                group_by(PLOT.ID, LONG_FIADB, LAT_FIADB, COMMON_NAME) %>% spread(M, ntrees) %>% mutate(fraction.dead = ifelse(!is.na(`1`) & !is.na(`0`),`1`/(`0`+`1`), 
                                                                                                                    ifelse(is.na(`1`), 0, 
                                                                                                                                     ifelse(is.na(`0`), 1 ))))
    plt.spp$pct.mort <- ifelse(plt.spp$fraction.dead <= 0.1, "< 0.1",
                                      ifelse(plt.spp$fraction.dead >0.1 & plt.spp$fraction.dead <= 0.49, "0.1 - 0.49", 
                               ifelse(plt.spp$fraction.dead >0.49 & plt.spp$fraction.dead <0.9, "0.5 - 0.89", "> 0.9")))
    
    plt.spp$pct.mort <-factor(plt.spp$pct.mort, levels = c("< 0.1", "0.1 - 0.49", "0.5 - 0.89", "> 0.9"))
    # plot distribution 
   spp.map <-  spp.distribution %>% mutate(Distribution = ifelse(CODE == 1, "Species Distribution", "Outside Distribution"))%>%
      ggplot() +
      geom_polygon(data = canada, 
                   aes(x=long, y=lat, group = group), 
                   color = "black", fill = "white") +
      geom_polygon(data = state_sub, 
                   aes(x=long, y=lat, group = group), 
                   color = "black", fill = "white") +
      geom_sf(alpha = 0.75, aes(fill = Distribution))+
      scale_fill_manual(values = c("Species Distribution" = "forestgreen", "Outside Distribution" = "white"))+
      geom_point(data = plt.spp, aes(x = LONG_FIADB, y = LAT_FIADB, color = pct.mort), size = 0.5)+
      scale_color_manual(values = c("< 0.1" = "black", "0.1 - 0.49" = "#fdcc8a", "0.5 - 0.89" = "#fc8d59", "> 0.9" = "#d7301f"))+
      theme_bw()+ggtitle(unique(plt.spp$COMMON_NAME))+
      coord_sf(xlim = c(-85, -68), ylim = c(35, 48))+theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'lightblue'))
    spp.map 
    ggsave(height = 6, width = 8, units = "in", here("images/species_comp/", paste0(spp, "mortality_pct_distribution_maps.png")))
  spp.map
    }
  else{
    cat(paste0("no range map exists for ", spp))
  }
}

################################################################################
# Apply to all the species except hickory
################################################################################
# map out with presence only 
lapply(unique(plt.shp$shp.code), plot.cored.spp.distn)

# map out colored by mortality
species.maps <- lapply(unique(plt.shp$shp.code), plot.mort.spp.distn)
species.maps.legend <- cowplot::get_legend(species.maps[[1]])

png(height = 12, width = 14, units = "in", res = 150, "images/species_comp/all_spp_mortality.png")
cowplot::plot_grid(cowplot::plot_grid(species.maps[[1]] + theme(legend.position = "none"), 
                   species.maps[[2]]+ theme(legend.position = "none"), 
                   species.maps[[3]]+ theme(legend.position = "none"), 
                   species.maps[[4]]+ theme(legend.position = "none"), 
                   species.maps[[5]]+ theme(legend.position = "none"), 
                   species.maps[[7]]+ theme(legend.position = "none"), 
                   species.maps[[8]]+ theme(legend.position = "none"), 
                   species.maps[[9]]+ theme(legend.position = "none"), 
                   species.maps[[10]]+ theme(legend.position = "none"), 
                   species.maps[[11]]+ theme(legend.position = "none"), 
                   species.maps[[12]]+ theme(legend.position = "none"), 
                   species.maps[[13]]+ theme(legend.position = "none"), 
                   species.maps[[14]]+ theme(legend.position = "none"), 
                   species.maps[[15]]+ theme(legend.position = "none"), 
                   species.maps[[16]]+ theme(legend.position = "none"), 
                   #species.maps[[17]]+ theme(legend.position = "none"), 
                   ncol = 4), 
                  species.maps.legend, rel_widths = c(0.85, 0.1), ncol = 2)
dev.off()
################################################################################
# Make a unique hickory map with multiple species distributions
################################################################################
# get all the hickory shapefiles:
dir.path = "C:/Users/KellyHeilman/OneDrive - USDA/Documents/USTreeAtlas/shp/" 
all.hickory <- list.files(path = dir.path, pattern = "cary")

# list of species in NE based on the TREEs of NA book
# "caryaqua" # not in range
# "carycord" # yes
# "caryflor" # not in range
# "caryglab" # yes 
# "caryilli" # borderline--maybe in ohio
# "carylaci" # yes
# "carymyri" # out of range
# "caryovat" # yes
# "carypall" # borderline -- maybe in WV
# "carytexa" # out of range
# "carytome" # yes

in.range.hickory <- all.hickory[all.hickory %in% c("carycord", "caryglab","carylaci", "caryovat", "carytome", "carypall", "caryilli")]

spp.distribution <- list()
for(i in 1:length(in.range.hickory)){
  spp.distribution[[i]] <- st_read(paste0(dir.path, in.range.hickory[i], "/"))
  spp.distribution[[i]]$`Hickory Species` <- ifelse(spp.distribution[[i]]$CODE == 1, in.range.hickory[i], "out of range")
}

# 1,3 is not in the NE
plt.spp <- plt.shp %>% filter(shp.code %in% "caryspp.")
# Carya in range

#spp.distribution[[4]] %>% 
 hickory.presence.map <-  ggplot() +
  geom_polygon(data = canada, 
               aes(x=long, y=lat, group = group), 
               color = "black", fill = "white") +
  geom_polygon(data = state_sub, 
               aes(x=long, y=lat, group = group), 
               color = "black", fill = "white") +
 
  geom_sf(data = spp.distribution[[2]], alpha = 0.5, aes(fill = `Hickory Species`))+
  geom_sf(data = spp.distribution[[3]], alpha = 0.5, aes(fill = `Hickory Species`))+
  geom_sf(data = spp.distribution[[4]], alpha = 0.5, aes(fill = `Hickory Species`))+
   geom_sf(data = spp.distribution[[5]], alpha = 0.5, aes(fill = `Hickory Species`))+
   geom_sf(data = spp.distribution[[6]], alpha = 0.5, aes(fill = `Hickory Species`))+
   geom_sf(data = spp.distribution[[7]], alpha = 0.5, aes(fill = `Hickory Species`)) +
  geom_sf(data = spp.distribution[[1]], alpha = 0.5, aes(fill = `Hickory Species`))+
  scale_fill_manual(values = c("carycord" = "forestgreen", 
                               "caryglab" = "#beaed4",
                               "caryilli" = "#fdc086",
                               "carylaci" = "#ffff99", 
                               "caryovat" = "#386cb0", 
                               "carypall" = "#f0027f",
                               "carytome" = "#bf5b17",
                               "out of range" = "white"))+
    geom_point(data = plt.spp, aes(x = LONG_FIADB, y = LAT_FIADB), size = 0.75)+
    #scale_color_manual(values = c("< 0.1" = "black", "0.1 - 0.49" = "#fdcc8a", "0.5 - 0.89" = "#fc8d59", "> 0.9" = "#d7301f"))+
    theme_bw()+ggtitle(unique(plt.spp$COMMON_NAME))+
    coord_sf(xlim = c(-85, -68), ylim = c(35, 48))+theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'lightblue'))
 hickory.presence.map
  ggsave(height = 6, width = 8, units = "in", here("images/species_comp/", paste0(spp, "_distribution_maps.png")))
  
# make the map with percent mortality
  plt.spp <-  plt.spp %>% group_by(PLOT.ID, LONG_FIADB, LAT_FIADB, COMMON_NAME, M) %>% summarise(ntrees = n()) %>% ungroup()%>%
    group_by(PLOT.ID, LONG_FIADB, LAT_FIADB, COMMON_NAME) %>% spread(M, ntrees) %>% mutate(fraction.dead = ifelse(!is.na(`1`) & !is.na(`0`),`1`/(`0`+`1`), 
                                                                                                                  ifelse(is.na(`1`), 0, 
                                                                                                                         ifelse(is.na(`0`), 1 ))))
  plt.spp$pct.mort <- ifelse(plt.spp$fraction.dead <= 0.1, "< 0.1",
                             ifelse(plt.spp$fraction.dead >0.1 & plt.spp$fraction.dead <= 0.49, "0.1 - 0.49", 
                                    ifelse(plt.spp$fraction.dead >0.49 & plt.spp$fraction.dead <0.9, "0.5 - 0.89", "> 0.9")))
  
  plt.spp$pct.mort <-factor(plt.spp$pct.mort, levels = c("< 0.1", "0.1 - 0.49", "0.5 - 0.89", "> 0.9"))
  
 all.hickory.spp.map <-  ggplot() +
    geom_polygon(data = canada, 
                 aes(x=long, y=lat, group = group), 
                 color = "black", fill = "white") +
    geom_polygon(data = state_sub, 
                 aes(x=long, y=lat, group = group), 
                 color = "black", fill = "white") +
    geom_sf(data = spp.distribution[[1]], alpha = 0.5, aes(fill = `Hickory Species`))+
    geom_sf(data = spp.distribution[[2]], alpha = 0.5, aes(fill = `Hickory Species`))+
    geom_sf(data = spp.distribution[[3]], alpha = 0.5, aes(fill = `Hickory Species`))+
    geom_sf(data = spp.distribution[[4]], alpha = 0.5, aes(fill = `Hickory Species`))+
    geom_sf(data = spp.distribution[[5]], alpha = 0.5, aes(fill = `Hickory Species`))+
    geom_sf(data = spp.distribution[[6]], alpha = 0.5, aes(fill = `Hickory Species`))+
    geom_sf(data = spp.distribution[[7]], alpha = 0.5, aes(fill = `Hickory Species`)) + 
    scale_fill_manual(values = c( 
                                 "caryglab" = "#beaed4",
                                 "caryilli" = "#fdc086",
                                 "carylaci" = "#ffff99",
                                 "carycord" = "forestgreen",
                                 "caryovat" = "#386cb0", 
                                 "carypall" = "#f0027f",
                                 "carytome" = "#bf5b17",
                                 "out of range" = "white"))+
    geom_point(data = plt.spp, aes(x = LONG_FIADB, y = LAT_FIADB, color = pct.mort), size = 0.75)+
  scale_color_manual(values = c("< 0.1" = "black", "0.1 - 0.49" = "#fdcc8a", "0.5 - 0.89" = "#fc8d59", "> 0.9" = "#d7301f"))+
    theme_bw()+ggtitle(unique(plt.spp$COMMON_NAME))+
    coord_sf(xlim = c(-85, -68), ylim = c(35, 48))+theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'lightblue'))
 all.hickory.spp.map 
  ggsave(height = 6, width = 8, units = "in", here("images/species_comp/", paste0(spp, "pct_mort_distribution_maps.png")))
  

# To do: Make one big or two big figures of species distributions with overlay
  hickory.legend <- cowplot::get_legend(hickory.presence.map)
  png(height = 15, width = 12, units = "in", res = 150, "images/species_comp/all_spp_mortality2.png")
  cowplot::plot_grid(cowplot::plot_grid(species.maps[[15]]+ theme(legend.position = "none")+ theme(axis.title = element_blank()),
                                        species.maps[[16]]+ theme(legend.position = "none")+ theme(axis.title = element_blank()),
                                        #species.maps[[17]]+ theme(legend.position = "none")+ theme(axis.title = element_blank()),
                                        species.maps[[14]]+ theme(legend.position = "none")+ theme(axis.title = element_blank()),
                                        species.maps[[5]]+ theme(legend.position = "none")+ theme(axis.title = element_blank()),
                                        species.maps[[1]] + theme(legend.position = "none") + theme(axis.title  = element_blank()), 
                                        species.maps[[11]]+ theme(legend.position = "none")+ theme(axis.title = element_blank()),
                                        
                                        species.maps[[13]]+ theme(legend.position = "none")+ theme(axis.title = element_blank()),
                                        species.maps[[3]]+ theme(legend.position = "none") + theme(axis.title = element_blank()),
                                        species.maps[[7]]+ theme(legend.position = "none")+ theme(axis.title = element_blank()),
                                        species.maps[[8]]+ theme(legend.position = "none")+ theme(axis.title = element_blank()),
                                        species.maps[[10]]+ theme(legend.position = "none")+ theme(axis.title = element_blank()),
                                        
                                        species.maps[[2]]+ theme(legend.position = "none") + theme(axis.title = element_blank()),
                                         species.maps[[4]]+ theme(legend.position = "none")+ theme(axis.title = element_blank()),
                                         species.maps[[9]]+ theme(legend.position = "none")+ theme(axis.title = element_blank()),
                                        
                                        species.maps[[12]]+ theme(legend.position = "none")+ theme(axis.title = element_blank()),
                                        
                                        
                                      all.hickory.spp.map + theme(legend.position = "none")+ theme(axis.title = element_blank()),
                                      hickory.legend, ncol = 3), 
                     species.maps.legend, rel_widths = c(0.85, 0.15), ncol = 2)
  dev.off()
  
#########################################################################################
# Generate covariate data figures for Figure 1 ----
#########################################################################################

# plot PRISM climate data for each location----
#all.monthly <- readRDS("C:/Users/KellyHeilman/Box/01. kelly.heilman Workspace/mortality/Eastern-Mortality/mortality_models/data/PRISM_monthly_NE_periodic.RDS" )

all.seas <- readRDS("C:/Users/KellyHeilman/Box/01. kelly.heilman Workspace/mortality/Eastern-Mortality/mortality_models/data/seasonal_climate_periodic_NE_plots.RDS" )
# 
# all.seas
# join up to get the statecode information
trimmed.data <- cleaned.data %>% dplyr::select(PLOT.ID, state, remper, date) %>% distinct()
unique(trimmed.data$state)
# state.id <- 9
# state.trim <- trimmed.data %>% filter(state %in% state.id )

# just get the plots
climate.plots <- all.seas %>% filter( PLOT.ID %in% trimmed.data$PLOT.ID )#

state.summary <- climate.plots %>% group_by(lat, lon, PLOT.ID, year) %>% 
  select(-ID, -cn, -plot, -invyr, -plot_status_cd, -cycle, -elev) %>% distinct() %>%
  group_by(statecd, year) %>%
  summarise(mean_PPT = mean(yr_PPT, na.rm = TRUE), 
            mean_Tmax = mean(yr_MeanTmax, na.rm = TRUE))
state.summary$State <- ref_statecd[match(state.summary$statecd, ref_statecd$VALUE),]$MEANING

region.summary <- climate.plots %>% group_by(lat, lon, PLOT.ID, year) %>% 
  select(-ID, -cn, -plot, -invyr, -plot_status_cd, -cycle, -elev) %>% distinct() %>%
  group_by(year) %>%
  summarise(mean_PPT_all = mean(yr_PPT, na.rm = TRUE), 
            mean_Tmax_all = mean(yr_MeanTmax, na.rm = TRUE))

ggplot()+geom_line(data = state.summary, aes(x = year, y = mean_PPT, group = State), color = "lightblue")+
  geom_line(data = region.summary, aes(x = year, y = mean_PPT_all), color = "blue", linewidth = 1.1)+
  theme_minimal(base_size = 14)+ylab("Mean Precipitation (mm)")+xlab("Year")

ggplot()+geom_line(data = state.summary, aes(x = year, y = mean_Tmax, group = State, color = State), alpha = 0.5)+
  geom_line(data = region.summary, aes(x = year, y = mean_Tmax_all), color = "red", linewidth = 1.1)+
  theme_minimal(base_size = 14)+ylab("Average Annual Tmax (C)")+xlab("Year")+
  geom_label(data = state.summary %>% filter(year == last(year)), aes(label = State, 
                                                               x = year + 6, 
                                                               y = mean_Tmax, 
                                                               color = State))+
  theme(legend.position = "none")


plot.rempers <- cleaned.data %>% select(state, remper, date, PLOT.ID) %>% distinct() %>%
  group_by(state, remper, date) %>%
  summarise(n()) %>%
  mutate(T1 = date-remper, 
         T2 = date)
plotcommon.remper <- plot.rempers %>% group_by(state) %>% mutate(n.max = max(`n()`)) %>% 
  filter(`n()` == n.max) %>%
  rename("statecd" = "state")

state.summary.remper <- left_join(state.summary, plotcommon.remper) %>%
  mutate(`State Remeasurement` = ifelse(year >=T1 & year <= T2, 0.2, 0.1))
remper.vals <- state.summary.remper %>% filter(year %in% T1:T2) %>% group_by(statecd, State, T1, T2) %>% 
  summarise(PPT.val = mean(mean_PPT), 
            Tmax.val = mean(mean_Tmax))

# set up a general state color scale:
state.scales <- c("#FFAA00", 
                         "#d94801", 
                         "goldenrod", 
                         "#98D851",
                         "darkgreen",
                         "#01665e",
                         "#35978f", 
                         "#1d91c0",
                         "#225ea8" ,
                         "#253494",
                         "#081d58", 
                         
                         "#810f7c", 
                         "#4d004b")
names(state.scales) <- c("Ohio","Delaware","Vermont","Maryland","West Virginia",
                                "Maine","Rhode Island","New Jersey","Connecticut","New Hampshire",
                                "New York","Massachusetts","Pennsylvania" )




# get a general remper for the span of our periodic data:
T1.T2periodic <- state.summary.remper %>% select(State, T1, T2) %>% distinct() %>%
  ungroup()%>%
  summarise(T1.all = min(T1), 
            T2.all = max(T2))


state.locations <- state.summary.remper  %>% filter(year == last(year)) %>% arrange(mean_Tmax)
state.locations$Yloc <- c(11.9, 12.4, 13, 13.5, 15.0, 15.9, 16.5, 17.0, 17.5, 18.1)

ggplot()+
  geom_rect(data = T1.T2periodic, aes(xmin = T1.all, 
                                      xmax = T2.all, 
                                      ymin = -Inf, 
                                      ymax = Inf), alpha = 0.25)+
  geom_line(data = state.summary.remper , aes(x = year, y = mean_Tmax, group = State, color = State))+
  scale_linewidth(range = c(0.1, 1))+
  geom_line(data = region.summary, aes(x = year, y = mean_Tmax_all), color = "black", linewidth = 1.1)+
  theme_bw()+ylab("Average Annual Tmax (C)")+xlab("Year")+
  geom_label_repel(data = state.locations, aes(label = State, 
                                                                      x = year, 
                                                                      y = mean_Tmax, 
                                                                      color = State, 
                                               x_nudge = year + 6), 
                   direction = "y", box.padding = 0.5,  nudge_x = 2, min.segment.length = 1,
                   
                   #vjust = 0.5,
                   hjust = 0.5)+
  theme(legend.position = "none")+scale_color_manual(values = state.scales)+
  
  xlim(1900, 2038)+ ylim (8,19)
 
ggsave(height = 5, width = 8, units = "in", dpi = 300, "images/all_state_Tmax_time_series.png")


tave.ts.plt <- ggplot()+
  geom_rect(data = T1.T2periodic, aes(xmin = T1.all, 
                                      xmax = T2.all, 
                                      ymin = -Inf, 
                                      ymax = Inf), alpha = 0.25)+
  geom_line(data = state.summary.remper , aes(x = year, y = mean_Tmax, group = State, color = State))+
  scale_linewidth(range = c(0.1, 1))+
  geom_line(data = region.summary, aes(x = year, y = mean_Tmax_all), color = "black", linewidth = 1.1)+
  theme_bw()+ylab("Average Annual Tmax (C)")+xlab("Year")+
  geom_label_repel(data = state.locations, aes(label = State, 
                                               x = year, 
                                               y = mean_Tmax, 
                                               color = State, 
                                               x_nudge = year + 6), 
                   direction = "y", box.padding = 0.5,  nudge_x = 2, min.segment.length = 1,
                   
                   #vjust = 0.5,
                   hjust = 0.5)+
  theme(legend.position = "none")+scale_color_manual(values = state.scales)+
  
  xlim(1925, 2038)+ ylim (9,19)

ggsave(plot = tave.ts.plt, height = 5, width = 5, units = "in", dpi = 300, "images/all_state_Tmax_time_series_1925_2030.png")


# do the same with precipitation


state.locations.ppt <- state.summary.remper  %>%  filter(year == last(year)) %>% arrange(mean_PPT)
state.locations.ppt$yloc <- c(900, 950, 1000, 1050, 1100, 1150, 1200, 1260, 1350, 1400)

ggplot()+
  geom_rect(data = T1.T2periodic, aes(xmin = T1.all, 
                                      xmax = T2.all, 
                                      ymin = -Inf, 
                                      ymax = Inf), alpha = 0.25)+
  geom_line(data = state.summary.remper, aes(x = year, y = mean_PPT, group = State, color = State), alpha = 0.9)+
  scale_linewidth(range = c(0.1, 1))+
  geom_line(data = region.summary, aes(x = year, y = mean_PPT_all), color = "black", linewidth = 1.1)+
  theme_bw(base_size = 12)+ylab("Mean Precipitation (mm)")+xlab("Year")+
xlim(1900, 2038)+
  geom_label(data = state.locations.ppt, aes(label = State, 
                                         x = year + 8, 
                                         y = yloc, 
                                         color = State), 
             direction = "y", box.padding = 100, min.segment.length = 2)+
  scale_color_manual(values = state.scales)+
  theme(legend.position = "none")
ggsave(height = 5, width = 8, units = "in", dpi = 300, "images/all_state_PPT_time_series.png")


ppt.ts.plt <- ggplot()+
  geom_rect(data = T1.T2periodic, aes(xmin = T1.all, 
                                      xmax = T2.all, 
                                      ymin = -Inf, 
                                      ymax = Inf), alpha = 0.25)+
  geom_line(data = state.summary.remper, aes(x = year, y = mean_PPT, group = State, color = State), alpha = 0.9)+
  scale_linewidth(range = c(0.1, 1))+
  geom_line(data = region.summary, aes(x = year, y = mean_PPT_all), color = "black", linewidth = 1.1)+
  theme_bw(base_size = 12)+ylab("Mean Precipitation (mm)")+xlab("Year")+
  xlim(1925, 2038)+
  geom_label(data = state.locations.ppt, aes(label = State, 
                                             x = year + 8, 
                                             y = yloc, 
                                             color = State), 
             direction = "y", box.padding = 100, min.segment.length = 2)+
  scale_color_manual(values = state.scales)+
  theme(legend.position = "none")
ggsave(plot = ppt.ts.plt, height = 5, width = 5, units = "in", dpi = 300, "images/all_state_PPT_time_series_1925_2030.png")


# get the N deposition data
boxdir <- "C:/Users/KellyHeilman/Box/01. kelly.heilman Workspace/mortality/Eastern-Mortality/mortality_models/"
N.oxidized <- read_delim(paste0(boxdir,"data/Ndep/Atmospheric_Oxidized.txt"))
N.oxidized$county <- as.numeric(N.oxidized$CountyFIPS)
N.oxidized$state <- as.numeric(N.oxidized$StateFIPS)

N.reduced <- read_delim(paste0(boxdir, "data/Ndep/Atmospheric_Reduced.txt"))
N.reduced$county <- as.numeric(N.reduced$CountyFIPS)
N.reduced$statecd <- as.numeric(N.reduced$StateFIPS)

Ntotal <- N.reduced
Ntotal[,5:92] <- N.reduced[,5:92]+N.oxidized[,5:92] # add together to get the total:
N.total.dep <- Ntotal %>% dplyr::select(statecd, county,`AREA (HA)`, y1930:y2017)
N.total.dep


Ndep.melt <- reshape2::melt(N.total.dep, id.vars = c( "county", "statecd", "AREA (HA)") )
Ndep.melt$year <- substring(Ndep.melt$variable, 2)

Ndep.totals <- Ndep.melt %>% filter(statecd %in% unique(plot.rempers$state)) %>% group_by(statecd, year) %>%
  mutate(areal.value = value/`AREA (HA)`) %>% 
  summarise(Avg.Ndep = mean(value), 
            Total.State.Ndep = sum(value)) %>%
  mutate(year = as.numeric(year))

Ndep.totals$State <- ref_statecd[match(Ndep.totals$statecd, ref_statecd$VALUE),]$MEANING

Ndep.total.remper <- left_join(Ndep.totals, plotcommon.remper) %>%
  mutate(`State Remeasurement` = ifelse(year >=T1 & year <= T2, 0.2, 0.1))


state.locations.ndep <- Ndep.total.remper %>% filter(year == T2) %>% arrange(Avg.Ndep)
state.locations.ndep$Ndep.loc <- c(9, 10, 11, 13, 14, 15, 16, 17, 19, 20)

state.scales.Ndep <- c("Maine"="#543005",
                      "Vermont"="#8c510a",
                      "New Hampshire"="#bf812d",
                      "New York" ="#b2abd2",
                      "Ohio"="#b2abd2",
                      "Maryland"="#8073ac",
                      "West Virginia"="#807dba",
                      "Connecticut" = "#6a51a3",
                      "Pennsylvania"="#542788",
                      "New Jersey"="#2d004b")



ggplot()+
  geom_rect(data = T1.T2periodic, aes(xmin = T1.all, 
                                      xmax = T2.all, 
                                      ymin = -Inf, 
                                      ymax = Inf), alpha = 0.25)+
  geom_line(data = Ndep.total.remper, aes(x = year, y = Avg.Ndep, group = State, color = State), alpha = 0.9)+
  scale_linewidth(range = c(0.1, 1))+
  #geom_line(data = region.summary, aes(x = year, y = mean_PPT_all), color = "black", linewidth = 1.1)+
  theme_bw(base_size = 12)+ylab("Reduced N deposition (kg/ha/year)")+xlab("Year")+
  xlim(1900, 2038)+
  geom_label(data = state.locations.ndep, aes(label = State,
                                             x = 1920,
                                             y = Ndep.loc,
                                             color = State))+
  scale_color_manual(values = state.scales)+
  theme(legend.position = "none")
ggsave(height = 5, width = 8, units = "in", dpi = 300, "images/all_state_Ndep_time_series.png")


ndep.ts.plt <- ggplot()+
  geom_rect(data = T1.T2periodic, aes(xmin = T1.all, 
                                      xmax = T2.all, 
                                      ymin = -Inf, 
                                      ymax = Inf), alpha = 0.25)+
  geom_line(data = Ndep.total.remper, aes(x = year, y = Avg.Ndep, group = State, color = State))+
  scale_linewidth(range = c(0.1, 1))+
  #geom_line(data = region.summary, aes(x = year, y = mean_PPT_all), color = "black", linewidth = 1.1)+
  theme_bw(base_size = 12)+ylab("Reduced N deposition (kg/ha/year)")+xlab("Year")+
  xlim(1935, 2038)+
  geom_label(data = state.locations.ndep, aes(label = State,
                                              x = 2025,
                                              y = Ndep.loc,
                                              color = State))+
  scale_color_manual(values = state.scales)+
  theme(legend.position = "none")
ggsave(plot = ndep.ts.plt, height = 5, width = 5, units = "in", dpi = 300, "images/all_state_Ndep_time_series_1925_2030.png")


#### Plot pests and disturbances over time ----
# read in the pest disturbance datafiles:
# spongy moth --defoliation records by state
# downloaded from here:https://apps.fs.usda.gov/nicportal/lddigest/cfm/dsp/dsplddigesthome.cfm
# spruce budworm data --defoliation recoreds by state
# beech bark disease: make table from: https://borealisdata.ca/dataset.xhtml?persistentId=doi:10.7939/DVN/10835
# additional beech bark information: https://www.nrs.fs.usda.gov/pubs/jrnl/2007/nrs_2007_morin_001.pdf
# Citation: By downloading these data the user agrees to cite the following references in any and all publications in which these data are used:
# Cale, Jonathan A; Morin, Randall S. "County-level detection dates for beech scale in Canada and the United States", http://dx.doi.org/10.7939/DVN/10835 V1 [Version]
# 
# Cale, Jonathan A; Garrison-Johnston; Mariann T; Teale, Stephen A; Castello, John D. Beech bark disease: Over a century of research revisited. 2017. Forest Ecology and Management. DOI: 10.1016/j.foreco.2017.03.031.
# 
# Morin RS, Liebhold AM, Tobin PC, Gottschalk KW, Luzader E. Spread of beech bark disease in the eastern United States and its relationship to regional forest composition. 2007. Canadian Journal of Forest Research. 37: 726-736. DOI: 10.1139/X06-281

# hemlock wooley adelgid infestations:
# hard to track down: using this data here: https://hemlock-woolly-adelgid-national-initiative-gmsts.hub.arcgis.com/pages/distribution-maps
# 
# get county level area information to plot disturbances by area:
# since this is just the time it was first observed, we want county area or something:
# Get county data for all the states
counties <- counties(state = c("DE","RI","ME", "MA","MD", "NH", "NJ", "NY", "OH", "PA", "VT", "WV"), cb = TRUE) %>% 
  data.frame()%>% 
  rename("State" = "STATE_NAME", 
         "CONAME" = "NAME") %>% select(State, STUSPS, CONAME, ALAND, AWATER)
# separate handling for CT, which changed thier county names around 2017
CT.counties <- counties(state = "CT", cb = FALSE, year = 2015) %>% 
  data.frame()%>% 
  rename( "CONAME" = "NAME") %>%
  mutate(State = "Connecticut", 
         STUSPS = "CT") %>% select(State, STUSPS, CONAME, ALAND, AWATER)

counties <- rbind(counties, CT.counties)



# spongy month defoliation records by state:
spongy <- read.csv("data/NE_spongy_moth_outbreaks.csv") %>% 
  rename( "statecd"="STATECD") %>% left_join(., plotcommon.remper) %>%
  group_by(statecd)%>%
  mutate(max.def = max(Acres.Defoliated))

# get total land and convert to hectares
# according to tigris census.gov all units are in square meters
state.area <- counties %>% mutate(ALAND.ha = ALAND/10000) %>%
  group_by(State, STUSPS)%>% summarise(total.land.ha = sum(ALAND.ha))
# join county area up to spongy moth and convert spongy moth defoliation to hectares too
spongy <- spongy %>% left_join(., state.area) %>% 
  mutate(HA.Defoliated = Acres.Defoliated/2.47105381)%>% #conversion
  mutate(fraction.Defoliated = HA.Defoliated/total.land.ha)


state.locations.spongy <- spongy %>% filter(Acres.Defoliated == max.def ) %>% arrange(Acres.Defoliated) %>%
  mutate(State.year = paste(State, year))
# state.locations.spongy$acre.loc <- c(10000,15000, 
#                                      120000, 130000,  
#                                      600000, 650000, 
#                                      700000, 1000000, 1500000, 
#                                      1900000, 2400000, 
#                                      30000000, 44000000)

spongy$State <- factor(spongy$State, levels = rev(c(state.locations.spongy$State)))


state.scales.spongy <- c("#FFAA00", 
                          "#d94801", 
                         "goldenrod", 
                          "#98D851",
                         "darkgreen",
                         "#01665e",
                          "#35978f", 
                          "#1d91c0",
                          "#225ea8" ,
                           "#253494",
                          "#081d58", 
                          
                          "#810f7c", 
                         "#4d004b")
names(state.scales.spongy) <- c(state.locations.spongy$State)


ggplot()+
  geom_rect(data = T1.T2periodic, aes(xmin = T1.all, 
                                      xmax = T2.all, 
                                      ymin = -Inf, 
                                      ymax = Inf), alpha = 0.25)+
  geom_line(data = spongy %>% filter(!is.na(State)), aes(x = year, y = fraction.Defoliated, group = State, color = State), alpha = 0.9, position = "identity", size = 1)+
  #scale_linewidth(range = c(0.1, 1))+
  #geom_line(data = region.summary, aes(x = year, y = mean_PPT_all), color = "black", linewidth = 1.1)+
  theme_bw(base_size = 12)+ylab("Fraction of Land Area Defoliated (Lymantria Dispar)")+xlab("Year")+
  xlim(1900, 2038)+
  geom_label_repel(data = state.locations.spongy, aes(label = State.year,
                                              x = year,
                                              y = fraction.Defoliated,
                                              color = State, 
                                              nudge_x = ifelse(year < 1990, year - 25, year + 15), 
                                              #nudge_y = ifelse(State %in% c("Ohio", "Maine"),1500000, fraction.Defoliated+200000)
                                              ), 
                   box.padding = 3, max.overlaps = Inf, min.segment.length = 0, segment.size = 0.5, 
                  
                 
                  direction = "y",
                  
                   segment.color = "black")+
  scale_fill_manual(values = state.scales.spongy)+
  scale_color_manual(values = state.scales.spongy)+
  theme(legend.position = "none")
ggsave(height = 5, width = 8, units = "in", dpi = 300, "images/all_state_Spongy_fraction_defoliated_time_series.png")

ggplot()+
  geom_rect(data = T1.T2periodic, aes(xmin = T1.all, 
                                      xmax = T2.all, 
                                      ymin = -Inf, 
                                      ymax = Inf), alpha = 0.25)+
  geom_line(data = spongy %>% filter(!is.na(State)), aes(x = year, y = fraction.Defoliated, group = State, color = State), alpha = 0.9, position = "identity", size = 1)+
  #scale_linewidth(range = c(0.1, 1))+
  #geom_line(data = region.summary, aes(x = year, y = mean_PPT_all), color = "black", linewidth = 1.1)+
  theme_bw(base_size = 12)+ylab("Fraction of Land Area Defoliated (Lymantria Dispar)")+xlab("Year")+
  xlim(1925, 2038)+
  geom_label_repel(data = state.locations.spongy, aes(label = State.year,
                                                      x = year,
                                                      y = fraction.Defoliated,
                                                      color = State, 
                                                      nudge_x = ifelse(year < 1990, year - 25, year + 15), 
                                                      #nudge_y = ifelse(State %in% c("Ohio", "Maine"),1500000, fraction.Defoliated+200000)
  ), 
  box.padding = 3, max.overlaps = Inf, min.segment.length = 0, segment.size = 0.5, 
  
  
  direction = "y",
  
  segment.color = "black")+
  scale_fill_manual(values = state.scales.spongy)+
  scale_color_manual(values = state.scales.spongy)+
  theme(legend.position = "none")
ggsave(height = 5, width = 5, units = "in", dpi = 300, 
       "images/all_state_Spongy_fraction_defoliated_time_series_1925_2030.png")



# plot up the total area defoliated:
ggplot()+
  geom_rect(data = T1.T2periodic, aes(xmin = T1.all, 
                                      xmax = T2.all, 
                                      ymin = -Inf, 
                                      ymax = Inf), alpha = 0.25)+
  geom_area(data = spongy %>% filter(!is.na(State)), aes(x = year, y = HA.Defoliated, group = State, fill = State), alpha = 0.9, position = "stack")+
  #scale_linewidth(range = c(0.1, 1))+
  #geom_line(data = region.summary, aes(x = year, y = mean_PPT_all), color = "black", linewidth = 1.1)+
  theme_bw(base_size = 12)+ylab("Hectares Defoliated (Lymantria Dispar)")+xlab("Year")+
  xlim(1900, 2038)+
  geom_label_repel(data = state.locations.spongy, aes(label = State.year,
                                                      x = year,
                                                      y = HA.Defoliated,
                                                      color = State, 
                                                      nudge_x = ifelse(year < 1990, year - 25, year + 15), 
                                                      nudge_y = ifelse(State %in% c("Ohio", "Maine"),1500000, HA.Defoliated+200000)
  ), 
  box.padding = 3, max.overlaps = Inf, min.segment.length = 0, segment.size = 0.5, 
  
  
  direction = "y",
  
  segment.color = "black")+
  scale_fill_manual(values = state.scales.spongy)+
  scale_color_manual(values = state.scales.spongy)+
  theme(legend.position = "none")
ggsave(height = 5, width = 8, units = "in", dpi = 300, "images/all_state_Spongy_time_series.png")
# plot up the total area defoliated:
spongy.ts.plt <- ggplot()+
  geom_rect(data = T1.T2periodic, aes(xmin = T1.all, 
                                      xmax = T2.all, 
                                      ymin = -Inf, 
                                      ymax = Inf), alpha = 0.25)+
  geom_area(data = spongy %>% filter(!is.na(State)), aes(x = year, y = HA.Defoliated, group = State, fill = State), alpha = 0.9, position = "stack")+
  #scale_linewidth(range = c(0.1, 1))+
  #geom_line(data = region.summary, aes(x = year, y = mean_PPT_all), color = "black", linewidth = 1.1)+
  theme_bw(base_size = 12)+ylab("Hectares Defoliated (Lymantria Dispar)")+xlab("Year")+
  xlim(1925, 2038)+
  geom_label_repel(data = state.locations.spongy, aes(label = State.year,
                                                      x = year,
                                                      y = HA.Defoliated,
                                                      color = State, 
                                                      nudge_x = ifelse(year < 1990, year - 25, year + 15), 
                                                      nudge_y = ifelse(State %in% c("Ohio", "Maine"),1500000, HA.Defoliated+200000)
  ), 
  box.padding = 3, max.overlaps = Inf, min.segment.length = 0, segment.size = 0.5, 
  
  
  direction = "y",
  
  segment.color = "black")+
  scale_fill_manual(values = state.scales.spongy)+
  scale_color_manual(values = state.scales.spongy)+
  theme(legend.position = "none")
ggsave(plot = spongy.ts.plt, height = 5, width = 5, units = "in", dpi = 300, 
       "images/all_state_Spongy_time_series_1925_2030.png")





# spruce budworm defoliation records by state:
budworm <- read.csv("data/NE_spruce_budworm_outbreaks.csv") %>% 
  rename( "statecd"="STATECD") %>% left_join(., plotcommon.remper) %>%
  mutate(K.Acres.Defoliated = Spruce.Budworm.Acres.Defoliated/1000)%>%
  group_by(statecd)%>%
  mutate(max.def = max(Spruce.Budworm.Acres.Defoliated))%>% 
  # join up with the state area
  left_join(., state.area) %>% 
  mutate(HA.Defoliated = Spruce.Budworm.Acres.Defoliated/2.47105381)%>% #conversion
  mutate(fraction.Defoliated = Spruce.Budworm.Acres.Defoliated/total.land.ha)



state.locations.budworm <- budworm %>% filter(Spruce.Budworm.Acres.Defoliated == max.def & !State %in% "Total" & !is.na(State)) %>% 
  arrange(Spruce.Budworm.Acres.Defoliated) %>% mutate(earlier.yr = min(Year))%>%
  filter(Year == earlier.yr)
#state.locations.budworm$acre.loc <- c(10000,120000, 130000,  600000, 650000, 700000, 1500000, 1900000, 2400000, 30000000, 44000000)

budworm$State <- factor(budworm$State, levels = rev(c(state.locations.budworm$State)))

# state.scales.budworm <- c( "Ohio" = "#fdae61", 
#                           "West Virginia" = "#d6604d", 
#                           "Maryland" = "#f4a582", 
#                           "Vermont" = "#d6604d", 
#                           "New Hampshire" = "#41b6c4", 
#                           "New York" = "#1d91c0",
#                           "Connecticut" = "#225ea8" ,
#                           "Maine" = "#253494",
#                           "Pennsylvania"="#081d58")


ggplot()+
  geom_rect(data = T1.T2periodic, aes(xmin = T1.all, 
                                      xmax = T2.all, 
                                      ymin = -Inf, 
                                      ymax = Inf), alpha = 0.25)+
  geom_area(data = budworm %>% filter(!is.na(State) & ! State %in% "Total" & !is.na(State)), aes(x = Year, y = fraction.Defoliated, group = State, color = State, fill = State), alpha = 0.5, position = "identity")+
  #scale_linewidth(range = c(0.1, 1))+
  #geom_line(data = region.summary, aes(x = year, y = mean_PPT_all), color = "black", linewidth = 1.1)+
  theme_bw(base_size = 12)+ylab("Spruce Budworm Defoliation (Fraction of State Land Area)")+xlab("Year")+
  xlim(1900, 2038)+
  geom_label_repel(data = state.locations.budworm, aes(label = State,
                                                      x = Year,
                                                      y = fraction.Defoliated,
                                                      color = State, 
                                                      nudge_x = ifelse(State == "Vermont", 1990, Year - 25)#, 
                                                      #nudge_y = Spruce.Budworm.Acres.Defoliated+95
                                                      ),
                   box.padding = 1, max.overlaps = Inf, min.segment.length = 0, segment.size = 0.5, 
                  
                   segment.color = "black")+
  scale_fill_manual(values = state.scales.spongy)+
  scale_color_manual(values = state.scales.spongy)+
  theme(legend.position = "none")
ggsave(height = 5, width = 8, units = "in", dpi = 300, "images/all_state_budworm_fraction_time_series_identity.png")


ggplot()+
  geom_rect(data = T1.T2periodic, aes(xmin = T1.all, 
                                      xmax = T2.all, 
                                      ymin = -Inf, 
                                      ymax = Inf), alpha = 0.25)+
  geom_line(data = budworm %>% filter(!is.na(State) & ! State %in% "Total" & !is.na(State)), aes(x = Year, y = fraction.Defoliated, group = State, color = State, fill = State), size = 1.2, position = "identity")+
  #scale_linewidth(range = c(0.1, 1))+
  #geom_line(data = region.summary, aes(x = year, y = mean_PPT_all), color = "black", linewidth = 1.1)+
  theme_bw(base_size = 12)+ylab("Spruce Budworm Defoliation (Fraction of State Land Area)")+xlab("Year")+
  xlim(1900, 2038)+
  geom_label_repel(data = state.locations.budworm, aes(label = State,
                                                       x = Year,
                                                       y = fraction.Defoliated,
                                                       color = State, 
                                                       nudge_x = ifelse(State == "Vermont", 1990, Year - 25)#, 
                                                       #nudge_y = Spruce.Budworm.Acres.Defoliated+95
  ),
  box.padding = 1, max.overlaps = Inf, min.segment.length = 0, segment.size = 0.5, 
  
  segment.color = "black")+
  scale_fill_manual(values = state.scales.spongy)+
  scale_color_manual(values = state.scales.spongy)+
  theme(legend.position = "none")
ggsave(height = 5, width = 8, units = "in", dpi = 300, "images/all_state_budworm_fraction_time_series_line.png")

# plot relativized acreage defoliated:
ggplot()+
  geom_rect(data = T1.T2periodic, aes(xmin = T1.all, 
                                      xmax = T2.all, 
                                      ymin = -Inf, 
                                      ymax = Inf), alpha = 0.25)+
  geom_area(data = budworm %>% filter(!is.na(State) & ! State %in% "Total" & !is.na(State)), aes(x = Year, y = HA.Defoliated, group = State, fill = State), alpha = 0.9, position = "stack")+
  #scale_linewidth(range = c(0.1, 1))+
  #geom_line(data = region.summary, aes(x = year, y = mean_PPT_all), color = "black", linewidth = 1.1)+
  theme_bw(base_size = 12)+ylab("Spruce Budworm Defoliation (Hectares)")+xlab("Year")+
  xlim(1900, 2038)+
  geom_label_repel(data = state.locations.budworm, aes(label = State,
                                                       x = Year,
                                                       y = HA.Defoliated,
                                                       color = State, 
                                                       nudge_x = ifelse(State == "Vermont", 1990, Year - 25)#, 
                                                       #nudge_y = Spruce.Budworm.Acres.Defoliated+95
  ),
  box.padding = 1, max.overlaps = Inf, min.segment.length = 0, segment.size = 0.5, 
  
  segment.color = "black")+
  scale_fill_manual(values = state.scales.spongy)+
  scale_color_manual(values = state.scales.spongy)+
  theme(legend.position = "none")
ggsave(height = 5, width = 8, units = "in", dpi = 300, "images/all_state_budworm_time_series.png")

# plot relativized acreage defoliated:
budworm.ts.plt <- ggplot()+
  geom_rect(data = T1.T2periodic, aes(xmin = T1.all, 
                                      xmax = T2.all, 
                                      ymin = -Inf, 
                                      ymax = Inf), alpha = 0.25)+
  geom_area(data = budworm %>% filter(!is.na(State) & ! State %in% "Total" & !is.na(State)), aes(x = Year, y = HA.Defoliated, group = State, fill = State), alpha = 0.9, position = "stack")+
  #scale_linewidth(range = c(0.1, 1))+
  #geom_line(data = region.summary, aes(x = year, y = mean_PPT_all), color = "black", linewidth = 1.1)+
  theme_bw(base_size = 12)+ylab("Spruce Budworm Defoliation (Hectares)")+xlab("Year")+
  xlim(1925, 2038)+
  geom_label_repel(data = state.locations.budworm, aes(label = State,
                                                       x = Year,
                                                       y = HA.Defoliated,
                                                       color = State, 
                                                       nudge_x = ifelse(State == "Vermont", 1990, Year - 25)#, 
                                                       #nudge_y = Spruce.Budworm.Acres.Defoliated+95
  ),
  box.padding = 1, max.overlaps = Inf, min.segment.length = 0, segment.size = 0.5, 
  
  segment.color = "black")+
  scale_fill_manual(values = state.scales.spongy)+
  scale_color_manual(values = state.scales.spongy)+
  theme(legend.position = "none")
ggsave(plot = budworm.ts.plt, height = 5, width = 5, units = "in", dpi = 300, 
       "images/all_state_budworm_time_series_1925_2030.png")



### Beech Bark disease spread ---
beech.bark <- read.csv("data/beech_bark_spread/Cale_Morin-BeechScaleDatesCanadaUS.csv") %>%
  filter(COUNTRY %in% "USA") %>%
  rename("State" = "PRVSTTNAME") %>% 
  filter(State %in% unique(state.summary.remper$State)) # just gets the remper states


# join to beech bark data
beech.bark.cos <- counties %>% left_join(.,beech.bark %>% select(-REFERENCE)) %>% select(-STUSPS) 
beech.bark.cos %>% filter(is.na(ALAND)) # make sure all have matches

# now need to fill out timeseries:
beech.scale <- beech.bark.cos  %>% distinct() %>%
  expand_grid(year = 1900:2020) %>% # expand the DF to include years 
  mutate(BB.county = ifelse(SCALEYR <= year, "Beech Bark Present", "No observed Beech Bark"))%>% # categorical variable that says if there is BB present in this year
  # calculate rolling total of land where BB is present in each year
  group_by(State, year) %>%
           group_by(State, year, BB.county) %>%
           summarise(Total.Land.ha = sum(ALAND/10000)) %>% spread(BB.county, Total.Land.ha, fill = 0) %>%
  mutate(Total.State.ha = `Beech Bark Present` + `No observed Beech Bark`)%>%
  mutate(Percent.infested = (`Beech Bark Present`/Total.State.ha)*100)

# get the labels of the first year:

state.locations.beech.scale <- beech.scale %>% left_join(., beech.bark.cos %>% group_by(State) %>% summarise(first.yr = min(SCALEYR, na.rm = TRUE))) %>% 
  filter(year == first.yr) %>% arrange(first.yr)%>%
  mutate(State.year = paste(State, year))



ggplot()+
  geom_rect(data = T1.T2periodic, aes(xmin = T1.all, 
                                      xmax = T2.all, 
                                      ymin = -Inf, 
                                      ymax = Inf), alpha = 0.25)+
  geom_line(data = beech.scale, aes(x = year, y = Percent.infested, group = State, color = State), linewidth = 2)+
  #scale_linewidth(range = c(0.1, 1))+
  #geom_line(data = region.summary, aes(x = year, y = mean_PPT_all), color = "black", linewidth = 1.1)+
  theme_bw(base_size = 12)+ylab("Beech Scale Presence (% land area)")+xlab("Year")+
  xlim(1900, 2020)+
  geom_label_repel(data = state.locations.beech.scale, aes(label = State.year,
                                                       x = year,
                                                       y = Percent.infested,
                                                       color = State,
                                                      nudge_x = ifelse(year < 1975, year -25, year + 25),
                                                      nudge_y = ifelse(State %in% c("New York", "Maine"), Percent.infested-10, 
                                                                       ifelse(State %in% c("Pennsylvania", "New Hampshire"), Percent.infested+7, Percent.infested))
                                                       ),
                   box.padding = 0.25, max.overlaps = Inf, min.segment.length = 0, segment.size = 0.5,
                   direction = "x", segment.color = "black")+
  #scale_fill_manual(values = state.scales.budworm)+
  scale_color_manual(values = state.scales.spongy)+
  theme(legend.position = "none")
ggsave(height = 5, width = 8, units = "in", dpi = 300, "images/all_state_Beech_scale_time_series.png")

ggplot()+
  geom_rect(data = T1.T2periodic, aes(xmin = T1.all, 
                                      xmax = T2.all, 
                                      ymin = -Inf, 
                                      ymax = Inf), alpha = 0.25)+
  geom_area(data = beech.scale, aes(x = year, y = Percent.infested, group = State, color = State, fill = State), alpha = 0.5, linewidth = 2, position = "identity")+
  #scale_linewidth(range = c(0.1, 1))+
  #geom_line(data = region.summary, aes(x = year, y = mean_PPT_all), color = "black", linewidth = 1.1)+
  theme_bw(base_size = 12)+ylab("Beech Scale Presence (% land area)")+xlab("Year")+
  xlim(1900, 2020)+
  geom_label_repel(data = state.locations.beech.scale, aes(label = State.year,
                                                           x = year,
                                                           y = Percent.infested,
                                                           color = State,
                                                           nudge_x = ifelse(year < 1975, year -25, year + 25),
                                                           nudge_y = ifelse(State %in% c("New York", "Maine"), Percent.infested-10, 
                                                                            ifelse(State %in% c("Pennsylvania", "New Hampshire"), Percent.infested+7, Percent.infested))
  ),
  box.padding = 0.25, max.overlaps = Inf, min.segment.length = 0, segment.size = 0.5,
  direction = "x", segment.color = "black")+
  scale_fill_manual(values = state.scales.spongy)+
  scale_color_manual(values = state.scales.spongy)+
  theme(legend.position = "none")
#ggsave(height = 5, width = 8, units = "in", dpi = 300, "images/all_state_Beech_scale_time_series.png")

# plot the total area with beech scale presence:
ggplot()+
  geom_rect(data = T1.T2periodic, aes(xmin = T1.all, 
                                      xmax = T2.all, 
                                      ymin = -Inf, 
                                      ymax = Inf), alpha = 0.25)+
  geom_area(data = beech.scale, aes(x = year, y = `Beech Bark Present`, group = State, color = State, fill = State), alpha = 0.9, linewidth = 1, position = "stack")+
  #scale_linewidth(range = c(0.1, 1))+
  #geom_line(data = region.summary, aes(x = year, y = mean_PPT_all), color = "black", linewidth = 1.1)+
  theme_bw(base_size = 12)+ylab("Total area with Beech Scale (hectares)")+xlab("Year")+
  xlim(1900, 2020)+
  geom_label_repel(data = state.locations.beech.scale, aes(label = State.year,
                                                           x = year,
                                                           y = `Beech Bark Present`,
                                                           color = State,
                                                           nudge_x = ifelse(year < 1975, year -25, year + 25) ,
                                                          nudge_y = ifelse(State %in% c("New York", "Maine"), `Beech Bark Present`-1, 
                                                                            ifelse(State %in% c("Pennsylvania", "New Hampshire"), `Beech Bark Present`+7, `Beech Bark Present`))
  ),
  box.padding = 1, max.overlaps = Inf, min.segment.length = 0, segment.size = 0.5,
  direction = "both", segment.color = "black")+
  scale_fill_manual(values = state.scales.spongy)+
  scale_color_manual(values = state.scales.spongy)+
  theme(legend.position = "none")
ggsave(height = 5, width = 8, units = "in", dpi = 300, "images/all_state_Beech_scale_time_series_Hectares_by_state.png")

beech.scale.ts.plt <- ggplot()+
  geom_rect(data = T1.T2periodic, aes(xmin = T1.all, 
                                      xmax = T2.all, 
                                      ymin = -Inf, 
                                      ymax = Inf), alpha = 0.25)+
  geom_area(data = beech.scale, aes(x = year, y = `Beech Bark Present`, group = State, color = State, fill = State), alpha = 0.9, linewidth = 1, position = "stack")+
  #scale_linewidth(range = c(0.1, 1))+
  #geom_line(data = region.summary, aes(x = year, y = mean_PPT_all), color = "black", linewidth = 1.1)+
  theme_bw(base_size = 12)+ylab("Total area with Beech Scale (hectares)")+xlab("Year")+
  xlim(1925, 2020)+
  geom_label_repel(data = state.locations.beech.scale, aes(label = State.year,
                                                           x = year,
                                                           y = `Beech Bark Present`,
                                                           color = State#,
                                                           #nudge_x = ifelse(year < 1975, year -25, year + 25) ,
                                                           #nudge_y = ifelse(State %in% c("New York", "Maine"), `Beech Bark Present`-1, 
                                                            #                ifelse(State %in% c("Pennsylvania", "New Hampshire"), `Beech Bark Present`+7, `Beech Bark Present`))
  ),
  box.padding = 1, max.overlaps = Inf, min.segment.length = 0, segment.size = 0.5,
  direction = "y", segment.color = "black")+
  scale_fill_manual(values = state.scales.spongy)+
  scale_color_manual(values = state.scales.spongy)+
  theme(legend.position = "none")
ggsave(plot = beech.scale.ts.plt, height = 5, width = 5, units = "in", dpi = 300, "images/all_state_Beech_scale_time_series_Hectares_by_state_1925_2030.png")


### Hemlock wooley adelgid disease spread ---
HWA <- read.csv("data/HWA_invested_counties_2023_3_6_24.csv")  # just gets the remper states

# since this is just the time it was first observed, we want county area or something:
# Get county data for all the states
HWA.cos <- counties %>% mutate(CNTYNAME = toupper(CONAME)) %>% 
  rename("ST" = "STUSPS")%>% left_join(.,HWA %>% select(-X) %>% distinct() ) # some duplicates in the data

HWA.cos

HWA.infest <- HWA.cos  %>% distinct() %>%
  expand_grid(year = 1900:2020) %>% # expand the DF to include years 
  mutate(HWA.county = ifelse(year_txt <= year, "HWA Present", "No observed HWA"))%>% # categorical variable that says if there is BB present in this year
  # calculate rolling total of land where BB is present in each year
  group_by(State, year) %>%
  group_by(State, year, HWA.county) %>%
  summarise(Total.Land.ha = sum(ALAND/10000)) %>% 
  spread(HWA.county, Total.Land.ha, fill = 0) %>%
  mutate(Total.State.ha = `HWA Present` + `No observed HWA`)%>%
  mutate(Percent.infested = (`HWA Present`/Total.State.ha)*100)

# get the labels of the first year:

state.locations.HWA.infest <- HWA.infest %>% 
  left_join(., HWA.cos %>% group_by(State)%>% summarise(first.yr = min(year_txt, na.rm =TRUE))) %>% 
  filter(year == first.yr) %>% arrange(first.yr)

# set color scale:
state.scales.budworm <- c("#b11e02",
                          "#ff355e", 
                          "#fd5b78",
                          "#ff6037", 
                          "#ffcc33",
                          "goldenrod",
                          "#4aa682", 
                          "darkgreen", 
                          "#16d0cb", 
                          "#50bfe6",
                          "#4067a8",
                          "#9c3ed6", 
                          "#ee34d2") 
names(state.scales.budworm) <- state.locations.HWA.infest$State


# alternate labelling:
state.locations.HWA.infest <- state.locations.HWA.infest %>%
  mutate(State.year = paste(State, year))

ggplot()+
  geom_rect(data = T1.T2periodic, aes(xmin = T1.all, 
                                      xmax = T2.all, 
                                      ymin = -Inf, 
                                      ymax = Inf), alpha = 0.25)+
  geom_line(data = HWA.infest, aes(x = year, y = Percent.infested, group = State, color = State), linewidth = 2)+
  #scale_linewidth(range = c(0.1, 1))+
  #geom_line(data = region.summary, aes(x = year, y = mean_PPT_all), color = "black", linewidth = 1.1)+
  theme_bw(base_size = 12)+ylab("Hemlock Wooley Adelgid Presence (% land area)")+xlab("Year")+
  xlim(1900, 2020)+
  geom_label_repel(data = state.locations.HWA.infest, aes(label = State.year,
                                                          x = year,
                                                          y = Percent.infested,
                                                          color = State,
                                                          nudge_x = ifelse(State %in% c("Vermont", "Ohio"),2035, 
                                                                           ifelse(State %in% c("Maine", "New Hampshire"), year+10, year - 35)),
                                                          nudge_y = ifelse(State == "Ohio",12, Percent.infested)
  ),
  box.padding =2, max.overlaps = Inf, min.segment.length = 0, segment.size = 0.5,
  direction = "y", 
  segment.color = "black")+
  #scale_fill_manual(values = state.scales.budworm)+
  scale_color_manual(values = state.scales.spongy)+
  theme(legend.position = "none")+xlim(1900, 2030)
ggsave(height = 5, width = 8, units = "in", dpi = 300, "images/all_state_HWA_time_series_alternate.png")

hwa.ts.plt <- ggplot()+
  geom_rect(data = T1.T2periodic, aes(xmin = T1.all, 
                                      xmax = T2.all, 
                                      ymin = -Inf, 
                                      ymax = Inf), alpha = 0.25)+
  geom_line(data = HWA.infest, aes(x = year, y = Percent.infested, group = State, color = State), linewidth = 2)+
  #scale_linewidth(range = c(0.1, 1))+
  #geom_line(data = region.summary, aes(x = year, y = mean_PPT_all), color = "black", linewidth = 1.1)+
  theme_bw(base_size = 12)+ylab("Hemlock Wooley Adelgid Presence (% land area)")+xlab("Year")+
  geom_label_repel(data = state.locations.HWA.infest, aes(label = State.year,
                                                          x = year,
                                                          y = Percent.infested,
                                                          color = State,
                                                          nudge_x = ifelse(State %in% c("Vermont", "Ohio"),2035, 
                                                                           ifelse(State %in% c("Maine", "New Hampshire"), year+10, year - 35)),
                                                          nudge_y = ifelse(State == "Ohio",12, Percent.infested)
  ),
  box.padding =2, max.overlaps = Inf, min.segment.length = 0, segment.size = 0.5,
  direction = "y", 
  segment.color = "black")+
  #scale_fill_manual(values = state.scales.budworm)+
  scale_color_manual(values = state.scales.spongy)+
  theme(legend.position = "none")+xlim(1925, 2030)
ggsave(plot = hwa.ts.plt, height = 5, width = 5, units = "in", dpi = 300, "images/all_state_HWA_time_series_alternate_1925_2030.png")

ggplot()+
  geom_rect(data = T1.T2periodic, aes(xmin = T1.all, 
                                      xmax = T2.all, 
                                      ymin = -Inf, 
                                      ymax = Inf), alpha = 0.25)+
  geom_area(data = HWA.infest, aes(x = year, y = `HWA Present`,  color = State, fill = State), linewidth= 0.5, position = "stack")+
  #scale_linewidth(range = c(0.1, 1))+
  #geom_line(data = region.summary, aes(x = year, y = mean_PPT_all), color = "black", linewidth = 1.1)+
  theme_bw(base_size = 12)+ylab("Hemlock Wooley Adelgid Presence \n (Hectares)")+xlab("Year")+
  geom_label_repel(data = state.locations.HWA.infest, aes(label = State.year,
                                                          x = year,
                                                          y = Percent.infested,
                                                          color = State,
                                                          # nudge_x = ifelse(State %in% c("Vermont", "Ohio"),2035, 
                                                          #                  ifelse(State %in% c("Maine", "New Hampshire"), year+10, year - 35)),
                                                          # nudge_y = ifelse(State == "Ohio",12, Percent.infested)
  ),
  box.padding =1, max.overlaps = Inf, min.segment.length = 0, segment.size = 0.5,
  direction = "both", 
  segment.color = "black")+
  scale_fill_manual(values = state.scales.spongy)+
  scale_color_manual(values = state.scales.spongy)+
  theme(legend.position = "none")+xlim(1925, 2030)
ggsave(height = 5, width = 5, units = "in", dpi = 300, "images/all_state_HWA_time_series_land_area_1925_2030.png")



# plot up the other drivers of veg change
# make up state-level summarise of variables, using all trees of focal species, including cut trees:
dmg.causes <- data.frame(Tree.status = rep(c("cut", "dead", "live"),each = 10),
                         damage = rep(c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90),3), 
                         agent = rep(c("No damage/unknown cause", 
                                   "Insect", 
                                   "Disease", 
                                   "Fire", 
                                   "Animal", 
                                   "Weather", 
                                   "Suppression", 
                                   "Miscellaneous", 
                                   "Logging", 
                                   "Form"), 3)) %>%
  mutate(`Damage or Mortality` = ifelse(Tree.status %in% "cut", "Logging", 
                                        ifelse(Tree.status %in% "dead" & agent %in% "No damage/unknown cause", "Unknown cause", 
                                               ifelse(Tree.status %in% "live" & agent %in% "No damage/unknown cause", "No Damage", 
                                                      ifelse(agent %in% "Miscellaneous", "Unknown cause", agent)))))


TREE.dmg.all <- TREE.data.all %>%
  mutate(Top.Species = ifelse(SPCD %in% c(nspp[1:17,]$SPCD, 372,132), Species, 
                              ifelse(SW_HW %in% "S", "other softwood","other hardwood"))) %>% # label for non-top species
  left_join(., dmg.causes)
TREE.dmg.all%>%
  group_by(Top.Species, Tree.status, `Damage or Mortality`) %>% summarise(n())

damage.fill <- scale_fill_manual(values = c(
   "Fire" = "#e41a1c",
 "Suppression" =  "#377eb8",
   "Insect" = "#4daf4a", 
 "Disease" =  "#984ea3" ,
  "Unknown cause" =  "#ff7f00" ,
   "No Damage" = "#ffff33" ,
  "Animal" =  "#a65628" ,
 "Weather" =  "#f781bf" ,
  "Logging" = "#999999" 
))

ggplot(TREE.dmg.all, aes(x = Tree.status, fill = `Damage or Mortality`))+
  geom_bar()+theme_bw(base_size = 14)+
  xlab("Tree Status")+ylab("# of trees") + damage.fill+
  theme(panel.grid = element_blank())
ggsave(height = 4, width = 4, dpi = 300, "images/all_live_dead_by_damage_cause.png")

dmg.cause.plt <- ggplot(TREE.dmg.all %>% filter(!`Damage or Mortality` %in% "No Damage"), aes(x = Tree.status, fill = `Damage or Mortality`))+
  geom_bar()+theme_bw(base_size = 14)+
  xlab("Tree Status")+ylab("# of trees") + damage.fill+
  theme(panel.grid = element_blank())
ggsave(plot = dmg.cause.plt, height = 4, width = 4, dpi = 300, "images/all_damaged_dead_by_cause.png")


TREE.dmg.all$STNAME <- factor(TREE.dmg.all$STNAME, levels = c("Maine", "New York", "Pennsylvania", "West Virginia", 
                                                              "Ohio", "New Hampshire","Vermont", "Maryland", "Connecticut", "New Jersey" ))

dmg.cause.plt.state <- ggplot(TREE.dmg.all %>% filter(!`Damage or Mortality` %in% "No Damage") %>% mutate(L.D.status = ifelse(Tree.status %in% "cut", "dead", 
                                                                                                       Tree.status)), 
       aes(x = STNAME, fill = `Damage or Mortality`))+
  geom_bar()+theme_bw(base_size = 14)+
  xlab("State")+ylab("# of trees") +  damage.fill+
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1))+facet_wrap(~L.D.status, scales = "free_y")
ggsave(plot = dmg.cause.plt.state, height = 4, width = 8, dpi = 300, paste0(output.dir, "images/all_damaged_dead_by_cause_state.png"))



ggplot(data = TREE.dmg.all %>% filter(!Tree.status %in% "cut" & dbhold > 5))+
  geom_density(aes(x = dbhold, fill = Tree.status), alpha = 0.5)

ggplot(data = TREE.dmg.all %>% filter(!Tree.status %in% "cut" & dbhold > 5))+
  geom_density(aes(x = elev, fill = Tree.status), alpha = 0.5)

ggplot(data = TREE.dmg.all %>% filter(!Tree.status %in% "cut" & dbhold > 5))+
  geom_density(aes(x = DIA_DIFF, fill = Tree.status), alpha = 0.5)

##########################################################################
# Combine map, composition, and disturbance records in Figure 1-------
##########################################################################
library(cowplot)

png(height = 12, 
    width = 16, 
    units = "in",
    res = 350, 
    paste0(output.dir, "images/Figure_1_option_1.png"))
# arrage in a grid
plot_grid(
  plot_grid( ggdraw(mortality_percent_map) + 
               draw_plot(legend.large, x = .75, y = .25, width = .25, height = .25, scale = 0.15),  labels = c("a)") , 
             plot_grid( st.comp.ba, dmg.cause.plt, ncol = 2, labels = c("b)", "c)")),
             ncol = 2,  rel_widths = c(0.9, 1), rel_heights = c(1, 0.8, 0.8)),# row 2
  
  
  plot_grid(tave.ts.plt, ppt.ts.plt, ndep.ts.plt, ncol = 3, align = "hv", labels = c("d)", "e)", "f)")),# row 1





plot_grid(spongy.ts.plt, hwa.ts.plt, budworm.ts.plt , ncol = 3, align = "hv", labels =c("h)", "i)", "j)")),# row 3
ncol = 1, rel_heights = c(1, 0.6, 0.6))
dev.off()


  
#####################################################################################################state######################################################################################################################
#Map out species-level model predicted probability of mortality for Figure 5 ----
#by species, by plot----
######################################################################################################################
  SPCD.df <- data.frame(SPCD = nspp[1:17, ]$SPCD, 
                        Species = nspp[1:17, ]$COMMON)


# loop to read in the predictions for pmort
  for(i in 1:17){# run for each of the 17 species
    common.name <- nspp[1:17, ] %>% filter(SPCD %in% SPCD.df[i,]$SPCD) %>% dplyr::select(COMMON)
     
    SPCD.id <- SPCD.df[i,]$SPCD
      model.number <- 6
      
      remper.correction <- 0.5
      model.name <- paste0("mort_model_",model.number,"_SPCD_", SPCD.id, "_remper_correction_", remper.correction)
      remp.cor <- 0.5
      
      SPCD.id <-  SPCD.df[i,]$SPCD
        
      cat(paste("mapping out put from stan mortality model ",model.number, " for SPCD", SPCD.df[i,]$SPCD, common.name$COMMON, " remper correction", remper.correction))
    
      fit.1 <- readRDS( paste0(output.dir, "SPCD_stanoutput_full_standardized_v3/samples/model_",model.number,"_SPCD_",SPCD.id, "_remper_correction_", remper.correction, ".RDS"))
      
      # tload the data used to fit the model 
      load(paste0("SPCD_standata_general_full_standardized_v3/SPCD_",SPCD.id, "remper_correction_", remper.correction,"model_",model.number, ".Rdata")) # load the species code data
      
      
      species.table <- unique(train.data[,c("SPCD","SPP")])
      species.table$COMMON <- FIESTA::ref_species[match(species.table$SPCD, FIESTA::ref_species$SPCD),]$COMMON_NAME
      species.table$SPP <- as.character(species.table$SPP)
      
      
      names(fit.1) <- c("alpha_SPP", colnames(mod.data$xM),
                        ## in sample predicted status
                        paste0("yrep[",1:mod.data$N, "]"),
                        paste0("psurv[",1:mod.data$N, "]"),
                        ## out of  sample predicted status
                        paste0("yhat[",1:mod.data$Nrep, "]"),
                        ## out of sample predicted prob mor
                        paste0("psurv.hat[",1:mod.data$Nrep, "]"),
                        ## in sample predicted status
                        #paste0("log_lik[",1:mod.data$N, "]"),
                        "lp__")
      par.names = c("alpha_SPP", colnames(mod.data$xM)) #,
      
      nvariables <- length(par.names)
      
      species.table <- unique(train.data[,c("SPCD","SPP")])
      species.table$COMMON <- FIESTA::ref_species[match(species.table$SPCD, FIESTA::ref_species$SPCD),]$COMMON_NAME
      species.table$SPP <- as.character(species.table$SPP)
      
      fit_ssm_df <- as.data.frame(fit.1) 
     
     # get the probability of survival estimates for in-sample data  
      psurv.estimates <- fit_ssm_df %>% dplyr::select( paste0("psurv[",1:mod.data$N, "]")) 
      psurv.m <- reshape2::melt(psurv.estimates)
      
      psurv.quant <- psurv.m %>% group_by(variable) %>% summarise(median = quantile(value, 0.5, na.rm =TRUE),
                                                                  ci.lo = quantile(value, 0.005, na.rm =TRUE),
                                                                  ci.hi = quantile(value, 0.975, na.rm =TRUE))
      
      
      
      # add data from the training dataset
      psurv.quant$Mobs <- as.character(mod.data$y)
      psurv.quant$COMMON <- unique(species.table$COMMON)
      
      train.data.attributes <- train.data %>% select(state, county, pltnum, point, tree, PLOT.ID, spp, dbhcur, dbhold, Species, LAT_FIADB, LONG_FIADB)
      psurvival.train <- cbind(psurv.quant, train.data.attributes)
      
    
      # get the probability of survival estimates for out of sample data 
      psurv.hat.estimates <- fit_ssm_df %>% dplyr::select( paste0("psurv.hat[",1:mod.data$Nrep, "]")) 
      psurv.hat.m <- reshape2::melt(psurv.hat.estimates)
      
      psurv.hat.quant <- psurv.hat.m %>% group_by(variable) %>% summarise(median = quantile(value, 0.5, na.rm =TRUE),
                                                                  ci.lo = quantile(value, 0.005, na.rm =TRUE),
                                                                  ci.hi = quantile(value, 0.975, na.rm =TRUE))
      
      
      
      # add data from the training dataset
      psurv.hat.quant$Mobs <- as.character(mod.data$ytest)
      psurv.hat.quant$COMMON <- unique(species.table$COMMON)
      
      test.data.attributes <- test.data %>% select(state, county, pltnum, point, tree, PLOT.ID, spp, dbhcur, dbhold, Species, LAT_FIADB, LONG_FIADB)
      psurvival.test <- cbind(psurv.hat.quant, test.data.attributes)
      psurvival.test$dataset <- "held-out"
      psurvival.train$dataset <- "in-sample"
      psurvival.all <-  rbind(psurvival.train, psurvival.test)
      
      
      
      # calculate mortality
     pmort.summaries <-  psurvival.all %>% mutate(pmort.med = 1-median) %>% group_by(PLOT.ID,LAT_FIADB, LONG_FIADB, spp) %>% 
        summarise(mean.pmort = median(pmort.med), 
                  pct.mortality = median(pmort.med)*100,
                  sum.pmort = sum(pmort.med))
     # generate bins
     pmort.summaries <- pmort.summaries %>% mutate(Mort.quantiles = cut(pct.mortality, 
                                                                            breaks = c(0, 1, 2.5, 5, 10, 20, 30, 40, 50, 100), 
                                                                            include.lowest=TRUE))
    
     large.size.font <- ggplot() +
       #geom_sf(alpha = 0.75, aes(fill = as.character(Distribution)))+
       #scale_fill_manual(values = c("Species Distribution" = "forestgreen", "Outside Distribution" = "white"))+
       geom_jitter(data = pmort.summaries, aes(x = LONG_FIADB, y = LAT_FIADB, color = Mort.quantiles), size = 5)+theme_bw()+
       coord_sf(xlim = c(-85, -68), ylim = c(35, 48))+theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'lightblue'), 
                                                            #legend.position = c(0.87, 0.25),
                                                            legend.title = element_blank(),
                                                            legend.background = element_rect(fill = "white", color = "black"))+
       scale_color_manual(values = c(
         "[0,1]" ="lightgrey",
         "(1,2.5]" = "#fcc5c0", 
         "(2.5,5]" = "#fa9fb5", 
         "(5,10]" = "#feb24c", 
         "(10,20]" = "#fd8d3c", 
         "(20,30]" = "#fc4e2a", 
         "(30,40]" = "#e31a1c", 
         "(40,50]" =  "#bd0026", 
         "(50,100]" = "#49006a" 
         
         
       ))
     
     
     legend.large <- cowplot::get_legend(large.size.font)
  
     mortality_percent_map_spp <- ggplot() +
       geom_polygon(data = canada, 
                    aes(x=long, y=lat, group = group), 
                    color = "black", fill = "white") +
       geom_polygon(data = state_sub, 
                    aes(x=long, y=lat, group = group), 
                    color = "black", fill = "white") +
       #geom_sf(alpha = 0.75, aes(fill = as.character(Distribution)))+
       #scale_fill_manual(values = c("Species Distribution" = "forestgreen", "Outside Distribution" = "white"))+
       geom_jitter(data = pmort.summaries, aes(x = LONG_FIADB, y = LAT_FIADB, color = Mort.quantiles), size = 0.5)+theme_bw()+
       coord_sf(xlim = c(-84, -65), ylim = c(37, 47.5))+theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'lightblue'), 
                                                              #legend.position = c(0.87, 0.25),
                                                              legend.title = element_blank(),
                                                              legend.background = element_rect(fill = "white", color = "black"), 
                                                              legend.position = "none")+
       scale_color_manual(values = c(
         "[0,1]" ="lightgrey",
         "(1,2.5]" = "#fcc5c0", 
         "(2.5,5]" = "#fa9fb5", 
         "(5,10]" = "#feb24c", 
         "(10,20]" = "#fd8d3c", 
         "(20,30]" = "#fc4e2a", 
         "(30,40]" = "#e31a1c", 
         "(40,50]" =  "#bd0026", 
         "(50,100]" = "#49006a" 
         
         
       ))+ggtitle(paste0("Plot averaged predicted mortality probability for ", SPCD.df[i,]$Species))
     
    # png(height = 5, width = 8.5, res = 300, units = "in", here(paste0("images/predicted_mortality/average_predicted_pmort_SPCD_",SPCD.df[i,]$SPCD,"_map.png") ))
    cowplot::plot_grid(legend.large, mortality_percent_map_spp,  align = "v", rel_widths = c(0.15, 0.89))+
       theme(panel.background = element_rect(fill = "white", colour = NA),
             plot.background = element_rect(fill = "white", colour = NA))
     ggsave(height = 5, width = 8.5, units = "in", paste0(output.dir, "SPCD_stanoutput_full_standardized_v3/images/predicted_mortality/average_predicted_pmort_SPCD_",SPCD.df[i,]$SPCD,"_map.png") )
     #dev.off()
     
     
    saveRDS(psurvival.all, paste0(output.dir, "SPCD_stanoutput_full_standardized_v3/pmort_LL_model_",model.number,"_SPCD_",SPCD.id, "_remper_correction_", remper.correction, ".RDS"))
     
  }
  
### now read in all the species-level data frames and make one big plot level summary map:
  
pmort.files <- paste0(output.dir, "SPCD_stanoutput_full_standardized_v3/", list.files(output.dir, "SPCD_stanoutput_full_standardized_v3/", "pmort_LL"))
pmort.list <- lapply(pmort.files, readRDS) 
pmort.df <- do.call(rbind, pmort.list)


# calculate mortality
pmort.summaries <- pmort.df %>% mutate(pmort.med = 1-median) %>% group_by(PLOT.ID, date, LAT_FIADB, LONG_FIADB, spp) %>% 
  summarise(mean.pmort = median(pmort.med), 
            pct.mortality = median(pmort.med)*100,
            sum.pmort = sum(pmort.med))
# generate bins
pmort.summaries <- pmort.summaries %>% mutate(Mort.quantiles = cut(pct.mortality, 
                                                                   breaks = c(0, 1, 2.5, 5, 10, 20, 30, 40, 50, 100), 
                                                                   include.lowest=TRUE))

large.size.font <- ggplot() +
  #geom_sf(alpha = 0.75, aes(fill = as.character(Distribution)))+
  #scale_fill_manual(values = c("Species Distribution" = "forestgreen", "Outside Distribution" = "white"))+
  geom_jitter(data = pmort.summaries, aes(x = LONG_FIADB, y = LAT_FIADB, color = Mort.quantiles), size = 5)+theme_bw()+
  coord_sf(xlim = c(-85, -68), ylim = c(35, 48))+theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'lightblue'), 
                                                       #legend.position = c(0.87, 0.25),
                                                       
                                                       legend.background = element_rect(fill = "white", color = "black"))+
  scale_color_manual(values = c(
    "[0,1]" ="lightgrey",
    "(1,2.5]" = "#fcc5c0", 
    "(2.5,5]" = "#fa9fb5", 
    "(5,10]" = "#feb24c", 
    "(10,20]" = "#fd8d3c", 
    "(20,30]" = "#fc4e2a", 
    "(30,40]" = "#e31a1c", 
    "(40,50]" =  "#bd0026", 
    "(50,100]" = "#49006a" 
    
    
  ))+ guides(color=guide_legend(title="Predicted\n mortality\n probability"))


legend.large <- cowplot::get_legend(large.size.font)

mortality_percent_map_all <- ggplot() +
  geom_polygon(data = canada, 
               aes(x=long, y=lat, group = group), 
               color = "black", fill = "white") +
  geom_polygon(data = state_sub, 
               aes(x=long, y=lat, group = group), 
               color = "black", fill = "white") +
  #geom_sf(alpha = 0.75, aes(fill = as.character(Distribution)))+
  #scale_fill_manual(values = c("Species Distribution" = "forestgreen", "Outside Distribution" = "white"))+
  geom_jitter(data = pmort.summaries, aes(x = LONG_FIADB, y = LAT_FIADB, color = Mort.quantiles), size = 0.5)+theme_bw()+
  coord_sf(xlim = c(-84, -65), ylim = c(37, 47.5))+theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'lightblue'), 
                                                         #legend.position = c(0.87, 0.25),
                                                         legend.title = element_blank(),
                                                         legend.background = element_rect(fill = "white", color = "black"), 
                                                         legend.position = "none")+
  scale_color_manual(values = c(
    "[0,1]" ="lightgrey",
    "(1,2.5]" = "#fcc5c0", 
    "(2.5,5]" = "#fa9fb5", 
    "(5,10]" = "#feb24c", 
    "(10,20]" = "#fd8d3c", 
    "(20,30]" = "#fc4e2a", 
    "(30,40]" = "#e31a1c", 
    "(40,50]" =  "#bd0026", 
    "(50,100]" = "#49006a" 
    
    
  ))#+ggtitle(paste0("Plot averaged predicted mortality probability predicted from species-level models"))

cowplot::plot_grid(legend.large, mortality_percent_map_all,  align = "v", rel_widths = c(0.15, 0.89))+
  theme(panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA))
ggsave(height = 5, width = 8.5, dpi = 350, paste0(output.dir, "SPCD_stanoutput_full_standardized_v3/images/predicted_mortality/average_predicted_pmort_alltrees_map.png") )


######################################################################
# Map out probability of mortlaity from the hierarchical models ---- 
# read in the psurvival predictions from the hierarchical models 
psurvival.train <- readRDS(paste0(output.dir, "SPCD_stanoutput_joint/ll.train.pmort.RDS"))
psurvival.test <- readRDS(paste0(output.dir, "SPCD_stanoutput_joint/ll.test.pmort.RDS"))
psurvival.test$dataset <- "held-out"
psurvival.train$dataset <- "in-sample"
psurvival.all <-  rbind(psurvival.train, psurvival.test)



# calculate mortality
pmort.summaries.joint <-  psurvival.all %>% mutate(pmort.med = 1-median) %>% group_by(PLOT.ID, date, LAT_FIADB, LONG_FIADB, spp) %>% 
  summarise(mean.pmort = median(pmort.med), 
            pct.mortality = median(pmort.med)*100,
            sum.pmort = sum(pmort.med))
# generate bins
pmort.summaries.joint <- pmort.summaries.joint %>% mutate(Mort.quantiles = cut(pct.mortality, 
                                                                   breaks = c(0, 1, 2.5, 5, 10, 20, 30, 40, 50, 100), 
                                                                   include.lowest=TRUE))

large.size.font <- ggplot() +
  #geom_sf(alpha = 0.75, aes(fill = as.character(Distribution)))+
  #scale_fill_manual(values = c("Species Distribution" = "forestgreen", "Outside Distribution" = "white"))+
  geom_jitter(data = pmort.summaries, aes(x = LONG_FIADB, y = LAT_FIADB, color = Mort.quantiles), size = 5)+theme_bw()+
  coord_sf(xlim = c(-85, -68), ylim = c(35, 48))+theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'lightblue'), 
                                                       #legend.position = c(0.87, 0.25),
                                                       legend.title = element_blank(),
                                                       legend.background = element_rect(fill = "white", color = "black"))+
  scale_color_manual(values = c(
    "[0,1]" ="lightgrey",
    "(1,2.5]" = "#fcc5c0", 
    "(2.5,5]" = "#fa9fb5", 
    "(5,10]" = "#feb24c", 
    "(10,20]" = "#fd8d3c", 
    "(20,30]" = "#fc4e2a", 
    "(30,40]" = "#e31a1c", 
    "(40,50]" =  "#bd0026", 
    "(50,100]" = "#49006a" 
    
    
  ))


legend.large <- cowplot::get_legend(large.size.font)

mortality_percent_map_spp_joint <- ggplot() +
  geom_polygon(data = canada, 
               aes(x=long, y=lat, group = group), 
               color = "black", fill = "white") +
  geom_polygon(data = state_sub, 
               aes(x=long, y=lat, group = group), 
               color = "black", fill = "white") +
  #geom_sf(alpha = 0.75, aes(fill = as.character(Distribution)))+
  #scale_fill_manual(values = c("Species Distribution" = "forestgreen", "Outside Distribution" = "white"))+
  geom_jitter(data = pmort.summaries.joint, aes(x = LONG_FIADB, y = LAT_FIADB, color = Mort.quantiles), size = 0.5)+theme_bw()+
  coord_sf(xlim = c(-84, -65), ylim = c(37, 47.5))+theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'lightblue'), 
                                                         #legend.position = c(0.87, 0.25),
                                                         legend.title = element_blank(),
                                                         legend.background = element_rect(fill = "white", color = "black"), 
                                                         legend.position = "none")+
  scale_color_manual(values = c(
    "[0,1]" ="lightgrey",
    "(1,2.5]" = "#fcc5c0", 
    "(2.5,5]" = "#fa9fb5", 
    "(5,10]" = "#feb24c", 
    "(10,20]" = "#fd8d3c", 
    "(20,30]" = "#fc4e2a", 
    "(30,40]" = "#e31a1c", 
    "(40,50]" =  "#bd0026", 
    "(50,100]" = "#49006a" 
    
    
  ))+ggtitle(paste0("Plot averaged predicted mortality probability - hierarchical model"))

# png(height = 5, width = 8.5, res = 300, units = "in", here(paste0("images/predicted_mortality/average_predicted_pmort_SPCD_",SPCD.df[i,]$SPCD,"_map.png") ))
cowplot::plot_grid(legend.large, mortality_percent_map_spp_joint,  align = "v", rel_widths = c(0.15, 0.89))+
  theme(panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA))
ggsave(height = 5, width = 8.5, units = "in", here(paste0("images/predicted_mortality/average_predicted_pmort_joint_SPCD_map.png") ))
#dev.off()

cowplot::plot_grid(legend.large, cowplot::plot_grid(mortality_percent_map_all, mortality_percent_map_spp_joint + ggtitle(""), align = 
                                             "hv", 
                                           nrow = 2, 
                                           labels = "AUTO"),
                     align = "v", rel_widths = c(0.15, 0.89))+
  theme(panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA))
ggsave(height = 10, width = 8.5, units = "in", here(paste0("images/predicted_mortality/average_predicted_pmort_alltrees_map.png") ))





#-------------------------------------------------------------------------------------
# generate predictions for additional species

cleaned.data <- cleaned.data %>% filter(!is.na(ba) & !is.na(slope) & ! is.na(physio) & !is.na(aspect))%>% dplyr::select(state, county, pltnum, cndtn, point, tree, PLOT.ID, cycle, spp, dbhcur, dbhold, damage, Species, SPCD,
                                                                                                                        remper, LAT_FIADB, LONG_FIADB, elev, DIA_DIFF, annual.growth, M, relative.growth, si, physio:RD) %>% distinct()
plot.medians <- unique(cleaned.data %>% ungroup()%>% dplyr::select(PLOT.ID, si, ba, slope, aspect, MAP, MATmin, MATmax, damage.total, elev, Ndep.remper.avg, physio, RD)) %>% ungroup() %>% summarise(si.median = median(si, na.rm =TRUE), 
                                                                                                                                                                                                      RD.median = median(RD, na.rm = TRUE),
                                                                                                                                                                                                      ba.median = median(ba, na.rm =TRUE), 
                                                                                                                                                                                                      slope.median = median(slope, na.rm = TRUE), 
                                                                                                                                                                                                      aspect.median = median(aspect, na.rm = TRUE),
                                                                                                                                                                                                      damage.median = median(damage.total, na.rm =TRUE),
                                                                                                                                                                                                      elev.median = median(elev, na.rm =TRUE),
                                                                                                                                                                                                      Ndep.median = median(Ndep.remper.avg, na.rm =TRUE),
                                                                                                                                                                                                      physio.median = median(physio, na.rm = TRUE),
                                                                                                                                                                                                      
                                                                                                                                                                                                      MAP.median = median(MAP, na.rm =TRUE), 
                                                                                                                                                                                                      MATmin.median = median(MATmin, na.rm =TRUE), 
                                                                                                                                                                                                      MATmax.median = median(MATmax, na.rm =TRUE), 
                                                                                                                                                                                                      
                                                                                                                                                                                                      RD.sd = sd(RD, na.rm = TRUE),
                                                                                                                                                                                                      ba.sd = sd(ba, na.rm =TRUE),
                                                                                                                                                                                                      si.sd = sd(si, na.rm =TRUE), 
                                                                                                                                                                                                      slope.sd = sd(slope, na.rm =TRUE),
                                                                                                                                                                                                      aspect.sd = sd(aspect, na.rm = TRUE),
                                                                                                                                                                                                      damage.sd = sd(damage.total, na.rm =TRUE),
                                                                                                                                                                                                      elev.sd = sd(elev, na.rm =TRUE),
                                                                                                                                                                                                      Ndep.sd = sd(Ndep.remper.avg, na.rm =TRUE),
                                                                                                                                                                                                      physio.sd = sd(physio, na.rm = TRUE),
                                                                                                                                                                                                      
                                                                                                                                                                                                      MAP.sd = sd(MAP, na.rm =TRUE), 
                                                                                                                                                                                                      MATmin.sd = sd(MATmin, na.rm =TRUE), 
                                                                                                                                                                                                      MATmax.sd = sd(MATmax, na.rm =TRUE)
)

#View(cleaned.data %>% group_by(SPGRPCD, SPCD) %>% summarise(n()))

cleaned.data.full <- cleaned.data

# get all the observations from non-focal trees
SPCD.id <- nspp[18:126,]$SPCD

cleaned.data <- cleaned.data.full %>% filter(SPCD %in% SPCD.id)
remper.correction <- 0.5

# scale the cleaned data tree-level diameters by species?
cleaned.data <- cleaned.data %>% ungroup()  %>% group_by(SPCD) %>% 
  mutate(rempercur = ifelse(M ==1, remper*remper.correction, remper), 
         annual.growth = DIA_DIFF/rempercur) %>% mutate(DIA.median = median(dbhcur, na.rm =TRUE), 
                                                        DIA.sd = sd(dbhcur, na.rm = TRUE),  
                                                        BAL.median = median(BAL, na.rm=TRUE),
                                                        BAL.sd = sd(BAL, na.rm = TRUE),
                                                        RD.median = median(RD, na.rm=TRUE), 
                                                        RD.sd = sd(RD, na.rm =TRUE),
                                                        annual.growth.median = median(annual.growth, na.rm = TRUE), 
                                                        annual.growth.sd = sd(annual.growth, na.rm = TRUE)) %>% 
  ungroup() %>% mutate(DIA_scaled = (dbhcur - DIA.median)/DIA.sd, 
                       annual.growth.scaled = (annual.growth - annual.growth.median)/annual.growth.sd, 
                       RD.scaled = (RD-RD.median)/RD.sd,
                       BAL.scaled = (BAL-BAL.median)/BAL.sd,
                       si.scaled = (si - plot.medians$si.median)/plot.medians$si.sd, 
                       ba.scaled = (ba - plot.medians$ba.median)/plot.medians$ba.sd,
                       aspect.scaled = (aspect - plot.medians$aspect.median)/plot.medians$aspect.sd,
                       slope.scaled = (slope - plot.medians$slope.median)/plot.medians$slope.sd,
                       damage.scaled = (damage.total - plot.medians$damage.median)/plot.medians$damage.sd,
                       MAP.scaled = (MAP-plot.medians$MAP.median)/plot.medians$MAP.sd, 
                       elev.scaled = (elev-plot.medians$elev.median)/plot.medians$elev.sd,
                       Ndep.scaled = (Ndep.remper.avg- plot.medians$Ndep.median)/plot.medians$Ndep.sd,
                       physio.scaled = (physio-plot.medians$physio.median)/plot.medians$physio.sd,
                       MATmin.scaled = (MATmin- plot.medians$MATmin.median)/plot.medians$MATmin.sd,
                       MATmax.scaled = (MATmax - plot.medians$MATmax.median)/plot.medians$MATmax.sd)

SPP.df <- data.frame(SPCD = unique(cleaned.data$SPCD), 
                     SPP = 1:length(unique(cleaned.data$SPCD)))

cleaned.data <- left_join(cleaned.data, SPP.df) 
cleaned.data <- cleaned.data %>%  filter(!is.na(si) & !is.na(dbhcur)& !is.na(M) & !is.na(annual.growth.scaled) & !is.na(ppt.anom))
#summary(cleaned.data$BAL.scaled)
cleaned.data$S <- ifelse(cleaned.data$M == 1, 0, 1)


train.data <- cleaned.data

# model.6 data
# 6. All Fixed effects and all growth + Diameter interactions
mod.data.6 <- list(N = nrow(train.data), 
                   y = train.data$S, 
                   xM = as.matrix(train.data %>% dplyr::select(annual.growth.scaled, 
                                                               DIA_scaled, 
                                                               RD.scaled, 
                                                               ba.scaled, 
                                                               BAL.scaled, 
                                                               damage.scaled,
                                                               MATmax.scaled, 
                                                               MATmin.scaled, 
                                                               MAP.scaled,
                                                               ppt.anom, 
                                                               tmin.anom, 
                                                               tmax.anom, 
                                                               slope.scaled, 
                                                               aspect.scaled,
                                                               elev.scaled, 
                                                               Ndep.scaled,
                                                               physio.scaled) %>%
                                    # generate growth interactions
                                    mutate_at(.funs = list(growth.int = ~.*annual.growth.scaled), 
                                              .vars = vars(DIA_scaled:physio.scaled)) %>% 
                                    mutate_at(.funs = list(DIA.int = ~.*DIA_scaled), 
                                              .vars = vars(RD.scaled:physio.scaled)) ))


# generate predictions 
# read in the alpha population 
mualpha <- readRDS("SPCD_stanoutput_joint/alpha.p_model_6_1000samples.rds")
# read in the betas populations
mubetas <- readRDS("SPCD_stanoutput_joint/beta_model_6_1000samples.rds")

posterior.predict <- function(x){
 exp.x <- exp(-1*as.data.frame( mualpha) %>% dplyr::select(alpha) + rowSums(as.data.frame(mubetas) %>% dplyr::select(paste0("mu_beta[", 1:48,"]"))*mod.data.6$xM[x,]))
  p.mort <- exp.x / (1+ exp.x)
  p.mort$alpha
  }

pred.list <- lapply(1:mod.data.6$N, FUN = posterior.predict)
pmort.hat.oos <- do.call(cbind, pred.list)
pmort.hat.oos.m <- reshape2::melt(pmort.hat.oos)
colnames(pmort.hat.oos.m) <- c("sample", "treeid", "pmort")
pmort.hat.oos.summary <- pmort.hat.oos.m %>% group_by(treeid) %>% summarise(median = quantile(pmort, 0.5, na.rm =TRUE), 
                                                        ci.lo = quantile(pmort, 0.025, na.rm =TRUE), 
                                                        ci.hi = quantile(pmort, 0.975, na.rm =TRUE))

pmort.hat.oos.summary$S <- train.data$S
pmort.hat.oos.summary$SPCD <- train.data$SPCD
pmort.hat.oos.summary$COMMON <-  FIESTA::ref_species[match(pmort.hat.oos.summary$SPCD, FIESTA::ref_species$SPCD),]$COMMON
pmort.hat.oos.summary$GENUS <-  FIESTA::ref_species[match(pmort.hat.oos.summary$SPCD, FIESTA::ref_species$SPCD),]$GENUS
pmort.hat.oos.summary$SPECIES <-  FIESTA::ref_species[match(pmort.hat.oos.summary$SPCD, FIESTA::ref_species$SPCD),]$SPECIES


# AUC using mltools auc_roc function
# for in sample data
actuals = train.data$S
preds = as.vector(pmort.hat.oos.summary$median)
auc.oos.nonfocal.spp <-auc_roc(preds, actuals)


auc.oos.spp <- pmort.hat.oos.summary %>% group_by(SPCD, COMMON) %>% summarise(`AUC` = auc_roc(median, S), 
                                                                              `# observations` = n())

hist(auc.oos.spp$AUC)

auc.oos.spp %>% filter(!is.na(AUC)) %>% filter(AUC > 0.5)
auc.oos.spp %>% filter(!is.na(AUC)) %>% filter(AUC > 0.5) %>% ungroup()%>% summarise(total.trees = sum(`# observations`))
auc.oos.spp %>% filter(!is.na(AUC)) %>% filter(AUC < 0.5)
auc.oos.spp %>% filter(!is.na(AUC)) %>% filter(AUC < 0.5) %>% ungroup()%>% summarise(total.trees = sum(`# observations`))


auc.oos.spp$GENUS <-  FIESTA::ref_species[match(auc.oos.spp$SPCD, FIESTA::ref_species$SPCD),]$GENUS
auc.oos.spp$SPECIES <-  FIESTA::ref_species[match(auc.oos.spp$SPCD, FIESTA::ref_species$SPCD),]$SPECIES

nspp$GENUS <- FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD),]$GENUS
nSPP.GENUS <- nspp[1:17,]$GENUS

auc.oos.spp %>% filter(!is.na(AUC)) %>% arrange(desc(AUC)) %>% ungroup()|> gt() %>%
  data_color(
    columns = GENUS,
    rows =  GENUS %in% nSPP.GENUS,
    method = "factor",
    #palette = c("red", "green"),
    #domain = c(0, 50)
  )%>%
  gtsave("tab_1.html", path = getwd())

auc.oos.spp %>% filter(!is.na(AUC)) %>% arrange(desc(AUC)) %>% group_by(AUC < 0.5, GENUS %in% nSPP.GENUS) %>%
  summarise(`total trees` = sum(`# observations`)) %>% ungroup()|> gt() %>%
  # data_color(
  #   columns = GENUS,
  #   rows =  GENUS %in% nSPP.GENUS,
  #   method = "factor",
  #   #palette = c("red", "green"),
  #   #domain = c(0, 50)
  # )%>%
  gtsave("summary_oos_species.html", path = getwd())

