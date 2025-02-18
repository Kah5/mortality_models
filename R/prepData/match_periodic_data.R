# script to match up TREEs between surveys
library(tidyverse)
library(gt)
library(odbc)
# load the aggregated tree file we are working with
boxdir <- "C:/Users/KellyHeilman/Box/01. kelly.heilman Workspace/mortality/Eastern-Mortality/mortality_models/"

TREE <- read_delim(paste0(boxdir, "data/formatted_older_matching_plts_TREE.txt"))
TREE$dbhcur <- as.numeric(TREE$dbhcur)
TREE$dbhold <- as.numeric(TREE$dbhold)


TREE$SPCD <- TREE$spp
TREE$Species <- FIESTA::ref_species[match(TREE$SPCD, FIESTA::ref_species$SPCD),]$COMMON

# for a state with multiple measurements, see if we can match up some trees:

# for Maine, we have multiple measurement years
TREE %>% filter(stname %in% "ME") %>% group_by(date) %>% summarise(n())
ME.1982 <- TREE %>% filter(stname %in% "ME" & date == 1982)
ME.1995 <- TREE %>% filter(stname %in% "ME" & date == 1995)

# check that tree level variables match up
unique(ME.1982$state)
unique(ME.1995$state)
# there is an extra county in maine 1982 (county 103)
unique(ME.1995$county)
# 1  3  5  7  9 11 13 15 17 19 21 23 25 27 29 31
unique(ME.1982$county)
# 29  19   3 103   9  11  13  15  27  21  25   7  17   1   5  23  31

# some plotnumbers match
unique(ME.1982$pltnum) %in% unique(ME.1995$pltnum)

# condition ids only exist in ME 1995
unique(ME.1982$cndtn)
unique(ME.1995$cndtn)

# all of 1982 has point == 0
unique(ME.1982$point)

# all of 1995 has a value for point
unique(ME.1995$point)

# tree numbers
unique(ME.1982$tree) # vary from 1-99
unique(ME.1995$tree) # varies from 0-202?


# concatenate the all of first tree-level variables
TREE <- TREE %>% mutate(TREE.ID = paste0(state,"_", county, "_", pltnum,"_", cndtn,"_", point, "_", tree))
ME.1982 <- TREE %>% filter(stname %in% "ME" & date == 1982)
ME.1995 <- TREE %>% filter(stname %in% "ME" & date == 1995)

# no matches when using all of the variables--NA in cndtn and point values are different
summary(ME.1982$TREE.ID %in% ME.1995$TREE.ID)

# what if we match over everything except for the point and condition?
TREE <- TREE %>% mutate(TREE.ID = paste0(state,"_", county, "_", pltnum, "_", tree))
ME.1982 <- TREE %>% filter(stname %in% "ME" & date == 1982)
ME.1995 <- TREE %>% filter(stname %in% "ME" & date == 1995)

# we get some matches
summary(ME.1982$TREE.ID %in% ME.1995$TREE.ID)
# Mode   FALSE    TRUE 
# logical   57497   38523 

# filter all trees from 1995 that match up to ME.1982:
ME.1995.m <- ME.1995 %>% filter(TREE.ID %in% ME.1982$TREE.ID)


# for all of those that match, check that the dbhold in 1995 matches dbhcur in 1982:
ME.1982.s <- ME.1982 %>% select(state, county, pltnum, tree, TREE.ID, status, SPCD, Species, dbhcur, dbhold) %>%
  rename("dbh_t2" = "dbhcur", 
         "dbh_t1" = "dbhold", 
         "status_t2" = "status")
ME.1995.s <- ME.1995 %>% select(state, county, pltnum, tree, TREE.ID, status, SPCD, Species, dbhcur, dbhold) %>%
  rename("dbh_t3" = "dbhcur", 
         "dbh_t2" = "dbhold", 
         "status_t3" = "status")
ME.time <- left_join(ME.1982.s, ME.1995.s)

# only 2422 trees have a 3rd matching tree and 93598 trees have a NA or non-matching tree
ME.time %>% group_by(is.na(dbh_t3)) %>% summarise(n())

# if the 1995 tree matches 1982, does the status code match or change
ME.time %>% filter(!is.na(dbh_t3)) %>% 
  mutate(status_change = ifelse(status_t2 == status_t3, "same status", "status change")) %>%
  group_by(status_change) %>% summarise(n())

# see if all the dead trees are a status code change

ME.time %>% filter(!is.na(dbh_t3)) %>% 
  mutate(status_change = ifelse(status_t2 > 1 & status_t3 > 1, "dead at both times", 
                                ifelse(status_t2 == 1 & status_t3 == 1, "live at both times", 
                                   ifelse(status_t2 == 1 & status_t3 > 1, "live to dead status change", "marked as dead to live")))) %>%
  group_by(status_change) %>% summarise(`# of trees` = n())|> gt() |> 
  gtsave(  "images/filtering_exploration/Maine_1982_1995_status_cd_change_ntrees_matching.png")

# only 37 trees show that they are dead at both times
View(ME.time %>% filter(status_t2 == 1 & status_t3 >1))

ME.time.diff <- ME.time %>% filter(!is.na(dbh_t3)) %>% 
  mutate(status_change = ifelse(status_t2 > 1 & status_t3 > 1, "dead at both times", 
                                ifelse(status_t2 == 1 & status_t3 == 1, "live at both times", 
                                       ifelse(status_t2 == 1 & status_t3 > 1, "live to dead status change",  "marked as dead to live"))))%>%
  mutate(DIA_diff = dbh_t3 - dbh_t2)


summary(ME.time.diff$DIA_diff)
ME.time.diff %>% ggplot(., aes(DIA_diff))+geom_histogram()+facet_wrap(~status_change, scales = "free")
ggsave("images/filtering_exploration/Maine_1982_1995_status_cd_change_dia_diff.png")

# group by species:

ME.time %>% filter(!is.na(dbh_t3)) %>% 
  mutate(status_change = ifelse(status_t2 > 1 & status_t3 > 1, "dead at both times", 
                                ifelse(status_t2 == 1 & status_t3 == 1, "live at both times", 
                                       ifelse(status_t2 == 1 & status_t3 > 1, "live to dead status change", "marked as dead to live")))) %>%
  group_by( Species, status_change) %>% summarise(ntrees = n()) %>%
  mutate(ntrees = ifelse(is.na(ntrees), 0, ntrees)) %>% 
  spread(status_change, ntrees) %>% ungroup()|> 
  gt()  |> 
  gtsave(  "images/filtering_exploration/Maine_1982_1995_status_cd_change_ntrees_matching_species.png")

#######################################################################################################
# lets check the other states that have multiple remeasurements
#######################################################################################################
TREE %>% group_by(stname) %>% select(stname, date) %>% distinct() %>% 
  group_by(stname) %>% summarise(n.measurements = n())

# four of the 13 states have two measurements:
stnames.remeasure = data.frame(stname = c("ME", "NH", "NY", "VT"), 
                               t2_date = c(1982, 1983, 1980, 1983), 
                               t3_date = c(1995, 1997, 1993, 1997))

# get the t2 earlier measurement 
TREE.t2 <- TREE %>% filter(stname %in% stnames.remeasure$stname ) %>% left_join(., stnames.remeasure) %>% 
  filter(date == t2_date) %>% select(stname, state, county, pltnum, tree, TREE.ID, status, SPCD, Species, dbhcur, dbhold, t2_date) %>%
                              rename("dbh_t2" = "dbhcur", 
                                     "dbh_t1" = "dbhold", 
                                     "status_t2" = "status")  

# get the later measurement
TREE.t3 <- TREE %>% filter(stname %in% stnames.remeasure$stname ) %>% left_join(., stnames.remeasure)%>% 
  filter(date == t3_date) %>% select(stname, state, county, pltnum, tree, TREE.ID, status, SPCD, Species, dbhcur, dbhold, t3_date) %>%
                          rename("dbh_t3" = "dbhcur", 
                                 "dbh_t2" = "dbhold", 
                                 "status_t3" = "status")

# join the two tables together:
TREE.time <- left_join(TREE.t2, TREE.t3)

# only 8652 trees have a 3rd matching tree and 261015 trees have a NA or non-matching tree
TREE.time %>% mutate(t2_t3_remeasurements = ifelse(is.na(dbh_t3), "no linked t2 - t3 measurements", "linked t2 - t3 measurements")) %>% 
  group_by(t2_t3_remeasurements, stname) %>% summarise(n()) %>% 
  spread(t2_t3_remeasurements , `n()`) |> gt()|>
  grand_summary_rows(
    columns = c("no linked t2 - t3 measurements","linked t2 - t3 measurements"),
    fns = list(
      total ~ sum(., na.rm =TRUE)
    ),
    fmt = ~ fmt_number(., use_seps = FALSE, decimals = 0)
  )|>
  gtsave(  "images/filtering_exploration/ALL_states_t2_t3_ntrees_matching.png")

# # A tibble: 2 × 2
# `is.na(dbh_t3)`  `n()`
# <lgl>            <int>
# 1 FALSE           8652
# 2 TRUE            261015

# check for status code changes:
TREE.time %>% filter(!is.na(dbh_t3)) %>% 
  mutate(status_change = ifelse(status_t2 == status_t3, "same status", "status change")) %>%
  group_by(status_change, stname) %>% summarise(n()) %>% spread(status_change, `n()`)

# see if all the dead trees are a status code change:

TREE.time %>% filter(!is.na(dbh_t3)) %>% 
  mutate(status_change = ifelse(status_t2 > 1 & status_t3 > 1, "dead at both times", 
                                ifelse(status_t2 == 1 & status_t3 == 1, "live at both times", 
                                       ifelse(status_t2 == 1 & status_t3 > 1, "live to dead status change", "marked as dead to live")))) %>%
  group_by(status_change, stname) %>% summarise(`# of trees` = n()) %>% spread(status_change, `# of trees`) %>% ungroup()|> gt() |> 
  grand_summary_rows(
    columns = c("dead at both times", "live at both times", "live to dead status change", "marked as dead to live"),
    fns = list(
      total ~ sum(., na.rm =TRUE)
    ),
    fmt = ~ fmt_number(., use_seps = FALSE, decimals = 0)
  )|>
  gtsave(  "images/filtering_exploration/ALL_states_t2_t3_status_cd_change_ntrees_matching.png")

# only 37 trees show that they are dead at both times
View(TREE.time %>% filter(status_t2 == 1 & status_t3 >1))

TREE.time.diff <- TREE.time %>% filter(!is.na(dbh_t3)) %>% 
  mutate(status_change = ifelse(status_t2 > 1 & status_t3 > 1, "dead at both times", 
                                ifelse(status_t2 == 1 & status_t3 == 1, "live at both times", 
                                       ifelse(status_t2 == 1 & status_t3 > 1, "live to dead status change",  "marked as dead to live"))))%>%
  mutate(DIA_diff = dbh_t3 - dbh_t2)


summary(TREE.time.diff$DIA_diff)
TREE.time.diff %>% ggplot(., aes(DIA_diff))+geom_histogram()+
  facet_grid(cols = vars(status_change), rows = vars(stname), scales = "free")
ggsave("images/filtering_exploration/All_states_t2_t3_status_cd_change_dia_diff.png", height = 6, width = 10)

# group by species:

TREE.time %>% filter(!is.na(dbh_t3)) %>% 
  mutate(status_change = ifelse(status_t2 > 1 & status_t3 > 1, "dead at both times", 
                                ifelse(status_t2 == 1 & status_t3 == 1, "live at both times", 
                                       ifelse(status_t2 == 1 & status_t3 > 1, "live to dead status change", "marked as dead to live")))) %>%
  group_by( Species, status_change) %>% summarise(ntrees = n()) %>%
  mutate(ntrees = ifelse(is.na(ntrees), 0, ntrees)) %>% 
  spread(status_change, ntrees) %>% ungroup()|> 
  gt()  |> 
  grand_summary_rows(
    columns = c("dead at both times", "live at both times", "live to dead status change", "marked as dead to live"),
    fns = list(
      total ~ sum(., na.rm =TRUE)
    ),
    fmt = ~ fmt_number(., use_seps = FALSE, decimals = 0)
  )|>
  gtsave(  "images/filtering_exploration/all_trees_t2_t3_status_cd_change_ntrees_matching_species.png")


######################################################################################################
# try matching up trees that dont directly match up the first way
######################################################################################################
# match TREE.ID.T2 = stname_county_pltnum_species_dbhold_T2date from the later inventory to:
# match TREE.ID.T2 = stname_county_pltnum_species_dbhcur_T2date from the earlier inventory
# lets call this new ids 


# four of the 13 states have two measurements:
stnames.remeasure = data.frame(stname = c("ME", "NH", "NY", "VT"), 
                               t2_date = c(1982, 1983, 1980, 1983), 
                               t3_date = c(1995, 1997, 1993, 1997))

# get the t2 earlier measurement and create TREE.ID that matches based on species and diameter, not tree id
TREE.t2 <- TREE %>% filter(stname %in% stnames.remeasure$stname ) %>% left_join(., stnames.remeasure) %>% 
  filter(date == t2_date) %>% select(stname, state, county, pltnum, tree, TREE.ID, status, SPCD, Species, dbhcur, dbhold, t2_date) %>%
  rename("dbh_t2" = "dbhcur", 
         "dbh_t1" = "dbhold", 
         "status_t2" = "status") %>% 
  mutate(TREE.ID.T2 = paste0(stname, "_", county, "_", pltnum, "_", SPCD, "_", dbh_t2))

# get the later measurement
TREE.t3 <- TREE %>% filter(stname %in% stnames.remeasure$stname ) %>% left_join(., stnames.remeasure)%>% 
  filter(date == t3_date) %>% select(stname, state, county, pltnum, tree, TREE.ID, status, SPCD, Species, dbhcur, dbhold, t3_date) %>%
  rename("dbh_t3" = "dbhcur", 
         "dbh_t2" = "dbhold", 
         "status_t3" = "status")%>% 
  mutate(TREE.ID.T2 = paste0(stname, "_", county, "_", pltnum, "_", SPCD, "_", dbh_t2))

# join the two tables together:
TREE.time.spcd.dbh <- left_join(TREE.t2, TREE.t3)

TREE.time.spcd.dbh %>% mutate(t2_t3_remeasurements = ifelse(is.na(dbh_t3), "no linked t2 - t3 measurements", "linked t2 - t3 measurements")) %>% 
  group_by(t2_t3_remeasurements, stname) %>% summarise(n()) %>% 
  spread(t2_t3_remeasurements , `n()`) |> gt()|>
  grand_summary_rows(
    columns = c("no linked t2 - t3 measurements","linked t2 - t3 measurements"),
    fns = list(
      total ~ sum(., na.rm =TRUE)
    ),
    fmt = ~ fmt_number(., use_seps = FALSE, decimals = 0)
  )|>
  gtsave(  "images/filtering_exploration/ALL_states_t2_t3_ntrees_matching_spcd_diameter.png")



TREE.time.spcd.dbh.diff <- TREE.time.spcd.dbh %>% filter(!is.na(dbh_t3)) %>% 
  mutate(status_change = ifelse(status_t2 > 1 & status_t3 > 1, "dead at both times", 
                                ifelse(status_t2 == 1 & status_t3 == 1, "live at both times", 
                                       ifelse(status_t2 == 1 & status_t3 > 1, "live to dead status change",  "marked as dead to live"))))%>%
  mutate(DIA_diff = dbh_t3 - dbh_t2)


summary(TREE.time.diff$DIA_diff)
TREE.time.spcd.dbh.diff %>% ggplot(., aes(DIA_diff))+geom_histogram()+
  facet_grid(cols = vars(status_change), rows = vars(stname), scales = "free")
ggsave("images/filtering_exploration/All_states_t2_t3_spcd_dbh_match_status_cd_change_dia_diff.png", height = 6, width = 10)

# group by species:

TREE.time.spcd.dbh %>% filter(!is.na(dbh_t3)) %>% 
  mutate(status_change = ifelse(status_t2 > 1 & status_t3 > 1, "dead at both times", 
                                ifelse(status_t2 == 1 & status_t3 == 1, "live at both times", 
                                       ifelse(status_t2 == 1 & status_t3 > 1, "live to dead status change", "marked as dead to live")))) %>%
  group_by( Species, status_change) %>% summarise(ntrees = n()) %>%
  mutate(ntrees = ifelse(is.na(ntrees), 0, ntrees)) %>% 
  spread(status_change, ntrees) %>% ungroup()|> 
  gt()  |> 
  grand_summary_rows(
    columns = c("dead at both times", "live at both times", "live to dead status change", "marked as dead to live"),
    fns = list(
      total ~ sum(., na.rm =TRUE)
    ),
    fmt = ~ fmt_number(., use_seps = FALSE, decimals = 0)
  )|>
  gtsave(  "images/filtering_exploration/all_trees_t2_t3_spcd_dbh_match_status_cd_change_ntrees_matching_species.png")


#################################################################################
# see if there are additional surveys in the fiadb
#################################################################################
ocon <- dbConnect(odbc(), "fiadb01p")


NE_plot <- dbGetQuery(ocon, "SELECT cn, statecd, unitcd, countycd, plot, invyr, plot_status_cd, cycle, lat, lon, elev, designcd
                      FROM fs_fiadb.plot
                      WHERE statecd =  ANY(09, 10, 23, 24, 25, 33, 34, 36, 39, 42, 44, 50, 54) 
                      and invyr < 2000") %>%
  as_tibble()%>%
  rename_with(tolower)
unique(NE_plot$statecd)
length(unique(NE_plot$cn))


unique(NE_plot$statecd)
length(unique(NE_plot$cn))

# when done, disconnect from ORACLE FIADB
dbDisconnect(ocon)
rm(ocon)

# see how many inventories wer have for each state
View(NE_plot %>% select(statecd, invyr)%>% distinct())
NE_plot %>% select(statecd, invyr, cn) %>% group_by(statecd, invyr) %>% summarise(n()) %>%
  spread(invyr, `n()`)

NE_plot <- NE_plot %>% mutate(PLOT_ID = paste0(statecd, "_", countycd, "_", plot))
matching.plots <- NE_plot %>% group_by(PLOT_ID) %>% summarise(n())
summary(matching.plots$`n()`)
hist(matching.plots$`n()`)
multi.meas.plot <- matching.plots %>% filter(`n()` > 1 )
NE_plot_remeas <- NE_plot %>% filter(PLOT_ID %in% unique(multi.meas.plot$PLOT_ID))

# get the trees from FIA_DB for the matching plots
ocon <- dbConnect(odbc(), "fiadb01p")


NE_tree <- dbGetQuery(ocon, "SELECT cn, prev_tre_cn, statecd, unitcd, countycd, plot, subp, tree, invyr, dia, statuscd, spcd
                      FROM fs_fiadb.tree
                      WHERE statecd =  ANY(09, 10, 23, 24, 25, 33, 34, 36, 39, 42, 44, 50, 54) 
                      and invyr < 2000") %>%
  as_tibble()%>%
  rename_with(tolower)

# get the trees in plots with two measurements
NE_tree_remeas <- NE_tree %>% mutate(PLOT_ID = paste0(statecd, "_", countycd, "_", plot)) %>% 
  filter(PLOT_ID %in% unique(multi.meas.plot$PLOT_ID))

unique(NE_tree_remeas$prev_tre_cn) # no prev_tre_cn numbers are included
# lets assign each an inventory year based on the state code. T2 == early inventory and T3 == later inventory
NE.survey.years <- unique(NE_tree_remeas %>% select(statecd, invyr))
NE.survey.years$INV_T <- rep(c("T2", "T3"), 9) # it looks like they are ordered properly
NE_tree_remeas <- left_join(NE_tree_remeas, NE.survey.years)

# make a unique tree id and see if we can match based on those
NE_tree_remeas <- NE_tree_remeas %>% mutate(TREE.ID = paste0(statecd, "_", countycd, "_", plot, "_", tree, "_", spcd))
NE_tree_T2 <- NE_tree_remeas %>% filter(INV_T %in% "T2") %>% rename("dbh_t2" = "dia",
                                                                    "status_t2" = "statuscd", 
                                                                    "invyr_t2" = "invyr")
NE_tree_T3 <- NE_tree_remeas %>% filter(INV_T %in% "T3")%>% rename("dbh_t3" = "dia",
                                                                   "status_t3" = "statuscd", 
                                                                   "invry_t3" = "invyr")

# are these tree ids unique?
NE_tree_T2[duplicated(NE_tree_T2),] # no duplicates
NE_tree_T3[duplicated(NE_tree_T3),] # no duplicates

# looks like there are some matches
unique(NE_tree_T2$TREE.ID)[1] %in% unique(NE_tree_T3$TREE.ID)

NE_tree_joined <- left_join(NE_tree_T2 %>% select(-cn, -subp, -prev_tre_cn, -INV_T), 
                            NE_tree_T3 %>% select(-cn, -subp, -prev_tre_cn, -INV_T))
NE_tree_joined_remeas <- NE_tree_joined %>% filter(!is.na(invry_t3))
length(unique(NE_tree_joined_remeas$TREE.ID))

NE_tree_T3 %>% filter(TREE.ID %in% unique(NE_tree_T2$TREE.ID)[1])
NE_tree_T2 %>% filter(TREE.ID %in% unique(NE_tree_T2$TREE.ID)[1])


NE_tree_T3 %>% filter(TREE.ID %in% unique(NE_tree_T2$TREE.ID)[2])
NE_tree_T2 %>% filter(TREE.ID %in% unique(NE_tree_T2$TREE.ID)[2])

# see if there are trees T2 that match to multiple trees at another T3
n.match.list <- lapply(1:length(unique(NE_tree_joined_remeas$TREE.ID)), function(x){
  n.match <- nrow(NE_tree_joined_remeas %>% filter(TREE.ID %in% unique(NE_tree_joined_remeas$TREE.ID)[x]))
  n.match.df <- data.frame(TREE.ID = unique(NE_tree_joined_remeas$TREE.ID)[x], 
                          num.match = n.match)  
  n.match.df  
})

n.match.df.all <- do.call(rbind, n.match.list)
n.match.df.all %>% group_by(num.match) %>% summarise (n()) |> gt()|>
  gtsave(  "images/filtering_exploration/FIADB_t2_t3_spcd_dbh_match_ntrees.png")

# of these ones that match only one other tree, how many have a status code change?
single.match <- n.match.df.all %>% filter(num.match == 1)
NE_tree_joined_remeas %>% filter(TREE.ID %in% unique(single.match$TREE.ID)) %>% 
  mutate(status.change = ifelse(status_t2 == status_t3, "live at both times", 
                                ifelse(status_t2 == 1 & status_t3 > 1, "live at t2, dead at t3", NA))) %>%
  group_by(status.change) %>% summarise(n())

double.match <- n.match.df.all %>% filter(num.match == 2)
NE_tree_joined_remeas %>% filter(TREE.ID %in% unique(double.match$TREE.ID)) %>% 
  mutate(status.change = ifelse(status_t2 == status_t3, "live at both times", 
                                ifelse(status_t2 == 1 & status_t3 > 1, "live at t2, dead at t3", NA))) %>%
  group_by(status.change) %>% summarise(n())

# to do:
# look at duplicated matches 
# decide which trees match