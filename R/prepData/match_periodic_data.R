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


ocon <- dbConnect(odbc(), "fiadb01p")

NE_designcd <- dbGetQuery(
  ocon,
  "SELECT s.ann_inventory, p.statecd, p.invyr, t.subp, p.designcd, COUNT(*)

FROM fs_fiadb.survey s

JOIN fs_fiadb.plot p

ON (s.cn = p.srv_cn)

JOIN fs_fiadb.tree t

ON (p.cn = t.plt_cn)

WHERE s.rscd = 24

AND s.p3_ozone_ind = 'N'

GROUP BY s.ann_inventory, p.statecd, p.invyr, t.subp, p.designcd

ORDER BY p.statecd, p.invyr, t.subp, p.designcd"
)

# notes state-by state
# state 9 connecticut
NE_designcd %>% filter(STATECD == 9)
#SUBP 101, designcd 101 match 1985 to 1998

# state 10 delaware
NE_designcd %>% filter(STATECD == 10)
# designcd 101 for 1986-- subps 100-105
# for 1999 designcd are 111, 112, 113, 1, 115--try matching up to 

# state 23 Maine
NE_designcd %>% filter(STATECD == 23)
# designcd 101 for 1995 and 1999-- subps 101
# for 1999 designcd are 1, 115, 116, 117,--116 and 17 for suplot 1 may be remeasurements? of 101

# state 24 Maryland
NE_designcd %>% filter(STATECD == 24)
# designcd 101 for 1986 -- subps 101-105
# for 1999 designcd are 111, 112, 113 for suplot 1 may be remeasurements? of 101

# state 24 Massachusetts
NE_designcd %>% filter(STATECD == 25)
# designcd 101 for 1985 -- subps 101-105
# designcd 111 for 1998 -- subs 101-104
# for 2003 designcd are 1 and 115 for suplot 1 may be remeasurements?? 115 is overlaid on FHM 4 subp design

# state 33 New Hampshire
NE_designcd %>% filter(STATECD == 33)
# designcd 101 for 1983 -- subps 100
# designcd 111 for 1997 -- subs 101-104
# for 2002 designcd are 1 and 115 for suplot 1 may be remeasurements?? 115 is overlaid on FHM 4 subp design

# state 34 New Jersey
NE_designcd %>% filter(STATECD == 34)
# designcd 101 for 1987 -- subps 101-114
# designcd 111, 112, 113 for 1999 -- subs 1-4
# for 2004 designcd are 1 and 115 for suplot 1 may be remeasurements?? 115 is overlaid on FHM 4 subp design

# state 36 New York
NE_designcd %>% filter(STATECD == 36)
# designcd 101 for 1993 -- subps 101-120
# for 2002 designcd are 1, 117and 115 for suplot 1 may be remeasurements?? 115 is overlaid on FHM 4 subp design

# state 39 Ohio
NE_designcd %>% filter(STATECD == 39)
# designcd 101 for 1991 -- subps 101-111
# for 2001 designcd are 1, 116, 117, 118 and 115 for suplot 1 may be remeasurements?? 115 is overlaid on FHM 4 subp design

# state 42 Pennsylvania
NE_designcd %>% filter(STATECD == 42)
# designcd 101 for 1989 -- subps 101-106
# for 2000 designcd are 1, 116, 117, and 115 for suplot 1 may be remeasurements?? 115 is overlaid on FHM 4 subp design

# state 50 Vermont
NE_designcd %>% filter(STATECD == 50)
# designcd 101 for 1983 -- subps 100
# designcd 111 for 1997 subp 101-104
# for 2003 designcd are 1,and 115 for suplot 1 may be remeasurements?? 115 is overlaid on FHM 4 subp design

# state 54 West Virginia
NE_designcd %>% filter(STATECD == 54)
# designcd 101 for 1989 -- subps 101-120
# for 2000 designcd are 111,112,and 113 for suplot 1 may be remeasurements?? 

##########################################
# lets start with # state 9 connecticut
NE_designcd %>% filter(STATECD == 9) 
#SUBP 101, designcd 101 match 1985 to 1998


ocon <- dbConnect(odbc(), "fiadb01p")

CT_plot <- dbGetQuery(ocon, 
                      "SELECT cn, statecd, unitcd, countycd, plot, invyr, designcd
FROM fs_fiadb.plot
WHERE statecd =  ANY(09) 
                      and invyr < 2000
") %>%
  as_tibble()%>%
  rename_with(tolower) %>% 
  rename("plt_cn" = "cn") %>%
  mutate(PLOT.ID = paste0(statecd, "_", unitcd, "_", countycd, "_", unitcd, "_", plot))
# get CT_plot
CN.plot <- CT_plot %>% mutate(PLOT.ID = paste0(statecd, "_", unitcd, "_", countycd, "_", unitcd, "_", plot))%>% group_by(PLOT.ID) %>% summarise(n())

CT_tree <- dbGetQuery(ocon, "SELECT cn, plt_cn, prev_tre_cn, statecd, unitcd, countycd, plot, subp, tree, invyr, dia, statuscd, spcd
                      FROM fs_fiadb.tree
                      WHERE statecd =  ANY(09) 
                      and invyr < 2000
                      ") %>%
  as_tibble()%>%
  rename_with(tolower)%>% left_join(., CT_plot)


# subplot 101 may be remeasured between 1985 and 1998
CT_tree %>% filter(PLOT.ID %in% unique(CT_tree$PLOT.ID)[1]) %>% group_by(subp, invyr) %>% summarise(n())
CT_tree %>% filter(PLOT.ID %in% unique(CT_tree$PLOT.ID)[1]) %>% filter(subp == 101) #%>% summarise(n())

# it looks like the tree ids match, so perhaps joining on PLOT, subp and tree in CT would work
CT_tree_T2 <- CT_tree %>% filter(invyr == 1985) %>% rename("dbh_t2" = "dia",
                                                           "status_t2" = "statuscd", 
                                                           "invyr_t2" = "invyr", 
                                                           "designcd_t2" = "designcd")
CT_tree_T3 <- CT_tree %>% filter(invyr == 1998)%>% rename("dbh_t3" = "dia",
                                                          "status_t3" = "statuscd", 
                                                          "invry_t3" = "invyr", 
                                                          "designcd_t3" = "designcd")

CT_tree_remeas <- left_join(CT_tree_T2 %>% select(-cn, -plt_cn, -prev_tre_cn), 
                            CT_tree_T3%>% select(-cn, -plt_cn, -prev_tre_cn))

# all the trees are live at both times in CT
CT_tree_remeas %>% filter(!is.na(invry_t3)) %>%  mutate(status.change = ifelse(status_t2 == status_t3, "live at both times", 
                                                                               ifelse(status_t2 == 1 & status_t3 > 1, "live at t2, dead at t3", NA))) %>%
  group_by(status.change) %>% summarise(n())



##########################################
# lets start with  state 10 delaware
NE_designcd %>% filter(STATECD == 10)
# designcd 101 for 1986-- subps 100-105
# for 1999 designcd are 111, 112, 113, 1, 115--try matching up to 


ocon <- dbConnect(odbc(), "fiadb01p")
DE_plot <- dbGetQuery(ocon, 
                      "SELECT cn, statecd, unitcd, countycd, plot, invyr, designcd
FROM fs_fiadb.plot
WHERE statecd =  ANY(10) 
                      and invyr < 2000
") %>%
  as_tibble()%>%
  rename_with(tolower) %>% 
  rename("plt_cn" = "cn") %>%
  mutate(PLOT.ID = paste0(statecd, "_", unitcd, "_", countycd, "_", unitcd, "_", plot))
# get CT_plot
CN.plot <- DE_plot %>% mutate(PLOT.ID = paste0(statecd, "_", unitcd, "_", countycd, "_", unitcd, "_", plot))%>% group_by(PLOT.ID) %>% summarise(n())

DE_tree <- dbGetQuery(ocon, "SELECT cn, plt_cn, prev_tre_cn, statecd, unitcd, countycd, plot, subp, tree, invyr, dia, statuscd, spcd
                      FROM fs_fiadb.tree
                      WHERE statecd =  ANY(10) 
                      and invyr < 2000
                      ") %>%
  as_tibble()%>%
  rename_with(tolower)%>% left_join(., DE_plot)


# subplot 101 may be remeasured between 1985 and 1998
DE_tree %>% filter(PLOT.ID %in% unique(DE_tree$PLOT.ID)[3]) %>% group_by(subp, invyr) %>% summarise(n())
View(DE_tree %>% filter(PLOT.ID %in% unique(DE_tree$PLOT.ID)[3]))  #%>% summarise(n())

# it the tree ids do not match, so we would need to change tree id 
DE_tree_T2 <- DE_tree %>% filter(invyr == 1986) %>% rename("dbh_t2" = "dia",
                                                           "status_t2" = "statuscd", 
                                                           "invyr_t2" = "invyr", 
                                                           "designcd_t2" = "designcd")%>% 
  mutate(subp.match = ifelse(subp == 101, 1, NA), 
         tree.match = tree - 10000)
DE_tree_T3 <- DE_tree %>% filter(invyr == 1999)%>% rename("dbh_t3" = "dia",
                                                          "status_t3" = "statuscd", 
                                                          "invry_t3" = "invyr", 
                                                          "designcd_t3" = "designcd") %>% 
  mutate(subp.match = ifelse(subp == 1, 1, NA), 
         tree.match = tree)

DE_tree_remeas <- left_join(DE_tree_T2 %>% select(-cn, -plt_cn, -prev_tre_cn, -subp, -tree), 
                            DE_tree_T3 %>% select(-cn, -plt_cn, -prev_tre_cn, -subp, -tree))

# no matches
DE_tree_remeas %>% filter(!is.na(invry_t3)) %>%  mutate(status.change = ifelse(status_t2 == status_t3, "live at both times", 
                                                                               ifelse(status_t2 == 1 & status_t3 > 1, "live at t2, dead at t3", NA))) %>%
  group_by(status.change) %>% summarise(n())

##########################################
# state 23 Maine
NE_designcd %>% filter(STATECD == 23)
# designcd 101 for 1995 and 1999-- subps 101
# for 1999 designcd are 1, 115, 116, 117,--116 and 17 for suplot 1 may be remeasurements? of 101


ocon <- dbConnect(odbc(), "fiadb01p")
ME_plot <- dbGetQuery(ocon, 
                      "SELECT cn, statecd, unitcd, countycd, plot, invyr, designcd
FROM fs_fiadb.plot
WHERE statecd =  ANY(23) 
                      and invyr < 2000
") %>%
  as_tibble()%>%
  rename_with(tolower) %>% 
  rename("plt_cn" = "cn") %>%
  mutate(PLOT.ID = paste0(statecd, "_", unitcd, "_", countycd, "_", unitcd, "_", plot))
# get CT_plot
CN.plot <- ME_plot %>% mutate(PLOT.ID = paste0(statecd, "_", unitcd, "_", countycd, "_", unitcd, "_", plot))%>% group_by(PLOT.ID) %>% summarise(n())
remeas.plots <- CN.plot %>% filter(`n()` ==2)

ME_tree <- dbGetQuery(ocon, "SELECT cn, plt_cn, prev_tre_cn, statecd, unitcd, countycd, plot, subp, tree, invyr, dia, statuscd, spcd
                      FROM fs_fiadb.tree
                      WHERE statecd =  ANY(23) 
                      and invyr < 2000
                      ") %>%
  as_tibble()%>%
  rename_with(tolower)%>% left_join(., ME_plot)


# subplot 101 may be remeasured between 1985 and 1998
ME_tree %>% filter(PLOT.ID %in% unique(remeas.plots$PLOT.ID)[3]) %>% group_by(subp, invyr) %>% summarise(n())
View(ME_tree %>% filter(PLOT.ID %in% unique(remeas.plots$PLOT.ID)[1]))  #%>% summarise(n())

# it the tree ids do not match, so we would need to change tree id 
ME_tree_T2 <- ME_tree %>% filter(invyr == 1995) %>% rename("dbh_t2" = "dia",
                                                           "status_t2" = "statuscd", 
                                                           "invyr_t2" = "invyr", 
                                                           "designcd_t2" = "designcd")%>% 
  mutate(subp.match = ifelse(subp == 101, 1, NA), 
         tree.match = tree - 10000)
ME_tree_T3 <- ME_tree %>% filter(invyr == 1999)%>% rename("dbh_t3" = "dia",
                                                          "status_t3" = "statuscd", 
                                                          "invry_t3" = "invyr", 
                                                          "designcd_t3" = "designcd") %>% 
  mutate(subp.match = ifelse(subp == 1, 1, NA), 
         tree.match = tree)

ME_tree_remeas <- left_join(ME_tree_T2 %>% select(-cn, -plt_cn, -prev_tre_cn, -subp, -tree), 
                            ME_tree_T3 %>% select(-cn, -plt_cn, -prev_tre_cn, -subp, -tree))

# no matches
ME_tree_remeas %>% filter(!is.na(invry_t3)) %>%  mutate(status.change = ifelse(status_t2 == status_t3, "live at both times", 
                                                                               ifelse(status_t2 == 1 & status_t3 > 1, "live at t2, dead at t3", NA))) %>%
  group_by(status.change) %>% summarise(n())

##########################################
# state 25 Massachusetts
NE_designcd %>% filter(STATECD == 25)
# designcd 101 for 1995 and 1999-- subps 101
# for 1999 designcd are 1, 115, 116, 117,--116 and 17 for suplot 1 may be remeasurements? of 101


ocon <- dbConnect(odbc(), "fiadb01p")
MA_plot <- dbGetQuery(ocon, 
                      "SELECT cn, statecd, unitcd, countycd, plot, invyr, designcd
FROM fs_fiadb.plot
WHERE statecd =  ANY(25) 
                      and invyr < 2000
") %>%
  as_tibble()%>%
  rename_with(tolower) %>% 
  rename("plt_cn" = "cn") %>%
  mutate(PLOT.ID = paste0(statecd, "_", unitcd, "_", countycd, "_", unitcd, "_", plot))
# get CT_plot
CN.plot <- MA_plot %>% mutate(PLOT.ID = paste0(statecd, "_", unitcd, "_", countycd, "_", unitcd, "_", plot))%>% group_by(PLOT.ID) %>% summarise(n())
remeas.plots <- CN.plot %>% filter(`n()` ==2)

MA_tree <- dbGetQuery(ocon, "SELECT cn, plt_cn, prev_tre_cn, statecd, unitcd, countycd, plot, subp, tree, invyr, dia, statuscd, spcd
                      FROM fs_fiadb.tree
                      WHERE statecd =  ANY(25) 
                      and invyr < 2000
                      ") %>%
  as_tibble()%>%
  rename_with(tolower)%>% left_join(., MA_plot)


# subplot 101 may be remeasured between 1985 and 1998
MA_tree %>% filter(PLOT.ID %in% unique(remeas.plots$PLOT.ID)[3]) %>% group_by(subp, invyr) %>% summarise(n())
View(MA_tree %>% filter(PLOT.ID %in% unique(remeas.plots$PLOT.ID)[1]))  #%>% summarise(n())

# it the tree ids do not match, so we would need to change tree id 
MA_tree_T2 <- MA_tree %>% filter(invyr == 1985) %>% rename("dbh_t2" = "dia",
                                                           "status_t2" = "statuscd", 
                                                           "invyr_t2" = "invyr", 
                                                           "designcd_t2" = "designcd")

MA_tree_T3 <- MA_tree %>% filter(invyr == 1998)%>% rename("dbh_t3" = "dia",
                                                          "status_t3" = "statuscd", 
                                                          "invry_t3" = "invyr", 
                                                          "designcd_t3" = "designcd")

MA_tree_remeas <- left_join(MA_tree_T2 %>% select(-cn, -plt_cn, -prev_tre_cn), 
                            MA_tree_T3 %>% select(-cn, -plt_cn, -prev_tre_cn))

# no matches
MA_tree_remeas %>% filter(!is.na(invry_t3)) %>%  mutate(status.change = ifelse(status_t2 == status_t3, "live at both times", 
                                                                               ifelse(status_t2 == 1 & status_t3 > 1, "live at t2, dead at t3", NA))) %>%
  group_by(status.change) %>% summarise(n())


######################################################################################################
# Matching fiadb trees that are in eastwide
NE_designcd %>% filter(INVYR < 2000) %>% select(STATECD, INVYR) %>% distinct() 



##########################################
# state 33 New Hampshire
NE_designcd %>% filter(STATECD == 33)
# designcd 101 for 1983 and 1997-- subps 100 for 1983 and 101, 102, 103, 104
# for 1997 designcd 111


ocon <- dbConnect(odbc(), "fiadb01p")
NH_plot <- dbGetQuery(ocon, 
                      "SELECT cn, statecd, unitcd, countycd, plot, invyr, designcd, kindcd, plot_status_cd, plot_nonsample_reasn_cd, designcd_p2a
FROM fs_fiadb.plot
WHERE statecd =  ANY(33) 
                      and invyr < 2000
") %>%
  as_tibble()%>%
  rename_with(tolower) %>% 
  rename("plt_cn" = "cn") %>%
  mutate(PLOT.ID = paste0(statecd,  "_", countycd,  "_", plot))
# get CT_plot
CN.plot <- NH_plot %>% mutate(PLOT.ID = paste0(statecd, "_", countycd, "_", plot))%>% group_by(PLOT.ID) %>% summarise(n())
remeas.plots <- CN.plot %>% filter(`n()` ==2)

NH_tree <- dbGetQuery(ocon, "SELECT cn, plt_cn, prev_tre_cn, statecd, unitcd, countycd, plot, subp, tree, invyr, dia, statuscd, spcd, tpa_unadj, prevdia, reconcilecd, prev_status_cd, P2A_GRM_FLG
                      FROM fs_fiadb.tree
                      WHERE statecd =  ANY(33) 
                      and invyr < 2000
                      ") %>%
  as_tibble()%>%
  rename_with(tolower)%>% left_join(., NH_plot)


# subplot 101 may be remeasured between 1985 and 1998
NH_tree %>% filter(PLOT.ID %in% unique(remeas.plots$PLOT.ID)[3]) %>% group_by(subp, invyr) %>% summarise(n())
View(NH_tree %>% filter(PLOT.ID %in% unique(remeas.plots$PLOT.ID)[1]))  #%>% summarise(n())
summary(NH_tree$prevdia)

# check if reconcilecd or prev_status_cd are filled in to help explain why status == 0
summary(NH_tree$reconcilecd)
summary(NH_tree$prev_status_cd)
unique(NH_tree$designcd_p2a) # all NA
# designcd_p2a:

unique(NH_tree$p2a_grm_flg) 

# p2a_grm_flag all "N" or NA:
# FROM FIA USERS MANUAL: Periodic to annual growth, removal, and mortality flag. A code indicating if this tree is 
# part of a periodic inventory that is only included for the purposes of computing growth, 
# removals and/or mortality estimates. The flag is set to 'Y' for those trees that are needed 
# for estimation and otherwise is left blank (null)
summary(NH_tree$plot_status_cd)
NH_plot %>% group_by(plot_status_cd, plot_nonsample_reasn_cd)%>% summarise(n())

# it the tree ids do not match, so we would need to change tree id 
NH_tree_T2 <- NH_tree %>% filter(invyr == 1983) %>% rename("dbh_t2" = "dia",
                                                           "status_t2" = "statuscd", 
                                                           "invyr_t2" = "invyr", 
                                                           "designcd_t2" = "designcd", 
                                                           "tpa_unadj_t2" = "tpa_unadj", 
                                                           "prevdia_t2" = "prevdia")  %>%
  mutate(subp.match = ifelse(subp == 101, 100, NA), 
         tree.match = tree)

NH_tree_T3 <- NH_tree %>% filter(invyr == 1997)%>% rename("dbh_t3" = "dia",
                                                          "status_t3" = "statuscd", 
                                                          "invry_t3" = "invyr", 
                                                          "designcd_t3" = "designcd", 
                                                          "tpa_unadj_t3" = "tpa_unadj", 
                                                          "prevdia_t3" = "prevdia")%>%
  mutate(subp.match = ifelse(subp == 101, 100, NA), 
         tree.match = ifelse(tree < 20000, tree, 
                             ifelse(tree > 30000, tree - 20000, tree -10000)))

NH_tree_remeas <- left_join(NH_tree_T2 %>% select(-cn, -plt_cn, -prev_tre_cn, -tree, -subp), 
                            NH_tree_T3 %>% select(-cn, -plt_cn, -prev_tre_cn, -tree, -subp))

unique(NH_tree_T2$PLOT.ID) %in% "33_2_3_2_156" 

# get the eastwide data that matches
eastwide.nh.match <- TREE %>% filter(stname %in% "NH" & date %in% 1997 & !dbhold == 0) 

# need to rename the plot ids to match those in the fiadb:
eastwide.nh.match <- eastwide.nh.match %>% ungroup() %>% 
  mutate(PLOT.ID =  paste0(state, "_", county, "_", pltnum))

# unique plots with matching eastwide data that we can search the fiadb for:
unique(eastwide.nh.match$PLOT.ID) %in% NH_tree_T3$PLOT.ID
unique(eastwide.nh.match$PLOT.ID) %in% NH_tree_T2$PLOT.ID

# plot i == 6 was all cut
# plot i == 7 has some cut trees and some live
unique(eastwide.nh.match$PLOT.ID) %in% Plots.with.mulit.match$PLOT.ID[1]


NH_tree_T2$dbh_t2 <- as.numeric(NH_tree_T2$dbh_t2)
NH_tree_T2$spcd <- as.numeric(NH_tree_T2$spcd)
NH_tree_T3$dbh_t3 <- as.numeric(NH_tree_T3$dbh_t3)
NH_tree_T3$spcd <- as.numeric(NH_tree_T3$spcd)

unique(eastwide.nh.match$PLOT.ID) %in% unique(dead.at.both.nh.potential.mismatch$PLOT.ID)

eastwide.nh.fiadb.match <- list()
for( i in 1:length(unique(eastwide.nh.match$PLOT.ID))){
  plot.i.T2 <- NH_tree_T2 %>% filter(PLOT.ID %in% unique(eastwide.nh.match$PLOT.ID)[i])
  plot.i.T3 <- NH_tree_T3 %>% filter(PLOT.ID %in% unique(eastwide.nh.match$PLOT.ID)[i])
  nrow(plot.i.T2)
  nrow(plot.i.T3)
  
  eastwide.plot.match <- eastwide.nh.match %>% filter(PLOT.ID %in% unique(eastwide.nh.match$PLOT.ID)[i])
  #View(eastwide.plot.match)
  
  # create columns to match up the current and old tree cns
  eastwide.plot.match$old_tree_cn <- NA
  eastwide.plot.match$cur_tree_cn <- NA
  
  # create columns to match up the current and old tree status
  eastwide.plot.match$old_tree_fiadb_status <- NA
  eastwide.plot.match$cur_tree_fiadb_status <- NA
  
  # create columns to match up the current and old tpa_unadj
  eastwide.plot.match$cur_tpa_unadj <- NA
  eastwide.plot.match$old_tpa_unadj <- NA
  
  # create columns to that describe how the match was done
  eastwide.plot.match$cur_match_note <- NA
  eastwide.plot.match$old_match_note <- NA
  
  # create columns to that have the current and old subp
  eastwide.plot.match$cur_subp <- NA
  eastwide.plot.match$old_subp <- NA
  
  #eastwide.plot.match$tree <- as.character(eastwide.plot.match$tree)

  for(j in 1:nrow(eastwide.plot.match)){
    potential.tree.ids <- as.character(c(10000, 20000, 30000, 40000, 50000, 60000) + eastwide.plot.match[j,]$tree)
    
    i.match.t3 <- plot.i.T3 %>% filter(dbh_t3 %in% as.character(eastwide.plot.match[j,]$dbhcur) & spcd %in% as.character(eastwide.plot.match[j,]$SPCD))
    i.match.t2 <-  plot.i.T2 %>% filter(dbh_t2 %in% as.character(eastwide.plot.match[j,]$dbhold) & spcd %in%  as.character(eastwide.plot.match[j,]$SPCD))
    
    if(length(i.match.t2$cn) == 0){
      if(length(plot.i.T2$cn)==0){
      eastwide.plot.match[j,]$old_match_note <- "no matching PLOT.ID"
      }else{
      eastwide.plot.match[j,]$old_match_note <- "no matching dbh and species"
    }
    }
    if(length(i.match.t3$cn) == 0){
      if(length(plot.i.T3$cn)==0){
        eastwide.plot.match[j,]$cur_match_note <- "no matching PLOT.ID"
      }else{
      eastwide.plot.match[j,]$cur_match_note <- "no matching dbh and species"
    }
    }
    # if there is a match for the old dataset, add it to the eastwid.plot.match
    if(length(i.match.t2$cn) == 1){
      eastwide.plot.match[j,]$old_tree_cn <- i.match.t2$cn
      eastwide.plot.match[j,]$old_tree_fiadb_status <- i.match.t2$status_t2
      eastwide.plot.match[j,]$old_tpa_unadj <- i.match.t2$tpa_unadj_t2
      eastwide.plot.match[j,]$old_subp <- i.match.t2$subp
      eastwide.plot.match[j,]$old_match_note <- "simple match, single tree"
    }
    # if there is a match in the current dataset, add it to the eastwide.plot.match
    if(length(i.match.t3$cn) == 1){
      eastwide.plot.match[j,]$cur_tree_cn <- i.match.t3$cn
      eastwide.plot.match[j,]$cur_tree_fiadb_status <- i.match.t3$status_t3
      eastwide.plot.match[j,]$cur_tpa_unadj <- i.match.t3$tpa_unadj_t3
      eastwide.plot.match[j,]$cur_subp <- i.match.t3$subp
      eastwide.plot.match[j,]$cur_match_note <- "simple match, single tree"
    }
    
    
    # if both time points have mulitple matching trees
    if(length(i.match.t3$cn) > 1 & length(i.match.t2$cn) > 1){
      
      
      # if the later inventory has mulitple matches, see if we can use tree.id or status to match up the tree 
      i.tree.match.t3 <- i.match.t3 %>% filter(tree.match %in% potential.tree.ids &
                                                 status_t3 == eastwide.plot.match[j,]$status) 
      
      # special case for snags and salvgable dead
      if(eastwide.plot.match[j,]$status == 5 | eastwide.plot.match[j,]$status == 4){
        i.tree.match.t3 <- i.match.t3 %>% filter(tree.match %in% potential.tree.ids & status_t3 == 2) 
      }
      
      
      # if the current FIADB status code is zero for the only matching tree id, we could match on tree ids alone?
     
      if(length(i.tree.match.t3$cn) == 0 ){ # if there are still no matches:
        tree.id.only.match <- i.match.t3 %>% filter(tree.match %in% potential.tree.ids)
        if(length(tree.id.only.match$status_t3) == 1){# and the fiadb says status cd == 0, then just match on potential ids
        i.tree.match.t3 <- i.match.t3 %>% filter(tree.match %in% potential.tree.ids)
      }
      }
      
      if(nrow(i.tree.match.t3)==1){
        eastwide.plot.match[j,]$cur_tree_cn <- i.tree.match.t3$cn
        eastwide.plot.match[j,]$cur_tree_fiadb_status <- i.tree.match.t3$status_t3
        eastwide.plot.match[j,]$cur_tpa_unadj <- i.tree.match.t3$tpa_unadj_t3
        eastwide.plot.match[j,]$cur_match_note <- "matched by potential tree id & status"
       
         if(eastwide.plot.match[j,]$status == 5 | eastwide.plot.match[j,]$status == 4){
          eastwide.plot.match[j,]$cur_match_note <- "matched by potential tree id & status; needed to change snag/salvegable status to 2"
         }
        
       
        # if the tree number matches up with the t3 tree, then select that one
        i.tree.match.t2 <- i.match.t2 %>% filter(tree.match %in% potential.tree.ids)
        
        if(nrow(i.tree.match.t2)==1){
          i.tree.match.t2 <- i.tree.match.t2
          eastwide.plot.match[j,]$old_match_note <- "matched by potential tree id"
          
        }else{ # if there are still more tree matches and all are 0 or live, just take the first one
          
          if(2 %in% unique(i.match.t2$status_t2)){
            only.live.t2 <- i.match.t2 %>% filter(!status_t2 == 2)
            if(length(only.live.t2$cn) == 1){
              i.tree.match.t2 <- only.live.t2
              eastwide.plot.match[j,]$old_match_note <- "treeids don't match; selected only non-dead matching tree"
            }else{
              i.tree.match.t2 <- only.live.t2[1,]
              eastwide.plot.match[j,]$old_match_note <- "treeids don't match; selected first non-dead matching tree"
            }
          }else{
          i.tree.match.t2 <- i.match.t2[1,]
          eastwide.plot.match[j,]$old_match_note <- "treeids don't match; selected first matching tree"
          }  
        }
        eastwide.plot.match[j,]$old_tree_cn <- i.tree.match.t2$cn
        eastwide.plot.match[j,]$old_tree_fiadb_status <- i.tree.match.t2$status_t2
        eastwide.plot.match[j,]$old_tpa_unadj <- i.tree.match.t2$tpa_unadj_t2
        
      }else{
        i.tree.subp.match.t3 <- i.tree.match.t3 %>% filter(subp == 101)
        
        if(length(i.tree.subp.match.t3$cn) == 1){
          eastwide.plot.match[j,]$cur_tree_cn <-  i.tree.subp.match.t3$cn
          eastwide.plot.match[j,]$cur_tree_fiadb_status <-  i.tree.subp.match.t3$status_t3
          eastwide.plot.match[j,]$cur_tpa_unadj <-  i.tree.subp.match.t3$tpa_unadj_t3
          eastwide.plot.match[j,]$cur_subp <- i.tree.subp.match.t3$subp
          eastwide.plot.match[j,]$cur_match_note <- "matched by potential tree id, status, and subp"
          
          

        }
        else{
        eastwide.plot.match[j,]$cur_tree_cn <- "multiple matches" #i.match.t3$cn
      }
      
      }
      
      # once we have a single current tree, see if the time 2 tre matches up
      i.tree.match.t2 <- i.match.t2 %>% filter(tree.match %in% potential.tree.ids)
      
      if(nrow(i.tree.match.t2)==1){
        i.tree.match.t2 <- i.tree.match.t2
        eastwide.plot.match[j,]$old_match_note <- "matched by potential tree id"
        
      }else{ # if there are still more tree matches, just take the first one
        i.tree.match.t2 <- i.match.t2[1,]
        eastwide.plot.match[j,]$old_match_note <- "treeids don't match; selected first matching tree"
        
      }
      eastwide.plot.match[j,]$old_tree_cn <- i.tree.match.t2$cn
      eastwide.plot.match[j,]$old_tree_fiadb_status <- i.tree.match.t2$status_t2
      eastwide.plot.match[j,]$old_tpa_unadj <- i.tree.match.t2$tpa_unadj_t2
      eastwide.plot.match[j,]$old_subp <- i.tree.match.t2$subp
      
      
      
    }
    
    # if there is only multiple matching trees in the t3 
    if(length(i.match.t3$cn) > 1  & length(i.match.t2$cn) <= 1){
      # if the later inventory has mulitple matches, see if we can use tree.id or status to match up the tree 
      
      i.tree.match.t3 <- i.match.t3 %>% filter(tree.match %in% potential.tree.ids &
                                                 status_t3 == eastwide.plot.match[j,]$status) 
      # create a special case for snags:
      if(eastwide.plot.match[j,]$status == 5 | eastwide.plot.match[j,]$status == 4){
        
        i.tree.match.t3 <- i.match.t3 %>% filter(tree.match %in% potential.tree.ids & status_t3 == 2)
        
      }
      
      # if the current FIADB status code is zero for the only matching tree id, we could match on tree ids alone?
      
      if(length(i.tree.match.t3$cn) == 0 ){ # if there are still no matches:
        tree.id.only.match <- i.match.t3 %>% filter(tree.match %in% potential.tree.ids)
        if(length(tree.id.only.match$status_t3) == 1){# and the fiadb says status cd == 0, then just match on potential ids
          i.tree.match.t3 <- i.match.t3 %>% filter(tree.match %in% potential.tree.ids)
        }else{
          # if there are still multiple matches, find the tree that matches statuscd (if it is marked dead) and subplot
          if(length(tree.id.only.match)> 1 & eastwide.plot.match[j,]$status > 1){
            i.tree.subp.match.t3 <- i.match.t3 %>% filter(subp %in% 101 & 
                                  status_t3 == 2)
            if(i.tree.subp.match.t3$cn == 0){ #if no trees match now, just select by statuscd
            i.tree.subp.match.t3 <- i.match.t3 %>% filter(status_t3 == 2)
            }
          }
          i.tree.match.t3 <- i.tree.subp.match.t3
        }
      }
      
      if(nrow(i.tree.match.t3)==1){
        eastwide.plot.match[j,]$cur_tree_cn <- i.tree.match.t3$cn
        eastwide.plot.match[j,]$cur_tree_fiadb_status <- i.tree.match.t3$status_t3
        eastwide.plot.match[j,]$cur_tpa_unadj <- i.tree.match.t3$tpa_unadj_t3
        eastwide.plot.match[j,]$cur_tpa_unadj <- i.tree.match.t3$subp
        eastwide.plot.match[j,]$cur_match_note <- "matched by potential tree id & status"
        
        if(eastwide.plot.match[j,]$status == 5 | eastwide.plot.match[j,]$status == 4 | eastwide.plot.match[j,]$status == 3){
          eastwide.plot.match[j,]$cur_match_note <- "matched by potential tree id & status; needed to change snag/salvegable status to 2"
        }
        # if(length(tree.id.only.match$status_t3) == 1 & tree.id.only.match$status_t3 == 0){
        #   eastwide.plot.match[j,]$cur_match_note <- paste0("matched by potential tree id only, FIADB statuscd == 0 but EWDB statcd ==", eastwide.plot.match[j,]$status)
        # }
        
      }else{
        # if there are still multiple matches, see if you can select the tree in suplot 101
        i.tree.subp.match.t3 <- i.tree.match.t3 %>% filter(subp == 101)
        if(length(i.tree.subp.match.t3$cn) == 1){
          eastwide.plot.match[j,]$cur_tree_cn <-  i.tree.subp.match.t3$cn
          eastwide.plot.match[j,]$cur_tree_fiadb_status <-  i.tree.subp.match.t3$status_t3
          eastwide.plot.match[j,]$cur_tpa_unadj <-  i.tree.subp.match.t3$tpa_unadj_t3
          eastwide.plot.match[j,]$cur_tpa_unadj <- i.tree.subp.match.t3$subp
          eastwide.plot.match[j,]$cur_match_note <- "matched by potential tree id, status, and subp"
        }
        else{
        eastwide.plot.match[j,]$cur_tree_cn <- "multiple matches"#i.match.t3$cn
        }
      }
      
    }
    # if there is more than one match for the old dataset, add a -999
    if(length(i.match.t2$cn) > 1 & length(i.match.t3$cn) <=1){
      # if the tree number matches up with the t3 tree, then select that one
      i.tree.match.t2 <- i.match.t2 %>% filter(tree.match %in% potential.tree.ids)
      
      if(nrow(i.tree.match.t2)==1){
        i.tree.match.t2 <- i.tree.match.t2
        eastwide.plot.match[j,]$old_match_note <- "matched by potential tree id"
      }else{ # if there are still more tree matches, just take the first one
        i.tree.match.t2 <- i.match.t2[1,]
        eastwide.plot.match[j,]$old_match_note <- "treeids don't match; selected first matching tree"
        
      }
      eastwide.plot.match[j,]$old_tree_cn <- i.tree.match.t2$cn
      eastwide.plot.match[j,]$old_tree_fiadb_status <- i.tree.match.t2$status_t2
      eastwide.plot.match[j,]$old_tpa_unadj <- i.tree.match.t2$tpa_unadj_t2
      eastwide.plot.match[j,]$old_tpa_unadj <- i.tree.match.t2$subp
      
      
    }
    

    
    
  }
  
  
  eastwide.nh.fiadb.match[[i]] <- eastwide.plot.match
  rm(eastwide.plot.match)
}
eastwide.nh.fiadb.updated.match.df <- do.call(rbind, eastwide.nh.fiadb.match)

unique(eastwide.nh.fiadb.updated.match.df$old_match_note)
# how many trees have matches now?
eastwide.nh.fiadb.updated.match.df %>% group_by(!is.na(cur_tree_cn) & !is.na(old_tree_cn) & !old_tree_cn %in% "multiple matches" & !cur_tree_cn %in% "multiple matches") %>% summarise(n())

# only 3 trees are unmatched for both time series
# 15 trees only have matches in the past
# 2215 trees are matched in the current inventory but not the past
eastwide.nh.fiadb.updated.match.df %>% group_by(cur_match_note %in% "no matching dbh and species", old_match_note %in% c("no matching dbh and species", "no matching PLOT.ID")) %>% summarise(n())

# all matches have a note now:
eastwide.nh.fiadb.updated.match.df %>% filter(cur_match_note %in% "matched by potential tree id, status, and subp" &
                                                is.na(old_match_note)) %>% select(PLOT.ID)

#lets see how most of the trees matched up:
eastwide.nh.fiadb.updated.match.df %>% group_by(!is.na(cur_tree_cn) & !is.na(old_tree_cn) & 
                                                  !old_tree_cn %in% "multiple matches" & 
                                                  !cur_tree_cn %in% "multiple matches", 
                                                old_match_note, cur_match_note
            ) %>% summarise(n())%>% ungroup()|>gt()

# lets see if all all the status code changes align:
unique(eastwide.nh.fiadb.updated.match.df$cur_tree_fiadb_status)
unique(eastwide.nh.fiadb.updated.match.df$old_tree_fiadb_status)


eastwide.nh.fiadb.updated.match.df %>% filter(!is.na(cur_tree_fiadb_status))%>% mutate(cur_status_consistent = ifelse(status == 1 & cur_tree_fiadb_status == 1, "yes, live", 
                                                                             ifelse(status > 1 & cur_tree_fiadb_status > 1, "yes, dead", 
                                                                                    ifelse(status > 1 & cur_tree_fiadb_status == 0, "no, EWDB = dead, current status = 0", 
                                                                                           ifelse(status == 1 & cur_tree_fiadb_status == 0, "no, EWDB = live, current status = 0", 
                                                                                                  ifelse(is.na(cur_tree_fiadb_status), "no matching tree in current FIADB", "no matching tree in current FIADB")))))) %>% 
  group_by(cur_status_consistent) %>% #filter(is.na(cur_status_consistent))%>% 
  #select(tree, status, cur_tree_fiadb_status)
summarise(`# of trees` = n())
# all the NA values don't have matching trees in the current FIADB, so I filtered these out

NA.consistent.status<- eastwide.nh.fiadb.updated.match.df %>% mutate(cur_status_consistent = ifelse(status == 1 & cur_tree_fiadb_status == 1, "yes, live", 
                                                                             ifelse(status > 1 & cur_tree_fiadb_status > 1, "yes, dead", 
                                                                                    ifelse(status > 1 & cur_tree_fiadb_status == 0, "no, EWDB = dead, current status = 0", 
                                                                                           ifelse(status == 1 & cur_tree_fiadb_status == 0, "no, EWDB = live, current status = 0", 
                                                                                                  ifelse(is.na(cur_tree_fiadb_status), "no matching tree in current FIADB", NA)))))) %>% 
  group_by(cur_status_consistent) %>% filter(is.na(cur_status_consistent))
# yes then mostly align, with the exception of current status == 0
# verify that these are status code changes:
eastwide.nh.fiadb.updated.match.df %>% mutate(status_change = ifelse(old_tree_fiadb_status == 1 & cur_tree_fiadb_status == 1, "live at both times", 
                                                                             ifelse(old_tree_fiadb_status > 1 & cur_tree_fiadb_status > 1, "dead at both times", 
                                                                                    ifelse(old_tree_fiadb_status == 1 & cur_tree_fiadb_status > 1, "live at t2, dead at t3", 
                                                                                    ifelse(old_tree_fiadb_status > 1 & cur_tree_fiadb_status == 0, "dead at t2, zero at t3", 
                                                                                           ifelse(old_tree_fiadb_status == 1 & cur_tree_fiadb_status == 0, "live at t2, status at t3 = 0", 
                                                                                                  ifelse(old_tree_fiadb_status == 0 & cur_tree_fiadb_status > 1, "status at t2 = 0, live at t3", 
                                                                                                         ifelse(old_tree_fiadb_status > 1 & cur_tree_fiadb_status == 1, "status at t2 = dead, but status at t3 == live",
                                                                                                                ifelse(is.na(old_tree_fiadb_status) | is.na(cur_tree_fiadb_status),"NA values for at least one status","no matching trees in FIADB"))))))))) %>% 
  group_by(status_change) %>% #filter(status_change %in% "no matching trees in FIADB") %>% select(old_tree_fiadb_status, cur_tree_fiadb_status, old_match_note , cur_match_note)#%>% 
  mutate(status_change = ifelse(is.na(status_change) == TRUE, "could not match tree to current or old FIADB", status_change))%>%
  #select(tree, status, cur_tree_fiadb_status)
  summarise(`# of trees` = n()) |> gt()|>
  gtsave("images/filtering_exploration/NH_fiadb_to_eastwide_matching_status_cd_change.png")

eastwide.nh.fiadb.updated.match.df %>%  mutate(status_change = ifelse(old_tree_fiadb_status == 1 & cur_tree_fiadb_status == 1, "live at both times", 
                                                                      ifelse(old_tree_fiadb_status > 1 & cur_tree_fiadb_status > 1, "dead at both times", 
                                                                             ifelse(old_tree_fiadb_status == 1 & cur_tree_fiadb_status > 1, "live at t2, dead at t3", 
                                                                                    ifelse(old_tree_fiadb_status > 1 & cur_tree_fiadb_status == 0, "dead at t2, zero at t3", 
                                                                                           ifelse(old_tree_fiadb_status == 1 & cur_tree_fiadb_status == 0, "live at t2, status at t3 = 0", 
                                                                                                  ifelse(old_tree_fiadb_status == 0 & cur_tree_fiadb_status > 1, "status at t2 = 0, live at t3", 
                                                                                                         ifelse(old_tree_fiadb_status > 1 & cur_tree_fiadb_status == 1, "status at t2 = dead, but status at t3 == live",
                                                                                                                ifelse(is.na(old_tree_fiadb_status) | is.na(cur_tree_fiadb_status),"NA values for at least one status","no matching trees in FIADB"))))))))) %>% 
  group_by(status_change, status) %>% #filter(status_change %in% "no matching trees in FIADB") %>% select(old_tree_fiadb_status, cur_tree_fiadb_status, old_match_note , cur_match_note)#%>% 
  mutate(status_change = ifelse(is.na(status_change) == TRUE & is.na(cur_tree_fiadb_status) == TRUE, "no matching tree in current FIADB", 
                                ifelse(is.na(status_change) == TRUE & is.na(old_tree_fiadb_status) == TRUE, "no matching tree in old FIADB", status_change)))%>%
  #select(tree, status, cur_tree_fiadb_status)
  summarise(`# of trees` = n()) %>% ungroup()|> gt()|>
  gtsave("images/filtering_exploration/NH_fiadb_to_eastwide_matching_status_cd_change_by_EW_status.png")

# is it possible that the dead at both times categories are a mismatch in t2?

EW.nh.fiadb.matches.status <- eastwide.nh.fiadb.updated.match.df %>%  mutate(status_change = ifelse(old_tree_fiadb_status == 1 & cur_tree_fiadb_status == 1, "live at both times", 
                                                                      ifelse(old_tree_fiadb_status > 1 & cur_tree_fiadb_status > 1, "dead at both times", 
                                                                             ifelse(old_tree_fiadb_status == 1 & cur_tree_fiadb_status > 1, "live at t2, dead at t3", 
                                                                                    ifelse(old_tree_fiadb_status > 1 & cur_tree_fiadb_status == 0, "dead at t2, zero at t3", 
                                                                                           ifelse(old_tree_fiadb_status == 1 & cur_tree_fiadb_status == 0, "live at t2, status at t3 = 0", 
                                                                                                  ifelse(old_tree_fiadb_status == 0 & cur_tree_fiadb_status > 1, "status at t2 = 0, dead at t3", 
                                                                                                         ifelse(old_tree_fiadb_status > 1 & cur_tree_fiadb_status == 1, "status at t2 = dead, but status at t3 == live",
                                                                                                                ifelse(is.na(old_tree_fiadb_status) | is.na(cur_tree_fiadb_status),"NA values for at least one status","no matching trees in FIADB"))))))))) %>% 
  group_by(status_change, status) %>% #filter(status_change %in% "no matching trees in FIADB") %>% select(old_tree_fiadb_status, cur_tree_fiadb_status, old_match_note , cur_match_note)#%>% 
  mutate(status_change = ifelse(is.na(status_change) == TRUE & is.na(cur_tree_fiadb_status) == TRUE, "no matching tree in current FIADB", 
                                ifelse(is.na(status_change) == TRUE & is.na(old_tree_fiadb_status) == TRUE, "no matching tree in old FIADB", status_change)))


# explore some of the different/non matching status changes
unique(EW.nh.fiadb.matches.status$status_change)
# 3 dead to live matches:
# two are simple matches for old fia with just one tree (possibly a tree misclassified as dead, or a mismatch of a tree)
View(EW.nh.fiadb.matches.status %>% filter(status_change %in% "status at t2 = dead, but status at t3 == live")%>% 
  select(old_subp, cur_subp, tree, status, cur_tree_fiadb_status, old_tree_fiadb_status, old_match_note, cur_match_note))

# a lot of the live to dead matches are cut trees
EW.nh.fiadb.matches.status %>% filter(status_change %in% "live at t2, dead at t3" )%>% 
       select(old_subp, cur_subp, tree, status, cur_tree_fiadb_status, old_tree_fiadb_status, old_match_note, cur_match_note) %>% 
  group_by(status, cur_tree_fiadb_status) %>%
  summarise(n())

# most of the status == 0, but dead in FIADB are status == 5 in EWDB
EW.nh.fiadb.matches.status %>% filter(status_change %in% "status at t2 = 0, dead at t3" )%>% 
  select(old_subp, cur_subp, tree, status, cur_tree_fiadb_status, old_tree_fiadb_status, old_match_note, cur_match_note) %>% 
  group_by(status, cur_tree_fiadb_status) %>%
  summarise(n())

# most of the status == 0, but dead in FIADB are status == 5 in EWDB
EW.nh.fiadb.matches.status %>% filter(status_change %in% "dead at both times" )%>% 
  select(old_subp, cur_subp, tree, status, cur_tree_fiadb_status, old_tree_fiadb_status, old_match_note, cur_match_note) %>% 
  group_by(old_match_note, cur_match_note) %>%
  summarise(n())

dead.at.both.nh <- EW.nh.fiadb.matches.status %>% filter(status_change %in% "dead at both times" ) 
unique(dead.at.both.nh$PLOT.ID)
dead.at.both.nh.potential.mismatch <- dead.at.both.nh %>% filter(old_match_note %in% c("matched by potential tree id", 
                                                 "treeids don't match; selected first matching tree"))

unique(dead.at.both.nh.potential.mismatch$PLOT.ID)
# summarise how all these trees were matched up:
eastwide.nh.fiadb.updated.match.df %>% 
  group_by(cur_match_note) %>% #filter(status_change %in% "dead at both times")%>% 
  #select(tree, status, cur_tree_fiadb_status)
  summarise(`# of trees` = n()) %>% ungroup()|> gt()|>
  gtsave("images/filtering_exploration/NH_fiadb_to_eastwide_matching_notes_T3.png")
# 

eastwide.nh.fiadb.updated.match.df %>% 
  group_by(old_match_note) %>% #filter(status_change %in% "dead at both times")%>% 
  #select(tree, status, cur_tree_fiadb_status)
  summarise(`# of trees` = n()) %>% ungroup()|> gt()|>
  gtsave("images/filtering_exploration/NH_fiadb_to_eastwide_matching_notes_T3.png")

hist(eastwide.nh.fiadb.updated.match.df$cur_tpa_unadj)
hist(eastwide.nh.fiadb.updated.match.df$old_tpa_unadj)
hist(eastwide.nh.fiadb.updated.match.df$volfac)

ggplot(eastwide.nh.fiadb.updated.match.df, aes(x = dbhcur, y = cur_tpa_unadj))+geom_point()
ggplot(eastwide.nh.fiadb.updated.match.df, aes(x = dbhcur, y = volfac))+geom_point()
ggplot(eastwide.nh.fiadb.updated.match.df, aes(x = dbhold, y = old_tpa_unadj))+geom_point()

# 
#####################################e###############################################
# do VT matches in FIADB match the trees in CT from eastewide?
# state 33 New Hampshire
NE_designcd %>% filter(STATECD == 9)
TREE %>% select(stname, date) %>% distinct()
# designcd 101 for 1983 and 1997-- subps 100 for 1983 and 101, 102, 103, 104
# for 1997 designcd 111


ocon <- dbConnect(odbc(), "fiadb01p")
VT_plot <- dbGetQuery(ocon, 
                      "SELECT cn, statecd, unitcd, countycd, plot, invyr, designcd, kindcd, plot_status_cd, plot_nonsample_reasn_cd, designcd_p2a
FROM fs_fiadb.plot
WHERE statecd =  ANY(50) 
                      and invyr < 2000
") %>%
  as_tibble()%>%
  rename_with(tolower) %>% 
  rename("plt_cn" = "cn") %>%
  mutate(PLOT.ID = paste0(statecd,  "_", countycd,  "_", plot))
# get CT_plot
CN.plot <- VT_plot %>% mutate(PLOT.ID = paste0(statecd, "_", countycd, "_", plot))%>% group_by(PLOT.ID) %>% summarise(n())
remeas.plots <- CN.plot %>% filter(`n()` ==2)

VT_tree <- dbGetQuery(ocon, "SELECT cn, plt_cn, prev_tre_cn, statecd, unitcd, countycd, plot, subp, tree, invyr, dia, statuscd, spcd, tpa_unadj, prevdia, reconcilecd, prev_status_cd, P2A_GRM_FLG
                      FROM fs_fiadb.tree
                      WHERE statecd =  ANY(50) 
                      and invyr < 2000
                      ") %>%
  as_tibble()%>%
  rename_with(tolower)%>% left_join(., VT_plot)

hist(VT_tree$tpa_unadj)

# check for any reconcilecd to indicate why we have status == 0
# also check for a prev_status_cd
summary(VT_tree$reconcilecd)
summary(VT_tree$prev_status_cd)
# all are NA values 

unique(VT_tree$p2a_grm_flg) # all "N" or NA

VT_tree %>% group_by(tpa_unadj == 0 , designcd, invyr, subp) %>%summarise(n()) %>% ungroup() |> gt()|> gtsave("images/filtering_exploration/VT_designcd_tpa_zeros.png")
# subplot 101 may be remeasured between 1985 and 1998
VT_tree %>% filter(PLOT.ID %in% unique(remeas.plots$PLOT.ID)[3]) %>% group_by(subp, invyr) %>% summarise(n())
View(VT_tree %>% filter(PLOT.ID %in% unique(remeas.plots$PLOT.ID)[1]))  #%>% summarise(n())

# it the tree ids do not match, so we would need to change tree id 
VT_tree_T2 <- VT_tree %>% filter(invyr == 1983) %>% rename("dbh_t2" = "dia",
                                                           "status_t2" = "statuscd", 
                                                           "invyr_t2" = "invyr", 
                                                           "designcd_t2" = "designcd", 
                                                           "tpa_unadj_t2" = "tpa_unadj", 
                                                           "prevdia_t2" = "prevdia")  %>%
  mutate(subp.match = ifelse(subp == 101, 100, subp), 
         tree.match = tree)

VT_tree_T3 <- VT_tree %>% filter(invyr == 1997)%>% rename("dbh_t3" = "dia",
                                                          "status_t3" = "statuscd", 
                                                          "invry_t3" = "invyr", 
                                                          "designcd_t3" = "designcd", 
                                                          "tpa_unadj_t3" = "tpa_unadj", 
                                                          "prevdia_t3" = "prevdia")%>%
  mutate(subp.match = ifelse(subp == 101, 100, NA), 
         tree.match = ifelse(tree < 20000, tree, 
                             ifelse(tree > 30000, tree - 20000, tree -10000)))

VT_tree_remeas <- left_join(VT_tree_T2 %>% select(-cn, -plt_cn, -prev_tre_cn, -tree, -subp), 
                            VT_tree_T3 %>% select(-cn, -plt_cn, -prev_tre_cn, -tree, -subp))

unique(VT_tree_T2$PLOT.ID) %in% "33_2_3_2_156" 

# get the eastwide data that matches
eastwide.VT.match <- TREE %>% filter(stname %in% "VT" & date %in% 1997 & !dbhold == 0) 

# need to rename the plot ids to match those in the fiadb:
eastwide.VT.match <- eastwide.VT.match %>% ungroup() %>% 
  mutate(PLOT.ID =  paste0(state, "_", county, "_", pltnum))

# unique plots with matching eastwide data that we can search the fiadb for:
unique(eastwide.VT.match$PLOT.ID) %in% VT_tree_T3$PLOT.ID
unique(eastwide.VT.match$PLOT.ID) %in% VT_tree_T2$PLOT.ID

# plot i == 6 was all cut
# plot i == 7 has some cut trees and some live

eastwide.VT.fiadb.match <- list()
# plot i == 6 was all cut
# plot i == 7 has some cut trees and some live
unique(eastwide.nh.match$PLOT.ID) %in% Plots.with.mulit.match$PLOT.ID[1]


VT_tree_T2$dbh_t2 <- as.numeric(VT_tree_T2$dbh_t2)
VT_tree_T2$spcd <- as.numeric(VT_tree_T2$spcd)
VT_tree_T3$dbh_t3 <- as.numeric(VT_tree_T3$dbh_t3)
VT_tree_T3$spcd <- as.numeric(VT_tree_T3$spcd)

unique(eastwide.VT.match$PLOT.ID) %in% "50_19_144" 

eastwide.VT.fiadb.match <- list()
for( i in 1:length(unique(eastwide.VT.match$PLOT.ID))){
  plot.i.T2 <- VT_tree_T2 %>% filter(PLOT.ID %in% unique(eastwide.VT.match$PLOT.ID)[i])
  plot.i.T3 <- VT_tree_T3 %>% filter(PLOT.ID %in% unique(eastwide.VT.match$PLOT.ID)[i])
  nrow(plot.i.T2)
  nrow(plot.i.T3)
  
  eastwide.plot.match <- eastwide.VT.match %>% filter(PLOT.ID %in% unique(eastwide.VT.match$PLOT.ID)[i])
  #View(eastwide.plot.match)
  
  # create columns to match up the current and old tree cns
  eastwide.plot.match$old_tree_cn <- NA
  eastwide.plot.match$cur_tree_cn <- NA
  
  # create columns to match up the current and old tree status
  eastwide.plot.match$old_tree_fiadb_status <- NA
  eastwide.plot.match$cur_tree_fiadb_status <- NA
  
  # create columns to match up the current and old tpa_unadj
  eastwide.plot.match$cur_tpa_unadj <- NA
  eastwide.plot.match$old_tpa_unadj <- NA
  
  # create columns to that describe how the match was done
  eastwide.plot.match$cur_match_note <- NA
  eastwide.plot.match$old_match_note <- NA
  
  # create columns to that have the current and old subp
  eastwide.plot.match$cur_subp <- NA
  eastwide.plot.match$old_subp <- NA
  
  #eastwide.plot.match$tree <- as.character(eastwide.plot.match$tree)
  
  for(j in 1:nrow(eastwide.plot.match)){
    potential.tree.ids <- as.character(c(10000, 20000, 30000, 40000, 50000, 60000) + eastwide.plot.match[j,]$tree)
    
    i.match.t3 <- plot.i.T3 %>% filter(dbh_t3 %in% as.character(eastwide.plot.match[j,]$dbhcur) & spcd %in% as.character(eastwide.plot.match[j,]$SPCD))
    i.match.t2 <-  plot.i.T2 %>% filter(dbh_t2 %in% as.character(eastwide.plot.match[j,]$dbhold) & spcd %in%  as.character(eastwide.plot.match[j,]$SPCD))
    
    if(length(i.match.t2$cn) == 0){
      if(length(plot.i.T2$cn)==0){
        eastwide.plot.match[j,]$old_match_note <- "no matching PLOT.ID"
      }else{
        eastwide.plot.match[j,]$old_match_note <- "no matching dbh and species"
      }
    }
    if(length(i.match.t3$cn) == 0){
      if(length(plot.i.T3$cn)==0){
        eastwide.plot.match[j,]$cur_match_note <- "no matching PLOT.ID"
      }else{
        eastwide.plot.match[j,]$cur_match_note <- "no matching dbh and species"
      }
    }
    # if there is a match for the old dataset, add it to the eastwid.plot.match
    if(length(i.match.t2$cn) == 1){
      eastwide.plot.match[j,]$old_tree_cn <- i.match.t2$cn
      eastwide.plot.match[j,]$old_tree_fiadb_status <- i.match.t2$status_t2
      eastwide.plot.match[j,]$old_tpa_unadj <- i.match.t2$tpa_unadj_t2
      eastwide.plot.match[j,]$old_subp <- i.match.t2$subp
      eastwide.plot.match[j,]$old_match_note <- "simple match, single tree"
    }
    # if there is a match in the current dataset, add it to the eastwide.plot.match
    if(length(i.match.t3$cn) == 1){
      eastwide.plot.match[j,]$cur_tree_cn <- i.match.t3$cn
      eastwide.plot.match[j,]$cur_tree_fiadb_status <- i.match.t3$status_t3
      eastwide.plot.match[j,]$cur_tpa_unadj <- i.match.t3$tpa_unadj_t3
      eastwide.plot.match[j,]$cur_subp <- i.match.t3$subp
      eastwide.plot.match[j,]$cur_match_note <- "simple match, single tree"
    }
    
    
    # if both time points have mulitple matching trees
    if(length(i.match.t3$cn) > 1 & length(i.match.t2$cn) > 1){
      
      
      # if the later inventory has mulitple matches, see if we can use tree.id or status to match up the tree 
      i.tree.match.t3 <- i.match.t3 %>% filter(tree.match %in% potential.tree.ids &
                                                 status_t3 == eastwide.plot.match[j,]$status) 
      
      # special case for snags and salvgable dead
      if(eastwide.plot.match[j,]$status == 5 | eastwide.plot.match[j,]$status == 4){
        i.tree.match.t3 <- i.match.t3 %>% filter(tree.match %in% potential.tree.ids & status_t3 == 2) 
      }
      
      
      # if the current FIADB status code is zero for the only matching tree id, we could match on tree ids alone?
      
      if(length(i.tree.match.t3$cn) == 0 ){ # if there are still no matches:
        tree.id.only.match <- i.match.t3 %>% filter(tree.match %in% potential.tree.ids)
        if(length(tree.id.only.match$status_t3) == 1){# and the fiadb says status cd == 0, then just match on potential ids
          i.tree.match.t3 <- i.match.t3 %>% filter(tree.match %in% potential.tree.ids)
        }
      }
      
      if(nrow(i.tree.match.t3)==1){
        eastwide.plot.match[j,]$cur_tree_cn <- i.tree.match.t3$cn
        eastwide.plot.match[j,]$cur_tree_fiadb_status <- i.tree.match.t3$status_t3
        eastwide.plot.match[j,]$cur_tpa_unadj <- i.tree.match.t3$tpa_unadj_t3
        eastwide.plot.match[j,]$cur_match_note <- "matched by potential tree id & status"
        
        if(eastwide.plot.match[j,]$status == 5 | eastwide.plot.match[j,]$status == 4){
          eastwide.plot.match[j,]$cur_match_note <- "matched by potential tree id & status; needed to change snag/salvegable status to 2"
        }
        
        
        # if the tree number matches up with the t3 tree, then select that one
        i.tree.match.t2 <- i.match.t2 %>% filter(tree.match %in% potential.tree.ids)
        
        if(nrow(i.tree.match.t2)==1){
          i.tree.match.t2 <- i.tree.match.t2
          eastwide.plot.match[j,]$old_match_note <- "matched by potential tree id"
          
        }else{ # if there are still more tree matches and all are 0 or live, just take the first one
          
          if(2 %in% unique(i.match.t2$status_t2)){
            only.live.t2 <- i.match.t2 %>% filter(!status_t2 == 2)
            if(length(only.live.t2$cn) == 1){
              i.tree.match.t2 <- only.live.t2
              eastwide.plot.match[j,]$old_match_note <- "treeids don't match; selected only non-dead matching tree"
            }else{
              i.tree.match.t2 <- only.live.t2[1,]
              eastwide.plot.match[j,]$old_match_note <- "treeids don't match; selected first non-dead matching tree"
            }
          }else{
            i.tree.match.t2 <- i.match.t2[1,]
            eastwide.plot.match[j,]$old_match_note <- "treeids don't match; selected first matching tree"
          }  
        }
        eastwide.plot.match[j,]$old_tree_cn <- i.tree.match.t2$cn
        eastwide.plot.match[j,]$old_tree_fiadb_status <- i.tree.match.t2$status_t2
        eastwide.plot.match[j,]$old_tpa_unadj <- i.tree.match.t2$tpa_unadj_t2
        
      }else{
        i.tree.subp.match.t3 <- i.tree.match.t3 %>% filter(subp == 101)
        
        if(length(i.tree.subp.match.t3$cn) == 1){
          eastwide.plot.match[j,]$cur_tree_cn <-  i.tree.subp.match.t3$cn
          eastwide.plot.match[j,]$cur_tree_fiadb_status <-  i.tree.subp.match.t3$status_t3
          eastwide.plot.match[j,]$cur_tpa_unadj <-  i.tree.subp.match.t3$tpa_unadj_t3
          eastwide.plot.match[j,]$cur_subp <- i.tree.subp.match.t3$subp
          eastwide.plot.match[j,]$cur_match_note <- "matched by potential tree id, status, and subp"
          
          
          
        }
        else{
          
          # at least one case where prev_dia is useful
          i.tree.match.t3$prevdia_t3 <- as.character(i.tree.match.t3$prevdia_t3)
          i.tree.subp.match.t3 <- i.tree.match.t3 %>% filter(prevdia_t3 == eastwide.plot.match[j,]$dbhold)
          
          if(length(i.tree.subp.match.t3$cn) == 1){
            eastwide.plot.match[j,]$cur_tree_cn <-  i.tree.subp.match.t3$cn
            eastwide.plot.match[j,]$cur_tree_fiadb_status <-  i.tree.subp.match.t3$status_t3
            eastwide.plot.match[j,]$cur_tpa_unadj <-  i.tree.subp.match.t3$tpa_unadj_t3
            eastwide.plot.match[j,]$cur_subp <- i.tree.subp.match.t3$subp
            eastwide.plot.match[j,]$cur_match_note <- "matched by potential tree id and prevdia column"
          }
          
          eastwide.plot.match[j,]$cur_tree_cn <- "multiple matches" #i.match.t3$cn
        }
        
      }
      
      # once we have a single current tree, see if the time 2 tre matches up
      i.tree.match.t2 <- i.match.t2 %>% filter(tree.match %in% potential.tree.ids)
      
      if(nrow(i.tree.match.t2)==1){
        i.tree.match.t2 <- i.tree.match.t2
        eastwide.plot.match[j,]$old_match_note <- "matched by potential tree id"
        
      }else{ # if there are still more tree matches, just take the first one
        i.tree.match.t2 <- i.match.t2[1,]
        eastwide.plot.match[j,]$old_match_note <- "treeids don't match; selected first matching tree"
        
      }
      eastwide.plot.match[j,]$old_tree_cn <- i.tree.match.t2$cn
      eastwide.plot.match[j,]$old_tree_fiadb_status <- i.tree.match.t2$status_t2
      eastwide.plot.match[j,]$old_tpa_unadj <- i.tree.match.t2$tpa_unadj_t2
      eastwide.plot.match[j,]$old_subp <- i.tree.match.t2$subp
      
      
      
    }
    
    # if there is only multiple matching trees in the t3 
    if(length(i.match.t3$cn) > 1  & length(i.match.t2$cn) <= 1){
      # if the later inventory has mulitple matches, see if we can use tree.id or status to match up the tree 
      
      i.tree.match.t3 <- i.match.t3 %>% filter(tree.match %in% potential.tree.ids &
                                                 status_t3 == eastwide.plot.match[j,]$status) 
      # create a special case for snags:
      if(eastwide.plot.match[j,]$status == 5 | eastwide.plot.match[j,]$status == 4){
        
        i.tree.match.t3 <- i.match.t3 %>% filter(tree.match %in% potential.tree.ids & status_t3 == 2)
        
      }
      
      # if the current FIADB status code is zero for the only matching tree id, we could match on tree ids alone?
      
      if(length(i.tree.match.t3$cn) == 0 ){ # if there are still no matches:
        tree.id.only.match <- i.match.t3 %>% filter(tree.match %in% potential.tree.ids)
        if(length(tree.id.only.match$status_t3) == 1){# and the fiadb says status cd == 0, then just match on potential ids
          i.tree.match.t3 <- i.match.t3 %>% filter(tree.match %in% potential.tree.ids)
        }else{
          # if there are still multiple matches, find the tree that matches statuscd (if it is marked dead) and subplot
          if(length(tree.id.only.match)> 1 & eastwide.plot.match[j,]$status > 1){
            i.tree.subp.match.t3 <- i.match.t3 %>% filter(subp %in% 101 & 
                                                            status_t3 == 2)
            if(i.tree.subp.match.t3$cn == 0){ #if no trees match now, just select by statuscd
              i.tree.subp.match.t3 <- i.match.t3 %>% filter(status_t3 == 2)
            }
          }
          i.tree.match.t3 <- i.tree.subp.match.t3
        }
      }
      
      if(nrow(i.tree.match.t3)==1){
        eastwide.plot.match[j,]$cur_tree_cn <- i.tree.match.t3$cn
        eastwide.plot.match[j,]$cur_tree_fiadb_status <- i.tree.match.t3$status_t3
        eastwide.plot.match[j,]$cur_tpa_unadj <- i.tree.match.t3$tpa_unadj_t3
        eastwide.plot.match[j,]$cur_tpa_unadj <- i.tree.match.t3$subp
        eastwide.plot.match[j,]$cur_match_note <- "matched by potential tree id & status"
        
        if(eastwide.plot.match[j,]$status == 5 | eastwide.plot.match[j,]$status == 4 | eastwide.plot.match[j,]$status == 3){
          eastwide.plot.match[j,]$cur_match_note <- "matched by potential tree id & status; needed to change snag/salvegable status to 2"
        }
        # if(length(tree.id.only.match$status_t3) == 1 & tree.id.only.match$status_t3 == 0){
        #   eastwide.plot.match[j,]$cur_match_note <- paste0("matched by potential tree id only, FIADB statuscd == 0 but EWDB statcd ==", eastwide.plot.match[j,]$status)
        # }
        
      }else{
        # if there are still multiple matches, see if you can select the tree in suplot 101
        i.tree.subp.match.t3 <- i.tree.match.t3 %>% filter(subp == 101)
        if(length(i.tree.subp.match.t3$cn) == 1){
          eastwide.plot.match[j,]$cur_tree_cn <-  i.tree.subp.match.t3$cn
          eastwide.plot.match[j,]$cur_tree_fiadb_status <-  i.tree.subp.match.t3$status_t3
          eastwide.plot.match[j,]$cur_tpa_unadj <-  i.tree.subp.match.t3$tpa_unadj_t3
          eastwide.plot.match[j,]$cur_tpa_unadj <- i.tree.subp.match.t3$subp
          eastwide.plot.match[j,]$cur_match_note <- "matched by potential tree id, status, and subp"
        }
        else{
          eastwide.plot.match[j,]$cur_tree_cn <- "multiple matches"#i.match.t3$cn
        }
      }
      
    }
    # if there is more than one match for the old dataset, add a -999
    if(length(i.match.t2$cn) > 1 & length(i.match.t3$cn) <=1){
      # if the tree number matches up with the t3 tree, then select that one
      i.tree.match.t2 <- i.match.t2 %>% filter(tree.match %in% potential.tree.ids)
      
      if(nrow(i.tree.match.t2)==1){
        i.tree.match.t2 <- i.tree.match.t2
        eastwide.plot.match[j,]$old_match_note <- "matched by potential tree id"
      }else{ # if there are still more tree matches, just take the first one
        i.tree.match.t2 <- i.match.t2[1,]
        eastwide.plot.match[j,]$old_match_note <- "treeids don't match; selected first matching tree"
        
      }
      eastwide.plot.match[j,]$old_tree_cn <- i.tree.match.t2$cn
      eastwide.plot.match[j,]$old_tree_fiadb_status <- i.tree.match.t2$status_t2
      eastwide.plot.match[j,]$old_tpa_unadj <- i.tree.match.t2$tpa_unadj_t2
      eastwide.plot.match[j,]$old_tpa_unadj <- i.tree.match.t2$subp
      
      
    }
    
    
    
    
  }
  
  
  eastwide.VT.fiadb.match[[i]] <- eastwide.plot.match
  rm(eastwide.plot.match)
}
eastwide.VT.fiadb.updated.match.df <- do.call(rbind, eastwide.VT.fiadb.match)

eastwide.VT.fiadb.match.df <- do.call(rbind, eastwide.VT.fiadb.match)
eastwide.VT.fiadb.match.df %>% group_by(!is.na(cur_tree_cn) & !is.na(old_tree_cn) & !old_tree_cn %in% "multiple matches" & !cur_tree_cn %in% "multiple matches") %>% summarise(n())
# 3881 tree records have matching values in FIADB

single.tree.matches <- eastwide.VT.fiadb.match.df %>% filter(!is.na(cur_tree_cn) & !is.na(old_tree_cn) & !old_tree_cn %in% "multiple matches" & !cur_tree_cn %in% "multiple matches")
single.tree.matches %>% mutate(fiadb.status.change = ifelse(old_tree_fiadb_status == 1 & cur_tree_fiadb_status == 1, "live to live", 
                                                            ifelse(old_tree_fiadb_status == 1 & cur_tree_fiadb_status ==2, "live to dead", 
                                                                   ifelse(old_tree_fiadb_status == 1 & cur_tree_fiadb_status == 3, "live to cut/removed", 
                                                                          ifelse(old_tree_fiadb_status == 2 & cur_tree_fiadb_status == 2, "dead to dead", "no status at t3"))))) %>%
  group_by(fiadb.status.change) %>% summarise(n()) |> gt() |> gtsave("images/filtering_exploration/VT_fiadb_to_eastwide_matching_status_summary.png")

single.tree.matches %>% mutate(fiadb.status.change = ifelse(old_tree_fiadb_status == 1 & cur_tree_fiadb_status == 1, "live to live", 
                                                            ifelse(old_tree_fiadb_status == 1 & cur_tree_fiadb_status ==2, "live to dead", 
                                                                   ifelse(old_tree_fiadb_status == 1 & cur_tree_fiadb_status == 3, "live to cut/removed", 
                                                                          ifelse(old_tree_fiadb_status == 2 & cur_tree_fiadb_status == 2, "dead to dead", "no status at t3"))))) %>%
  rename("eastwide status" = "status")%>% group_by(fiadb.status.change, `eastwide status`) %>% summarise(n()) %>%
  ungroup()|> gt() |> gtsave("images/filtering_exploration/VT_fiadb_to_eastwide_matching_status_summary_ew.png")

# how many trees have matches now?
eastwide.VT.fiadb.updated.match.df %>% group_by(!is.na(cur_tree_cn) & !is.na(old_tree_cn) & !old_tree_cn %in% "multiple matches" & !cur_tree_cn %in% "multiple matches") %>% summarise(n())

# only 3 trees are unmatched for both time series
# 15 trees only have matches in the past
# 2215 trees are matched in the current inventory but not the past
eastwide.VT.fiadb.updated.match.df %>% group_by(cur_match_note %in% "no matching dbh and species", old_match_note %in% "no matching dbh and species") %>% summarise(n())

eastwide.VT.fiadb.updated.match.df %>% filter(cur_match_note %in% "matched by potential tree id, status, and subp" &
                                                is.na(old_match_note)) %>% select(PLOT.ID)
unique(eastwide.VT.match$PLOT.ID) %in% "33_9_382" 
# compare to the past matches?
eastwide.VT.fiadb.match.df %>% group_by(!is.na(cur_tree_cn) & !is.na(old_tree_cn) & !old_tree_cn %in% "multiple matches" & !cur_tree_cn %in% "multiple matches") %>% summarise(n())
#lets see how most of the trees matched up:
eastwide.VT.fiadb.updated.match.df %>% group_by(!is.na(cur_tree_cn) & !is.na(old_tree_cn) & 
                                                  !old_tree_cn %in% "multiple matches" & 
                                                  !cur_tree_cn %in% "multiple matches", 
                                                old_match_note, cur_match_note
) %>% summarise(n())

NA.consistent.status <- eastwide.VT.fiadb.updated.match.df %>% mutate(cur_status_consistent = ifelse(status == 1 & cur_tree_fiadb_status == 1, "yes, live", 
                                                                                                    ifelse(status > 1 & cur_tree_fiadb_status > 1, "yes, dead", 
                                                                                                           ifelse(status > 1 & cur_tree_fiadb_status == 0, "no, EWDB = dead, current status = 0", 
                                                                                                                  ifelse(status == 1 & cur_tree_fiadb_status == 0, "no, EWDB = live, current status = 0", 
                                                                                                                         ifelse(is.na(cur_tree_fiadb_status), "no matching tree in current FIADB", NA)))))) %>% 
  group_by(cur_status_consistent) %>% filter(is.na(cur_status_consistent))
# yes then mostly align, with the exception of current status == 0
# verify that these are status code changes:
eastwide.VT.fiadb.updated.match.df %>% mutate(status_change = ifelse(old_tree_fiadb_status == 1 & cur_tree_fiadb_status == 1, "live at both times", 
                                                                     ifelse(old_tree_fiadb_status > 1 & cur_tree_fiadb_status > 1, "dead at both times", 
                                                                            ifelse(old_tree_fiadb_status == 1 & cur_tree_fiadb_status > 1, "live at t2, dead at t3", 
                                                                                   ifelse(old_tree_fiadb_status > 1 & cur_tree_fiadb_status == 0, "dead at t2, zero at t3", 
                                                                                          ifelse(old_tree_fiadb_status == 1 & cur_tree_fiadb_status == 0, "live at t2, status at t3 = 0", 
                                                                                                 ifelse(old_tree_fiadb_status == 0 & cur_tree_fiadb_status > 1, "status at t2 = 0, live at t3", 
                                                                                                        ifelse(old_tree_fiadb_status > 1 & cur_tree_fiadb_status == 1, "status at t2 = dead, but status at t3 == live",
                                                                                                               ifelse(is.na(old_tree_fiadb_status) | is.na(cur_tree_fiadb_status),"NA values for at least one status","no matching trees in FIADB"))))))))) %>% 
  group_by(status_change) %>% #filter(status_change %in% "no matching trees in FIADB") %>% select(old_tree_fiadb_status, cur_tree_fiadb_status, old_match_note , cur_match_note)#%>% 
  mutate(status_change = ifelse(is.na(status_change) == TRUE, "could not match tree to current or old FIADB", status_change))%>%
  #select(tree, status, cur_tree_fiadb_status)
  summarise(`# of trees` = n()) |> gt()|>
  gtsave("images/filtering_exploration/VT_fiadb_to_eastwide_matching_status_cd_change.png")

eastwide.VT.fiadb.updated.match.df %>%  mutate(status_change = ifelse(old_tree_fiadb_status == 1 & cur_tree_fiadb_status == 1, "live at both times", 
                                                                      ifelse(old_tree_fiadb_status > 1 & cur_tree_fiadb_status > 1, "dead at both times", 
                                                                             ifelse(old_tree_fiadb_status == 1 & cur_tree_fiadb_status > 1, "live at t2, dead at t3", 
                                                                                    ifelse(old_tree_fiadb_status > 1 & cur_tree_fiadb_status == 0, "dead at t2, zero at t3", 
                                                                                           ifelse(old_tree_fiadb_status == 1 & cur_tree_fiadb_status == 0, "live at t2, status at t3 = 0", 
                                                                                                  ifelse(old_tree_fiadb_status == 0 & cur_tree_fiadb_status > 1, "status at t2 = 0, live at t3", 
                                                                                                         ifelse(old_tree_fiadb_status > 1 & cur_tree_fiadb_status == 1, "status at t2 = dead, but status at t3 == live",
                                                                                                                ifelse(is.na(old_tree_fiadb_status) | is.na(cur_tree_fiadb_status),"NA values for at least one status","no matching trees in FIADB"))))))))) %>% 
  group_by(status_change, status) %>% #filter(status_change %in% "no matching trees in FIADB") %>% select(old_tree_fiadb_status, cur_tree_fiadb_status, old_match_note , cur_match_note)#%>% 
  mutate(status_change = ifelse(is.na(status_change) == TRUE & is.na(cur_tree_fiadb_status) == TRUE, "no matching tree in current FIADB", 
                                ifelse(is.na(status_change) == TRUE & is.na(old_tree_fiadb_status) == TRUE, "no matching tree in old FIADB", status_change)))%>%
  #select(tree, status, cur_tree_fiadb_status)
  summarise(`# of trees` = n()) %>% ungroup()|> gt()|>
  gtsave("images/filtering_exploration/VT_fiadb_to_eastwide_matching_status_cd_change_by_EW_status.png")

# is it possible that the dead at both times categories are a mismatch in t2?

EW.VT.fiadb.matches.status <- eastwide.VT.fiadb.updated.match.df %>%  mutate(status_change = ifelse(old_tree_fiadb_status == 1 & cur_tree_fiadb_status == 1, "live at both times", 
                                                                                                    ifelse(old_tree_fiadb_status > 1 & cur_tree_fiadb_status > 1, "dead at both times", 
                                                                                                           ifelse(old_tree_fiadb_status == 1 & cur_tree_fiadb_status > 1, "live at t2, dead at t3", 
                                                                                                                  ifelse(old_tree_fiadb_status > 1 & cur_tree_fiadb_status == 0, "dead at t2, zero at t3", 
                                                                                                                         ifelse(old_tree_fiadb_status == 1 & cur_tree_fiadb_status == 0, "live at t2, status at t3 = 0", 
                                                                                                                                ifelse(old_tree_fiadb_status == 0 & cur_tree_fiadb_status > 1, "status at t2 = 0, dead at t3", 
                                                                                                                                       ifelse(old_tree_fiadb_status > 1 & cur_tree_fiadb_status == 1, "status at t2 = dead, but status at t3 == live",
                                                                                                                                              ifelse(is.na(old_tree_fiadb_status) | is.na(cur_tree_fiadb_status),"NA values for at least one status","no matching trees in FIADB"))))))))) %>% 
  group_by(status_change, status) %>% #filter(status_change %in% "no matching trees in FIADB") %>% select(old_tree_fiadb_status, cur_tree_fiadb_status, old_match_note , cur_match_note)#%>% 
  mutate(status_change = ifelse(is.na(status_change) == TRUE & is.na(cur_tree_fiadb_status) == TRUE, "no matching tree in current FIADB", 
                                ifelse(is.na(status_change) == TRUE & is.na(old_tree_fiadb_status) == TRUE, "no matching tree in old FIADB", status_change)))


# explore some of the different/non matching status changes
unique(EW.VT.fiadb.matches.status$status_change)
# 3 dead to live matches:
# two are simple matches for old fia with just one tree (possibly a tree misclassified as dead, or a mismatch of a tree)
View(EW.VT.fiadb.matches.status %>% filter(status_change %in% "status at t2 = dead, but status at t3 == live")%>% 
       select(old_subp, cur_subp, tree, status, cur_tree_fiadb_status, old_tree_fiadb_status, old_match_note, cur_match_note))

# a lot of the live to dead matches are cut trees
EW.VT.fiadb.matches.status %>% filter(status_change %in% "live at t2, dead at t3" )%>% 
  select(old_subp, cur_subp, tree, status, cur_tree_fiadb_status, old_tree_fiadb_status, old_match_note, cur_match_note) %>% 
  group_by(status, cur_tree_fiadb_status) %>%
  summarise(n())

# most of the status == 0, but dead in FIADB are status == 5 in EWDB
EW.VT.fiadb.matches.status %>% filter(status_change %in% "status at t2 = 0, dead at t3" )%>% 
  select(old_subp, cur_subp, tree, status, cur_tree_fiadb_status, old_tree_fiadb_status, old_match_note, cur_match_note) %>% 
  group_by(status, cur_tree_fiadb_status) %>%
  summarise(n())

# most of the status == 0, but dead in FIADB are status == 5 in EWDB
EW.VT.fiadb.matches.status %>% filter(status_change %in% "dead at both times" )%>% 
  select(old_subp, cur_subp, tree, status, cur_tree_fiadb_status, old_tree_fiadb_status, old_match_note, cur_match_note) %>% 
  group_by(old_match_note, cur_match_note) %>%
  summarise(n())

dead.at.both.VT <- EW.VT.fiadb.matches.status %>% filter(status_change %in% "dead at both times" ) 
unique(dead.at.both.VT$PLOT.ID)
dead.at.both.VT.potential.mismatch <- dead.at.both.VT %>% filter(old_match_note %in% c("matched by potential tree id", 
                                                                                       "treeids don't match; selected first matching tree"))

unique(dead.at.both.VT.potential.mismatch$PLOT.ID)

eastwide.VT.fiadb.updated.match.df %>% filter(is.na(cur_match_note)) %>% select(PLOT.ID)

# summarise how all these trees were matched up:
eastwide.VT.fiadb.updated.match.df %>% 
  group_by(cur_match_note) %>% #filter(status_change %in% "dead at both times")%>% 
  #select(tree, status, cur_tree_fiadb_status)
  summarise(`# of trees` = n()) %>% ungroup()|> gt()|>
  gtsave("images/filtering_exploration/VT_fiadb_to_eastwide_matching_notes_T3.png")
# 

eastwide.VT.fiadb.updated.match.df %>% 
  group_by(old_match_note) %>% #filter(status_change %in% "dead at both times")%>% 
  #select(tree, status, cur_tree_fiadb_status)
  summarise(`# of trees` = n()) %>% ungroup()|> gt()|>
  gtsave("images/filtering_exploration/VT_fiadb_to_eastwide_matching_notes_T3.png")


ggplot(eastwide.VT.fiadb.updated.match.df, aes(x = dbhcur, y = cur_tpa_unadj))+geom_point()
ggplot(eastwide.VT.fiadb.updated.match.df, aes(x = dbhcur, y = volfac))+geom_point()
ggplot(eastwide.VT.fiadb.updated.match.df, aes(x = dbhold, y = old_tpa_unadj))+geom_point()

#####################################e###############################################
# do NY matches in FIADB match the trees in CT from eastewide?
# state 33 New Hampshire
NE_designcd %>% filter(STATECD == 36)

NE_designcd %>% filter(INVYR < 2000) %>% select(DESIGNCD) %>% distinct()

TREE %>% select(stname, date) %>% distinct()
# designcd 101 for 1983 and 1997-- subps 100 for 1983 and 101, 102, 103, 104
# for 1997 designcd 111


ocon <- dbConnect(odbc(), "fiadb01p")
NY_plot <- dbGetQuery(ocon, 
                      "select cn, statecd, unitcd, countycd, plot, invyr, designcd
FROM fs_fiadb.plot
WHERE statecd =  ANY(36) 
                      and invyr < 2000
") %>%
  as_tibble()%>%
  rename_with(tolower) %>% 
  rename("plt_cn" = "cn") %>%
  mutate(PLOT.ID = paste0(statecd,  "_", countycd,  "_", plot))
# get NY_plot
CN.plot <- NY_plot %>% mutate(PLOT.ID = paste0(statecd, "_", countycd, "_", plot))%>% group_by(PLOT.ID) %>% summarise(n())
remeas.plots <- CN.plot %>% filter(`n()` ==2)

NY_tree <- dbGetQuery(ocon, "select cn, plt_cn, prev_tre_cn, statecd, unitcd, countycd, plot, subp, tree, invyr, dia, statuscd, spcd, prevdia
                      FROM fs_fiadb.tree
                      WHERE statecd =  ANY(36) 
                      and invyr < 2000
                      ") %>%
  as_tibble()%>%
  rename_with(tolower)#%>% left_join(., NY_plot)
summary(NY_tree$prevdia)
NY_tree %>% filter(!is.na(prevdia)) %>% select(invyr, statecd) %>% distinct()
# there is no 1980 survey in FIADB


#####################################e###############################################
# do ME matches in FIADB match the trees in CT from eastewide?
# state 33 New Hampshire
NE_designcd %>% filter(STATECD == 23)
TREE %>% select(stname, date) %>% distinct()
# designcd 101 for 1983 and 1997-- subps 100 for 1983 and 101, 102, 103, 104
# for 1997 designcd 111


ocon <- dbConnect(odbc(), "fiadb01p")
ME_plot <- dbGetQuery(ocon, 
                      "select cn, statecd, unitcd, countycd, plot, invyr, designcd
FROM fs_fiadb.plot
WHERE statecd =  ANY(23) 
                      and invyr < 2000
") %>%
  as_tibble()%>%
  rename_with(tolower) %>% 
  rename("plt_cn" = "cn") %>%
  mutate(PLOT.ID = paste0(statecd,  "_", countycd,  "_", plot))
# get ME_plot
CN.plot <- ME_plot %>% mutate(PLOT.ID = paste0(statecd, "_", countycd, "_", plot))%>% group_by(PLOT.ID) %>% summarise(n())
remeas.plots <- CN.plot %>% filter(`n()` ==2)

ME_tree <- dbGetQuery(ocon, "select cn, plt_cn, prev_tre_cn, statecd, unitcd, countycd, plot, subp, tree, invyr, dia, statuscd, spcd, prevdia
                      FROM fs_fiadb.tree
                      WHERE statecd =  ANY(23) 
                      and invyr < 2000
                      ") %>%
  as_tibble()%>%
  rename_with(tolower)%>% left_join(., ME_plot)


# subplot 101 may be remeasured between 1985 and 1998
ME_tree %>% filter(PLOT.ID %in% unique(remeas.plots$PLOT.ID)[3]) %>% group_by(subp, invyr) %>% summarise(n())
View(ME_tree %>% filter(PLOT.ID %in% unique(remeas.plots$PLOT.ID)[1]))  #%>% summarise(n())
ME_tree %>% filter(!is.na(prevdia)) %>% select(invyr, statecd) %>% distinct()
