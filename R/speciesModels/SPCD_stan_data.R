SPCD.stan.data <- function(SPCD.id, remper.correction, cleaned.data.full){
  cleaned.data <- cleaned.data.full %>% filter(SPCD %in% SPCD.id) %>% 
                                        filter(dbhold > 5 & ! dbhcur-dbhold == 0)
  
  cleaned.data_old <- cleaned.data.full %>% filter(SPCD %in% SPCD.id)
  
  cleaned.data %>% group_by(M) %>% summarise(n())
  
  cleaned.data_old %>% group_by(M, dbhcur-dbhold == 0) %>% summarise(n())
  
  
  View(cleaned.data_old %>% select(M))
 
if(remper.correction == "random"){
  # uniform sample across the board for the remper year correction
 
  # scale the cleaned data tree-level diameters by species
  cleaned.data <- cleaned.data %>% ungroup()  %>% group_by(SPCD) %>% 
    mutate(remper.sample = runif(length(cleaned.data$state)))%>%
    mutate(rempercur = ifelse(M ==1, remper*remper.sample, remper), 
           annual.growth = DIA_DIFF/rempercur) %>% mutate(DIA.median = median(dbhcur, na.rm =TRUE), 
                                                          DIA.sd = sd(dbhcur, na.rm = TRUE),  
                                                          DIA.DIFF.median = median(DIA_DIFF, na.rm =TRUE), 
                                                          DIA.DIFF.sd = sd(DIA_DIFF, na.rm =TRUE),
                                                          BAL.median = median(BAL, na.rm=TRUE),
                                                          BAL.sd = sd(BAL, na.rm = TRUE),
                                                          RD.median = median(RD, na.rm=TRUE), 
                                                          RD.sd = sd(RD, na.rm =TRUE),
                                                          nonSPCD_BA_tot.sd = sd(non_SPCD_BA, na.rm = TRUE),
                                                          SPCD_BA.sd = sd(SPCD_BA, na.rm =
                                                                            TRUE),
                                                          prop.focal.ba.median = median(SPCD_BA/BA_total, na.rm =TRUE), 
                                                          prop.focal.ba.sd = sd(SPCD_BA/BA_total, na.rm =TRUE), 
                                                          BA_tot.median = median(BA_total, na.rm =
                                                                                   TRUE),
                                                          nonSPCD_BA_tot.median = median(non_SPCD_BA, na.rm = TRUE),
                                                          SPCD_BA.median = median(SPCD_BA, na.rm =
                                                                                    TRUE),
                                                          annual.growth.median = median(annual.growth, na.rm = TRUE), 
                                                          annual.growth.sd = sd(annual.growth, na.rm = TRUE)) %>% 
    ungroup() %>% mutate(DIA_scaled = (dbhcur - DIA.median)/DIA.sd,
                         DIA_DIFF_scaled = (DIA_DIFF - DIA.DIFF.median)/DIA.DIFF.sd,
                         annual.growth.scaled = (annual.growth - annual.growth.median)/annual.growth.sd,
                         RD.scaled = (RD-RD.median)/RD.sd,
                         BAL.scaled = (BAL-BAL.median)/BAL.sd,
                         SPCD.BA.scaled = (SPCD_BA - SPCD_BA.median)/SPCD_BA.sd,
                         non_SPCD.BA.scaled = (non_SPCD_BA - nonSPCD_BA_tot.median)/nonSPCD_BA_tot.sd,
                         prop.focal.ba.scaled = ((SPCD_BA/BA_total) - prop.focal.ba.median)/prop.focal.ba.sd, 
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
    }else{
    cleaned.data <- cleaned.data %>% ungroup()  %>% group_by(SPCD) %>% 
      mutate(rempercur = ifelse(M ==1, remper*remper.correction, remper), 
             annual.growth = DIA_DIFF/rempercur) %>% mutate(DIA.median = median(dbhcur, na.rm =TRUE), 
                                                            DIA.sd = sd(dbhcur, na.rm = TRUE),  
                                                            DIA.DIFF.median = median(DIA_DIFF, na.rm =TRUE), 
                                                            DIA.DIFF.sd = sd(DIA_DIFF, na.rm =TRUE),
                                                            BAL.median = median(BAL, na.rm=TRUE),
                                                            BAL.sd = sd(BAL, na.rm = TRUE),
                                                            RD.median = median(RD, na.rm=TRUE), 
                                                            RD.sd = sd(RD, na.rm =TRUE),
                                                            nonSPCD_BA_tot.sd = sd(non_SPCD_BA, na.rm = TRUE),
                                                            SPCD_BA.sd = sd(SPCD_BA, na.rm =
                                                                              TRUE),
                                                            prop.focal.ba.median = median(SPCD_BA/BA_total, na.rm =TRUE), 
                                                            prop.focal.ba.sd = sd(SPCD_BA/BA_total, na.rm =TRUE), 
                                                            BA_tot.median = median(BA_total, na.rm =
                                                                                     TRUE),
                                                            nonSPCD_BA_tot.median = median(non_SPCD_BA, na.rm = TRUE),
                                                            SPCD_BA.median = median(SPCD_BA, na.rm =
                                                                                      TRUE),
                                                            annual.growth.median = median(annual.growth, na.rm = TRUE), 
                                                            annual.growth.sd = sd(annual.growth, na.rm = TRUE)) %>% 
      # rescale to 
      ungroup() %>% mutate(DIA_scaled = (dbhcur - DIA.median)/DIA.sd,
                           DIA_DIFF_scaled = (DIA_DIFF - DIA.DIFF.median)/DIA.DIFF.sd,
                           annual.growth.scaled = (annual.growth - annual.growth.median)/annual.growth.sd,
                           
                           RD.scaled = (RD-RD.median)/RD.sd,
                           BAL.scaled = (BAL-BAL.median)/BAL.sd,
                           SPCD.BA.scaled = (SPCD_BA - SPCD_BA.median)/SPCD_BA.sd,
                           non_SPCD.BA.scaled = (non_SPCD_BA - nonSPCD_BA_tot.median)/nonSPCD_BA_tot.sd,
                           prop.focal.ba.scaled = ((SPCD_BA/BA_total) - prop.focal.ba.median)/prop.focal.ba.sd, 
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
   }
  # rescale to values of 0 to 1
     # ungroup()%>% mutate(DIA_scaled = rescale(dbhold,, 
     #                    annual.growth.scaled = rescale(annual.growth,,
     #                    RD.scaled = rescale(RD,,
     #                    BAL.scaled = rescale(BAL,,
     #                    
     #                    SPCD.BA.scaled = rescale(SPCD_BA,,
     #                    non_SPCD.BA.scaled = rescale(non_SPCD_BA,,
     #                    prop.focal.ba = rescale(SPCD_BA/BA_total,, 
     #                    
     #                    density.scaled = rescale(density_total,,
     #                    SPCD.density.scaled = rescale(SPCD_density,,
     #                    non.SPCD.density.scaled = rescale(non_SPCD_density,,
     #                    prop.focal.density = rescale(SPCD_density/density_total,, 
     #                    
     #                    si.scaled = rescale(si,,
     #                    ba.scaled = rescale(BA_total,,
     #                    aspect.scaled = rescale(aspect,,
     #                    slope.scaled = rescale(slope,,
     #                    damage.scaled = rescale(damage,,
     #                    MAP.scaled = rescale(MAP,,
     #                    elev.scaled = rescale(elev,,
     #                    Ndep.scaled = rescale(Ndep.remper.avg,,
     #                    physio.scaled = rescale(physio,,
     #                    MATmin.scaled = rescale(MATmin,,
     #                    MATmax.scaled = rescale(MATmax,, 
     #                    ppt.anom = rescale(ppt.anom,, 
     #                    tmax.anom = rescale(tmax.anom,, 
     #                    tmin.anom = rescale(tmin.anom,)
  # old method of scaling                      
   

  SPP.df <- data.frame(SPCD = unique(cleaned.data$SPCD), 
                       SPP = 1:length(unique(cleaned.data$SPCD)))
  
  cleaned.data<- left_join(cleaned.data, SPP.df) 
  cleaned.data <- cleaned.data %>%  filter(!is.na(si) & !is.na(dbhcur)& !is.na(M) & 
                                             !is.na(annual.growth.scaled) & !is.na(ppt.anom))# & !is.na(prop.focal.ba) & !is.na(prop.focal.density))
  #summary(cleaned.data$BAL.scaled)
  cleaned.data$S <- ifelse(cleaned.data$M == 1, 0, 1)
  
  N_train <- nrow(cleaned.data)*0.7
  N_test <- nrow(cleaned.data)*0.3
  train_ind <- sample(c(1:nrow(cleaned.data)), size = N_train, replace = FALSE)
  
  train.data <- cleaned.data[train_ind,]
  test.data <- cleaned.data[-train_ind, ]
  
  # Option B:
  # 1. Annual growth 
  
  
  ggplot(test.data, aes(x= as.character(M), y = annual.growth))+geom_violin()+facet_wrap(~SPCD, scales = "free_y")
  ggplot(test.data, aes(x= as.character(M), y = prop.focal.ba.scaled))+geom_violin()+facet_wrap(~SPCD, scales = "free_y")
  ggplot(test.data, aes(x= as.character(M), y = SPCD.BA.scaled))+geom_violin()+facet_wrap(~SPCD, scales = "free_y")
  #ggplot(test.data, aes(x= as.character(M), y = SPCD.density.scaled))+geom_violin()+facet_wrap(~SPCD, scales = "free_y")
  ggplot(test.data, aes(x= as.character(M), y = non_SPCD.BA.scaled))+geom_violin()+facet_wrap(~SPCD, scales = "free_y")
  #ggplot(test.data, aes(x= as.character(M), y = non.SPCD.density.scaled))+geom_violin()+facet_wrap(~SPCD, scales = "free_y")
  
  
  
  # model.1 data
  mod.data.1 <- list(N = nrow(train.data), 
                     y = train.data$S,                       
                     Remper = train.data$remper, 
                     xM = as.matrix(train.data %>% dplyr::select(DIA_DIFF_scaled)))
  # model.2 data
  # 2. Diameter + Annual growth
  mod.data.2 <- list(N = nrow(train.data), 
                     y = train.data$S,                       
                     Remper = train.data$remper, 
                     xM = as.matrix(train.data %>% dplyr::select(DIA_DIFF_scaled, 
                                                                 DIA_scaled)))
  
  # model.3 data
  # 3. Diameter + Annual growth + competition variables (RD.scaled, BAL, damage)
  mod.data.3 <- list(N = nrow(train.data), 
                     y = train.data$S,                       
                     Remper = train.data$remper, 
                     xM = as.matrix(train.data %>% dplyr::select(DIA_DIFF_scaled, 
                                                                 DIA_scaled, 
                                                                 #RD.scaled, 
                                                                 ba.scaled, 
                                                                 BAL.scaled,
                                                                 #non_SPCD.BA.scaled,
                                                                 damage.scaled 
                                                                 
                                                                 )))
  # model.4 data
  # 4. Diameter + Annual growth + competition variables (RD.scaled, BAL, damage) + Climate variables
  
  mod.data.4 <- list(N = nrow(train.data), 
                     y = train.data$S,                       
                     Remper = train.data$remper, 
                     xM = as.matrix(train.data %>% dplyr::select(DIA_DIFF_scaled, 
                                                                 DIA_scaled, 
                                                                 #RD.scaled, 
                                                                 ba.scaled, 
                                                                 BAL.scaled, 
                                                                 #non_SPCD.BA.scaled,
                                                                 damage.scaled,
                                                                 
                                                                 MATmax.scaled, 
                                                                 #MATmin.scaled, 
                                                                 MAP.scaled,
                                                                 ppt.anom, 
                                                                 #tmin.anom, 
                                                                 tmax.anom)))
  # model.5 data
  # 5. Diameter + Annual growth + competition variables (RD.scaled, BAL, damage) + 
  #Climate variables + site/soil effects + ndep
  mod.data.5 <- list(N = nrow(train.data), 
                     y = train.data$S,                       
                     Remper = train.data$remper, 
                     xM = as.matrix(train.data %>% dplyr::select(DIA_DIFF_scaled, 
                                                                 DIA_scaled, 
                                                                 #RD.scaled, 
                                                                 ba.scaled, 
                                                                 BAL.scaled,
                                                                 #non_SPCD.BA.scaled,
                                                                 damage.scaled,
                                                                
                                                                 MATmax.scaled, 
                                                                 #MATmin.scaled, 
                                                                 MAP.scaled,
                                                                 ppt.anom, 
                                                                 #tmin.anom, 
                                                                 tmax.anom, 
                                                                 slope.scaled, 
                                                                 aspect.scaled,
                                                                 #elev.scaled, 
                                                                 Ndep.scaled))) #,
                                                                 #physio.scaled)))
  # model.6 data
  # 6. All Fixed effects and all growth + Diameter interactions
  mod.data.6 <- list(N = nrow(train.data), 
                     y = train.data$S,                       
                     Remper = train.data$remper, 
                     xM = as.matrix(train.data %>% dplyr::select(DIA_DIFF_scaled, 
                                                                 DIA_scaled, 
                                                                 #RD.scaled, 
                                                                 ba.scaled, 
                                                                 BAL.scaled, 
                                                                 #non_SPCD.BA.scaled,
                                                                 damage.scaled,
                                                                
                                                                 MATmax.scaled, 
                                                                 #MATmin.scaled, 
                                                                 MAP.scaled,
                                                                 ppt.anom, 
                                                                 #tmin.anom, 
                                                                 tmax.anom, 
                                                                 slope.scaled, 
                                                                 aspect.scaled,
                                                                 #elev.scaled, 
                                                                 Ndep.scaled)%>% #,
                                                                 #physio.scaled) %>%
                                      # generate growth interactions
                                      mutate_at(.funs = list(growth.int = ~.*DIA_DIFF_scaled), 
                                                .vars = vars(DIA_scaled:Ndep.scaled)) %>% 
                                      mutate_at(.funs = list(DIA.int = ~.*DIA_scaled), 
                                                .vars = vars(ba.scaled:Ndep.scaled)) )) #%>% 
                                     # mutate_at(scale, .vars = vars(DIA_scaled_growth.int:physio.scaled_DIA.int))))
  
  
  
  
  xM <- mod.data.6$xM
  summary(xM[,1]*xM[,2])
  summary(xM[,"DIA_scaled_growth.int"])
  # model.7 data
  # 7. model 5 + competition interactions
  mod.data.7 <- list(N = nrow(train.data), 
                     y = train.data$S,                       
                     Remper = train.data$remper, 
                     xM = as.matrix(train.data %>% dplyr::select(DIA_DIFF_scaled, 
                                                                 DIA_scaled, 
                                                                 #RD.scaled, 
                                                                 ba.scaled, 
                                                                 BAL.scaled, 
                                                                 #non_SPCD.BA.scaled,
                                                                 damage.scaled,
                                                                
                                                                 MATmax.scaled, 
                                                                 #MATmin.scaled, 
                                                                 MAP.scaled,
                                                                 ppt.anom, 
                                                                 #tmin.anom, 
                                                                 tmax.anom, 
                                                                 slope.scaled, 
                                                                 aspect.scaled,
                                                                 #elev.scaled, 
                                                                 Ndep.scaled) %>% #,
                                                                   #physio.scaled) %>%,
                                                                 #physio.scaled) %>%
                                      # generate growth interactions
                                      mutate_at(.funs = list(growth.int = ~.*DIA_DIFF_scaled), 
                                                .vars = vars(DIA_scaled:Ndep.scaled)) %>% 
                                      # generate diameter interactions
                                      mutate_at(.funs = list(DIA.int = ~.*DIA_scaled), 
                                                .vars = vars(ba.scaled:Ndep.scaled)) %>%
                                      
                                      
                                      # # generate RD.scaled interactions
                                      # mutate_at(.funs = list(RD.scaled.int = ~.*RD.scaled), 
                                      #           .vars = vars(ba.scaled:physio.scaled))%>%
                                      
                                      # generate ba interactions
                                      mutate_at(.funs = list(ba.int = ~.*ba.scaled), 
                                                .vars = vars(BAL.scaled:Ndep.scaled))%>%
                                      # generate BAL interactions
                                      mutate_at(.funs = list(BAL.int = ~.*BAL.scaled), 
                                                .vars = vars(damage.scaled:Ndep.scaled))%>%
                                      
                                      # generate damage interactions
                                      mutate_at(.funs = list(damage.int = ~.*damage.scaled), 
                                                .vars = vars(MATmax.scaled:Ndep.scaled))#%>%
                                      # make sure interactions terms are scaled to be closer to values
                                      #mutate_at(.funs = function(x)(x/10), .vars = vars(physio.scaled_damage.int))
                                    # mutate_at(.funs = function(x)(x/10), .vars = vars( RD.scaled_growth.int:physio.scaled_damage.int))
                     ))
  
  
  
  # model.8 data
  # 8. model 6 + climate interactions
  mod.data.8 <- list(N = nrow(train.data), 
                     y = train.data$S,                       
                     Remper = train.data$remper, 
                     xM = as.matrix(train.data %>% dplyr::select(DIA_DIFF_scaled, 
                                                                 DIA_scaled, 
                                                                 #RD.scaled, 
                                                                 ba.scaled, 
                                                                 BAL.scaled,
                                                                 #non_SPCD.BA.scaled,
                                                                 damage.scaled,
                                                                 
                                                                 MATmax.scaled, 
                                                                 #MATmin.scaled, 
                                                                 MAP.scaled,
                                                                 ppt.anom, 
                                                                 #tmin.anom, 
                                                                 tmax.anom, 
                                                                 slope.scaled, 
                                                                 aspect.scaled,
                                                                 #elev.scaled, 
                                                                Ndep.scaled) %>% #,
                                      #physio.scaled) %>%,
                                      #physio.scaled) %>%
                                      # generate growth interactions
                                      mutate_at(.funs = list(growth.int = ~.*DIA_DIFF_scaled), 
                                                .vars = vars(DIA_scaled:Ndep.scaled)) %>% 
                                      # generate diameter interactions
                                      mutate_at(.funs = list(DIA.int = ~.*DIA_scaled), 
                                                .vars = vars(ba.scaled:Ndep.scaled)) %>%
                                      
                                      
                                      # # generate RD.scaled interactions
                                      # mutate_at(.funs = list(RD.scaled.int = ~.*RD.scaled), 
                                      #           .vars = vars(ba.scaled:physio.scaled))%>%
                                      
                                      # generate ba interactions
                                      mutate_at(.funs = list(ba.int = ~.*ba.scaled), 
                                                .vars = vars(BAL.scaled:Ndep.scaled))%>%
                                      # generate BAL interactions
                                      mutate_at(.funs = list(BAL.int = ~.*BAL.scaled), 
                                                .vars = vars(damage.scaled:Ndep.scaled))%>%
                                      
                                      # generate damage interactions
                                      mutate_at(.funs = list(damage.int = ~.*damage.scaled), 
                                                .vars = vars(MATmax.scaled:Ndep.scaled)) %>%
                                      
                                      # generate MATmax interactions
                                      mutate_at(.funs = list(MATmax.scaled.int = ~.*MATmax.scaled), 
                                                .vars = vars(MAP.scaled:Ndep.scaled))%>%
                                      # # generate MATmin.scaled interactions
                                      # mutate_at(.funs = list(MATmin.scaled.int = ~.*MATmin.scaled), 
                                      #           .vars = vars(MAP.scaled:physio.scaled))%>%
                                      
                                      # generate MAP.scaled interactions
                                      mutate_at(.funs = list(MAP.scaled.int = ~.*MAP.scaled), 
                                                .vars = vars(ppt.anom:Ndep.scaled))%>%
                                      
                                      # generate ppt.anom interactions
                                      mutate_at(.funs = list(ppt.anom.int = ~.*ppt.anom), 
                                                .vars = vars(tmax.anom:Ndep.scaled))%>%
                                      # # generate tmin.anom interactions
                                      # mutate_at(.funs = list(tmin.anom.int = ~.*tmin.anom), 
                                      #           .vars = vars(tmax.anom:Ndep.scaled))%>%
                                      
                                      # generate tmax.anom interactions
                                      mutate_at(.funs = list(tmax.anom.int = ~.*tmax.anom), 
                                                .vars = vars(slope.scaled:Ndep.scaled))#%>%
                                      # make sure interactions terms are scaled to be closer to values
                                      #mutate_at(.funs = function(x)(x/10), .vars = vars(ba.scaled_RD.scaled.int:physio.scaled_damage.int))
                                      #mutate_at(.funs = function(x)(x/10), .vars = vars( RD.scaled_growth.int:physio.scaled_tmax.anom.int))
                     ))
  
  # model.9 data
  
  # 9. All Fixed effects and all interactions
  mod.data.9 <- list(N = nrow(train.data), 
                     y = train.data$S,                       
                     Remper = train.data$remper, 
                     xM = as.matrix(train.data %>% dplyr::select(DIA_DIFF_scaled, 
                                                                 DIA_scaled, 
                                                                 #RD.scaled, 
                                                                 ba.scaled, 
                                                                 BAL.scaled, 
                                                                 #non_SPCD.BA.scaled,
                                                                 damage.scaled,
                                                                 MATmax.scaled, 
                                                                 #MATmin.scaled, 
                                                                 MAP.scaled,
                                                                 ppt.anom, 
                                                                 #tmin.anom, 
                                                                 tmax.anom, 
                                                                 slope.scaled, 
                                                                 aspect.scaled,
                                                                 #elev.scaled, 
                                                                 Ndep.scaled) %>% #,
                                      #physio.scaled) %>%,
                                      #physio.scaled) %>%
                                      # generate growth interactions
                                      mutate_at(.funs = list(growth.int = ~.*DIA_DIFF_scaled), 
                                                .vars = vars(DIA_scaled:Ndep.scaled)) %>% 
                                      # generate diameter interactions
                                      mutate_at(.funs = list(DIA.int = ~.*DIA_scaled), 
                                                .vars = vars(ba.scaled:Ndep.scaled)) %>%
                                      
                                      
                                      # # generate RD.scaled interactions
                                      # mutate_at(.funs = list(RD.scaled.int = ~.*RD.scaled), 
                                      #           .vars = vars(ba.scaled:physio.scaled))%>%
                                      
                                      # generate ba interactions
                                      mutate_at(.funs = list(ba.int = ~.*ba.scaled), 
                                                .vars = vars(BAL.scaled:Ndep.scaled))%>%
                                      # generate BAL interactions
                                      mutate_at(.funs = list(BAL.int = ~.*BAL.scaled), 
                                                .vars = vars(damage.scaled:Ndep.scaled))%>%
                                      
                                      # generate damage interactions
                                      mutate_at(.funs = list(damage.int = ~.*damage.scaled), 
                                                .vars = vars(MATmax.scaled:Ndep.scaled)) %>%
                                      
                                      # generate MATmax interactions
                                      mutate_at(.funs = list(MATmax.scaled.int = ~.*MATmax.scaled), 
                                                .vars = vars(MAP.scaled:Ndep.scaled))%>%
                                      # # generate MATmin.scaled interactions
                                      # mutate_at(.funs = list(MATmin.scaled.int = ~.*MATmin.scaled), 
                                      #           .vars = vars(MAP.scaled:physio.scaled))%>%
                                      
                                      # generate MAP.scaled interactions
                                      mutate_at(.funs = list(MAP.scaled.int = ~.*MAP.scaled), 
                                                .vars = vars(ppt.anom:Ndep.scaled))%>%
                                      
                                      mutate_at(.funs = list(ppt.anom.int = ~.*ppt.anom), 
                                                .vars = vars(tmax.anom:Ndep.scaled))%>%
                                      # # generate tmin.anom interactions
                                      # mutate_at(.funs = list(tmin.anom.int = ~.*tmin.anom), 
                                      #           .vars = vars(tmax.anom:Ndep.scaled))%>%
                                      
                                      # generate tmax.anom interactions
                                      mutate_at(.funs = list(tmax.anom.int = ~.*tmax.anom), 
                                                .vars = vars(slope.scaled:Ndep.scaled)) %>%
                                      # generate slope.scaled interactions
                                      mutate_at(.funs = list(slope.int = ~.*slope.scaled), 
                                                .vars = vars(aspect.scaled:Ndep.scaled))%>%
                                      # generate aspect.scaled interactions
                                      mutate_at(.funs = list(aspect.int = ~.*aspect.scaled), 
                                                .vars = vars(Ndep.scaled))#%>%
                                      
                                      # # generate elev interactions
                                      # mutate_at(.funs = list(elev.int = ~.*elev.scaled), 
                                      #           .vars = vars(Ndep.scaled:physio.scaled)) %>%
                                      # generate Ndep interactions
                                      #mutate(Ndep.physio.int = physio.scaled*Ndep.scaled) #%>%
                                      #mutate_at(.funs = function(x)(x/10), .vars = vars( RD.scaled_growth.int:Ndep.physio.int))
                                    
                     ))
  
  
  ###############################################################################
  # do the same for the test data:
  
  # model.1 data
  mod.data.1.test <- list(N = nrow(test.data), 
                          y = test.data$S,                       
                          Remper = test.data$remper, 
                          xM = as.matrix(test.data %>% dplyr::select(DIA_DIFF_scaled)))
  # model.2 data
  # 2. Diameter + Annual growth
  mod.data.2.test <- list(N = nrow(test.data), 
                          y = test.data$S,                       
                          Remper = test.data$remper, 
                          xM = as.matrix(test.data %>% dplyr::select(DIA_DIFF_scaled, 
                                                                     DIA_scaled)))
  
  # model.3 data
  # 3. Diameter + Annual growth + competition variables (RD.scaled, BAL, damage)
  mod.data.3.test <- list(N = nrow(test.data), 
                          y = test.data$S,                       
                          Remper = test.data$remper, 
                          xM = as.matrix(test.data %>% dplyr::select(DIA_DIFF_scaled, 
                                                                     DIA_scaled, 
                                                                     #RD.scaled, 
                                                                     ba.scaled, 
                                                                     BAL.scaled, 
                                                                     #non_SPCD.BA.scaled,
                                                                     damage.scaled)))
  # model.4 data
  # 4. Diameter + Annual growth + competition variables (RD.scaled, BAL, damage) + Climate variables
  
  mod.data.4.test <- list(N = nrow(test.data), 
                          y = test.data$S,                       
                          Remper = test.data$remper, 
                          xM = as.matrix(test.data %>% dplyr::select(DIA_DIFF_scaled, 
                                                                     DIA_scaled, 
                                                                     #RD.scaled, 
                                                                     ba.scaled, 
                                                                     BAL.scaled, 
                                                                     #non_SPCD.BA.scaled,
                                                                     damage.scaled,
                                                                     MATmax.scaled, 
                                                                     #MATmin.scaled, 
                                                                     MAP.scaled,
                                                                     ppt.anom, 
                                                                     #tmin.anom, 
                                                                     tmax.anom)))
  # model.5 data
  # 5. Diameter + Annual growth + competition variables (RD.scaled, BAL, damage) + 
  #Climate variables + site/soil effects + ndep
  mod.data.5.test <- list(N = nrow(test.data), 
                          y = test.data$S,                       
                          Remper = test.data$remper, 
                          xM = as.matrix(test.data %>% dplyr::select(DIA_DIFF_scaled, 
                                                                     DIA_scaled, 
                                                                     #RD.scaled, 
                                                                     ba.scaled, 
                                                                     BAL.scaled, 
                                                                     #non_SPCD.BA.scaled,
                                                                     damage.scaled,
                                                                     MATmax.scaled, 
                                                                     #MATmin.scaled, 
                                                                     MAP.scaled,
                                                                     ppt.anom, 
                                                                     #tmin.anom, 
                                                                     tmax.anom, 
                                                                     slope.scaled, 
                                                                     aspect.scaled,
                                                                     #elev.scaled, 
                                                                     Ndep.scaled )))#,
                                                                     #physio.scaled)))
  # model.6 data
  # 6. All Fixed effects and all growth + Diameter interactions
  mod.data.6.test <- list(N = nrow(test.data), 
                          y = test.data$S,                       
                          Remper = test.data$remper, 
                          xM = as.matrix(test.data %>% dplyr::select(DIA_DIFF_scaled, 
                                                                     DIA_scaled, 
                                                                     #RD.scaled, 
                                                                     ba.scaled, 
                                                                     BAL.scaled, 
                                                                     #non_SPCD.BA.scaled,
                                                                     damage.scaled,
                                                                     
                                                                     MATmax.scaled, 
                                                                     #MATmin.scaled, 
                                                                     MAP.scaled,
                                                                     ppt.anom, 
                                                                     #tmin.anom, 
                                                                     tmax.anom, 
                                                                     slope.scaled, 
                                                                     aspect.scaled,
                                                                     #elev.scaled, 
                                                                     Ndep.scaled)%>% #,
                                           #physio.scaled) %>%
                                           # generate growth interactions
                                           mutate_at(.funs = list(growth.int = ~.*DIA_DIFF_scaled), 
                                                     .vars = vars(DIA_scaled:Ndep.scaled)) %>% 
                                           mutate_at(.funs = list(DIA.int = ~.*DIA_scaled), 
                                                     .vars = vars(ba.scaled:Ndep.scaled)) )) 
  
  
  
  
  
  
  # model.7 data
  # 7. model 5 + competition interactions
  mod.data.7.test <- list(N = nrow(test.data), 
                          y = test.data$S,                       
                          Remper = test.data$remper, 
                          xM = as.matrix(test.data %>% dplyr::select(DIA_DIFF_scaled, 
                                                                     DIA_scaled, 
                                                                     #RD.scaled, 
                                                                     ba.scaled, 
                                                                     BAL.scaled, 
                                                                     #non_SPCD.BA.scaled,
                                                                     damage.scaled,
                                                                     
                                                                     MATmax.scaled, 
                                                                     #MATmin.scaled, 
                                                                     MAP.scaled,
                                                                     ppt.anom, 
                                                                     #tmin.anom, 
                                                                     tmax.anom, 
                                                                     slope.scaled, 
                                                                     aspect.scaled,
                                                                     #elev.scaled, 
                                                                     Ndep.scaled) %>% #,
                                           #physio.scaled) %>%,
                                           #physio.scaled) %>%
                                           # generate growth interactions
                                           mutate_at(.funs = list(growth.int = ~.*DIA_DIFF_scaled), 
                                                     .vars = vars(DIA_scaled:Ndep.scaled)) %>% 
                                           # generate diameter interactions
                                           mutate_at(.funs = list(DIA.int = ~.*DIA_scaled), 
                                                     .vars = vars(ba.scaled:Ndep.scaled)) %>%
                                           
                                           
                                           # # generate RD.scaled interactions
                                           # mutate_at(.funs = list(RD.scaled.int = ~.*RD.scaled), 
                                           #           .vars = vars(ba.scaled:physio.scaled))%>%
                                           
                                           # generate ba interactions
                                           mutate_at(.funs = list(ba.int = ~.*ba.scaled), 
                                                     .vars = vars(BAL.scaled:Ndep.scaled))%>%
                                           # generate BAL interactions
                                           mutate_at(.funs = list(BAL.int = ~.*BAL.scaled), 
                                                     .vars = vars(damage.scaled:Ndep.scaled))%>%
                                           
                                           # generate damage interactions
                                           mutate_at(.funs = list(damage.int = ~.*damage.scaled), 
                                                     .vars = vars(MATmax.scaled:Ndep.scaled))#%>%
                                         # make sure interactions terms are scaled to be closer to values
                                         #mutate_at(.funs = function(x)(x/10), .vars = vars(physio.scaled_damage.int))
                                         # mutate_at(.funs = function(x)(x/10), .vars = vars( RD.scaled_growth.int:physio.scaled_damage.int))
                          ))
  
  
  
  
  # model.8 data
  # 8. model 6 + climate interactions
  mod.data.8.test <- list(N = nrow(test.data), 
                          y = test.data$S,                       
                          Remper = test.data$remper, 
                          xM = as.matrix(test.data %>% dplyr::select(DIA_DIFF_scaled, 
                                                                     DIA_scaled, 
                                                                     #RD.scaled, 
                                                                     ba.scaled, 
                                                                     BAL.scaled,
                                                                     #non_SPCD.BA.scaled,
                                                                     damage.scaled,
                                                                     
                                                                     MATmax.scaled, 
                                                                     #MATmin.scaled, 
                                                                     MAP.scaled,
                                                                     ppt.anom, 
                                                                     #tmin.anom, 
                                                                     tmax.anom, 
                                                                     slope.scaled, 
                                                                     aspect.scaled,
                                                                     #elev.scaled, 
                                                                     Ndep.scaled) %>% #,
                                           #physio.scaled) %>%,
                                           #physio.scaled) %>%
                                           # generate growth interactions
                                           mutate_at(.funs = list(growth.int = ~.*DIA_DIFF_scaled), 
                                                     .vars = vars(DIA_scaled:Ndep.scaled)) %>% 
                                           # generate diameter interactions
                                           mutate_at(.funs = list(DIA.int = ~.*DIA_scaled), 
                                                     .vars = vars(ba.scaled:Ndep.scaled)) %>%
                                           
                                           
                                           # # generate RD.scaled interactions
                                           # mutate_at(.funs = list(RD.scaled.int = ~.*RD.scaled), 
                                           #           .vars = vars(ba.scaled:physio.scaled))%>%
                                           
                                           # generate ba interactions
                                           mutate_at(.funs = list(ba.int = ~.*ba.scaled), 
                                                     .vars = vars(BAL.scaled:Ndep.scaled))%>%
                                           # generate BAL interactions
                                           mutate_at(.funs = list(BAL.int = ~.*BAL.scaled), 
                                                     .vars = vars(damage.scaled:Ndep.scaled))%>%
                                           
                                           # generate damage interactions
                                           mutate_at(.funs = list(damage.int = ~.*damage.scaled), 
                                                     .vars = vars(MATmax.scaled:Ndep.scaled)) %>%
                                           
                                           # generate MATmax interactions
                                           mutate_at(.funs = list(MATmax.scaled.int = ~.*MATmax.scaled), 
                                                     .vars = vars(MAP.scaled:Ndep.scaled))%>%
                                           # # generate MATmin.scaled interactions
                                           # mutate_at(.funs = list(MATmin.scaled.int = ~.*MATmin.scaled), 
                                           #           .vars = vars(MAP.scaled:physio.scaled))%>%
                                           
                                           # generate MAP.scaled interactions
                                           mutate_at(.funs = list(MAP.scaled.int = ~.*MAP.scaled), 
                                                     .vars = vars(ppt.anom:Ndep.scaled))%>%
                                           
                                           # generate ppt.anom interactions
                                           mutate_at(.funs = list(ppt.anom.int = ~.*ppt.anom), 
                                                     .vars = vars(tmax.anom:Ndep.scaled))%>%
                                           # # generate tmin.anom interactions
                                           # mutate_at(.funs = list(tmin.anom.int = ~.*tmin.anom), 
                                           #           .vars = vars(tmax.anom:Ndep.scaled))%>%
                                           
                                           # generate tmax.anom interactions
                                           mutate_at(.funs = list(tmax.anom.int = ~.*tmax.anom), 
                                                     .vars = vars(slope.scaled:Ndep.scaled))#%>%
                                         # make sure interactions terms are scaled to be closer to values
                                         #mutate_at(.funs = function(x)(x/10), .vars = vars(ba.scaled_RD.scaled.int:physio.scaled_damage.int))
                                         #mutate_at(.funs = function(x)(x/10), .vars = vars( RD.scaled_growth.int:physio.scaled_tmax.anom.int))
                          ))
  
  
  # model.9 data
  
  # 9. All Fixed effects and all interactions
  mod.data.9.test <- list(N = nrow(test.data), 
                          y = test.data$S,                       
                          Remper = test.data$remper, 
                          xM = as.matrix(test.data %>% dplyr::select(DIA_DIFF_scaled, 
                                                                     DIA_scaled, 
                                                                     #RD.scaled, 
                                                                     ba.scaled, 
                                                                     BAL.scaled, 
                                                                     #non_SPCD.BA.scaled,
                                                                     damage.scaled,
                                                                     MATmax.scaled, 
                                                                     #MATmin.scaled, 
                                                                     MAP.scaled,
                                                                     ppt.anom, 
                                                                     #tmin.anom, 
                                                                     tmax.anom, 
                                                                     slope.scaled, 
                                                                     aspect.scaled,
                                                                     #elev.scaled, 
                                                                     Ndep.scaled) %>% #,
                                           #physio.scaled) %>%,
                                           #physio.scaled) %>%
                                           # generate growth interactions
                                           mutate_at(.funs = list(growth.int = ~.*DIA_DIFF_scaled), 
                                                     .vars = vars(DIA_scaled:Ndep.scaled)) %>% 
                                           # generate diameter interactions
                                           mutate_at(.funs = list(DIA.int = ~.*DIA_scaled), 
                                                     .vars = vars(ba.scaled:Ndep.scaled)) %>%
                                           
                                           
                                           # # generate RD.scaled interactions
                                           # mutate_at(.funs = list(RD.scaled.int = ~.*RD.scaled), 
                                           #           .vars = vars(ba.scaled:physio.scaled))%>%
                                           
                                           # generate ba interactions
                                           mutate_at(.funs = list(ba.int = ~.*ba.scaled), 
                                                     .vars = vars(BAL.scaled:Ndep.scaled))%>%
                                           # generate BAL interactions
                                           mutate_at(.funs = list(BAL.int = ~.*BAL.scaled), 
                                                     .vars = vars(damage.scaled:Ndep.scaled))%>%
                                           
                                           # generate damage interactions
                                           mutate_at(.funs = list(damage.int = ~.*damage.scaled), 
                                                     .vars = vars(MATmax.scaled:Ndep.scaled)) %>%
                                           
                                           # generate MATmax interactions
                                           mutate_at(.funs = list(MATmax.scaled.int = ~.*MATmax.scaled), 
                                                     .vars = vars(MAP.scaled:Ndep.scaled))%>%
                                           # # generate MATmin.scaled interactions
                                           # mutate_at(.funs = list(MATmin.scaled.int = ~.*MATmin.scaled), 
                                           #           .vars = vars(MAP.scaled:physio.scaled))%>%
                                           
                                           # generate MAP.scaled interactions
                                           mutate_at(.funs = list(MAP.scaled.int = ~.*MAP.scaled), 
                                                     .vars = vars(ppt.anom:Ndep.scaled))%>%
                                           
                                           mutate_at(.funs = list(ppt.anom.int = ~.*ppt.anom), 
                                                     .vars = vars(tmax.anom:Ndep.scaled))%>%
                                           # # generate tmin.anom interactions
                                           # mutate_at(.funs = list(tmin.anom.int = ~.*tmin.anom), 
                                           #           .vars = vars(tmax.anom:Ndep.scaled))%>%
                                           
                                           # generate tmax.anom interactions
                                           mutate_at(.funs = list(tmax.anom.int = ~.*tmax.anom), 
                                                     .vars = vars(slope.scaled:Ndep.scaled)) %>%
                                           # generate slope.scaled interactions
                                           mutate_at(.funs = list(slope.int = ~.*slope.scaled), 
                                                     .vars = vars(aspect.scaled:Ndep.scaled))%>%
                                           # generate aspect.scaled interactions
                                           mutate_at(.funs = list(aspect.int = ~.*aspect.scaled), 
                                                     .vars = vars(Ndep.scaled))#%>%
                                         
                                         # # generate elev interactions
                                         # mutate_at(.funs = list(elev.int = ~.*elev.scaled), 
                                         #           .vars = vars(Ndep.scaled:physio.scaled)) %>%
                                         # generate Ndep interactions
                                         #mutate(Ndep.physio.int = physio.scaled*Ndep.scaled) #%>%
                                         #mutate_at(.funs = function(x)(x/10), .vars = vars( RD.scaled_growth.int:Ndep.physio.int))
                                         
                          ))
  
  
  
  
  
  
  model.name <- paste0("simple_logistic_SPCD_", SPCD.id, "remper_",remper.correction)
  
  mod.data <- mod.data.1
  mod.data.test <- mod.data.1.test
  mod.data$Nrep <- mod.data.test$N
  mod.data$xMrep <- mod.data.test$xM
  mod.data$ytest <- mod.data.test$y
  mod.data$Remperoos <- mod.data.test$Remper
  save(train.data, 
       test.data, 
       mod.data,
       mod.data.test, 
       model.name, 
       file = paste0("SPCD_standata_general_full_standardized_v3/SPCD_",SPCD.id,"remper_correction_",remper.correction,"model_1.Rdata"))
  
  mod.data <- mod.data.2
  mod.data.test <- mod.data.2.test
  mod.data$Nrep <- mod.data.test$N
  mod.data$xMrep <- mod.data.test$xM
  mod.data$ytest <- mod.data.test$y
  mod.data$Remperoos <- mod.data.test$Remper
  save(train.data, 
       test.data, 
       mod.data,
       mod.data.test, 
       model.name, 
       file = paste0("SPCD_standata_general_full_standardized_v3/SPCD_",SPCD.id,"remper_correction_",remper.correction,"model_2.Rdata"))
  
  mod.data <- mod.data.3
  mod.data.test <- mod.data.3.test
  mod.data$Nrep <- mod.data.test$N
  mod.data$xMrep <- mod.data.test$xM
  mod.data$ytest <- mod.data.test$y
  mod.data$Remperoos <- mod.data.test$Remper
  save(train.data, 
       test.data, 
       mod.data,
       mod.data.test, 
       model.name, 
       file = paste0("SPCD_standata_general_full_standardized_v3/SPCD_",SPCD.id,"remper_correction_",remper.correction,"model_3.Rdata"))
  
  mod.data <- mod.data.4
  mod.data.test <- mod.data.4.test
  mod.data$Nrep <- mod.data.test$N
  mod.data$xMrep <- mod.data.test$xM
  mod.data$ytest <- mod.data.test$y
  mod.data$Remperoos <- mod.data.test$Remper
  save(train.data, 
       test.data, 
       mod.data,
       mod.data.test, 
       model.name, 
       file = paste0("SPCD_standata_general_full_standardized_v3/SPCD_",SPCD.id,"remper_correction_",remper.correction,"model_4.Rdata"))
  
  mod.data <- mod.data.5
  mod.data.test <- mod.data.5.test
  mod.data$Nrep <- mod.data.test$N
  mod.data$xMrep <- mod.data.test$xM
  mod.data$ytest <- mod.data.test$y
  mod.data$Remperoos <- mod.data.test$Remper
  save(train.data, 
       test.data, 
       mod.data,
       mod.data.test, 
       model.name, 
       file = paste0("SPCD_standata_general_full_standardized_v3/SPCD_",SPCD.id,"remper_correction_",remper.correction,"model_5.Rdata"))
  
  mod.data <- mod.data.6
  mod.data.test <- mod.data.6.test
  mod.data$Nrep <- mod.data.test$N
  mod.data$xMrep <- mod.data.test$xM
  mod.data$ytest <- mod.data.test$y
  mod.data$Remperoos <- mod.data.test$Remper
  save(train.data, 
       test.data, 
       mod.data,
       mod.data.test, 
       model.name, 
       file = paste0("SPCD_standata_general_full_standardized_v3/SPCD_",SPCD.id,"remper_correction_",remper.correction,"model_6.Rdata"))
  
  mod.data <- mod.data.7
  mod.data.test <- mod.data.7.test
  mod.data$Nrep <- mod.data.test$N
  mod.data$xMrep <- mod.data.test$xM
  mod.data$ytest <- mod.data.test$y
  mod.data$Remperoos <- mod.data.test$Remper
  save(train.data, 
       test.data, 
       mod.data,
       mod.data.test, 
       model.name, 
       file = paste0("SPCD_standata_general_full_standardized_v3/SPCD_",SPCD.id,"remper_correction_",remper.correction,"model_7.Rdata"))
  
  mod.data <- mod.data.8
  mod.data.test <- mod.data.8.test
  mod.data$Nrep <- mod.data.test$N
  mod.data$xMrep <- mod.data.test$xM
  mod.data$ytest <- mod.data.test$y
  mod.data$Remperoos <- mod.data.test$Remper
  save(train.data, 
       test.data, 
       mod.data,
       mod.data.test, 
       model.name, 
       file = paste0("SPCD_standata_general_full_standardized_v3/SPCD_",SPCD.id,"remper_correction_",remper.correction,"model_8.Rdata"))
  
  mod.data <- mod.data.9
  mod.data.test <- mod.data.9.test
  mod.data$Nrep <- mod.data.test$N
  mod.data$xMrep <- mod.data.test$xM
  mod.data$ytest <- mod.data.test$y
  mod.data$Remperoos <- mod.data.test$Remper
  save(train.data, 
       test.data, 
       mod.data,
       mod.data.test, 
       model.name, 
       file = paste0("SPCD_standata_general_full_standardized_v3/SPCD_",SPCD.id,"remper_correction_",remper.correction,"model_9.Rdata"))
  
}
