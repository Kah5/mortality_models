# analyze_hierarchical_marginals.R
library(rstan)
library(dplyr)
library(tidyr)
library(ggplot2)
library(FIESTA)
library(posterior)

output.folder <- "/home/rstudio/"
#output.folder = "C:/Users/KellyHeilman/Box/01. kelly.heilman Workspace/mortality/Eastern-Mortality/mortality_models"

nspp <- data.frame(SPCD = c(316, 318, 833, 832, 261, 531, 802, 129, 762,  12, 541,  97, 621, 400, 371, 241, 375))
nspp$Species <- paste(FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD),]$GENUS, FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD),]$SPECIES)
nspp$COMMON <- paste(FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD),]$GENUS, FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD),]$SPECIES)


# Load input data ----
# mod.data.full contains all the infomration used to fit the model
#mod.data.full <- readRDS("SPCD_stanoutput_joint_v3/all_SPCD_model_6.RDS")


model.no <- 6
# set up species table
spp.table <- data.frame(SPCD = nspp[1:17,]$SPCD, 
                        SPP.id = 1:17, 
                        spp = 1:17,
                        Species = nspp[1:17,]$COMMON)
# but we also need the state id from the data, so add here from species data


cat(paste0("running model number ", model.no))

xM.list <-
  xMrep.list <-
  y.list <-
  y.test.list <-
  nSPP.list <-
  nSPP.rep.list <- 
  state.list<- state.rep.list <- 
  plot.list <- plot.rep.list <- remper.list <- remper.rep.list <- list()

for (i in 1:17) {
  # SPCD.id
  SPCD.id <- spp.table[i, ]$SPCD
  load(
    paste0(
      getwd(), "/SPCD_standata_general_full_standardized_v3/SPCD_",
      SPCD.id,
      "remper_correction_0.5model_",
      model.no,
      ".Rdata"
    )
  ) # load the species code data
  
  
  
  xM.list[[i]] <- mod.data$xM
  xMrep.list[[i]] <- mod.data$xMrep
  y.list[[i]] <- mod.data$y
  y.test.list[[i]] <- mod.data$ytest
  nSPP.list[[i]] <- rep(i, length(mod.data$y))
  nSPP.rep.list[[i]] <- rep(i, length(mod.data$ytest))
  
  remper.list[[i]] <- train.data$remper
  remper.rep.list[[i]] <- test.data$remper
  
  # get state information:
  state.list[[i]] <- train.data$state
  state.rep.list[[i]] <- test.data$state
  
  # get plot information:
  plot.list[[i]] <- train.data
  plot.rep.list[[i]] <- test.data
}




mod.data.full <-  list(
  xM = do.call(rbind, xM.list),
  xMrep = do.call(rbind, xMrep.list),
  y = unlist(y.list),
  ytest = unlist(y.test.list),
  SPP = unlist(nSPP.list),
  SPPrep = unlist(nSPP.rep.list),
  Remper = unlist(remper.list),
  Remperoos = unlist(remper.rep.list)
)

plot.data.train <- do.call(rbind, plot.list)
plot.data.test <- do.call(rbind, plot.rep.list)
plot.data.train$SPP <- mod.data.full$SPP


mod.data.full$K <- ncol(mod.data.full$xM)
mod.data.full$N <- length(mod.data.full$y)

mod.data.full$Nspp <- 17
mod.data.full$Nrep <- length(mod.data.full$ytest)

# save:
saveRDS (
  mod.data.full,
  paste0(
    output.folder,
    "/SPCD_stanoutput_joint_v3/all_SPCD_model_",
    model.no,
    ".RDS"
  )
)


# save:
saveRDS (
  plot.data.train,
  paste0(
    output.folder,
    "/SPCD_stanoutput_joint_v3/train_plot_data_SPCD_model_",
    model.no,
    ".RDS"
  )
)

# save:
saveRDS (
  plot.data.test,
  paste0(
    output.folder,
    "/SPCD_stanoutput_joint_v3/test_plot_data_SPCD_model_",
    model.no,
    ".RDS"
  )
)



saveRDS(spp.table,
        paste0(output.folder, "/SPCD_stanoutput_joint_v3/spp.table.rds"))
spp.table <-
  readRDS(paste0(output.folder, "/SPCD_stanoutput_joint_v3/spp.table.rds"))


# Load posterior estimates ----

# read in alphas by species
alpha_species <- readRDS(paste0(output.folder,"/SPCD_stanoutput_joint_v3/samples/alpha.spp_model_6_1000samples.rds"))

# get random sample of each parameter:
rand.draws <- sample(1:max(alpha_species$.iteration), 100)
alpha_species <- alpha_species %>% posterior::subset_draws(iteration = rand.draws)


# read in betas by species
beta_species <- readRDS(paste0(output.folder,"/SPCD_stanoutput_joint_v3/samples/u_betas_model_6_1000samples.rds"))%>% 
  posterior::subset_draws(iteration = rand.draws)

# read in posterior probability by individual (mMhat)
p_surv_train <- readRDS(paste0(output.folder,"/SPCD_stanoutput_joint_v3/samples/mMhat_model_6_1000samples.rds"))%>% 
  posterior::subset_draws(iteration = rand.draws)


# Set up groups and marginal covs ----
X <- data.frame(mod.data.full$xM)                # [N x K]
N <- nrow(X)
species <- mod.data.full$SPP     # [N]
state <- plot.data.train$state            # [N]
tree_ids <- 1:nrow(X)

# Predictors
predictor_names <- colnames(X)
X_matrix <- as.matrix(X)

X$tree =  tree_ids# set up tree id

N <- nrow(X) # number of individual trees

S <- dim(beta_species)[1] # number of iterations
K <- ncol(X_matrix) # number of predictors
n_species <- length(unique(species))

# --- PREP POSTERIOR COMPONENTS ---

# Species-level intercepts [S × n_species]
colnames(alpha_species)
alpha_species <- alpha_species %>% select(-.chain, -.iteration, -.draw) 

# Random slopes [S × n_species × K]
beta_species <- beta_species %>% select(-.chain, -.iteration, -.draw) 
colnames(beta_species)
#dimnames(beta_species) <- list(NULL, paste0("sp", 1:n_species), predictor_names)

# Fixed effects (unused here, but accessible)
#beta_fixed <- post$beta

# --- COMPUTE TREE-LEVEL POSTERIOR SURVIVAL PROBABILITIES----
# 
# alpha_i <- matrix(NA, N, S)
# for (s in 1:S) {
#   alpha_i[, s] <- alpha_species[s, species]
# }
# 
# eta <- matrix(0, N, S)
# for (s in 1:S) {
#   for (i in 1:N) {
#     species_i <- species[i]
#     slopes_i <- beta_species[s, species_i, ]
#     eta[i, s] <- alpha_i[i, s] + sum(X_matrix[i, ] * slopes_i)
#   }
# }
# p_surv <- plogis(eta)  # [N × S]

# just use p_surv_train in from posterior estimates:
colnames(p_surv_train)

# --- LONG FORMAT PREP ---
unique(plot.data.train$state)
unique(plot.data.train$SPP)

st <- 33

df_meta <- plot.data.train
df_meta$tree <- 1:N
df_meta$SPP <- mod.data.full$SPP

# p_long <- as.data.frame(p_surv_train) %>% select(-.chain, -.iteration, -.draw) 
# colnames(p_long) <- paste0("tree_", 1:N)
# p_long$draw <- 1:S

# filter only trees in a given state:
#st <- 36
st.names.by.num <- plot.data.train %>% group_by(state) %>% summarise(n())%>% arrange(`n()`)
get_statewide_marginal_variances <- function(st) {
  
  # create the state directory
  dir.create(paste0("SPCD_stanoutput_joint_v3/predicted_mort/state_",st))
  
  # get the metadata for the state only
  df_meta_st <- df_meta %>% filter(state %in% st)
  
  trees.per.spp <- df_meta_st %>% filter(state %in% st) %>%
    group_by(Species, spp, SPCD)%>%
    summarise(ntree_spp = n())
  
  # if there is less than 1 tree for that species/state, remove from calculation:
  
  
  # # get all of the draws for p(survival of) each tree 
  # p_long_st <-
  #   p_long %>% select(c(paste0("tree_", df_meta_st$tree), "draw"))
  
  # filter out the covariate data for each tree in the state
  X_st <- X %>% filter(tree %in% df_meta_st$tree)
  
  # get the species ID indexing
  species.st <- mod.data.full$SPP[df_meta$state %in% st]
  
  
  # set up a species indexing dataframe with the tree id, and species information
  Species.indexing <- data.frame(tree = X$tree,
                                 spp = mod.data.full$SPP) %>%
    left_join(., spp.table)
  
  # # reformat p_long_st (actually a wide format) to be a long format dataframe
  # p_long_long <- pivot_longer(
  #   p_long_st,
  #   cols = starts_with("tree_"),
  #   names_to = "tree",
  #   names_prefix = "tree_",
  #   values_to = "p_surv"
  # ) %>%
  #   mutate(tree = as.integer(tree)) #%>%
  #left_join(., Species.indexing)
  
  # join up the the long format data up to the state information
  # runing intoe some memory sisue here
  # df_k_st <- df_meta_st %>%
  #   left_join(p_long_long)# %>%
  # # # pivot_longer(cols = starts_with("draw_"), names_to = "draw", values_to = "p_surv") %>%
  # # mutate(draw = as.integer(gsub("draw_", "", draw))) %>%
  # bind_cols(X[rep(1:N, S), ])
  
  # for each tree get the beta samples
  
  
  
  # get the marginal effects including random slopes ----
  
  #marginal_all <- list()
  
  
  # get_marginal_pred_val <- function (k){
  #   #for (k in 1:length(predictor_names)) {
  #   cat(paste(
  #     "getting marginal response of",
  #     predictor_names[k],
  #     "in state",
  #     st, "\n"
  #   ))
  #   slope_k <- matrix(NA, nrow = length(species.st), ncol = S)
  #   
  #   
  #   species.betas <- paste0("u_beta[", species.st, ",", k, "]")
  #   
  #   df_k <- df_meta_st %>% #df_k_st %>%
  #     
  #     # get the name of beta of interest for each species
  #     mutate(beta_k_name = paste0("u_beta[", SPP, ",", k, "]")) #%>% #select(beta_k_name)
  #   # get the beta sample of interest
  #   
  #   # for each tree get the beta samples
  #   #marginal_tree_list <- list()
  #   
  #   
  #   
  #   species.level.marg.effect <- function(sp.id){
  #     #tree.no <- tre #unique(df_k$tree)[tre]
  #     df_tree <- df_k %>% filter(SPP %in% sp.id)
  #     
  #     beta_tree <- beta_species[, unique(df_tree$beta_k_name)] %>% 
  #       mutate(draw = 1:nrow(beta_species[,1]), 
  #              beta_k_name = unique(df_tree$beta_k_name)) %>%
  #       rename("beta_sample" = unique(df_tree$beta_k_name))
  #     
  #     marg_tree <- df_tree %>% select(tree, draw, state, county, SPCD, p_surv) %>%
  #       left_join(., beta_tree) %>% #, relationship = "many-to-many") %>%
  #       ungroup() %>%
  #       mutate(marginal_effect = p_surv * (1 - p_surv) * beta_sample)
  #     
  #     marg_tree
  #   }
  #   #species.level.marg.effect(sp.id = 2)
  #   
  #   # for (tre in 1:length(unique(df_k$tree))) {
  #   marg_tree_st <- lapply(unique(species.st), species.level.marg.effect) %>%
  #     do.call(rbind,.)
  #   
  #   
  #   marg_tree_st$species <-
  #     FIESTA::ref_species[match(marg_tree_st$SPCD, FIESTA::ref_species$SPCD), ]$COMMON
  #   
  #   # save the samples
  #   saveRDS(
  #     marg_tree_st,
  #     paste0(
  #       "SPCD_stanoutput_joint_v3/predicted_mort/marginal_state/marginal_trees_",
  #       unique(predictor_names)[k],
  #       "_state_",
  #       st ,
  #       ".RDS"
  #     )
  #   )
  #   
  #   df_k_summary <- marg_tree_st %>%
  #     group_by(species, state, draw) %>%
  #     summarize(marginal = mean(marginal_effect),
  #               .groups = "drop") %>%
  #     group_by(species, state) %>%
  #     summarize(
  #       mean = mean(marginal),
  #       lower = quantile(marginal, 0.05),
  #       upper = quantile(marginal, 0.95),
  #       .groups = "drop"
  #     ) %>%
  #     mutate(predictor = unique(predictor_names)[k])
  #   
  #   df_k_summary
  # }
  # marginal_all <- lapply(1:length(predictor_names), get_marginal_pred_val )
  # marginal_summary_all <- do.call(rbind, marginal_all)
  # #marginal_summary_all <- bind_rows(marginal_all)
  # 
  # write.csv(
  #   marginal_summary_all,
  #   paste0(
  #     "SPCD_stanoutput_joint_v3/predicted_mort/marginal_effects_by_species_state_",
  #     st,
  #     ".csv"
  #   ),
  #   row.names = FALSE
  # )
  
  # rm(marginal_all)
  # # --- VARIANCE PARTITIONING BY PREDICTOR ---
  # 
  # # Total group-level variance
  # group_mean_total <- df_k_st %>%
  #   group_by(SPCD, state, tree, draw) %>%
  #   summarize(mean_p = mean(p_surv), .groups = "drop") %>%
  #   group_by(SPCD, state, draw) %>%
  #   summarize(var_total = var(mean_p, na.rm = TRUE),
  #             .groups = "drop")
  # 
  # ggplot(group_mean_total, aes(var_total)) + geom_histogram() +
  #   facet_wrap( ~ SPCD) +
  #   xlab("Total Variance in tree-level predicted probability of surival")
  # 
  # saveRDS(
  #   df_k,
  #   paste0(
  #     "SPCD_stanoutput_joint_v3/predicted_mort/df_total_preds_spp_state_",
  #     st,
  #     ".RDS"
  #   )
  # )
  
  
  # Predictor-wise variance explained for this state
  cat("computing variance explained by each predictor")
  compute_var_explained <- function(k) {
    # if there is more than one tree per species
    spp.in.state <- data.frame(spp = species.st)%>% left_join(.,spp.table)%>%
      group_by(spp)%>% summarise(n()) %>% 
      filter(`n()` > 2)
    
    # calculate partial logits by species
    logits.partial <- lapply(spp.in.state$spp, function(sp_i) {
      
      
      #sp_i <- i #species.st[i]
      #df_tree <- df_k %>% filter(SPP %in% sp.id)
      species.index <- df_meta_st$SPP %in% sp_i
      slope_i <- beta_species[, paste0("u_beta[", sp_i, ",", k, "]")]
      
      slope.alpha = data.frame(alpha = alpha_species[paste0("alpha_SPP[", sp_i, "]")], 
                               slope_i = as.vector(slope_i[,1]), 
                               draw = 1:nrow(slope_i), 
                               predictor = predictor_names[k])
      colnames(slope.alpha) <- c("alpha_spp", "slope_spp", "draw", "predictor")
      
      df.out <- data.frame(
        state = st,
        tree = X_st[species.index, ]$tree,
        SPP.id = sp_i,
        SPCD = spp.table[sp_i, ]$SPCD,
        predictor = predictor_names[k],
        predictor_val = X_st[species.index, predictor_names[k]]) %>% 
        # link up alphas and betas for each species
        left_join(., slope.alpha, relationship = "many-to-many") %>% 
        # calculate partial logit
        mutate(logit.part = alpha_spp + predictor_val * slope_spp)%>%
        mutate(p_partial = plogis(logit.part))
      
      saveRDS(df.out ,paste0(output.folder,
                             "/SPCD_stanoutput_joint_v3/predicted_mort/state_",st,"/df_SPCD_",spp.table[sp_i, ]$SPCD,"_variance_state_",
                             st,"_predictor_", predictor_names[k],
                             ".RDS"
      ))
      
      df.out #%>% left_join(.,spp.table)
    })
    
    
    partial_logit_df <- do.call(rbind, logits.partial)
    
    # total_mean_by_k <- partial_logit_df %>%
    #   group_by(tree, SPCD, state, draw) %>%
    #   summarize(mean_p_k = mean(p_partial), .groups = "drop")
    
    
    mean_by_group_k <- partial_logit_df %>%
      group_by(tree, SPCD, state, draw) %>%
      summarize(mean_p_k = mean(p_partial, na.rm = TRUE), .groups = "drop") %>%
      group_by(draw, SPCD) %>%
      summarize(var_k = var(mean_p_k, na.rm =TRUE), .groups = "drop") %>%
      mutate(predictor = predictor_names[k])
    
    return(mean_by_group_k)
  }
  
  var_parts <-
    lapply(1:length(predictor_names), compute_var_explained) %>%
    do.call(rbind,.)
  
  # save the full variance explained table by state
  saveRDS(
    var_parts,
    paste0(
      "SPCD_stanoutput_joint_v3/predicted_mort/state_",st,"/variance_explained_state_",
      st,
      ".RDS"
    )
  )
  
  
  tot.var <- var_parts %>% group_by(SPCD, draw) %>%
    summarise(var_total = sum(var_k, na.rm =TRUE))
  
  var_summary <-
    left_join(var_parts, tot.var, by = c("SPCD", "draw")) %>%
    
    #group_mean_total, by = c("draw", "SPCD")) %>%
    mutate(rel_var = (var_k / var_total)) %>%
    group_by(predictor, SPCD) %>%
    summarize(
      mean = mean(rel_var),
      lower = quantile(rel_var, 0.05),
      upper = quantile(rel_var, 0.95)
    )
  var_summary$predictor <-
    factor(var_summary$predictor, levels = predictor_names)
  
  # ggplot(data = var_summary, aes(x = predictor, y = mean)) + geom_point() +
  #   geom_errorbar(aes(x = predictor, ymin = lower, ymax = upper)) +
  #   theme(axis.text = element_text(angle = 45, hjust = 1),
  #         panel.grid = element_blank()) + facet_wrap( ~ SPCD, scales = "free_y")
  # 
  # ggplot(data = var_summary %>% filter(!predictor %in% "DIA_DIFF_scaled"),
  #        aes(x = predictor, y = mean)) + geom_point() +
  #   geom_errorbar(aes(x = predictor, ymin = lower, ymax = upper)) +
  #   theme(axis.text = element_text(angle = 45, hjust = 1),
  #         panel.grid = element_blank()) +
  #   ylab ("relative variance in annaul prob(survival) explained by predictor") +
  #   facet_wrap( ~ SPCD, scales = "free_y")
  # 
  # var_summary %>% filter(SPCD %in% 129)
  # 
  # ggplot(
  #   data = var_summary %>% filter(!predictor %in% "DIA_DIFF_scaled" &
  #                                   SPCD %in% 129),
  #   aes(x = predictor, y = mean)
  # ) + geom_point() +
  #   geom_errorbar(aes(x = predictor, ymin = lower, ymax = upper)) +
  #   theme(axis.text = element_text(angle = 45, hjust = 1),
  #         panel.grid = element_blank()) +
  #   ylab ("relative variance in annaul prob(survival) explained by predictor") +
  #   facet_wrap( ~ SPCD, scales = "free_y")
  # 
  main.preds <- unique(predictor_names)[1:12]
  inter.preds <- unique(predictor_names)[13:33]
  var_summary$Species <-
    FIESTA::ref_species[match(var_summary$SPCD, FIESTA::ref_species$SPCD), ]$COMMON
  
  write.csv(
    var_summary,
    paste0(
      "SPCD_stanoutput_joint_v3/predicted_mort/state_",st,"/variance_partitioning_summary_by_predictor_state_",
      st,
      ".csv"
    ),
    row.names = FALSE
  )
  
  rm(X_st, df_k, df_k, df_k, df_meta_st, group_mean_total, marg_tree_st, 
     marginal_summary_all, mean_by_group_k, 
     p_long_long, p_long_st, partial_logit_df)
  
  var_parts <- readRDS(paste0("SPCD_stanoutput_joint_v3/predicted_mort/state_",st,"/variance_explained_state_",st,".RDS"))
  
  covariate.names <- read.csv("model_covariate_types.csv")
  
  
  # species - wide summary total variance
  tot.var <- var_parts %>%  group_by(SPCD, draw) %>%
    summarise(var_total = sum(var_k))
  
  # species -wide summary total variance but without dia_diff
  tot.var.no.dia.diff <- var_parts %>% filter(!predictor %in% "DIA_DIFF_scaled")%>% group_by(SPCD, draw) %>%
    summarise(var_total = sum(var_k))
  
  # state - wide summary total variance
  tot.var.st <- tot.var  %>% group_by(draw) %>% 
    summarise(var_total_st = sum(var_total))
  
  var_summary <-
    left_join(var_parts, tot.var, by = c("SPCD", "draw")) %>%
    
    #group_mean_total, by = c("draw", "SPCD")) %>%
    mutate(rel_var = (var_k / var_total)) %>%
    group_by(predictor, SPCD) %>%
    summarize(
      mean = mean(rel_var),
      lower = quantile(rel_var, 0.05),
      upper = quantile(rel_var, 0.95)
    ) %>% 
    rename("Covariate" = "predictor") %>% left_join(covariate.names)
  
  var_summary$Predictor <-
    factor(var_summary$Predictor, levels = covariate.names$Predictor)
  var_summary$COMMON <- FIESTA::ref_species[match(var_summary$SPCD, FIESTA::ref_species$SPCD),]$COMMON
  
  
  # ggplot(data = var_summary, aes(x = Predictor, y = mean)) + geom_point() +
  #   geom_errorbar(aes(x = Predictor, ymin = lower, ymax = upper)) +
  #   theme(axis.text = element_text(angle = 45, hjust = 1),
  #         panel.grid = element_blank()) + facet_wrap( ~ SPCD, scales = "free_y")
  # 
  # ggplot(data = var_summary %>% filter(!Predictor %in% "DIA_DIFF_scaled"),
  #        aes(x = Predictor, y = mean)) + geom_point() +
  #   geom_errorbar(aes(x = Predictor, ymin = lower, ymax = upper)) +
  #   theme(axis.text = element_text(angle = 45, hjust = 1),
  #         panel.grid = element_blank()) +
  #   ylab ("relative variance in annaul prob(survival) explained by predictor") +
  #   facet_wrap( ~ SPCD, scales = "free_y")
  
  
  
  
  
  # ggplot(data = var_summary)+
  #   geom_bar(aes(x = COMMON, y = mean, fill = predictor.class), stat = "identity", position = "stack")+
  #   theme_minimal()+ theme(axis.text.x = element_text(angle = 60, hjust = 1), 
  #                          panel.background = element_blank())+
  #   ylab("Average proportion of across-tree \n p(survival) variance explained")
  # 
  # ggplot(data = var_summary)+
  #   geom_bar(aes(x = COMMON, y = mean, fill = Predictor), stat = "identity", position = "stack")+
  #   theme_minimal()+ theme(axis.text.x = element_text(angle = 60, hjust = 1), 
  #                          panel.background = element_blank())+
  #   ylab("Average proportion of across-tree \n p(survival) variance explained")+
  #   xlab("Species")
  
  
  color.pred.class <- c(
    "Site x G & S"= "#1b9e77",
    "Competition x G & S"="#d95f02" ,
    "Climate x G & S"="#7570b3" ,
    "Climate"= "#e7298a" ,
    "Growth & Size" = "#66a61e" ,
    "Site Conditions" = "#e6ab02",
    "Competition" = "#a6761d"
  )
  
  var_summary$predictor.class <- factor(var_summary$predictor.class, 
                                        levels = c(
                                          "Growth & Size", 
                                          "Competition",
                                          "Climate",
                                          "Site Conditions",
                                          "Competition x G & S",
                                          "Climate x G & S", 
                                          "Site x G & S"
                                          
                                          
                                        ))
  ggplot(data = var_summary )+
    geom_bar(aes(x = COMMON, y = mean, fill = predictor.class), stat = "identity", position = "stack")+
    theme_minimal()+ theme(axis.text.x = element_text(angle = 60, hjust = 1), 
                           panel.background = element_blank())+
    ylab("Average proportion of across-tree \n p(survival) variance explained")+
    xlab("Species")+scale_fill_manual(values = color.pred.class, name = "")
  ggsave(height = 4, width = 6, dpi = 350, paste0("SPCD_stanoutput_joint_v3/predicted_mort/state_",st,"/Prop_variance_state_",st,".png"))
  
  # ggplot(data = var_summary )+
  #   geom_col(aes(x = "", y = mean, fill = predictor.class), color = "black",stat = "identity", position = "stack")+
  #   theme_minimal()+ theme(axis.text.x = element_text(angle = 60, hjust = 1), 
  #                          panel.background = element_blank())+
  #   ylab("Average proportion of across-tree \n p(survival) variance explained")+
  #   xlab("Species")+scale_fill_manual(values = color.pred.class, name = "") +
  #   coord_polar(theta = "y")+facet_wrap(~COMMON)
  # 
  # ggplot(data = var_summary %>% filter(!Covariate %in% "DIA_DIFF_scaled" ))+
  #   geom_col(aes(x = "", y = mean, fill = predictor.class), color = "black")+ #,stat = "identity", position = "stack")+
  #   theme_minimal()+ theme(axis.text.x = element_text(angle = 60, hjust = 1), 
  #                          panel.background = element_blank())+
  #   ylab("Average proportion of across-tree \n p(survival) variance explained")+
  #   xlab("Species")+scale_fill_manual(values = color.pred.class, name = "") +
  #   coord_polar(theta = "y")+facet_wrap(~COMMON)
  # 
  
  # non-DIA_diff variables
  ggplot(data = var_summary %>% filter(!Covariate %in% "DIA_DIFF_scaled"))+
    geom_bar(aes(x = COMMON, y = mean, fill = predictor.class), stat = "identity", position = "stack")+
    theme_minimal()+ theme(axis.text.x = element_text(angle = 60, hjust = 1), 
                           panel.background = element_blank())+
    ylab("Average proportion of across-tree \n p(survival) variance explained")+
    xlab("Species")+scale_fill_manual(values = color.pred.class, name = "")
  ggsave(height = 4, width = 6, dpi = 350, paste0("SPCD_stanoutput_joint_v3/predicted_mort/state_",st,"/Prop_variance_state_",st,"_no_diadiff.png"))
  
  # ggplot(data = var_summary %>% filter(!Covariate %in% "DIA_DIFF_scaled"))+
  #   geom_bar(aes(x = COMMON, y = mean, fill = Predictor), stat = "identity", position = "stack")+
  #   theme_minimal()+ theme(axis.text.x = element_text(angle = 60, hjust = 1), 
  #                          panel.background = element_blank())+
  #   ylab("Average proportion of across-tree \n p(survival) variance explained")+
  #   xlab("Species")#+scale_fill_manual(values = color.pred.class)
  # 
  
  # do the summary for all non-DIA_diff variables:
  var_summary_no_dia_diff <-
    left_join(var_parts %>% filter(!predictor %in% "DIA_DIFF_scaled"), tot.var.no.dia.diff, by = c("SPCD", "draw")) %>%
    
    #group_mean_total, by = c("draw", "SPCD")) %>%
    mutate(rel_var = (var_k / var_total)) %>%
    group_by(predictor, SPCD) %>%
    summarize(
      mean = mean(rel_var),
      lower = quantile(rel_var, 0.05),
      upper = quantile(rel_var, 0.95)
    ) %>% 
    rename("Covariate" = "predictor") %>% left_join(covariate.names)
  
  var_summary_no_dia_diff$Predictor <-
    factor(var_summary_no_dia_diff$Predictor, levels = covariate.names$Predictor)
  var_summary_no_dia_diff$COMMON <- FIESTA::ref_species[match(var_summary_no_dia_diff$SPCD, FIESTA::ref_species$SPCD),]$COMMON
  
  ggplot(data = var_summary_no_dia_diff )+
    geom_col(aes(x = "", y = mean, fill = predictor.class), color = "black") + #,stat = "identity", position = "stack")+
    theme_minimal()+ theme(axis.text.x = element_text(angle = 60, hjust = 1), 
                           panel.background = element_blank())+
    ylab("Average proportion of across-tree \n p(survival) variance explained")+
    xlab("")+scale_fill_manual(values = color.pred.class, name = "") +
    coord_polar(theta = "y")+facet_wrap(~COMMON)
  ggsave(height = 6, width = 6, dpi = 350, paste0("SPCD_stanoutput_joint_v3/predicted_mort/state_",st,"/Prop_variance_state_",st,"_no_dia_diff_pie_spp.png"))
  
  cat(paste("\n copying state", st, " outputs to cyverse output folder"))
  # copy to the data-store output
  system(paste(
    "cp -r",
    paste0("SPCD_stanoutput_joint_v3/predicted_mort/state_",st,"/"),
    "data-store/data/output/"
  ))
  
  
  
}

get_statewide_marginal_variances(st = st.names.by.num$state[1])#state_34 
get_statewide_marginal_variances(st = 42)

lapply(unique(plot.data.train$state), FUN = function(x){get_statewide_marginal_variances(st = x)} )



state.files <- list.files(paste0(output.dir, "SPCD_stanoutput_joint_v3/predicted_mort/"), pattern = "*.csv")
# Split strings and keep the last part
last_parts <- sapply(strsplit(state.files, "_"), function(x) tail(x, 1))
state.numbers <- sapply(strsplit(last_parts, ".csv"), function(x) tail(x, 1))

######################################################################################
# read in the total variance outputs for each state and species to get a regional summary for the species
var.all.files <- list.files(paste0(output.dir, "SPCD_stanoutput_joint_v3/predicted_mort/"), pattern = "*.RDS", full.names = TRUE)
var.all.list <- list()
for(i in 1:length(state.files)){
  
  var.all.list[[i]] <- readRDS(var.all.files[i])%>% mutate(state = state.numbers[i])
}
var.all.df <- do.call(rbind, var.all.list)

# get variance summaries by species:

tot.var <- var.all.df %>% group_by(SPCD, draw, state) %>%
  summarise(var_total = sum(var_k, na.rm =TRUE))

var_summary <-
  left_join(var.all.df, tot.var, by = c("SPCD", "draw", "state")) %>%
  
  #group_mean_total, by = c("draw", "SPCD")) %>%
  mutate(rel_var = (var_k / var_total)) %>%
  group_by(predictor, SPCD) %>%
  summarize(
    mean = mean(rel_var),
    lower = quantile(rel_var, 0.05),
    upper = quantile(rel_var, 0.95)
  )


covariate.names <- read.csv(paste0(output.dir, "data/covariate_names/model_covariate_types_v2.csv"))
state.dataframe <- data.frame(state = as.character(c(9, 23, 24, 33, 34, 36, 39, 42, 50, 54)), 
                              STNAME = c("Connecticut", "Maine", "Maryland", "New Hampshire", "New Jersey", "New York", 
                                         "Ohio", "Pennsylvania", "Vermont", "West Virginia"))

var_summary <- var_summary %>% 
  rename("Covariate" = "predictor") %>% left_join(., covariate.names)

var_summary$Predictor <-
  factor(var_summary$Predictor, levels = covariate.names$Predictor)
var_summary$COMMON <- FIESTA::ref_species[match(var_summary$SPCD, FIESTA::ref_species$SPCD),]$COMMON

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
  # main effects 
  "Size" = "darkgreen",
  "Change in Size" = "#66a61e",
  "Climate"= "#e7298a" ,
  "Site Conditions" = "#e6ab02",
  "Competition" = "#a6761d",
  "N deposition" = "red4", 
  "% Damage" = "purple", 
  "Growth x Size" =  "#1b9e77",
  
  "Site x G & S"= "yellow3",
  "Competition x G & S"="#d95f02" ,
  
  "Ndep x G & S" = "red",
  "Climate x G & S"="magenta" ,
  "Damage x G & S" = "#7570b3" 
 
  
  
)

var.summary.region <- var_summary
var.summary.region$COMMON <- factor(var.summary.region$COMMON, 
                                    levels = rev(c("balsam fir", "red spruce", "northern white-cedar", 
                                               "eastern hemlock", "American beech", 
                                               "black oak", "chestnut oak", "northern red oak", "white oak", "yellow birch", "paper birch", 
                                               "hickory spp.", "eastern white pine", "red maple", "sugar maple", 
                                               "black cherry", "white ash", "yellow-poplar")))
var.summary.region$predictor.class <- factor(var.summary.region$predictor.class, 
                                      levels = c(
                                        
                                        "Change in Size",
                                        "Size", 
                                        "Growth x Size", 
                                        "Competition",
                                        "Climate",
                                        "Site Conditions",
                                        "Competition x G & S",
                                        "Climate x G & S", 
                                        "Site x G & S"
                                        
                                        
                                      ))
regional.var.plt <- ggplot(data = var.summary.region )+
  geom_bar(aes(x = COMMON, y = mean, fill = predictor.class), stat = "identity", position = "stack")+
  theme_minimal(base_size = 16)+ theme(axis.text.x = element_text(angle = 60, hjust = 1), 
                    panel.background = element_blank())+
  ylab("Across-tree variance in  \n p(survival) explained")+
  scale_fill_manual(values = color.pred.class, name = "", 
                    guide = guide_legend(direction = "horizontal", 
                                         ncol = 9,nrow = 1,reverse = TRUE,
                                         label.position="top", label.hjust = 0.5, 
                                         label.vjust = 0.5,
                                         label.theme = element_text(angle = 90)))+
  coord_flip()+
  xlab("")+
  theme(panel.grid = element_blank(), 
        legend.position = "top")

ggsave(plot = regional.var.plt, height = 8, width = 6, units = "in", dpi = 300, 
       paste0(output.dir, "SPCD_stanoutput_joint_v3/predicted_mort/regional_species_var_partitioning.png"))

ggplot(data = var.summary.region)+
  geom_bar(aes(x = COMMON, y = mean, fill = predictor.class2), stat = "identity", position = "stack", color = "grey")+
  theme_minimal(base_size = 16)+ theme(axis.text.x = element_text(angle = 60, hjust = 1), 
                                       panel.background = element_blank())+
  ylab("Across-tree variance in  \n p(survival) explained")+
  scale_fill_manual(values = color.pred.class.2, name = "",
                    guide = guide_legend(direction = "horizontal",
                                         ncol = 14,nrow = 1,reverse = TRUE,
                                         label.position="top", label.hjust = 0.5,
                                         label.vjust = 0.5,
                                         label.theme = element_text(angle = 90)))+
  coord_flip()+
  xlab("")+
  theme(panel.grid = element_blank(), 
        legend.position = "top")

#########################################################################################
# read in the pre-summarised outputs for each state -----
var.part.files <- list.files(paste0(output.dir, "SPCD_stanoutput_joint_v3/predicted_mort/"), pattern = "*.csv", full.names = TRUE)
var.part.list <- list()
for(i in 1:length(state.files)){
 var.part.list[[i]] <- read.csv(var.part.files[i])%>% mutate(state = state.numbers[i])
}
var.part.df <- do.call(rbind, var.part.list)

covariate.names <- read.csv(paste0(output.dir, "data/covariate_names/model_covariate_types_v2.csv"))
state.dataframe <- data.frame(state = as.character(c(9, 23, 24, 33, 34, 36, 39, 42, 50, 54)), 
                              STNAME = c("Connecticut", "Maine", "Maryland", "New Hampshire", "New Jersey", "New York", 
                                         "Ohio", "Pennsylvania", "Vermont", "West Virginia"))

var_summary <- var.part.df %>% 
  rename("Covariate" = "predictor") %>% left_join(., covariate.names)%>%
  left_join(., state.dataframe)

var_summary$Predictor <-
  factor(var_summary$Predictor, levels = covariate.names$Predictor)
var_summary$COMMON <- FIESTA::ref_species[match(var_summary$SPCD, FIESTA::ref_species$SPCD),]$COMMON




var_summary$predictor.class <- factor(var_summary$predictor.class, 
                                      levels = c(
                                       
                                        "Change in Size",
                                        "Size", 
                                        "Growth x Size", 
                                        "Competition",
                                        "Climate",
                                        "Site Conditions",
                                        "Competition x G & S",
                                        "Climate x G & S", 
                                        "Site x G & S"
                                        
                                        
                                      ))


ggplot(data = var_summary )+
  geom_bar(aes(x = COMMON, y = mean, fill = predictor.class), stat = "identity", position = "stack")+
  theme_minimal()+ theme(axis.text.x = element_text(angle = 60, hjust = 1), 
                         panel.background = element_blank())+
  ylab("Average proportion of across-tree \n p(survival) variance explained")+
  xlab("Species")+scale_fill_manual(values = color.pred.class, name = "")+facet_wrap(~STNAME)




# balsam fir and red spruce
boreal.var.plt <- ggplot(data = var_summary %>% filter(COMMON %in% c("balsam fir", "red spruce", "northern white-cedar", "paper birch")))+
  geom_bar(aes(x = STNAME, y = mean, fill = predictor.class), stat = "identity", position = "stack")+
  theme_bw()+ theme(axis.text.x = element_text(angle = 60, hjust = 1), 
                         panel.background = element_blank())+
  ylab("Across-tree variance in  \n p(survival) explained")+
  scale_fill_manual(values = color.pred.class, name = "")+
coord_flip()+facet_wrap(~COMMON, scales = "free_y", ncol = 1)+
  xlab("")
ggsave(plot = boreal.var.plt, height = 6, width = 5, units = "in", dpi = 300, 
       paste0(output.dir, "SPCD_stanoutput_joint_v3/predicted_mort/subboreal_species_var_partitioning.png"))


# white pine, sugar maple, red maple
nmix.var.plt <- ggplot(data = var_summary %>% filter(COMMON %in% c("eastern white pine", "sugar maple", "red maple")))+
  geom_bar(aes(x = STNAME, y = mean, fill = predictor.class), stat = "identity", position = "stack")+
  theme_bw()+ theme(axis.text.x = element_text(angle = 60, hjust = 1), 
                    panel.background = element_blank())+
  ylab("Across-tree variance in  \n p(survival) explained")+
  scale_fill_manual(values = color.pred.class, name = "")+
  coord_flip()+facet_wrap(~COMMON, scales = "free_y", ncol = 1)+
  xlab("")
ggsave(height = 6, width = 5, units = "in", dpi = 300, 
       paste0(output.dir, "SPCD_stanoutput_joint_v3/predicted_mort/North_mixed_species_var_partitioning.png"))


# Beech, hemlock, YB
lsuccess.var.plt <- ggplot(data = var_summary %>% filter(COMMON %in% c("American beech", "eastern hemlock", "yellow birch")))+
  geom_bar(aes(x = STNAME, y = mean, fill = predictor.class), stat = "identity", position = "stack")+
  theme_bw()+ theme(axis.text.x = element_text(angle = 60, hjust = 1), 
                    panel.background = element_blank())+
  ylab("Across-tree variance in  \n p(survival) explained")+
  scale_fill_manual(values = color.pred.class, name = "")+
  coord_flip()+facet_wrap(~COMMON, scales = "free_y", ncol = 1)+
  xlab("")
ggsave(plot = lsuccess.var.plt, height = 6, width = 6, units = "in", dpi = 300, 
       paste0(output.dir, "SPCD_stanoutput_joint_v3/predicted_mort/North_late_successional_species_var_partitioning.png"))

#oaks and hickories
oak.hickory.var.plt <- ggplot(data = var_summary %>% filter(COMMON %in% c("chestnut oak", "white oak","northern red oak", "hickory spp.")))+
  geom_bar(aes(x = STNAME, y = mean, fill = predictor.class), stat = "identity", position = "stack")+
  theme_bw()+ theme(axis.text.x = element_text(angle = 60, hjust = 1), 
                    panel.background = element_blank())+
  ylab("Across-tree variance in  \n p(survival) explained")+
  scale_fill_manual(values = color.pred.class, name = "")+
  coord_flip()+facet_wrap(~COMMON, scales = "free_y", ncol = 1)+
  xlab("")
ggsave(plot = oak.hickory.var.plt, height = 8, width = 6, units = "in", dpi = 300, 
       paste0(output.dir, "SPCD_stanoutput_joint_v3/predicted_mort/Oak_hickory_species_var_partitioning.png"))

# other taxa
other.spp.var.plt <- ggplot(data = var_summary %>% filter(COMMON %in% c("black cherry", "white ash", "yellow-poplar")))+
  geom_bar(aes(x = STNAME, y = mean, fill = predictor.class), stat = "identity", position = "stack")+
  theme_bw()+ theme(axis.text.x = element_text(angle = 60, hjust = 1), 
                    panel.background = element_blank())+
  ylab("Across-tree variance in  \n p(survival) explained")+
  scale_fill_manual(values = color.pred.class, name = "")+
  coord_flip()+facet_wrap(~COMMON, scales = "free_y", ncol = 1)
ggsave(plot = other.spp.var.plt , height = 6, width = 6, units = "in", dpi = 300, 
       paste0(output.dir, "SPCD_stanoutput_joint_v3/predicted_mort/Other_species_var_partitioning.png"))

variance.legend <- get_legend(lsuccess.var.plt)


forest.type.variance <- plot_grid(variance.legend, plot_grid(
  nmix.var.plt+theme(legend.position = "none"), 
  lsuccess.var.plt+theme(legend.position = "none"), 
          oak.hickory.var.plt+theme(legend.position = "none"), 
          other.spp.var.plt+theme(legend.position = "none"),
          boreal.var.plt+theme(legend.position = "none"), ncol = 5, align = "hv" ),
          rel_widths = c(0.25, 1))
png(height = 7, width = 15, units = "in", res = 350, paste0(output.dir, "SPCD_stanoutput_joint_v3/predicted_mort/Forest_type_var_partitioning.png"))
forest.type.variance
dev.off()

ggplot(data = var_summary)+ #%>% filter(COMMON %in% c("eastern white pine", "sugar maple", "red maple")))+
  geom_bar(aes(x = STNAME, y = mean, fill = predictor.class), stat = "identity", position = "stack")+
  theme_bw()+ theme(axis.text.x = element_text(angle = 60, hjust = 1), 
                    panel.background = element_blank())+
  ylab("Average proportion of across-tree \n p(survival) variance explained")+
  scale_fill_manual(values = color.pred.class, name = "")+
  coord_flip()+facet_wrap(~COMMON)+
  xlab("")



# get the top species for each of the states:
cleaned.data <- readRDS( "data/cleaned.data.mortality.TRplots.RDS")

cleaned.data <- cleaned.data %>% filter(!is.na(ba) & !is.na(slope) & ! is.na(physio) & !is.na(aspect))%>% 
  dplyr::select(state, county, pltnum, cndtn, point, tree, PLOT.ID, date,cycle, spp, dbhcur, dbhold, status, damage, Species, SPCD,
                remper, LAT_FIADB, LONG_FIADB, elev, DIA_DIFF, annual.growth, M, relative.growth, si, physio:RD) %>% distinct()

nspp <- cleaned.data %>% group_by(SPCD) %>% summarise(n = n(), 
                                                      pct = n/nrow(cleaned.data)) %>% arrange (desc(`pct`))

nspp$cumulative.pct <- cumsum(nspp$pct)



# link up to the species table:
nspp$COMMON <- FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD),]$COMMON
nspp$Species <- paste(FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD),]$GENUS, FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD),]$SPECIES)

#View(nspp)

nspp[1:17,]$COMMON

cleaned.data$SPGRPCD <- FIESTA::ref_species[match(cleaned.data$SPCD, FIESTA::ref_species$SPCD),]$E_SPGRPCD

SPGRP.df <- FIESTA::ref_codes %>% filter(VARIABLE %in% "SPGRPCD") %>% filter(VALUE %in% unique(cleaned.data$SPGRPCD))
cleaned.data$SPGRPNAME <- SPGRP.df[match(cleaned.data$SPGRPCD, SPGRP.df$VALUE),]$MEANING
cleaned.data$STNAME <- ref_statecd[match(cleaned.data$state, ref_statecd$VALUE),]$MEANING


# set up custom colors for species
# set the species order using the factors:
SP.TRAITS <- read.csv("data/NinemetsSpeciesTraits.csv") %>% filter(COMMON_NAME %in% c(unique(nspp[1:17,]$COMMON), "sweet birch", "Virginia pine"))
# order the trait db by softwood-hardwood, then shade tolerance, then name (this puts all the oaks together b/c hickory and red oak have the same tolerance values)
SP.TRAITS <- SP.TRAITS %>% group_by(SFTWD_HRDWD) %>% arrange(desc(SFTWD_HRDWD), desc(ShadeTol), desc(COMMON_NAME))
n.species.state <- cleaned.data %>% rename("COMMON" = "Species")%>%
  group_by(STNAME, SPCD, COMMON) %>% summarise(n())%>% 
  arrange(STNAME, desc(`n()`)) 

state.groups <- data.frame(state.regions = c(rep("Ohio Valley", 3), rep("Atlantic", 4), 
                                             rep("North", 4)), 
                          STNAME = c("Ohio", "Pennsylvania", "West Virginia", 
                                     "New Jersey", "Maryland", "Massachusetts", "Connecticut", 
                                     "Maine", "New Hampshire", "Vermont", "New York") )

top.5.spp <- n.species.state %>% left_join(.,state.groups) %>% ungroup() %>% 
  group_by(STNAME)%>% slice_head(n=5)%>%
  mutate(species.ranking.st = "Top 5")

top.5.region <- cleaned.data %>% rename("COMMON" = "Species")%>%
  left_join(.,state.groups) %>%
  group_by(state.regions, SPCD, COMMON) %>% summarise(n())%>% 
  arrange(state.regions, desc(`n()`)) %>% ungroup() %>% 
  group_by(state.regions) %>% slice_head(n=5)%>%
  mutate(species.ranking.region = "Top 5")

var_summary <- var_summary %>% left_join(.,top.5.spp %>% select(STNAME, COMMON, species.ranking.st))%>%
  left_join(.,top.5.region %>% select(state.regions, COMMON, species.ranking.region))


# group by far north states
spp.north <- unique(var_summary %>% filter(STNAME %in% c("Maine", "New Hampshire", "Vermont")& !is.na(species.ranking.st))%>%select(COMMON))
spp.ovr <- unique(var_summary %>% filter(STNAME %in% c("Ohio", "Pennsylvania", "West Virginia")& !is.na(species.ranking.st))%>%select(COMMON))
spp.atl <- unique(var_summary %>% filter(STNAME %in% c("New Jersey", "Maryland", "Massachusetts", "Connecticut")& !is.na(species.ranking.st))%>%select(COMMON))

var_summary$COMMON <- factor(var_summary$COMMON, levels = rev(SP.TRAITS$COMMON_NAME))

north.st.plt <- ggplot(data = var_summary %>% filter(STNAME %in% c("Maine", "New Hampshire", "Vermont") & COMMON %in% spp.north$COMMON))+
  geom_bar(aes(x = COMMON, y = mean, fill = predictor.class), stat = "identity", position = "stack")+
  theme_bw()+ theme(axis.text.x = element_text(angle = 60, hjust = 1), 
                    panel.background = element_blank())+
  ylab("Across-tree variance in \n p(survival) explained")+
  scale_fill_manual(values = color.pred.class, name = "")+
  coord_flip()+facet_wrap(~STNAME)+
  xlab("")
ggsave(plot = north.st.plt, height = 4, width = 8, units = "in", dpi = 300, 
       paste0(output.dir, "SPCD_stanoutput_joint_v3/predicted_mort/Northern_state_var_partitioning.png"))


# ohio valley
ovr.var.plt <- ggplot(data = var_summary %>% filter(STNAME %in% c("Ohio", "Pennsylvania", "West Virginia")& COMMON %in% spp.ovr$COMMON))+
  geom_bar(aes(x = COMMON, y = mean, fill = predictor.class), stat = "identity", position = "stack")+
  theme_bw()+ theme(axis.text.x = element_text(angle = 60, hjust = 1), 
                    panel.background = element_blank())+
  ylab("Across-tree variance in \n p(survival) explained")+
  scale_fill_manual(values = color.pred.class, name = "")+
  coord_flip()+facet_wrap(~STNAME)+
  xlab("")
ggsave(plot = ovr.var.plt, height = 4, width = 8, units = "in", dpi = 300, 
       paste0(output.dir, "SPCD_stanoutput_joint_v3/predicted_mort/Ohio_valley_state_var_partitioning.png"))

# coastal
atlantic.var.plt <- ggplot(data = var_summary %>% filter(STNAME %in% c("New Jersey", "Maryland", "Massachusetts", "Connecticut")& COMMON %in% spp.atl$COMMON))+
  geom_bar(aes(x = COMMON, y = mean, fill = predictor.class), stat = "identity", position = "stack")+
  theme_bw()+ theme(axis.text.x = element_text(angle = 60, hjust = 1), 
                    panel.background = element_blank())+
  ylab("Across-tree variance in \n p(survival) explained")+
  scale_fill_manual(values = color.pred.class, name = "")+
  coord_flip()+facet_wrap(~STNAME)+
  xlab("")
ggsave(height = 4, width = 8, units = "in", dpi = 300, 
       paste0(output.dir, "SPCD_stanoutput_joint_v3/predicted_mort/Atlantic_state_var_partitioning.png"))



state.region.variance <- plot_grid(variance.legend, plot_grid(
 ovr.var.plt+theme(legend.position = "none"), 
  atlantic.var.plt+theme(legend.position = "none"), 
  north.st.plt+theme(legend.position = "none"), ncol = 3, align = "hv" ),
  rel_widths = c(0.25, 2))

png(height = 3.5, width = 15, units = "in", res = 350, paste0(output.dir, "SPCD_stanoutput_joint_v3/predicted_mort/region_var_partitioning.png"))
state.region.variance
dev.off()


SP.Color <- c(# softwoods
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



#ggplot(data = var_summary %>% filter(predictor %in% main.preds))+
# geom_bar(aes(x = as.factor(SPCD), y = mean, fill = predictor), stat = "identity", position = "stack")

ggplot(data = var_summary %>% filter(!predictor %in% "DIA_DIFF_scaled"))+
  geom_bar(aes(x = Species, y = mean, fill = predictor), stat = "identity", position = "stack")+
  theme_minimal()+theme(axis.text = element_text(hjust = 1, angle = 60))+
  ylab("Average proportion of variance explained")





