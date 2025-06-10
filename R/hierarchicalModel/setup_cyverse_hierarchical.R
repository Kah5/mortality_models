# install some extra packages
install.packages("FIESTA")
install.packages("mltools")

# copy the data-store files
system(paste("cp -r", "data-store/data/iplant/home/kellyheilman/SPCD_standata_general_v4/", 
             "SPCD_standata_general_full_standardized_v3"))


# code to create output folders and copy input folders:

dir.create("SPCD_stanoutput_joint_v3")
dir.create("SPCD_stanoutput_joint_v3/images")
dir.create("SPCD_stanoutput_joint_v3/computational_resources")
dir.create("SPCD_stanoutput_joint_v3/samples")
dir.create("SPCD_stanoutput_joint_v3/predicted_mort")


# create the R and stan code folders
dir.create("R/hierarchicalModels", recursive = T)
dir.create("modelcode")
# need to manually upload code right now to these folders



nspp <- data.frame(SPCD = c(316, 318, 833, 832, 261, 531, 802, 129, 762,  12, 541,  97, 621, 400, 371, 241, 375))
nspp$Species <- paste(FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD),]$GENUS, FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD),]$SPECIES)
nspp$COMMON <- paste(FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD),]$GENUS, FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD),]$SPECIES)


# copy the posterior estimates files (for variance partitioning, if the models have been run)
file.list <- system("find data-store/data/iplant/home/kellyheilman/analyses/Mortality_hierarchical_models-2025-04-07-18-56-47.7/SPCD_stanoutput_joint_v3_model_6/ -name '*1000samples.rds'", intern = TRUE ) 

system(paste("cp", paste(file.list, collapse = " "), "SPCD_stanoutput_joint_v3/samples"))
