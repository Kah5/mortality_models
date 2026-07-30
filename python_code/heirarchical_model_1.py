import cmdstanpy
import json
import os
from cmdstanpy import CmdStanModel
#import arviz as az
#import matplotlib.pyplot as plt

#az.style.use("arviz-darkgrid")


# 1. Compile the Stan model
my_stanfile = os.path.join('.', 'modelcode','QR_reparam_hierarchical.stan')
hierarchical_model = CmdStanModel(stan_file=my_stanfile)

# checking the model names
print(hierarchical_model.exe_file)
print(hierarchical_model.name)
print(hierarchical_model.stan_file)


# 2. provide path to the data and set up output directories
#for i in range(1, 10):
mod_data = os.path.join('.','modeldata/hierarchical_data_model_1.json')
os.makedirs('data/output/hierarchical_model_1', exist_ok=True)
output_dir = os.path.join(".",'data/output/hierarchical_model_1')


# 3. Run the Hamiltonian Monte Carlo (HMC) sampler
fit = hierarchical_model.sample(data=mod_data, chains=4, parallel_chains=4, inits=0.1)

# 4. Save the stan csv files to the output directory

fit.save_csvfiles(dir=output_dir)

# 5. Check diagnostics to ensure no convergence issues
#print(fit.diagnose())

# 5. View summary results
#mod_samples = az.from_cmdstanpy(fit) 

# plot traceplots
#az.plot_trace(mod_samples, var_names=[ "mu_beta", "mu_alpha", "sigma_beta", "sigma_alpha"]) 
# 6. View summary results
#print(fit.summary())

