import cmdstanpy
import json
import os
from cmdstanpy import CmdStanModel
print("CmdStan path:", cmdstanpy.cmdstan_path())


# 1. Compile the Stan model
my_stanfile = os.path.join('.', 'modelcode','QR_reparam_hierarchical.stan')
hierarchical_model = CmdStanModel(stan_file=my_stanfile)

# checking the model names
print(hierarchical_model.exe_file)
print(hierarchical_model.name)
print(hierarchical_model.stan_file)


# 2. provide path to the data and set up output directories
mod_data = os.path.join('.','modeldata/hierarchical_data_model_6.json')
os.makedirs('data/output/hierarchical_model_6', exist_ok=True)
output_dir = os.path.join('.','data/output/hierarchical_model_6')


# 3. Run the Hamiltonian Monte Carlo (HMC) sampler--note max_treedepth = 10
fit = hierarchical_model.sample(data=mod_data, 
                                chains=4, 
                                parallel_chains=4,
                                inits=0.1,  
                                max_treedepth=10)

# 4. Save the stan csv files to the output directory

fit.save_csvfiles(dir=output_dir)

# 5. Check diagnostics to ensure no convergence issues
#print(fit.diagnose())

# 5. View summary results
#mod_samples = az.from_cmdstanpy(fit) 


