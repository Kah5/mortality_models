import cmdstanpy
cmdstanpy.install_cmdstan()
import cmdstanpy
import json

# Check where CmdStan is installed
print("CmdStan path:", cmdstanpy.cmdstan_path())

# Test model compilation
model_code = """
data {
  int<lower=0> N;
  array[N] real y;
}
parameters {
  real mu;
  real<lower=0> sigma;
}
model {
  mu ~ normal(0, 10);
  sigma ~ cauchy(0, 5);
  y ~ normal(mu, sigma);
}
"""
with open("bernoulli.stan", "w") as f:
    f.write(model_code)

# Compile
model = cmdstanpy.CmdStanModel(stan_file="bernoulli.stan")


import os
from cmdstanpy import CmdStanModel, cmdstan_path

# 1. Point to the built-in Bernoulli model
bernoulli_stan = os.path.join(cmdstan_path(), 'examples', 'bernoulli', 'bernoulli.stan')
bernoulli_data = os.path.join(cmdstan_path(), 'examples', 'bernoulli', 'bernoulli.data.json')

# 2. Compile the Stan model
model = CmdStanModel(stan_file=bernoulli_stan)

# 3. Run the Hamiltonian Monte Carlo (HMC) sampler
fit = model.sample(data=bernoulli_data)

# 4. View summary results
print(fit.summary())

# 5. Check diagnostics to ensure no convergence issues
print(fit.diagnose())

model = CmdStanModel(stan_file= "mort_model_general.stan")


