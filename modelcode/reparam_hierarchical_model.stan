data {
  int<lower=0> N; // Number of tree level observations
  int<lower=1> Nspp; // number of unique species
  array[N] int<lower=1, upper=Nspp> SPP; // indexing describing which observations belong to each species
  int<lower=1> K; // Number of covariate predictors (varies by model number)
  array[N] row_vector[K] xM; // predictor matrix (note we already standardized our predictiosr)
  array[N] int<lower=0, upper=1> y; // observations of survival
  vector<lower=1>[N] Remper; // list of remeasurement periods (lower bound here is 1, but observations are higher)
  
  // Out-of-sample data for generated quantities
  // int<lower=0> Nrep; // Number of trees in held-out observations
  // array[Nrep] int<lower=1, upper=Nspp> SPPrep; // index describing which out-of-sample observations belong to each species
  // array[Nrep] row_vector[K] xMrep; // out-of-sample predictor covariate matrix
  // vector<lower=1>[Nrep] Remperoos; // out-of-sample remeasurement period by tree
}
parameters {
  vector[K] mu_beta; // population-level means for u_betas
  matrix[Nspp, K] z_beta; // non-centered prior information for species u_betas
  //array[Nspp] vector[K] z_beta;
  vector<lower=0>[K] sigma_beta; // scaling parameters for non-centered species u_betas
  real mu_alpha; // global mean for the species intercept
  real<lower=0> sigma_alpha; // scaling parameter for species-level non-centered intercept
  vector[Nspp] z_alpha_SPP; // non-centered prior information for species-level random intercepts
  
}
transformed parameters {
  vector[N] mM;//mean survival probability over remeasurement for bernoulli logit
  matrix[Nspp, K] u_beta; // actual species RE for the beta effects
 
  // calculate the non-centered u_betas from the sigm_beta and z_beta matrix
  for (k in 1:K) {
    u_beta[, k] = mu_beta[k] + sigma_beta[k] * z_beta[,k]; // multiply each column of z_beta by its sigma_beta[k]
  }

  vector[Nspp] alpha_SPP = mu_alpha + sigma_alpha * z_alpha_SPP; // non-centered parameterization for species random effect on intercept
 
  for (n in 1:N) {
       mM[n] = pow(inv_logit(alpha_SPP[SPP[n]] + dot_product(xM[n], u_beta[SPP[n]])), Remper[n]); // cumulative survival over remeasurement period = pSannual^remper
  }
 
}
model {
 
  // Priors
  //for (k in 1:K) {
   // mu_beta[k] ~ normal(0, 2); // regularized prior for beta effect population means
    //sigma_beta[k] ~ normal(0, 1); // regularized prior for species-level effect scales
    
    mu_beta ~ normal(0, 2); // normal prior for beta effect population means
    sigma_beta ~ normal(0, 2); // normal prior for species-level effect scales
  //}
 
  mu_alpha ~ normal(0, 1); // population mean species intercept
  sigma_alpha ~ normal(0, 1); // Prior for scaling of species-level random intercepts
 
  
  to_vector(z_beta) ~ normal(0, 1); // normal for non-
  //to_vector(z_beta) ~ normal(0, 1); // normal for non-centered parameterization of beta
  z_alpha_SPP ~ normal(0, 1); // Standard normal for latent species-level intercepts
 
  // Likelihood
 
  y ~ bernoulli(mM); // bernoulli liklihood for survival over remper
}
generated quantities {
  //log liklihood
  vector[N] log_lik;
  for (n in 1:N) {
   // //get point-wise log liklihood
  log_lik[n] = bernoulli_lpmf(y[n] | mM[n]);
  }
  
}
