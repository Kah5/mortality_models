functions{
  vector survival_log_lik(array[] int y, vector eta, vector Remper){
    int N = num_elements(eta);
    
    vector[N] log_p = Remper .* log_inv_logit(eta);
    
    vector[N] log1m_p = log1m_exp(log_p);
    vector[N] yv = to_vector(y);
    
    return(yv .* log_p + (1-yv) .* log1m_p);
    
  }
}
data {
  int<lower=0> N; // Number of tree level observations
  int<lower=1> Nspp; // number of unique species
  array[N] int<lower=1, upper=Nspp> SPP; // indexing describing which observations belong to each species
  int<lower=1> K; // Number of covariate predictors (varies by model number)
  matrix[N,K] xM; // predictor matrix (note we already standardized our predictiosr)
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
  vector<lower=1e-6>[K] sigma_beta; // scaling parameters for non-centered species u_betas
  real mu_alpha; // global mean for the species intercept
  real<lower=1e-6> sigma_alpha; // scaling parameter for species-level non-centered intercept
  vector[Nspp] z_alpha_SPP; // non-centered prior information for species-level random intercepts
  
}
transformed parameters {
  //vector[N] mM;//mean survival probability over remeasurement for bernoulli logit
  matrix[Nspp, K] u_beta; // actual species RE for the beta effects
 
  // calculate the non-centered u_betas from the sigm_beta and z_beta matrix
  profile("u_beta_alpahs_estimation"){
    for (k in 1:K) {
      u_beta[, k] = mu_beta[k] + sigma_beta[k] * z_beta[,k]; // multiply each column of z_beta by its sigma_beta[k]
    }
  }
    vector[Nspp] alpha_SPP = mu_alpha + sigma_alpha * z_alpha_SPP; // non-centered parameterization for species random effect on intercept
  

 // //calculate eta and use log space for calculating mM and log liklihood
  vector[N] eta = (alpha_SPP[SPP] + rows_dot_product(xM, u_beta[SPP]));
 
 
}
model {
 
 profile("priors"){
  // Priors
  
    mu_beta ~ normal(0, 2); // normal prior for beta effect population means
    sigma_beta ~ normal(0, 1); // normal prior for species-level effect scales
  
 
  mu_alpha ~ normal(0, 1); // population mean species intercept
  sigma_alpha ~ normal(0, 1); // Prior for scaling of species-level random intercepts
 
  
  to_vector(z_beta) ~ normal(0, 1); // normal for non-
 
  z_alpha_SPP ~ normal(0, 1); // Standard normal for latent species-level intercepts
 }
  // Likelihood
  profile("liklihood"){
    //sum point wise log-liklihoods implies bernoulli liklihood function y ~ bernoulli(mM); 
    target += sum(survival_log_lik(y, eta, Remper));
  }
 
}
generated quantities {
  vector[N] log_lik;
  //log liklihood
  profile("log_lik_generated"){
    log_lik = survival_log_lik(y, eta, Remper);
  }
}
