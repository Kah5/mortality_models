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
  vector<lower=0>[K] sigma_beta; // scaling parameters for non-centered species u_betas
  real mu_alpha; // global mean for the species intercept
  real<lower=0> sigma_alpha; // scaling parameter for species-level non-centered intercept
  vector[Nspp] z_alpha_SPP; // non-centered prior information for species-level random intercepts
}
transformed parameters {
  matrix[Nspp, K] u_beta; // actual species RE for the beta effects
  // calculate the non-centered u_betas from the sigm_beta and z_beta matrix
  for (k in 1:K) {
    u_beta[, k] = mu_beta[k] + sigma_beta[k] * z_beta[, k]; // multiply each column of z_beta by its sigma_beta[k]
  }

  vector[Nspp] alpha_SPP = mu_alpha + sigma_alpha * z_alpha_SPP; // non-centered parameterization for species random effect on intercept
}
model {
  // Priors
  for (k in 1:K) {
    mu_beta[k] ~ normal(0, 2); // regularized prior for beta effect population means
    sigma_beta[k] ~ normal(0, 1); // regularized prior for species-level effect scales
  }
  to_vector(z_beta) ~ normal(0, 1); // normal for non-centered parameterization of beta
  mu_alpha ~ normal(0, 1); // population mean species intercept
  sigma_alpha ~ normal(0, 1); // Prior for scaling of species-level random intercepts
  z_alpha_SPP ~ normal(0, 1); // Standard normal for latent species-level intercepts

  // Likelihood
  for (n in 1:N) {
    //real logit_p_annual = ; // annual survival probability is f(species intercept + u_betas*covariates)
    // {
    // real pSannual = inv_logit(alpha_SPP[SPP[n]] + dot_product(xM[n], u_beta[SPP[n]])); // annual survival probability on the right scale
    // }
    real mM = pow(inv_logit(alpha_SPP[SPP[n]] + dot_product(xM[n], u_beta[SPP[n]])), Remper[n]); // cumulative survival over remeasurement period = pSannual^remper
    y[n] ~ bernoulli(mM); // bernoulli liklihood for survival over remper
  }
}
generated quantities {
  //in sample predictions
  // array[N] int<lower=0, upper=1> y_hat; // predicted survival for each tree 
  // vector<lower=0, upper=1>[N] mMhat; // predicted survival probability for each tree
  // array[N] real<lower=0, upper=1> pSannual_hat; // Annual survival probabilities for in-sample predictions
  vector[N] log_lik;
  for (n in 1:N) {
    // //real logit_p_annual_hat =; // regression equations
    // pSannual_hat[n] = inv_logit( alpha_SPP[SPP[n]] + dot_product(xM[n], u_beta[SPP[n]])); // annual survival probability
    // mMhat[n] = pow(pSannual_hat[n], Remper[n]); // cumulative survival over remeasurement period
    // y_hat[n] = bernoulli_rng(mMhat[n]); // predicted in sample survival statues
    // 
          // //get point-wise log liklihood
  log_lik[n] = bernoulli_lpmf(y[n] | mM[n]);
  }
  
  //Out of sample predictions
  // array[Nrep] int<lower=0, upper=1> y_rep; //oos predicted survival
  // vector<lower=0, upper=1>[Nrep] mMrep; //oos predicted survival probability
  // array[Nrep] real<lower=0, upper=1> pSannual_rep; // annual survival probabilities for out-of-sample predictions
  // for (n in 1:Nrep) {
  //   //real logit_p_annual_rep = ; // regression predictions with posteriors
  //   pSannual_rep[n] = inv_logit(alpha_SPP[SPPrep[n]] + dot_product(xMrep[n], u_beta[SPPrep[n]])); // annual survival probability
  //   mMrep[n] = pow(pSannual_rep[n], Remperoos[n]); // cumulative survival over remeasurement period
  //   y_rep[n] = bernoulli_rng(mMrep[n]); // prediction oos survival
  // }
}
