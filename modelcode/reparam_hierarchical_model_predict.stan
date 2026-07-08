data {
  int<lower=0> N; // Number of tree level observations
  int<lower=1> Nspp; // number of unique species
  array[N] int<lower=1, upper=Nspp> SPP; // indexing describing which observations belong to each species
  int<lower=1> K; // Number of covariate predictors (varies by model number)
  array[N] row_vector[K] xM; // predictor matrix (note we already standardized our predictiosr)
  array[N] int<lower=0, upper=1> y; // observations of survival
  vector<lower=1>[N] Remper; // list of remeasurement periods (lower bound here is 1, but observations are higher)
  
  // Out-of-sample data for generated quantities
  int<lower=0> Nrep; // Number of trees in held-out observations
  array[Nrep] int<lower=1, upper=Nspp> SPPrep; // index describing which out-of-sample observations belong to each species
  array[Nrep] row_vector[K] xMrep; // out-of-sample predictor covariate matrix
  vector<lower=1>[Nrep] Remperoos; // out-of-sample remeasurement period by tree
}
parameters {
 
   
   // Hyperparameters for species-level slopes
  vector[K] mu_beta;
  vector<lower=1e-6>[K] sigma_beta;

  // Hyperparameters for species-level intercepts
  real mu_alpha;
  real<lower=1e-6> sigma_alpha;
  
  matrix[Nspp, K] u_beta; // actual species RE for the beta effects
  vector[Nspp] alpha_SPP;

}
model{
  
}
generated quantities {
  //in sample predictions
  array[N] int<lower=0, upper=1> y_hat; // predicted survival for each tree 
  vector<lower=0, upper=1>[N] mMhat; // predicted survival probability for each tree
  // array[N] real<lower=0, upper=1> pSannual_hat; // Annual survival probabilities for in-sample predictions
  //vector[N] log_lik;
  for (n in 1:N) {
    //real logit_p_annual_hat =; // regression equations
    //pSannual_hat[n] = ; // annual survival probability
    mMhat[n] = pow(inv_logit( alpha_SPP[SPP[n]] + dot_product(xM[n], u_beta[SPP[n]])), Remper[n]); // cumulative survival over remeasurement period
    y_hat[n] = bernoulli_rng(mMhat[n]); // predicted in sample survival statues

          // //get point-wise log liklihood
  //log_lik[n] = bernoulli_lpmf(y[n] | mM[n]);
  }
  
  //Out of sample predictions, ut in 1 of 17 species groups
  array[Nrep] int<lower=0, upper=1> y_rep; //oos predicted survival
  vector<lower=0, upper=1>[Nrep] mMrep; //oos predicted survival probability
  // array[Nrep] real<lower=0, upper=1> pSannual_rep; // annual survival probabilities for out-of-sample predictions
  for (n in 1:Nrep) {
    //real logit_p_annual_rep = ; // regression predictions with posteriors
    //pSannual_rep[n] = inv_logit(alpha_SPP[SPPrep[n]] + dot_product(xMrep[n], u_beta[SPPrep[n]])); // annual survival probability
    mMrep[n] = pow(inv_logit(alpha_SPP[SPPrep[n]] + dot_product(xMrep[n], u_beta[SPPrep[n]])), Remperoos[n]); // cumulative survival over remeasurement period
    y_rep[n] = bernoulli_rng(mMrep[n]); // prediction oos survival
  }
}
