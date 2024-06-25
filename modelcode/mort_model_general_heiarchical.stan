data {
  int<lower=0> N;// N. observations
  int<lower=1> Nspp;//number of unique species
  array[N] int<lower=1, upper=Nspp> SPP;//indexing describing which observations belong to each species
  // covariate data is in one big matrix
  int<lower=1> K;       // N. predictors 
  array[N] row_vector[K] xM;      // Predictor matrix
  array[N] int<lower=0, upper=1> y; //observations of mortality
  
  ///out of sample data for generated quantities
 int<lower=0> Nrep;// N. held out observations
 array[Nrep] int<lower=1, upper=Nspp> SPPrep;//
 //int SPPrep[N]; //index describing which out of sample observations belong to each species
 //covariate data is in one big matrix
 array[Nrep] row_vector[K] xMrep;      // Predictor matrix
 //matrix[Nrep,K] xMrep;        // Predictor matrix
 
}
parameters {
  array[K] real mu_beta; // population-level means for betas
 // array[K] real<lower=0> sigma;
  array[Nspp] vector[K] u_beta; //species by K array of ubetas
  real alpha; // population intercept
  array[Nspp] real alpha_SPP; //species-specific intercepts

}
model {
//priors for population-level parameters
 mu_beta ~ normal(0, 5);
 alpha ~ normal(0,5);
 
 //priors for species=level means, centered on population level means
 //maybe we want species-specific variances??
for (s in 1:Nspp) {
  alpha_SPP[s] ~ normal(alpha, 1);
  u_beta[s] ~ normal(mu_beta, 5);
}

// estimated mean for each observation 
vector[N] mM; 
  for (n in 1:N) {
    mM[n] = alpha_SPP[SPP[n]] + xM[n] * u_beta[SPP[n]];
  }
//liklihood function
  y ~ bernoulli_logit(mM);

}
generated quantities{
   // simulate data from the posterior

  vector[N] y_rep;//in sample predictions
  vector[Nrep] y_hat;//out of sample predictions
  array[N] real log_lik; //observations of mortality
  // log-likelihood posterior
  //vector[N] log_lik; //calculate log likilhoods

  // calculate the probabilities
  vector[N]  mMrep;//in sample probs
  vector[Nrep]  mMhat;//in sample probs
 //generate out of sample predictions ad yhat
  for (i in 1:Nrep) {
    mMhat[i] = inv_logit(alpha_SPP[SPPrep[i]] + xMrep[i]*u_beta[SPPrep[i]]);
   //inv_logit(alpha + beta*x_test[i])
    y_hat[i] = bernoulli_rng(mMhat[i]);
  }

  // individual log-likelihoods for use in loo
   for (n in 1:N) {
    log_lik[n] = bernoulli_logit_lpmf(y[n] | xM[n] * u_beta[SPP[n]]);

    //generate in sample predictions as yrep
    mMrep[n] = inv_logit(alpha_SPP[SPP[n]] + xM[n]*u_beta[SPP[n]]);
    y_rep[n] = bernoulli_rng(mMrep[n]);
  }

  }

