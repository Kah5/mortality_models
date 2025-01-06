data {
  int<lower=0> N;// N. observations

  // covariate data is in one big matrix
  int<lower=0> K;         // N. predictors 
  matrix[N,K] xM;        // Predictor matrix
  vector[N] Remper; //list of remper observations
  int<lower=0,upper=1> y[N]; //observations of mortality
  
  ///out of sample data for generated quantities
  int<lower=0> Nrep;// N. held out observations
  matrix[Nrep,K] xMrep;        // Predictor matrix
  int<lower = 0> Remperoos[Nrep]; //list of out of sample remper observations
  
  
  //basis function data spline terms
  int<lower=0> S; //S number of basis functions
  matrix[N, S] X_basis; // matrix of basis functions for the spline
  matrix[Nrep, S] X_basisMrep; // matrix of basis functions for the spline for out of sample data
  
 
  
}
parameters {
  
  vector[K] u_beta;   //vector of length K for each coeff mean
 // real alpha; //intercept mean
  real alpha_SPP; //species level intercept
  vector[S] eta; // coefficients for the basis functions
 real<lower = 0> sigma_beta[S];
}
model {

  //priors
  alpha_SPP ~ normal(0, 1);
  u_beta ~ normal(0, 5);
  //eta ~ normal(0, 1); // weakly informative priors for basis function coefficient
  sigma_beta ~ cauchy(0,2.5);
  eta[1] ~ normal(0,sigma_beta[1]);
  
  for(s in 2:S){
  eta[s] ~ normal(eta[s-1], sigma_beta[s]);
  }
  vector[N] pSannual;//mean annual surivival for bernoulli logit
  vector[N] mM;//mean survival probability over remeasurement for bernoulli logit

  //Liklihood function
  // annual survival probability
    pSannual[1:N] = alpha_SPP + xM[1:N]*u_beta + X_basis[1:N]*eta;
    
    // convert to remeasurement period survival rate
    mM[1:N] = pSannual[1:N]^Remper[1:N];
    
    y ~ bernoulli_logit(mM);

    
}
generated quantities{
   // simulate data from the posterior

  vector[N] y_rep;//in sample predictions
  vector[Nrep] y_hat;//out of sample predictions
  // log-likelihood posterior
  vector[N] log_lik; //calculate log likilhoods
  
   // calculate the probabilities
  vector[N]  pSannualrep;//annual survival probability for in sample data
  vector[N]  mMrep;//remeasurement probability of survival for in sample data
  vector[Nrep]  pSannualhat;//annual survival proability for out of sample data
  vector[Nrep]  mMhat;//remper survival probability for out of sample data
  
 //generate out of sample predictions ad yhat
  for (i in 1:Nrep) {
    
   //generate annual survival predictions
     pSannualhat[i] = inv_logit(alpha_SPP + xMrep[i]*u_beta  + X_basisMrep[i]*eta);
   
    // convert to remeasurement period survival rate
    mMhat[i] = pSannualhat[i]^Remperoos[i];
   
    y_hat[i] = bernoulli_rng(mMhat[i]);
  }
  
  // individual log-likelihoods for use in loo
   for (n in 1:N) {
     
    log_lik[n] = bernoulli_logit_lpmf(y[n] | xM[n] * u_beta + X_basis[n]*eta);
    
   //generate annual survival predictions
    pSannualrep[n] = inv_logit(alpha_SPP + xM[n]*u_beta + X_basis[n]*eta);
    
    // convert to remeasurement period survival rate
    mMrep[n] = pSannualrep[n]^Remper[n];
    
    //generate in sample predictions as yrep
   
    y_rep[n] = bernoulli_rng(mMrep[n]);
  }
  
  }

