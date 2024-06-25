data {
  int<lower=0> N;// N. observations
  // covariate data is in one big matrix
  int<lower=0> K;         // N. predictors 
  matrix[N,K] xM;        // Predictor matrix
  int<lower=0,upper=1> y[N]; //observations of mortality
  
  ///out of sample data for generated quantities
  int<lower=0> Nrep;// N. held out observations
  // covariate data is in one big matrix
 // int<lower=0> K;         // N. predictors 
  matrix[Nrep,K] xMrep;        // Predictor matrix
 
  
}
parameters {
  
  vector[K] u_beta;   //vector of length K for each coeff mean
 // real alpha; //intercept mean
  real alpha_SPP; //species level intercept
  

}
model {

  //priors
  alpha_SPP ~ normal(0, 1);
  u_beta ~ normal(0, 5);
  
  vector[N] mM;//mean for bernoulli logit

  //Liklihood function
 // y ~ bernoulli_logit(alpha_SPP[SPP] + beta_si * si + beta_growth * annual_growth + beta_DIA*DIA + beta_MAP*MAP + beta_MATmin*MATmin + beta_MATmax*MATmax);
//for(n in 1:N){
  //  mM[n] = alpha_SPP + xM[n]*u_beta;
    //}
    

    mM[1:N] = alpha_SPP + xM[1:N]*u_beta;
    y ~ bernoulli_logit(mM);

}
generated quantities{
   // simulate data from the posterior

  vector[N] y_rep;//in sample predictions
   vector[Nrep] y_hat;//out of sample predictions
  // log-likelihood posterior
  vector[N] log_lik; //calculate log likilhoods
  
  // calculate the probabilities
  vector[N]  mMrep;//in sample probs
  vector[Nrep]  mMhat;//in sample probs
 //generate out of sample predictions ad yhat
  for (i in 1:Nrep) {
    mMhat[i] = inv_logit(alpha_SPP + xMrep[i]*u_beta);
   //inv_logit(alpha + beta*x_test[i])
    y_hat[i] = bernoulli_rng(mMhat[i]);
  }
  
  // individual log-likelihoods for use in loo
   for (n in 1:N) {
    log_lik[n] = bernoulli_logit_lpmf(y[n] | xM[n] * u_beta);
  
    //generate in sample predictions as yrep
    mMrep[n] = inv_logit(alpha_SPP + xM[n]*u_beta);
    y_rep[n] = bernoulli_rng(mMrep[n]);
  }
  
  }

