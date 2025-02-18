data {
  int<lower=0> N;// N. observations
  int<lower=1> Nspp;//number of unique species
  array[N] int<lower=1, upper=Nspp> SPP;//indexing describing which observations belong to each species
  // covariate data is in one big matrix
  int<lower=1> K;       // N. predictors 
  array[N] row_vector[K] xM;      // Predictor matrix
  array[N] int<lower=0, upper=1> y; //observations of mortality
  vector<lower = 8, upper = 20>[N] Remper; //list of remper observations
  
  ///out of sample data for generated quantities
 int<lower=0> Nrep;// N. held out observations
 array[Nrep] int<lower=1, upper=Nspp> SPPrep;//index describing which out of sample observations belong to each species

 array[Nrep] row_vector[K] xMrep;      // out of sample Predictor matrix
 int<lower = 0> Remperoos[Nrep]; //list of out of sample remper observations
 
}
parameters {
  vector[K] mu_beta; // population-level means for betas
  matrix[Nspp, K] beta_tilde; //non-centered priors for species effects
  vector<lower = 0>[K]  sigma_s; // K array of across-species variance for each parameter
  
  real alpha; // population intercept
  real<lower = 0> sigma_alpha; //species shared variance on intercepts
  array[Nspp] real alpha_tilde; //non-centered prior for alpha
 
}

transformed parameters{
  vector[N] logit_p_annual;
  array[Nspp] vector[K] u_beta; //species by K array of ubetas
  array[Nspp] real alpha_SPP; //species-specific intercepts
  
for(s in 1:Nspp) { 
  alpha_SPP[s] = alpha + sigma_alpha * alpha_tilde[s];
  u_beta[s,] =  mu_beta +  beta_tilde[s,]*sigma_s ;
}
  
  for (n in 1:N) {
    // annual survival probability
  
    logit_p_annual[n] =  alpha_SPP[SPP[n]] + xM[n,] * to_vector(u_beta[SPP[n],]); //annual suvival probability
}    

   vector<lower = 0, upper  = 1>[N] pSannual;//mean annual surivival for bernoulli logit
   vector<lower = 0, upper  = 1> [N] mM;//mean survival probability over remeasurement for bernoulli logit
 
    pSannual = inv_logit(logit_p_annual);
    // convert to remeasurement period survival rate
    mM = pSannual^Remper;

}

model {
//priors for population-level parameters

 alpha ~ normal(0,5);
 alpha_tilde ~ normal(0,1);//prior for non-centered alpha
 sigma_alpha ~ cauchy(0,2.5); //prior for the across-species variance for alphas
 mu_beta[1:K] ~ normal(0, 5);
 
for(s in 1:Nspp){
 beta_tilde[s,1:K] ~ normal(0,1);//prior for non-centered betas
}
 sigma_s ~ cauchy(0,2.5);// prior for the K array of across-species variances for each parameter


 //priors for species=level means, centered on population level means
 
// for (s in 1:Nspp) {
//   alpha_SPP[s] ~ normal(alpha, sigma_alpha);//sigma_s[1]);
//   u_beta[s, 1:K] ~ normal(mu_beta[1:K], sigma_s[1:K]);
// }
  

  //Liklihood function

  y ~ bernoulli(mM);

}
 generated quantities{
   // simulate data from the posterior

  vector[N] y_rep;//in sample predictions
  vector[Nrep] y_hat;//out of sample predictions
  array[N] real log_lik; //observations of mortality
  // log-likelihood posterior
  //vector[N] log_lik; //calculate log likilhoods

   // calculate the probabilities
  vector[N]  pSannualrep;//annual survival probability for in sample data
  vector[N]  mMrep;//remeasurement probability of survival for in sample data
  vector[Nrep]  pSannualhat;//annual survival proability for out of sample data
  vector[Nrep]  mMhat;//remper survival probability for out of sample data

 //generate out of sample predictions ad yhat
  for (i in 1:Nrep) {
    //generate annual survival predictions
    pSannualhat[i] = inv_logit(alpha_SPP[SPPrep[i]] + xMrep[i]*u_beta[SPPrep[i]]);


    // convert to remeasurement period survival rate
    mMhat[i] = pSannualhat[i]^Remperoos[i];

    y_hat[i] = bernoulli_rng(mMhat[i]);
  }

  // individual log-likelihoods for use in loo
   for (n in 1:N) {
    log_lik[n] = bernoulli_logit_lpmf(y[n] | xM[n] * u_beta[SPP[n]]);

  //generate annual survival predictions
    pSannualrep[n] = inv_logit(alpha_SPP[SPP[n]] + xM[n]*u_beta[SPP[n]]);

    // convert to remeasurement period survival rate
    mMrep[n] = pSannualrep[n]^Remper[n];

    //generate in sample predictions as yrep

    y_rep[n] = bernoulli_rng(mMrep[n]);
  }

  }

