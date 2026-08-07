data {
  int<lower=1> Nspp; // number of unique species
  int<lower=1> K;    // number of covariate predictors
  int<lower=0> N;
  array[N] int<lower=1, upper=Nspp> SPP;
  matrix[N, K] xM;
  vector<lower=1>[N] Remper;

  // held-out prediction data -> y_rep, mMrep
  int<lower=0> Nrep;
  array[Nrep] int<lower=1, upper=Nspp> SPPrep;
  matrix[Nrep, K] xMrep;
  vector<lower=1>[Nrep] Remperoos;
}

parameters {
  // MUST match QR_reparam_hierarchical.stan's parameters block exactly --
  // same names, same types, same declared order -- since CmdStan validates
  // the fitted_params csv against this block.
  
  real mu_alpha;
  real<lower=1e-6> sigma_alpha;
  vector[Nspp] alpha_SPP;
  
  matrix[Nspp, K] u_beta;
  vector[K] mu_beta;
 
  
}

generated quantities {
 

  // in-sample predictions
  array[N] int<lower=0, upper=1> y_hat;
  vector<lower=0, upper=1>[N] mMhat;
  //for (n in 1:N) {
    mMhat = pow(inv_logit(alpha_SPP[SPP] +  rows_dot_product(xM, u_beta[SPP])), Remper);
    y_hat = bernoulli_rng(mMhat);
  //}

  // held-out / out-of-sample predictions
  array[Nrep] int<lower=0, upper=1> y_rep;
  vector<lower=0, upper=1>[Nrep] mMrep;
  //for (j in 1:Nrep) {
    mMrep = pow(inv_logit(alpha_SPP[SPPrep] + rows_dot_product(xMrep, u_beta[SPPrep])), Remperoos);
    y_rep = bernoulli_rng(mMrep);
  //}
}
