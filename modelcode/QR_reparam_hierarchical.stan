data {
  int<lower=0> N; // Number of tree level observations
  int<lower=1> Nspp; // number of unique species
  array[N] int<lower=1, upper=Nspp> SPP; // indexing describing which observations belong to each species
  int<lower=1> K; // Number of covariate predictors (varies by model number)
  matrix[N, K] xM; // predictor matrix (note we already standardized our predictiosr)
  array[N] int<lower=0, upper=1> y; // observations of survival
  vector<lower=1>[N] Remper; // list of remeasurement periods (lower bound here is 1, but observations are higher)
  
  // Out-of-sample data for generated quantities
  // int<lower=0> Nrep; // Number of trees in held-out observations
  // array[Nrep] int<lower=1, upper=Nspp> SPPrep; // index describing which out-of-sample observations belong to each species
  // array[Nrep] row_vector[K] xMrep; // out-of-sample predictor covariate matrix
  // vector<lower=1>[Nrep] Remperoos; // out-of-sample remeasurement period by tree
}
transformed data {
  // thin QR decomposition on the betas
 
  matrix[N, K] Q_ast = qr_thin_Q(xM) * sqrt(N - 1);
  matrix[K, K] R_ast = qr_thin_R(xM) / sqrt(N - 1);
  matrix[K, K] R_ast_inverse = inverse(R_ast);
}
parameters {
  vector[K] mu_theta; // population-level means for u_betas, rotated in Q space
  matrix[Nspp, K] z_theta; // non-centered prior information for species u_betas, in Q space
  //array[Nspp] vector[K] z_beta;
  vector<lower=1e-6>[K] sigma_theta; // scaling parameters for non-centered species u_betas, in Q space
  real mu_alpha; // global mean for the species intercept
  real<lower=1e-6> sigma_alpha; // scaling parameter for species-level non-centered intercept
  vector[Nspp] z_alpha_SPP; // non-centered prior information for species-level random intercepts
  
}
transformed parameters {
  vector[N] mM;//mean survival probability over remeasurement for bernoulli logit
  matrix[Nspp, K] theta_beta; // actual species RE for the beta effects, rotated inin Q space
  profile("u_beta_alpahs_estimation"){
  // calculate the non-centered u_betas from the sigm_beta and z_beta matrix
  for (k in 1:K) {
    theta_beta[, k] = mu_theta[k] + sigma_theta[k] * z_theta[,k]; // multiply each column of z_beta by its sigma_beta[k]
  }
}
  vector[Nspp] alpha_SPP = mu_alpha + sigma_alpha * z_alpha_SPP; // non-centered parameterization for species random effect on intercept
 
  // for (n in 1:N) {
  //      mM[n] = pow(inv_logit(alpha_SPP[SPP[n]] + dot_product(xM[n], u_beta[SPP[n]])), Remper[n]); // cumulative survival over remeasurement period = pSannual^remper
  // }
 // in
 // //calculate eta and use log space for calculating mM
 vector[N] eta = (alpha_SPP[SPP] + rows_dot_product(Q_ast, theta_beta[SPP]));
 mM = pow(inv_logit(eta), Remper); // cumulative survival over remeasurement period = pSannual^remper

 // loop is implied here
  // for (n in 1:N) {
 //mM = pow(inv_logit(alpha_SPP[SPP] + rows_dot_product(Q_ast, theta_beta[SPP])), Remper); // cumulative survival over remeasurement period = pSannual^remper
  // }
 
}
model {
 
  // Priors
   profile("priors"){
    mu_theta ~ normal(0, 2); // normal prior for beta effect population means
    sigma_theta ~ normal(0, 1); // normal prior for species-level effect scales
  //}
 
  mu_alpha ~ normal(0, 1); // population mean species intercept
  sigma_alpha ~ normal(0, 1); // Prior for scaling of species-level random intercepts
 
  
  to_vector(z_theta) ~ normal(0, 1); // normal for non-
  //to_vector(z_beta) ~ normal(0, 1); // normal for non-centered parameterization of beta
  z_alpha_SPP ~ normal(0, 1); // Standard normal for latent species-level intercepts
   }
   
  // Likelihood
 profile("liklihood"){
  y ~ bernoulli(mM); // bernoulli liklihood for survival over remper
 }
}
generated quantities {
  
  
    matrix[Nspp, K] u_beta;
    vector[K] mu_beta;
profile("recovering betas from QR"){
  //need to recover actual betas
  u_beta = theta_beta * R_ast_inverse';
  mu_beta = R_ast_inverse * mu_theta;
}  
  //log liklihood
  vector[N] log_lik;
  
profile("log_lik_generated"){
  for (n in 1:N) {
   // //get point-wise log liklihood
    log_lik[n] = bernoulli_lpmf(y[n] | mM[n]);
 
   }
  
}
}
