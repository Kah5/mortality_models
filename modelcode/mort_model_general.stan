data {
  int<lower=0> N;// N. observations
  // covariate data is in one big matrix
  int<lower=0> K;         // N. predictors 
  matrix[N,K] xM;        // Predictor matrix
  vector <lower = 1>[N] Remper; //list of remper observations
  int<lower=0, upper=1> y[N]; //observations of mortality
  
  ///out of sample data for generated quantities
  int<lower=0> Nrep;// N. held out observations
  // covariate data is in one big matrix
  //int<lower=0> K;         // N. predictors 
  matrix[Nrep,K] xMrep;        // Predictor matrix
  int<lower = 0> Remperoos[Nrep]; //list of out of sample remper observations
  
}
parameters {
  
  vector[K] u_beta;   //vector of length K for each coeff mean
 // real alpha; //intercept mean
  real alpha_SPP; //species level intercept
   

}
transformed parameters{
  vector<lower=0, upper=1>[N] p_annual; 
  vector[N] logit_p_annual;
  vector<lower = 0, upper = 1>[N] mM;//mean survival probability over remeasurement for bernoulli logit
 
  
  logit_p_annual = alpha_SPP + xM[1:N]*u_beta; //regression
  p_annual = inv_logit(logit_p_annual);//annual suvival probability
  
  
  mM[1:N] = p_annual[1:N]^Remper[1:N]; //remper survival probability
}
model {

  //priors
  alpha_SPP ~ normal(0, 5);
  u_beta ~  normal(0, 5);//student_t(4,0,1);
 // print(u_beta);
 //print(alpha_SPP);
  //vector[N] pSannual;//mean annual surivival for bernoulli logit
  //vector[N] mM;//mean survival probability over remeasurement for bernoulli logit
  
 
  //Liklihood function
  
   
    // remeasurement survival rate

    // convert to remeasurement period survival rate
   // mM[1:N] = pow(pSannual[1:N], Remper[1:N]);
    
    y ~ bernoulli(mM);

}
generated quantities{
//    // simulate data from the posterior
// 
   vector[N] y_rep;//in sample predictions
   vector[Nrep] y_hat;//out of sample predictions
//   // log-likelihood posterior
 //  vector[N] log_lik; //calculate log likilhoods
//   
//   // calculate the probabilities
   vector[N]  pSannualrep;//annual survival probability for in sample data
   vector[N]  mMrep;//remeasurement probability of survival for in sample data
   vector[Nrep]  pSannualhat;//annual survival proability for out of sample data
   vector[Nrep]  mMhat;//remper survival probability for out of sample data
//  // //generate out of sample predictions ad yhat
 for (i in 1:Nrep) {
//  //    
//  //    
       //pSannual_logit[i] = alpha_SPP + xMrep[i]*u_beta;
//  //    
//  //    //generate annual survival predictions
        pSannualhat[i] = inv_logit(alpha_SPP + xMrep[i]*u_beta);
//  //   
//  //    // convert to remeasurement period survival rate
         mMhat[i] = pSannualhat[i]^Remperoos[i];
//  // 
        y_hat[i] = bernoulli_rng(mMhat[i]);
 }
//  //  
//   vector[N] logit_p_annual_pred = xM * u_beta;
//   // individual log-likelihoods for use in loo
//   
 for (n in 1:N) {
//      
 //real p_annual_pred = inv_logit(logit_p_annual_pred[n]);
// p_total = p_annual_pred^Remper[n];
//      
//log_lik[n] = bernoulli_logit_lpmf(y[n] | p_total);
//   
//     //generate in sample predictions as yrep
//     
//     //generate annual survival predictions
//    // mMrep[n] = inv_logit(alpha_SPP + xM[n]*u_beta)^Remper[n];
pSannualrep[n] = inv_logit(alpha_SPP + xM[n]*u_beta);     
//     // convert to remeasurement period survival rate
mMrep[n] = pSannualrep[n]^Remper[n];
// 
 y_rep[n] = bernoulli_rng(mMrep[n]);
}
//   
   }
// 
