data {
  int<lower=0> N;// N. observations
  // covariate data is in one big matrix
  int<lower=0> K;         // N. predictors 
  matrix[N,K] xM;        // Predictor matrix
  vector <lower = 1>[N] Remper; //list of remper observations
  array[N] int<lower=0, upper=1> y;//observations of mortality
  
  ///out of sample data for generated quantities
  int<lower=0> Nrep;// N. held out observations
  // covariate data is in one big matrix
  //int<lower=0> K;         // N. predictors 
  matrix[Nrep,K] xMrep;        // Predictor matrix
  array[Nrep] int<lower=0> Remperoos;
  
}
parameters {
  
  vector[K] u_beta;   //vector of length K for each coeff mean
  real alpha_SPP; //species level intercept
   

}
generated quantities{
//    // simulate data from the posterior
// 
   vector<lower=0, upper=1>[N] y_rep;//in sample predictions
   vector<lower=0, upper=1>[Nrep] y_hat;//out of sample predictions

//   // calculate the probabilities
  // vector<lower=0, upper=1>[N]  pSannualrep;//annual survival probability for in sample data
   vector<lower=0, upper=1>[N]  mMrep;//remeasurement probability of survival for in sample data
  // vector<lower=0, upper=1>[Nrep]  pSannualhat;//annual survival proability for out of sample data
   vector<lower=0, upper=1>[Nrep]  mMhat;//remper survival probability for out of sample data
   //vector[N] log_lik;
//  // //generate out of sample predictions ad yhat
 for (i in 1:Nrep) {
//  //    
//  //    
       //pSannual_logit[i] = alpha_SPP + xMrep[i]*u_beta;
//  //    
//  //    //generate annual survival predictions
        //pSannualhat[i] = inv_logit(alpha_SPP + xMrep[i]*u_beta);
//  //   
//  //    // convert to remeasurement period survival rate
         mMhat[i] = inv_logit(alpha_SPP + xMrep[i]*u_beta)^Remperoos[i];
//  // 
        y_hat[i] = bernoulli_rng(mMhat[i]);
 
 }
//  //  
//   vector[N] logit_p_annual_pred = xM * u_beta;
//   // individual log-likelihoods for use in loo
//   
 for (n in 1:N) {
//     //generate in sample predictions as yrep
//     
//     //generate annual survival predictions
        //pSannualrep[n] = inv_logit(alpha_SPP + xM[n]*u_beta);     
//     // convert to remeasurement period survival rate
        mMrep[n] = inv_logit(alpha_SPP + xM[n]*u_beta)^Remper[n];
        // get pointwise predictions for y from probabilities 
        y_rep[n] = bernoulli_rng(mMrep[n]);
        
        // //get point-wise log liklihood
        // log_lik[n] = bernoulli_lpmf(y[n] | mMrep[n]);
 
}
//   
   }
// 
