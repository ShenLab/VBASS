data {
    int<lower=1> NN; //Number of genes
    int<lower=1> K; //Number of classes
    int<lower=1> NCdn; //Number of de novo classes
    int Ndn[NCdn]; //Number of trios

    int dataDN[NN, NCdn]; //denovo data: Kdn classes
    real mutRate[NN, NCdn]; //mutation rates: Kdn classes
    real betaPars[3]; //These parameters are used to adjust beta values as a non-linear function of gamma means

    //These below parameters should be default
    real<lower=0> upperPi0;
    real<lower=0> lowerPi0;
    real<lower=0> lowerHyperGamma; //Low limit for mean of relative risks
    real<lower=0> lowerGamma; //Low limit for relative risks
    real<lower=0> lowerBeta;
    real<lower=0> hyperBetaDN0[NCdn];
    int<lower=0> adjustHyperBeta; //Option to adjust betas or not; if this option is used (!=0) then beta is a non-linear function of gamma means
    real<lower=0> hyper2GammaMeanDN[NCdn]; //hyperGammaMeanDN ~ gamma(hyper2GammaMeanDN, hyper2BetaDN)
    real<lower=0> hyper2BetaDN[NCdn]; //hyperGammaMeanDN ~ gamma(hyper2GammaMeanDN, hyper2BetaDN)

    }

parameters {
    real<lower=lowerPi0,upper=upperPi0> pi0; //Proportion of risk genes
    real<lower=lowerHyperGamma> hyperGammaMeanDN[NCdn]; //Hyper parameters for de novo relative risks
    real<lower=lowerGamma> gammaMeanDN[NCdn]; //parameters (in the sampling process) for de novo relative risks

}

transformed parameters {
    real hyperBetaDN[NCdn];
    if (adjustHyperBeta != 0) {
      for (i2i in 1:NCdn){
            hyperBetaDN[i2i] = exp(betaPars[1]*hyperGammaMeanDN[i2i]^(betaPars[2]) + betaPars[3]);

       }
   }
    else {
        hyperBetaDN = hyperBetaDN0;
        }
    }
 model {
     int newIndex;
     real ps[K];
     real sumDN[2];
     pi0 ~ beta(1, 5); //prior for the proportion of risk genes

     //Both CC + DN


  //De novo data: sample for hyper priors (NPdn populations and Kdn categories)
    for (ip in 1:NCdn){
          hyperGammaMeanDN[ip]	~ gamma(hyper2GammaMeanDN[ip], hyper2BetaDN[ip]); //gamma(1, 0.1); //normal(15, 10);
    }

    for (ip in 1:NCdn){
          gammaMeanDN[ip] ~ gamma(hyperGammaMeanDN[ip]*hyperBetaDN[ip], hyperBetaDN[ip]);
           }

////Main program
//Loop through data points
////
     for (ii in 1:NN){
         ps[1] = log1m(pi0);
         ps[2] = log(pi0);
   //For de novo data
         for (jj in 1:NCdn){
             ps[1] = ps[1] + poisson_lpmf(dataDN[ii, jj] | Ndn[jj]*2*mutRate[ii, jj]); //Add Null hypothesis
             ps[2] = ps[2] + poisson_lpmf(dataDN[ii, jj] | Ndn[jj]*2*mutRate[ii, jj]*gammaMeanDN[jj]); //Add alternative hypothesis
             }
         target += log_sum_exp(ps);
         //increment_log_prob(log_sum_exp(ps));
         }
}
generated quantities {
  vector[NN] log_lik;
  real ps[2];
  for (ii in 1:NN){
         ps[1] = log1m(pi0);
         ps[2] = log(pi0);
   //For de novo data
         for (jj in 1:NCdn){
             ps[1] = ps[1] + poisson_lpmf(dataDN[ii, jj] | Ndn[jj]*2*mutRate[ii, jj]); //Add Null hypothesis
             ps[2] = ps[2] + poisson_lpmf(dataDN[ii, jj] | Ndn[jj]*2*mutRate[ii, jj]*gammaMeanDN[jj]); //Add alternative hypothesis
             }
         log_lik[ii] = log_sum_exp(ps);
  }
}