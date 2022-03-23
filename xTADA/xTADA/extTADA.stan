data {
  int<lower=1> NN; //Number of genes
  int<lower=1> K; //Number of classes
  int<lower=1> NCdn; //Number of de novo classes
  int<lower=1> Ncov; //Number of covariates
  int Ndn[NCdn]; //Number of trios

  int dataDN[NN, NCdn]; //denovo data: Kdn classes
  real mutRate[NN, NCdn]; //mutation rates: Kdn classes
  real covariates[NN, Ncov]; //gene expression rank percentile
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
  real<lower=0> A[Ncov];
  real<lower=0, upper=1> B[Ncov];
  real<lower=0, upper=1> C[Ncov];
  real<lower=lowerHyperGamma> hyperGammaMeanDN[NCdn]; //Hyper parameters for de novo relative risks
  real<lower=lowerGamma> gammaMeanDN[NCdn]; //parameters (in the sampling process) for de novo relative risks
}

transformed parameters {
  real L[Ncov];
  real hyperBetaDN[NCdn];
  if (adjustHyperBeta != 0) {
    for (i2i in 1:NCdn){
      hyperBetaDN[i2i] = exp(betaPars[1]*hyperGammaMeanDN[i2i]^(betaPars[2]) + betaPars[3]);
    }
  } else {
    hyperBetaDN = hyperBetaDN0;
  }
  for (i2i in 1:Ncov) {
    L[i2i] = (1-C[i2i]) * A[i2i] / (log(exp(A[i2i])+exp(B[i2i]*A[i2i])) - log(exp(B[i2i]*A[i2i])+1)); // scale factor for the distribution
  }
}

model {
  real ps[K];
  for (ip in 1:Ncov) {
    A[ip] ~ gamma(1, 0.01); // very week prior
    B[ip] ~ uniform(0, 1);
    C[ip] ~ uniform(0, 1);
  }
  pi0 ~ beta(1, 5); //prior for the proportion of risk genes

  //De novo data: sample for hyper priors (NCdn categories)
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
    //For covariates
    for (jj in 1:Ncov) {
      ps[1] = ps[1] + log(1 + pi0/(1-pi0)*(1-C[jj]-L[jj]/(1+exp(-A[jj]*(covariates[ii, jj]-B[jj])))));
      ps[2] = ps[2] + log(C[jj] + L[jj] / (1+exp(-A[jj]*(covariates[ii, jj]-B[jj]))));
    }
    target += log_sum_exp(ps);
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
    //For covariates
    for (jj in 1:Ncov) {
      ps[1] = ps[1] + log(1 + pi0/(1-pi0)*(1-C[jj]-L[jj]/(1+exp(-A[jj]*(covariates[ii, jj]-B[jj])))));
      ps[2] = ps[2] + log(C[jj] + L[jj] / (1+exp(-A[jj]*(covariates[ii, jj]-B[jj]))));
    }
    log_lik[ii] = log_sum_exp(ps);
  }
}