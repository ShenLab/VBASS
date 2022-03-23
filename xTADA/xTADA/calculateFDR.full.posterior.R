Bayesian.FDR.full.posterior <- function(PP, alpha=0.05) {
  # convert PP to FDR
  q0 <- 1 - PP # posterior probability of null model
  
  # the FDR at each PPA cutoff
  n <- length(PP)
  FDR <- numeric(n)
  for (i in 1:n) FDR[i] <- sum(q0[1:i]) / i 
  
  # the cutoff
  t <- 1
  while (t <= length(q0) & mean(q0[1:t]) <= alpha) { t <- t+1 }
  return (list(FDR=FDR, ND=t))
}


calculateFDR.full.posterior <- function(pars,
                                        caseData = NULL,
                                        controlData = NULL,
                                        dnData = NULL,
                                        mutData = NULL,
                                        geneName){
  # utilizing full posterior distribution from MCMC, support only DN data for now.
  mcmc.samples <- pars$mcmc.samples
  
  outData <- data.frame(geneName)
  if (!is.null(dnData))
    outData <- cbind(outData, dnData)
  if (!is.null(mutData))
    outData <- cbind(outData, mutData)
  if (!is.null(caseData))
    outData <- cbind(outData, caseData)
  if (!is.null(controlData))
    outData <- cbind(outData, controlData)
  
  pp.samples <- matrix(1, dim(outData)[1], dim(mcmc.samples$pi0)[1])
  
  for (i in 1:dim(mcmc.samples$pi0)[1]) {
    bfAll <- rep(1, dim(outData)[1])
    A = mcmc.samples$A[i,]
    B = mcmc.samples$B[i,]
    C = mcmc.samples$C[i,]
    if (!is.na(dim(mcmc.samples$pi0)[2])) {
      pi0 = mcmc.samples$pi0[i,]
    } else {
      pi0 = mcmc.samples$pi0[i]
    }
    
    covariates = pars$covariates
    L = (1-C) * A / (log(exp(A)+exp(B*A)) - log(exp(B*A)+1)) # scale factor for the distribution
    
    bfDN <- matrix(1, nrow = dim(dnData)[1], ncol = dim(dnData)[2]) ##De novo bayes factors
    for (j in 1:dim(bfDN)[2]) {
      e.hyperGammaMeanDN <- mcmc.samples$hyperGammaMeanDN[i, j]
      e.hyperBetaDN <- mcmc.samples$hyperBetaDN[i, j]
      e.bf <- bayes.factor.denovo(x =  dnData[, j],
                                  N = pars$nfamily[j],
                                  mu =  mutData[, j],
                                  gamma.mean = e.hyperGammaMeanDN,
                                  beta = e.hyperBetaDN)
      bfDN[, j] <- e.bf
    }
    tmpBF <- apply(bfDN, 1, prod)
    
    bfcovariates <- matrix(1, nrow = dim(dnData)[1], ncol = 1)
    for (ii in 1:dim(dnData)[1]) {
      for (jj in 1:dim(covariates)[2]) {
        bfcovariates[ii, 1] = bfcovariates[ii, 1] * 
          (C[jj] + L[jj] / (1+exp(-A[jj]*(covariates[ii, jj]-B[jj])))) / 
          (1 + pi0/(1-pi0)*(1-C[jj]-L[jj]/(1+exp(-A[jj]*(covariates[ii, jj]-B[jj])))))
      }
    }
    tmpBF <- tmpBF * as.vector(bfcovariates)
    
    pp.samples[, i] <- (tmpBF*pi0)/(1 - pi0 + pi0*tmpBF)
  }
  
  outData$PP <- apply(pp.samples, 1, mean)
  outData <- outData[order(-outData$PP),]
  outData$qvalue <- Bayesian.FDR.full.posterior(outData$PP)$FDR
  
  colnames(outData)[1] <- "Gene"
  return(outData)
}

