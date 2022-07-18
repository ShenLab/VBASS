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
                                        dnData = NULL,
                                        mutData = NULL,
                                        geneName){
  # a function to calculate FDR by utilizing the full posterior distribution from MCMC
  # only support de novo for now.
  # pars should consist mcmc.samples and nfamily
  mcmc.samples <- pars$mcmc.samples
  
  outData <- data.frame(geneName)
  if (!is.null(dnData))
    outData <- cbind(outData, dnData)
  if (!is.null(mutData))
    outData <- cbind(outData, mutData)

  pp.samples <- matrix(1, dim(outData)[1], dim(mcmc.samples$pi0)[1])
  
  for (i in 1:dim(mcmc.samples$pi0)[1]) {
    bfDN <- matrix(1, nrow = dim(dnData)[1], ncol = dim(dnData)[2]) ##De novo bayes factors
    for (j in 1:dim(bfDN)[2]) {
      e.hyperGammaMeanDN <- mcmc.samples$hyperGammaMeanDN[i, j]
      e.hyperBetaDN <- mcmc.samples$hyperBetaDN[i, j]
      e.bf <- bayes.factor.denovo(x = dnData[, j],
                                  N = pars$nfamily[j],
                                  mu =mutData[, j],
                                  gamma.mean = e.hyperGammaMeanDN,
                                  beta = e.hyperBetaDN)
      bfDN[, j] <- e.bf
    }
    tmpBF <- apply(bfDN, 1, prod)
    pp.samples[, i] <- (tmpBF*mcmc.samples$pi0[i])/(1 - mcmc.samples$pi0[i] + mcmc.samples$pi0[i]*tmpBF)
  }
  
  outData$PP <- apply(pp.samples, 1, mean)
  outData <- outData[order(-outData$PP),]
  outData$qvalue <- Bayesian.FDR.full.posterior(outData$PP)$FDR
  
  colnames(outData)[1] <- "Gene"
  return(outData)
}

