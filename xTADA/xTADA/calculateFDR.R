calculateFDR <- function(pars,
                         caseData = NULL,
                         controlData = NULL,
                         dnData = NULL,
                         mutData = NULL,
                         geneName){
  
  outData <- data.frame(geneName)
  if (!is.null(dnData))
    outData <- cbind(outData, dnData)
  if (!is.null(mutData))
    outData <- cbind(outData, mutData)
  if (!is.null(caseData))
    outData <- cbind(outData, caseData)
  if (!is.null(controlData))
    outData <- cbind(outData, controlData)
  
  
  bfAll <- rep(1, dim(outData)[1])
  A = pars$A
  B = pars$B
  C = pars$C
  covariates = pars$covariates
  L = (1-C) * A / (log(exp(A)+exp(B*A)) - log(exp(B*A)+1)) # scale factor for the distribution
  
  if ( length(pars$gammaMeanDN) == 0) {
    message("No parameters for de novo data; therefore, these categories are not calculated in this step.\n")
  }  else {
    bfDN <- matrix(1, nrow = dim(dnData)[1], ncol = dim(dnData)[2]) ##De novo bayes factors
    for (j2 in 1:dim(bfDN)[2]) {
      e.hyperGammaMeanDN <- pars$gammaMeanDN[j2]
      e.hyperBetaDN <- pars$betaDN[j2]
      e.bf <- bayes.factor.denovo(x =  dnData[, j2],
                                  N = pars$nfamily[j2],
                                  mu =  mutData[, j2],
                                  gamma.mean = e.hyperGammaMeanDN,
                                  beta = e.hyperBetaDN)
      bfDN[, j2] <- e.bf
    }
    bfAll <- bfAll*apply(bfDN, 1, prod)
  }
  
  if (length(pars$gammaMeanCC) == 0) {
    message("No parameters for case-control data;  therefore, these categories are not calculated in this step.\n")
  } else {
    
    bfCC <- matrix(1, ncol = dim(caseData)[2], nrow = dim(caseData)[1])
    for (cc3 in 1:dim(bfCC)[2]){
      e.hyperGammaMeanCC <- pars$gammaMeanCC[cc3]
      e.hyperBetaCC <- pars$betaCC[cc3]
      e.nu <- 200
      t.case <- caseData[, cc3]
      t.control <- controlData[, cc3]
      e.rho <- e.nu*mean(t.case + t.control)/(pars$ncase[cc3] + pars$ncontrol[cc3])
      e.bf <- BayesFactorCC3(Nsample = list(ca = pars$ncase[cc3], cn = pars$ncontrol[cc3]),
                             x.case = t.case, x.control = t.control,
                             gamma.meanCC = e.hyperGammaMeanCC, betaCC = e.hyperBetaCC,
                             rhoCC = e.rho, nuCC = e.nu)
      bfCC[, cc3] <- e.bf
    }
    
    bfAll <- bfAll*apply(bfCC, 1, prod)
  }
  bfAlls = bfAll
  bfcovariates <- matrix(1, nrow = dim(dnData)[1], ncol = 1)
  for (ii in 1:dim(dnData)[1]) {
    for (jj in 1:dim(covariates)[2]) {
      bfcovariates[ii, 1] = bfcovariates[ii, 1] * 
        (C[jj] + L[jj] / (1+exp(-A[jj]*(covariates[ii, jj]-B[jj])))) / 
        (1 + pars$pi0/(1-pars$pi0)*(1-C[jj]-L[jj]/(1+exp(-A[jj]*(covariates[ii, jj]-B[jj])))))
    }
  }
  
  bfAlls <- bfAlls * as.vector(bfcovariates)
  
  outData$BF <- bfAlls
  outData$BFAlls <- bfAlls
  #    outData$BFdn <- bfDN #cbind(outData, bfDN)
  
  tempPP <- pars$pi0*bfAlls
  outData$PP <- tempPP/(1-pars$pi0 + tempPP)
  
  outData <- outData[order(-outData$PP),]
  
  outData$qvalue <- Bayesian.FDR(outData$PP)$FDR
  
  # pvalues <- matrix(1, nrow = dim(dnData)[1], ncol = 2^dim(dnData)[2]-1)
  # row.names(pvalues) <- geneName
  # for (i in 1:(2^dim(dnData)[2]-1)) {
  #   pvalues[,i] <- tempPP[,i]/(pars$pi0[2^dim(dnData)[2]] + apply(tempPP, 1, sum))
  # }
  # qvalues <- pvalues
  # for (i in 1:(2^dim(dnData)[2]-1)) {
  #   qvalues <- qvalues[order(-qvalues[,i]),]
  #   qvalues[,i] <- Bayesian.FDR(qvalues[,i])$FDR
  # }
  # finalqvalues <- qvalues[match(outData$geneName, row.names(qvalues)),]
  colnames(outData)[1] <- "Gene"
  # outData$qvalues <- finalqvalues
  return(outData)
}

