sim_gene_set_exttada <- function(dnv_table, mutRate, gene_set, samplenumber) {
  # use extTADA to sample values of parameters
  geneset = gene_set
  dataDN1 = dnv_table
  colnames(dataDN1) = c('dn_LGD','dn_Dmis')
  rownames(dataDN1) <- NULL
  
  mutRate1 = mutRate
  colnames(mutRate1) = c('mut_LGD','mut_Dmis')
  rownames(mutRate1) <- NULL
  mutRate1[mutRate1 == 0] <- .Machine$double.eps
  # covariatesData <- data.frame(geneset, covariates_array)
  # colnames(covariatesData) <- c('gene', 'covariate')
  # covariates = matrix(NA, nrow = length(geneset), ncol = dim(covariatesData)[2]-1)
  # colnames_covariates = rep(NA, dim(covariatesData)[2]-1)
  # 
  # k = 1
  # for (i in colnames(covariatesData)) {
  #   if (i != "gene") {
  #     covariates_tmp = covariatesData[match(geneset, covariatesData$gene), i]
  #     covariates_tmp = dplyr::percent_rank(covariates_tmp)
  #     covariates_tmp[is.na(covariates_tmp)] = 0.5
  #     covariates[, k] = covariates_tmp
  #     colnames_covariates[k] = i
  #     k = k + 1
  #   }
  # }
  # colnames(covariates) = colnames_covariates
  # samplenumber = 2645
  # covariates1 = matrix(covariates_array, nrow = length(covariates_array), ncol = 1)
  # covariates1 = covariates_array
  
  options(mc.cores = parallel::detectCores())
  
  mcmcDD <- extTADAmcmc(modelName = DNextTADA,
                        dataDN = dataDN1,
                        mutRate = mutRate1,
                        Ndn = rep(samplenumber, 2),
                        nCore = 4, nChain = 4,
                        nIteration = 2000)
  # mcmcDD <- readRDS('tmp.RDS')

  options(repr.plot.width = 4, repr.plot.height = 3)
  par(mfrow = c(1,2))
  # plotParHeatmap(c("pi0[1]", "hyperGammaMeanDN[1]"), mcmcResult = mcmcDD)
  # plotParHeatmap(c("pi0[2]", "hyperGammaMeanDN[2]"), mcmcResult = mcmcDD)
  
  pars0 = estimateParsExtTADA(pars = c('pi0',
                                       'hyperGammaMeanDN[1]', 'hyperGammaMeanDN[2]',
                                       'hyperBetaDN[1]', 'hyperBetaDN[2]'),
                              mcmcResult = mcmcDD)

  parsFDR <- list(gammaMeanDN = pars0[,1][2:3],
                  betaDN = pars0[,1][4:5],
                  pi0 = pars0[,1][1],
                  nfamily = rep(samplenumber, 2))
  
  dataFDR <- calculateFDR(pars = parsFDR,
                          dnData = dataDN1,
                          mutData = mutRate1,
                          geneName = geneset)
  
  dataFDR.full.posterior <- calculateFDR.full.posterior(pars = list(mcmc.samples=extract(mcmcDD),
                                                                    nfamily = rep(samplenumber, 2)),
                                                        dnData = dataDN1,
                                                        mutData = mutRate1,
                                                        geneName = geneset)
  
  HGNC = reference[match(dataFDR$Gene, reference$GeneID),"HGNC"]
  dataFDR <- cbind(dataFDR, HGNC)
  
  HGNC = reference[match(dataFDR.full.posterior$Gene, reference$GeneID),"HGNC"]
  dataFDR.full.posterior <- cbind(dataFDR.full.posterior, HGNC)
  
  result = list("mcmcDD"=mcmcDD, "pars0"=pars0,
                "dataFDR"=dataFDR,
                "dataFDR.full.posterior"=dataFDR.full.posterior)
  return(result)
}
