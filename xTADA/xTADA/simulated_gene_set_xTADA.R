sim_gene_set_exttada <- function(dnv_table, mutRate, covariates_array, gene_set, samplenumber) {
  geneset = gene_set
  dataDN1 = dnv_table # dnv_table
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
  covariates = matrix(covariates_array, nrow = length(covariates_array), ncol = 1)
  #covariates = covariates_array
  
  options(mc.cores = parallel::detectCores())
  mcmcDD <- extTADAmcmc(modelName = DNextTADA,
                        dataDN = dataDN1,
                        mutRate = mutRate1,
                        covariates = covariates,
                        Ndn = rep(samplenumber, 2),
                        nIteration = 2000,
                        nChain = 4,
                        nCore = 4)
  # mcmcDD <- readRDS('tmp.RDS')

  mcmcDataFrame <- summary(mcmcDD)$summary
  # use mean as estimation, not mode
  pars0 = mcmcDataFrame[grep("hyper|pi|A|B|C", row.names(mcmcDataFrame)), ]

  # pars0 = estimatePars(pars = mcmcDD@sim$fnames_oi,
  #                      mcmcResult = mcmcDD)
  parsFDR <- list(A = as.numeric(pars0[grep('A', rownames(pars0)),1]),
                  B = as.numeric(pars0[grep('B', rownames(pars0)),1]),
                  C = as.numeric(pars0[grep('C', rownames(pars0)),1]),
                  pi0 = as.numeric(pars0[grep('pi0', rownames(pars0)),1]),
                  gammaMeanDN = as.numeric(pars0[grep('hyperGammaMeanDN', rownames(pars0)),1]),
                  betaDN = as.numeric(pars0[grep('hyperBetaDN', rownames(pars0)),1]),
                  nfamily = rep(samplenumber, 2),
                  covariates = covariates)

  dataFDR <- calculateFDR(pars = parsFDR,
                          dnData = dataDN1,
                          mutData = mutRate1,
                          geneName = geneset)

  HGNC = reference[match(dataFDR$Gene, reference$GeneID),"HGNC"]
  covariates_reordered <- covariates[match(dataFDR$Gene, geneset),]
  dataFDR <- cbind(dataFDR, HGNC, covariates_reordered)

  dataFDR.full.posterior <- calculateFDR.full.posterior(pars = list(mcmc.samples=extract(mcmcDD),
                                                                    nfamily = rep(samplenumber, 2),
                                                                    covariates = covariates),
                                                        dnData = dataDN1,
                                                        mutData = mutRate1,
                                                        geneName = geneset)
  HGNC = reference[match(dataFDR.full.posterior$Gene, reference$GeneID),"HGNC"]
  covariates_reordered <- covariates[match(dataFDR.full.posterior$Gene, geneset),]
  dataFDR.full.posterior <- cbind(dataFDR.full.posterior, HGNC, covariates_reordered)

  result = list("mcmcDD"=mcmcDD, "pars0"=pars0,
                "dataFDR"=dataFDR, "dataFDR.full.posterior"=dataFDR.full.posterior)
  return(result)
}
