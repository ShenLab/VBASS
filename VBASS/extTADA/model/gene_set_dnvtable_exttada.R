gene_set_dnvtable_exttada <- function(geneset, dnv_table, samplenumber, reference, cal_pval = FALSE, outputpath=NA) {
  # form a proper input for TADA
  mytada.data1 = dnv_table[match(geneset, rownames(dnv_table)),]
  # mytada.data1 <- dnv_table
  dataDN1 = mytada.data1[,c('dn.cls1','dn.cls2')]
  colnames(dataDN1) = c('dn_LGD','dn_Dmis')
  mutRate1 = mytada.data1[,c('mut.cls1','mut.cls2')]
  colnames(mutRate1) = c('mut_LGD','mut_Dmis')
  
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
  mcmcDataFrame <- summary(mcmcDD)$summary
  # use mean as estimation, not mode
  pars0 = mcmcDataFrame[grep("hyper|pi", row.names(mcmcDataFrame)), ]
  # pars0 = estimateParsExtTADA(pars = c('pi0',
  #                                      'hyperGammaMeanDN[1]', 'hyperGammaMeanDN[2]',
  #                                      'hyperBetaDN[1]', 'hyperBetaDN[2]'),
  #                             mcmcResult = mcmcDD)
  
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