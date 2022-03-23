gene_set_xtada <- function(geneset, cases, samplenumber, reference, covariatesData = NA, cal_pval = FALSE, outputpath=NA) {
  # form a proper input for TADA
  mytada.data1 = data.frame(gene.id=geneset,
                            mut.cls0=numeric(length(geneset)),
                            mut.cls1=numeric(length(geneset)),
                            mut.cls2=numeric(length(geneset)),
                            dn.cls0=numeric(length(geneset)),
                            dn.cls1=numeric(length(geneset)),
                            dn.cls2=numeric(length(geneset)))
  cases1 = cases
  samplenumber1 = samplenumber
  for (i in 1:length(geneset)) {
    index = match(as.character(geneset[i]), reference$GeneID)
    # check whether we found the gene
    if (!is.na(index)) {
      mytada.data1$mut.cls0[i] = reference$Mu_Silent[index]
      mytada.data1$mut.cls1[i] = reference$Mu_LoF[index]
      mytada.data1$mut.cls2[i] = reference$Mu_Dmis_REVEL0.5[index]
      
    } else {
      mytada.data1$mut.cls0[i] = .Machine$double.eps
      mytada.data1$mut.cls1[i] = .Machine$double.eps
      mytada.data1$mut.cls2[i] = .Machine$double.eps
      
    }
    if (mytada.data1$mut.cls0[i] <= 0) {
      mytada.data1$mut.cls0[i] = .Machine$double.eps
    }
    if (mytada.data1$mut.cls1[i] <= 0) {
      mytada.data1$mut.cls1[i] = .Machine$double.eps
    }
    if (mytada.data1$mut.cls2[i] <= 0) {
      mytada.data1$mut.cls2[i] = .Machine$double.eps
    }
    gene_cases1 = cases1[cases1$GeneID==as.character(geneset[i]),]
    
    mytada.data1$dn.cls0[i] = sum(gene_cases1$vclass=="syn")
    mytada.data1$dn.cls1[i] = sum(gene_cases1$vclass=="LGD")
    mytada.data1$dn.cls2[i] = sum(gene_cases1$vclass=="Dmis")
    
  }
  # first do sanity check using syn mutation
  message("sanity check before correction")
  syn_oe1 = sum(mytada.data1$dn.cls0) / (sum(mytada.data1$mut.cls0) * samplenumber1 * 2)
  message(paste0("syn_oe1 = ", syn_oe1))
  message("adjusting mutation rate according to syn_oe")
  constant1 = sum(mytada.data1$dn.cls0) / (sum(mytada.data1$mut.cls0) * samplenumber1 * 2)
  mytada.data1[,2:4]=mytada.data1[,2:4]*constant1
  
  syn_oe1 = sum(mytada.data1$dn.cls0) / (sum(mytada.data1$mut.cls0) * samplenumber1 * 2)
  LGD_oe1 = sum(mytada.data1$dn.cls1) / (sum(mytada.data1$mut.cls1) * samplenumber1 * 2)
  Dmis_oe1 = sum(mytada.data1$dn.cls2) / (sum(mytada.data1$mut.cls2) * samplenumber1 * 2)
  
  message("sanity check after correction")
  message(paste0("syn_oe1 = ", syn_oe1))
  message(paste0("LGD_oe1 = ", LGD_oe1))
  message(paste0("Dmis_oe1 = ", Dmis_oe1))
  # message("output geneset datafile")
  # write.table(mytada.data1, file = "mytada.data1.csv", quote = FALSE, sep = ",", row.names = FALSE)
  
  # use extTADA to sample values of parameters
  dataDN1 = mytada.data1[,c('dn.cls1','dn.cls2')]
  colnames(dataDN1) = c('dn_LGD','dn_Dmis')
  mutRate1 = mytada.data1[,c('mut.cls1','mut.cls2')]
  colnames(mutRate1) = c('mut_LGD','mut_Dmis')
  
  covariates = matrix(NA, nrow = length(geneset), ncol = dim(covariatesData)[2]-1)
  colnames_covariates = rep(NA, dim(covariatesData)[2]-1)
  k = 1
  for (i in colnames(covariatesData)) {
    if (i != "gene") {
      covariates_tmp = covariatesData[match(geneset, covariatesData$gene), i]
      covariates_tmp = dplyr::percent_rank(covariates_tmp)
      covariates_tmp[is.na(covariates_tmp)] = 0.5
      covariates[, k] = covariates_tmp
      colnames_covariates[k] = i
      k = k + 1
    }
  }
  colnames(covariates) = colnames_covariates
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
