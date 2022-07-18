gene_set_scalarVBASS <- function(geneset, cases, samplenumber, reference, covariatesData = NA, cal_pval = FALSE, outputpath=NA) {
  # form a proper input for TADA
  my.data = data.frame(gene.id=geneset,
                       mut.cls0=numeric(length(geneset)),
                       mut.cls1=numeric(length(geneset)),
                       mut.cls2=numeric(length(geneset)),
                       dn.cls0=numeric(length(geneset)),
                       dn.cls1=numeric(length(geneset)),
                       dn.cls2=numeric(length(geneset)))
  for (i in 1:length(geneset)) {
    index = match(as.character(geneset[i]), reference$GeneID)
    # check whether we found the gene
    if (!is.na(index)) {
      my.data$mut.cls0[i] = reference$Mu_Silent[index]
      my.data$mut.cls1[i] = reference$Mu_LoF[index]
      my.data$mut.cls2[i] = reference$Mu_Dmis_REVEL0.5[index]
      
    } else {
      my.data$mut.cls0[i] = .Machine$double.eps
      my.data$mut.cls1[i] = .Machine$double.eps
      my.data$mut.cls2[i] = .Machine$double.eps
      
    }
    if (my.data$mut.cls0[i] <= 0) {
      my.data$mut.cls0[i] = .Machine$double.eps
    }
    if (my.data$mut.cls1[i] <= 0) {
      my.data$mut.cls1[i] = .Machine$double.eps
    }
    if (my.data$mut.cls2[i] <= 0) {
      my.data$mut.cls2[i] = .Machine$double.eps
    }
    gene_case = cases[cases$GeneID==as.character(geneset[i]),]
    
    my.data$dn.cls0[i] = sum(gene_case$vclass=="syn")
    my.data$dn.cls1[i] = sum(gene_case$vclass=="LGD")
    my.data$dn.cls2[i] = sum(gene_case$vclass=="Dmis")
    
  }
  # first do sanity check using syn mutation
  message("sanity check before correction")
  syn_oe1 = sum(my.data$dn.cls0) / (sum(my.data$mut.cls0) * samplenumber * 2)
  message(paste0("syn_oe1 = ", syn_oe1))
  message("adjusting mutation rate according to syn_oe")
  constant1 = sum(my.data$dn.cls0) / (sum(my.data$mut.cls0) * samplenumber * 2)
  my.data[,2:4]=my.data[,2:4]*constant1
  
  syn_oe1 = sum(my.data$dn.cls0) / (sum(my.data$mut.cls0) * samplenumber * 2)
  LGD_oe1 = sum(my.data$dn.cls1) / (sum(my.data$mut.cls1) * samplenumber * 2)
  Dmis_oe1 = sum(my.data$dn.cls2) / (sum(my.data$mut.cls2) * samplenumber * 2)
  
  message("sanity check after correction")
  message(paste0("syn_oe1 = ", syn_oe1))
  message(paste0("LGD_oe1 = ", LGD_oe1))
  message(paste0("Dmis_oe1 = ", Dmis_oe1))
  # message("output geneset datafile")
  # write.table(my.data, file = "my.data.csv", quote = FALSE, sep = ",", row.names = FALSE)
  
  # use extTADA to sample values of parameters
  dataDN = my.data[,c('dn.cls1','dn.cls2')]
  colnames(dataDN) = c('dn_LGD','dn_Dmis')
  mutRate = my.data[,c('mut.cls1','mut.cls2')]
  colnames(mutRate) = c('mut_LGD','mut_Dmis')
  
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
  mcmcDD <- scalar_VBASS_mcmc(modelName = DNscalarVBASS,
                              dataDN = dataDN,
                              mutRate = mutRate,
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
                          dnData = dataDN,
                          mutData = mutRate,
                          geneName = geneset)
  
  HGNC = reference[match(dataFDR$Gene, reference$GeneID),"HGNC"]
  covariates_reordered <- covariates[match(dataFDR$Gene, geneset),]
  dataFDR <- cbind(dataFDR, HGNC, covariates_reordered)
  
  dataFDR.full.posterior <- calculateFDR.full.posterior(pars = list(mcmc.samples=extract(mcmcDD),
                                                                    nfamily = rep(samplenumber, 2),
                                                                    covariates = covariates),
                                                        dnData = dataDN,
                                                        mutData = mutRate,
                                                        geneName = geneset)
  HGNC = reference[match(dataFDR.full.posterior$Gene, reference$GeneID),"HGNC"]
  covariates_reordered <- covariates[match(dataFDR.full.posterior$Gene, geneset),]
  dataFDR.full.posterior <- cbind(dataFDR.full.posterior, HGNC, covariates_reordered)
  
  result = list("mcmcDD"=mcmcDD, "pars0"=pars0,
                "dataFDR"=dataFDR, "dataFDR.full.posterior"=dataFDR.full.posterior)
  return(result)
}
