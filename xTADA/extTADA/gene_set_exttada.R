gene_set_exttada <- function(geneset, cases, samplenumber, reference, cal_pval = FALSE, outputpath=NA) {
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
  message("output geneset datafile")
  write.table(mytada.data1, file = "mytada.data1.csv", quote = FALSE, sep = ",", row.names = FALSE)
  
  # use extTADA to sample values of parameters
  dataDN1 = mytada.data1[,c('dn.cls1','dn.cls2')]
  colnames(dataDN1) = c('dn_LGD','dn_Dmis')
  mutRate1 = mytada.data1[,c('mut.cls1','mut.cls2')]
  colnames(mutRate1) = c('mut_LGD','mut_Dmis')
  print(class(dataDN1))
  print(head(dataDN1))
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