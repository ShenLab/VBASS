gene_set_exttada_single <- function(geneset, cases, muttype, samplenumber, reference, cal_pval = FALSE, outputpath=NA) {
  # colnames of cases should be c('genes', 'classes', ...)
  # colnames of reference should be c('genes', 'Mu_rate', 'Mu_rate_silent')
  mytada.data = data.frame(gene.id=geneset,
                           mut=numeric(length(geneset)),
                           dn=numeric(length(geneset)),
                           mut.silent=numeric(length(geneset)),
                           dn.silent=numeric(length(geneset)))
  for (i in 1:length(geneset)) {
    index = match(as.character(geneset[i]), reference$genes)
    # check whether we found the gene
    if (!is.na(index)) {
      mytada.data$mut[i] = reference$Mu_rate[index]
      mytada.data$mut.silent[i] = reference$Mu_rate_silent[index]
    } else {
      mytada.data$mut[i] = .Machine$double.eps
      mytada.data$mut.silent[i] = .Machine$double.eps
    }
    if (mytada.data$mut[i] <= 0) {
      mytada.data$mut[i] = .Machine$double.eps
    }
    if (mytada.data$mut.silent[i] <= 0) {
      mytada.data$mut.silent[i] = .Machine$double.eps
    }
    gene_cases = cases[cases$genes==as.character(geneset[i]),]
    
    if (muttype=="Dmis") {
      mytada.data$dn[i] = sum(gene_cases$classes=="mis"
                              & gene_cases$REVEL>=0.5
                              & !is.na(gene_cases$REVEL))
    } else {
      mytada.data$dn[i] = sum(gene_cases$classes==muttype)
    }
    mytada.data$dn.silent[i] = sum(gene_cases$classes=="syn")
  }
  # first do sanity check using syn mutation
  message("sanity check before correction")
  syn_oe = sum(mytada.data$dn.silent) / (sum(mytada.data$mut.silent) * samplenumber * 2)
  message(paste0("syn_oe = ", syn_oe))
  message("adjusting mutation rate according to syn_oe")
  constant = sum(mytada.data$dn.silent) / (sum(mytada.data$mut.silent) * samplenumber * 2)
  mytada.data[,c(2,4)]=mytada.data[,c(2,4)]*syn_oe
  
  syn_oe = sum(mytada.data$dn.silent) / (sum(mytada.data$mut.silent) * samplenumber * 2)
  muttype_oe = sum(mytada.data$dn) / (sum(mytada.data$mut) * samplenumber * 2)
  
  message("sanity check after correction")
  message(paste0("syn_oe = ", syn_oe))
  message(paste0(muttype, "_oe = ", muttype_oe))
  
  # use extTADA to sample values of parameters
  dataDN = data.frame(mytada.data[,c('dn')])
  colnames(dataDN) = c(paste0('dn_', muttype))
  mutRate = data.frame(mytada.data[,c('mut')])
  colnames(mutRate) = c(paste0('mut_', muttype))
  
  options(mc.cores = parallel::detectCores())
  mcmcDD <- extTADAmcmc(modelName = DNextTADA,
                        dataDN = dataDN,
                        mutRate = mutRate,
                        Ndn = rep(samplenumber, 1),
                        nIteration = 1000)
  
  # options(repr.plot.width = 4, repr.plot.height = 3)
  # par(mfrow = c(1,2))
  # plotParHeatmap(c("pi0[1]", "hyperGammaMeanDN[1]"), mcmcResult = mcmcDD)
  # plotParHeatmap(c("pi0[2]", "hyperGammaMeanDN[2]"), mcmcResult = mcmcDD)
  
  pars0 = estimateParsExtTADA(pars = c('pi0',
                                       'hyperGammaMeanDN[1]',
                                       'hyperBetaDN[1]'),
                              mcmcResult = mcmcDD)
  
  parsFDR <- list(gammaMeanDN = pars0[,1][2],
                  betaDN = pars0[,1][3],
                  pi0 = pars0[,1][1],
                  nfamily = rep(samplenumber, 1))
  
  dataFDR <- calculateFDR(pars = parsFDR,
                          dnData = dataDN,
                          mutData = mutRate,
                          geneName = geneset)
  result = list("mcmcDD"=mcmcDD, "pars0"=pars0, "dataFDR"=dataFDR)
  return(result)
}