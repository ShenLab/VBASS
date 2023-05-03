args <- commandArgs(trailingOnly = TRUE)
gs.id <- as.numeric(args[1])
print(gs.id)
if (!file.exists(paste0('ASD.GO/', gs.id, '.RDS'))) {
  xInputGS.raw <- read.csv("~/ld1/VBASS_xTADA_paper/VBASS/figs/table.S5.csv", header = TRUE, as.is = TRUE)
  xInputGS <- xInputGS.raw[,c("X", "dn_LGD", "dn_Dmis", "mut_LGD", "mut_Dmis")]
  xInputGS$mut_LGD <- xInputGS$mut_LGD / 16616
  xInputGS$mut_Dmis <- xInputGS$mut_Dmis / 16616
  colnames(xInputGS)[1] <- "Gene"
  go.terms <- as.list(GSA::GSA.read.gmt(paste("~/Data/GSEA/msigdb.v7.0.symbols.gmt", sep = "")))
  to.select <- c()
  for (i in 1:length(go.terms$genesets)) {
    if (length(go.terms$genesets[[i]]) > 5) {
      to.select <- c(to.select, i)
    }
  }
  xInputGS$GS <- as.numeric(xInputGS.raw$X %in% go.terms$genesets[[to.select[gs.id]]])

  ##SET PARAMETERS for DECO
  nIteration = 1000 ###This should be higher (e.g, >=5,000)
  nChain = 2
  nCore = nChain
  ntrio = c(16616, 16616) #Trio numbers: there are three categories of de novo mutations
  ncase = c(0, 0) #Case numbers: there are three population samples of rare case/control variants
  ncontrol = c(0, 0) #Control numbers: there are three population samples of rare case/control variants

  ###RUN DECO
  source('script/DECO.R') #Load the source code

  ###Obtain results
  out <- DECO(
    inputData = xInputGS, ## Input data should be formatted as above
    Ndn = array(c(ntrio)), #rep(ntrio, 1), ##Two de novo categories
    Ncase = array(ncase), #rep(N$ca, 1), ##Two case categories
    Ncontrol = array(ncontrol), #rep(N$cn, 1), ##Two control categories
    nIteration = nIteration,# nIteration ## Number of iterations: should be upto higher for real data
    nChain = nChain, #Number of MCMC chains
    nCore = nCore #Number of computer cores
  )

  saveRDS(out, file = paste0('ASD.GO/', gs.id, '.RDS'))
}
