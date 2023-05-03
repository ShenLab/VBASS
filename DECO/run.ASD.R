xInputGS <- read.table("data/extTADA_SCZ_constrainedGenes.csv", header = TRUE, as.is = TRUE)
head(xInputGS)

##SET PARAMETERS for DECO
nIteration = 1000 ###This should be higher (e.g, >=5,000)
nChain = 2
nCore = nChain
ntrio = c(1077, 1077, 1077) #Trio numbers: there are three categories of de novo mutations
ncase = c(3157, 1091, 1353) #Case numbers: there are three population samples of rare case/control variants
ncontrol = c(4672, 1193, 4769) #Control numbers: there are three population samples of rare case/control variants

###RUN DECO
source('script/DECO.R') #Load the source code

###Obtain results
outSCZ <- DECO(
  inputData = xInputGS, ## Input data should be formatted as above
  Ndn = array(c(ntrio)), #rep(ntrio, 1), ##Two de novo categories
  Ncase = array(ncase), #rep(N$ca, 1), ##Two case categories
  Ncontrol = array(ncontrol), #rep(N$cn, 1), ##Two control categories
  nIteration = nIteration,# nIteration ## Number of iterations: should be upto higher for real data
  nChain = nChain, #Number of MCMC chains
  nCore = nCore #Number of computer cores
)

saveRDS(out, file = 'ASD.GO.RDS')
