xInputGS.raw <- read.csv("~/ld1/VBASS_xTADA_paper/scalar_VBASS/result/table.S4.csv", header = TRUE, as.is = TRUE)
xInputGS <- xInputGS.raw[,c("Gene", "dn_LGD", "dn_Dmis", "mut_LGD", "mut_Dmis")]
xInputGS$GS <- 0
# select top 25% as the geneset
xInputGS$GS[xInputGS$covariates_reordered >= 0.75] <- 1
xInputGS$covariates_reordered <- NULL

head(xInputGS)

##SET PARAMETERS for DECO
nIteration = 1000 ###This should be higher (e.g, >=5,000)
nChain = 2
nCore = nChain
ntrio = c(2645, 2645) #Trio numbers: there are three categories of de novo mutations
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

saveRDS(out, file = 'CHD.RDS')
