# simulation strategy: divide genes into 20 sets by mutation rate. Let the first 
# and last be disease risk genes. This is based on the knowledge that extTADA has
# low power on low mutation rate genes.
simulate_xTADA <- function(reference, samplesize, 
                           effect_size_mean = c(15, 10),
                           effect_size_var = c(1, 1),
                           pi = 0.05,
                           A = 80, B = 0.75, C = 0.1, seed=1, repeat.size=1) {
  set.seed(1)
  genes <- reference$GeneID
  reference <- reference[rep(1:dim(reference)[1], repeat.size),]
  reference$GeneID <- paste0(genes, ".", rep(1:repeat.size, each=length(genes), times=1))
  # simulate dnv table
  sum_mut_rate <- reference$Mu_LoF + reference$Mu_Dmis_REVEL0.5
  reference <- reference[order(sum_mut_rate),]
  dnv_table <- data.frame(row.names = reference$GeneID,
                          dn.cls1 = rep(0, dim(reference)[1]),
                          dn.cls2 = rep(0, dim(reference)[1]))
  mutRate <- data.frame(row.names = reference$GeneID,
                        mut.cls1 = reference$Mu_LoF,
                        mut.cls2 = reference$Mu_Dmis_REVEL0.5)
  indexes <- seq(1,dim(reference)[1])
  risk.genes <- sample(indexes, floor(dim(reference)[1])*pi)
  null.genes <- indexes[!indexes %in% risk.genes]
  set.seed(seed)
  dn_risk_LGD <- rnbinom(length(risk.genes), effect_size_mean[1]*effect_size_var[1],
                         effect_size_var[1]/(effect_size_var[1]+2*samplesize*mutRate$mut.cls1[risk.genes]))
  dn_risk_Dmis <- rnbinom(length(risk.genes), effect_size_mean[2]*effect_size_var[2],
                          effect_size_var[2]/(effect_size_var[2]+2*samplesize*mutRate$mut.cls2[risk.genes]))
  
  dn_null_LGD <- rpois(length(null.genes), mutRate$mut.cls1[null.genes] * samplesize * 2)
  dn_null_Dmis <- rpois(length(null.genes), mutRate$mut.cls2[null.genes] * samplesize * 2)
  
  dnv_table$dn.cls1 <- rep(NA, dim(reference)[1])
  dnv_table$dn.cls2 <- rep(NA, dim(reference)[1])
  
  dnv_table$dn.cls1[risk.genes] <- dn_risk_LGD
  dnv_table$dn.cls1[null.genes] <- dn_null_LGD
  dnv_table$dn.cls2[risk.genes] <- dn_risk_Dmis
  dnv_table$dn.cls2[null.genes] <- dn_null_Dmis
  
  # simulate covariates data
  L = (1-C)*A/(log(exp(A)+exp(A*B))-log(exp(A*B)+1))
  covariates_risk_gene <- c()
  while (length(covariates_risk_gene) < length(risk.genes)) {
    x <- runif(length(risk.genes)-length(covariates_risk_gene), min = 0, max = 1)
    accept_array <- runif(length(x), min = 0, max = C+L/(1+exp(-A*(1-B))))
    y <- C + L / (1+exp(-A*(x-B)))
    covariates_risk_gene <- c(covariates_risk_gene, x[accept_array <= y])
  }
  covariates_null_gene <- c()
  while (length(covariates_null_gene) < length(null.genes)) {
    x <- runif(length(null.genes)-length(covariates_null_gene), min = 0, max = 1)
    accept_array <- runif(length(x), min = 0, max = 1+0.05/0.95*(1-C-L/(1+exp(A*B))))
    y <- 1 + 0.05/0.95 * (1 - C - L / (1+exp(-A*(x-B))))
    covariates_null_gene <- c(covariates_null_gene, x[accept_array <= y])
  }
  covariates_array <- rep(1, dim(reference)[1])
  covariates_array[risk.genes] <- covariates_risk_gene
  covariates_array[null.genes] <- covariates_null_gene
  risk.gene.IDs <- reference$GeneID[risk.genes]
  result <- list(dnv_table, mutRate, covariates_array, risk.gene.IDs)
}

