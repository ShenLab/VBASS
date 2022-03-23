# simulation strategy: divide genes into 20 sets by mutation rate. Let the first 
# and last be disease risk genes. This is based on the knowledge that extTADA has
# low power on low mutation rate genes.
simulate_AutoEncoder <- function(effect_size_mean = c(16.6, 6.8),
                                 seed = 8
) {
  # do random sample on real expression data
  set.seed(seed)
  output.folder <- paste0('AutoEncoderTADA.input.simulation.realexp/')
  dir.create(output.folder)
  # first simulate expression data
  genes.x <- read.csv('AutoEncoderTADA.input/input.x.csv', row.names = 1)
  genes.x_var <- read.csv('AutoEncoderTADA.input/input.x_var.csv', row.names = 1)
  reference <- genes.x_var[,c(81,82)]
  colnames(reference) <- c('mut.cls1', 'mut.cls2')
  cell_number <- genes.x_var[,1:80]
  
  exp_mat <- genes.x[,1:80]
  exp_mat.mean <- exp_mat / cell_number
  
  # get risk genes
  # first transform mean of exp_mat to non-linear
  # use real model output pi to simulate
  real.pi <- read.csv('AutoEncoderTADA.input/all.prediction.csv', row.names = 1)
  
  pi <- (real.pi$X0)^2
  pi <- (real.pi$X0)
  SPARK.RDS <- readRDS('SPARK_extTADA.WES1.RDS')
  SPARK.pars0 <- summary(SPARK.RDS$mcmcDD)
  pi <- pi * SPARK.pars0$summary[1, 1] / mean(pi)
  # next sum and rank
  # exp_mat.sum.rank <- dplyr::percent_rank(rowSums(exp_mat.mean.non.linear))
  # define risk genes based on that
  # repeat 10 times
  set.seed(seed)
  indexes <- seq(1,dim(reference)[1])
  risk.genes <- which(runif(length(pi))<=pi)
  null.genes <- indexes[!indexes %in% risk.genes]
  
  input.x.all <- c()
  input.x_var.all <- c()
  input.x.label.all <- c()
  test.x.label.all <- c()
  for (r in 1:50) {
    # indexes <- seq(1,dim(reference)[1])
    # indexes_pi <- indexes[order(pi, decreasing = T)]
    # indexes_pi_rev <- indexes[order(pi)]
    dnv_table <- data.frame(row.names = rownames(reference),
                            dn.cls1 = rep(0, dim(reference)[1]),
                            dn.cls2 = rep(0, dim(reference)[1]))
    
    dn_risk_LGD <- rnbinom(length(risk.genes), effect_size_mean[1],
                           1/(1+2*reference$mut.cls1[risk.genes]))
    dn_risk_Dmis <- rnbinom(length(risk.genes), effect_size_mean[2],
                            1/(1+2*reference$mut.cls2[risk.genes]))
    
    dn_null_LGD <- rpois(length(null.genes), 2*reference$mut.cls1[null.genes])
    dn_null_Dmis <- rpois(length(null.genes), 2*reference$mut.cls2[null.genes])
    
    dnv_table$dn.cls1 <- rep(NA, dim(reference)[1])
    dnv_table$dn.cls2 <- rep(NA, dim(reference)[1])
    
    dnv_table$dn.cls1[risk.genes] <- dn_risk_LGD
    dnv_table$dn.cls1[null.genes] <- dn_null_LGD
    dnv_table$dn.cls2[risk.genes] <- dn_risk_Dmis
    dnv_table$dn.cls2[null.genes] <- dn_null_Dmis
    # combine the final outputs
    input.x <- cbind(exp_mat, dnv_table)
    input.x_var <- cbind(cell_number, reference)
    input.x.label <- data.frame(row.names = rownames(reference),
                                label = rep(NA, dim(reference)[1]))
    test.x.label <- data.frame(row.names = rownames(reference),
                               label = rep(NA, dim(reference)[1]))
    # set training genes by mutation rate
    # risk.genes.mut.sum <- reference$mut.cls1[risk.genes] + reference$mut.cls2[risk.genes]
    # null.genes.mut.sum <- reference$mut.cls1[null.genes] + reference$mut.cls2[null.genes]
    risk.genes.pi <- pi[risk.genes]
    null.genes.pi <- pi[null.genes]
    # get top risk genes and null genes
    training.pos.genes <- risk.genes[order(-risk.genes.pi)]
    training.neg.genes <- null.genes[order(null.genes.pi)]
    training.pos.genes <- training.pos.genes[1:100]
    training.neg.genes <- training.neg.genes[1:300]
    input.x.label$label[training.pos.genes] <- 1
    input.x.label$label[training.neg.genes] <- 0
    test.x.label$label[risk.genes] <- 1
    test.x.label$label[null.genes] <- 0
    
    input.x.all <- rbind(input.x.all, input.x)
    input.x_var.all <- rbind(input.x_var.all, input.x_var)
    input.x.label.all <- rbind(input.x.label.all, input.x.label)
    test.x.label.all <- rbind(test.x.label.all, test.x.label)
    write.csv(input.x, file = paste0(output.folder, 'input.x.', r, '.csv'))
    write.csv(input.x_var, file = paste0(output.folder, 'input.x_var.', r, '.csv'))
    write.csv(input.x.label, file = paste0(output.folder, 'input.label.', r, '.csv'))
    write.csv(test.x.label, file = paste0(output.folder, 'test.label.', r, '.csv'))
    
    real.log.BF <- log(pi) + 
      dnbinom(dnv_table$dn.cls1, effect_size_mean[1], 1/(1+2*reference$mut.cls1), log = TRUE) +
      dnbinom(dnv_table$dn.cls2, effect_size_mean[2], 1/(1+2*reference$mut.cls2), log = TRUE) -
      dpois(dnv_table$dn.cls1, 2*reference$mut.cls1, log = TRUE) -
      dpois(dnv_table$dn.cls2, 2*reference$mut.cls2, log = TRUE) -
      log(1-pi)
    real.log.BF <- data.frame(real.log.BF, row.names = rownames(dnv_table))
    write.csv(real.log.BF, file = paste0(output.folder, 'real.log.BF.', r, '.csv'))
  }
  
  dir.create(output.folder)
  write.csv(input.x.all, file = paste0(output.folder, 'input.x.csv'))
  write.csv(input.x_var.all, file = paste0(output.folder, 'input.x_var.csv'))
  write.csv(input.x.label.all, file = paste0(output.folder, 'input.label.csv'))
  write.csv(test.x.label.all, file = paste0(output.folder, 'test.label.csv'))
  
  # write.csv(dnv_table, file = paste0(output.folder, 'dnv.table.csv'))
  # write.csv(mutRate, file = paste0(output.folder, 'mut.rate.csv'))
  # effect.size <- data.frame(LGD=effect_size_all_LGD, Dmis=effect_size_all_Dmis)
  # write.csv(effect.size, file = paste0(output.folder, "effect.size.csv"))
  
  result <- list(input.x, input.x_var, input.x.label, training.pos.genes, training.neg.genes)
}
