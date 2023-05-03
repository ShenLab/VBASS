library(ggrepel)
library(ggplot2)
Bayesian.FDR.fromBF <- function(log_BF, alpha=0.05) {
  # convert BFs to PPA (posterior probability of alternative model)
  #pi <- 1-pi0
  #q <- apply(pi*BF, 1, prod)/(apply(t(1 - pi), 1, prod)+apply(pi*BF, 1, prod)) # PPA
  q0 <- 1/(exp(log_BF) + 1) # posterior probability of null model
  
  # the FDR at each PPA cutoff
  n <- length(log_BF)
  FDR <- numeric(n)
  for (i in 1:n) FDR[i] <- sum(q0[1:i]) / i 
  
  # the cutoff
  t <- 1
  while (t <= length(q0) & mean(q0[1:t]) <= alpha) { t <- t+1 }
  return (list(FDR=FDR, ND=t))
}
Bayesian.FDR.fromPP <- function(PP, alpha=0.05) {
  # convert BFs to PPA (posterior probability of alternative model)
  #pi <- 1-pi0
  #q <- apply(pi*BF, 1, prod)/(apply(t(1 - pi), 1, prod)+apply(pi*BF, 1, prod)) # PPA
  q0 <- 1- PP # posterior probability of null model
  
  # the FDR at each PPA cutoff
  n <- length(PP)
  FDR <- numeric(n)
  for (i in 1:n) FDR[i] <- sum(q0[1:i]) / i 
  
  # the cutoff
  t <- 1
  while (t <= length(q0) & mean(q0[1:t]) <= alpha) { t <- t+1 }
  return (list(FDR=FDR, ND=t))
}

train.input <- read.csv('VBASS.input/input.x.csv', row.names = 1)
train.input.var <- read.csv('VBASS.input/input.x_var.csv', row.names = 1)
train.labels <- read.csv('VBASS.input/input.label.csv')
train.positive <- as.character(train.labels$X[train.labels$label==1 & !is.na(train.labels$label)])

input.dnv.table <- cbind(train.input[,c('dn_LGD', 'dn_Dmis')],
                   train.input.var[,c('mut_LGD', 'mut_Dmis')])
# model result

result.folder <- 'train.real/'

A_risk <- readxl::read_excel('data/A-risk.xlsx')
log_BF <- read.csv(paste0(result.folder, 'log_BFs.csv'), row.names = 1)
log_logits <- read.csv(paste0(result.folder, 'y_logits.csv'), row.names = 1)
reconstruction_means <- read.csv(paste0(result.folder, 'reconstruction_means.csv'),
                                 row.names = 1)
reconstruction_vars <- read.csv(paste0(result.folder, 'reconstruction_vars.csv'),
                                row.names = 1)
log_logits$log_logits <- log_logits$X1 - log_logits$X0
reconstruction_means <- reconstruction_means[match(rownames(reconstruction_means),
                                                   rownames(input.dnv.table)),]
colnames(reconstruction_means) <- c("LGD_mean", "Dmis_mean")
reconstruction_vars <- reconstruction_vars[match(rownames(reconstruction_vars),
                                                 rownames(input.dnv.table)),]
colnames(reconstruction_vars) <- c("LGD_var", "Dmis_var")
dnv_table <- cbind(input.dnv.table,
                   train.labels=train.labels$label,
                   log_BF=log_BF[match(rownames(input.dnv.table), rownames(log_BF)),],
                   log_logits=log_logits$log_logits[match(rownames(input.dnv.table), rownames(log_logits))],
                   reconstruction_means,
                   reconstruction_vars,
                   A_risk_raw_score=A_risk$`A-risk raw score`[match(rownames(input.dnv.table), A_risk$GeneName)],
                   A_risk_rank_score=A_risk$Rank[match(rownames(input.dnv.table), A_risk$GeneName)])

dnv_table$mut.rate <- dplyr::percent_rank(dnv_table$mut_LGD + dnv_table$mut_Dmis)
dnv_table$pi <- 1/(1+exp(-dnv_table$log_logits))

dnv_table <- dnv_table[order(dnv_table$log_BF, decreasing = T),]
dnv_table$FDR <- Bayesian.FDR.fromBF(dnv_table$log_BF)$FDR
dnv_table$PPA <- 1/(exp(-dnv_table$log_BF) + 1)

# TADA.result <- readRDS('data/SPARK_extTADA.WES1.RDS')
TADA.result <- readRDS('data/SPARK_extTADA.WES1.RDS')

TADA_dnv_table <- TADA.result$dataFDR.full.posterior
dnv_table$TADA.PPA <- TADA_dnv_table$PP[match(rownames(dnv_table),
                                              TADA_dnv_table$HGNC)]
dnv_table$TADA.FDR <- TADA_dnv_table$qvalue[match(rownames(dnv_table),
                                                     TADA_dnv_table$HGNC)]


# supplementary table
dnv_table$GeneID <- TADA_dnv_table$Gene[match(rownames(dnv_table),TADA_dnv_table$HGNC)]
write.csv(dnv_table, file = 'figs/table.S6.csv')

# remove train positive genes
dnv_table <- dnv_table[!rownames(dnv_table) %in% train.positive,]
dnv_table <- dnv_table[order(dnv_table$log_BF, decreasing = T),]
dnv_table$FDR <- Bayesian.FDR.fromBF(dnv_table$log_BF)$FDR
dnv_table <- dnv_table[order(dnv_table$TADA.PPA, decreasing = T),]
dnv_table$TADA.FDR <- Bayesian.FDR.fromPP(dnv_table$TADA.PPA)$FDR

# FDR of dnv_table
dnv_table$group <- rep(NA, dim(dnv_table)[1])
threshold_1 <- 0.05
threshold_2 <- 0.1
dnv_table$group[dnv_table$FDR<=threshold_2 & dnv_table$TADA.FDR<=threshold_2] = 'both'
dnv_table$group[dnv_table$FDR>threshold_1 & dnv_table$TADA.FDR<=threshold_1] = 'FDR_0.05_extTADA_only'
dnv_table$group[dnv_table$FDR<=threshold_1 & dnv_table$TADA.FDR>threshold_1] = 'FDR_0.05_VBASS_only'
dnv_table$group[dnv_table$FDR>threshold_2 
                & dnv_table$TADA.FDR<=threshold_2
                & dnv_table$TADA.FDR>threshold_1] = 'FDR_0.1_extTADA_only'
dnv_table$group[dnv_table$FDR<=threshold_2
                & dnv_table$FDR>threshold_1
                & dnv_table$TADA.FDR>threshold_2] = 'FDR_0.1_VBASS_only'

dnv_table$group[dnv_table$FDR>threshold_2 & dnv_table$TADA.FDR>threshold_2] = 'NA'
dnv_table$label <- rownames(dnv_table)
dnv_table$label[dnv_table$group %in% c('NA', 'both')] <- ''
write.csv(dnv_table, file = 'figs/table.S7.csv')
# plot comparison of log_BF and TADA.FDR

p <- ggplot(dnv_table, aes(x=FDR, y=TADA.FDR, col=group, label=label)) +
  geom_point() +
  geom_vline(xintercept=c(threshold_1), col="red", alpha=0.5, linetype='dotdash') +
  geom_hline(yintercept=c(threshold_1), col="red", alpha=0.5, linetype='dotdash') + 
  geom_vline(xintercept=c(threshold_2), col="blue", alpha=0.5, linetype='dotdash') +
  geom_hline(yintercept=c(threshold_2), col="blue", alpha=0.5, linetype='dotdash') + 
  theme_light() +
  xlim(0, 0.25) +
  ylim(0, 0.25) +
  ylab('extTADA.FDR') +
  xlab('VBASS.FDR') +
  geom_text_repel(size=2.5, colour='black')
ggsave(plot = p, filename = paste0('figs/fig.5.A.pdf'),
       width = 6, height = 6)


