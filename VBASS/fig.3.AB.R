summary.folder <- 'train.simulation.out/'

library(ggrepel)
Bayesian.FDR <- function(log_BF, alpha=0.05) {
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

# model result
# args = commandArgs(trailingOnly=TRUE)
all_dnv_table <- c()
for (i in 1:50) {
  train.input <- read.csv(paste0('VBASS.simulation.with.RDS/input.x.',
                                 i, '.csv'),
                          row.names = 1)
  train.input.var <- read.csv(paste0('VBASS.simulation.with.RDS/input.x_var.',
                                     i, '.csv'),
                              row.names = 1)
  
  train.point.est <- train.input / train.input.var
  
  input.dnv.table <- cbind(train.input[,c('dn.cls1', 'dn.cls2')],
                           train.input.var[,c('mut.cls1', 'mut.cls2')])
  
  result.folder <- paste0(summary.folder, i, '.out/')
  
  # epoch_number <- args[2]
  # A_risk <- readxl::read_excel('~/fetal_brain/genetics/A-risk/A-risk.xlsx')
  
  log_BF <- read.csv(paste0(result.folder, 'log_BFs.csv'), row.names = 1)
  log_logits <- read.csv(paste0(result.folder, 'y_logits.csv'), row.names = 1)
  reconstruction_means <- read.csv(paste0(result.folder, 'reconstruction_means.csv'), row.names = 1)
  reconstruction_vars <- read.csv(paste0(result.folder, 'reconstruction_vars.csv'), row.names = 1)
  
  log_logits$log_logits <- log_logits$X1 - log_logits$X0
  reconstruction_means <- reconstruction_means[match(rownames(reconstruction_means),
                                                     rownames(input.dnv.table)),]
  colnames(reconstruction_means) <- c("LGD_mean", "Dmis_mean")
  reconstruction_vars <- reconstruction_vars[match(rownames(reconstruction_vars),
                                                   rownames(input.dnv.table)),]
  colnames(reconstruction_vars) <- c("LGD_var", "Dmis_var")
  
  real_BF <- read.csv(paste0('VBASS.simulation.with.RDS/real.log.BF.',
                             i, '.csv'),
                      row.names = 1)
  
  dnv_table <- cbind(input.dnv.table,
                     log_BF=log_BF[match(rownames(input.dnv.table), rownames(log_BF)),],
                     log_logits=log_logits$log_logits[match(rownames(input.dnv.table), rownames(log_logits))],
                     real_BF=real_BF$real.log.BF[match(rownames(input.dnv.table), rownames(real_BF))],
                     reconstruction_means,
                     reconstruction_vars)
  
  dnv_table$mut.rate <- dplyr::percent_rank(dnv_table$mut.cls1 + dnv_table$mut.cls2)
  dnv_table$pi <- 1/(1+exp(-dnv_table$log_logits))
  
  dnv_table <- dnv_table[order(dnv_table$log_BF, decreasing = T),]
  dnv_table$qvalue <- Bayesian.FDR(dnv_table$log_BF)$FDR
  dnv_table$PPA <- 1/(exp(dnv_table$log_BF) + 1)
  
  dnv_table <- dnv_table[order(dnv_table$real_BF, decreasing = T),]
  dnv_table$real.qvalue <- Bayesian.FDR(dnv_table$real_BF)$FDR
  dnv_table$real.PPA <- 1/(exp(dnv_table$real_BF) + 1)
  
  TADA.result <- readRDS(paste0('VBASS.simulation.with.RDS/extTADA.simulation.realexp.',
                                i, '.RDS'))
  TADA_dnv_table <- TADA.result$dataFDR.full.posterior
  dnv_table$TADA_qvalue <- TADA_dnv_table$qvalue[match(rownames(dnv_table),
                                                       TADA_dnv_table$Gene)]
  # PPA of dnv_table
  true_labels <- read.csv(paste0('VBASS.simulation.with.RDS/test.label.',
                                 i, '.csv'))
  training_true_labels <- read.csv(paste0('VBASS.simulation.with.RDS/input.label.',
                                          i, '.csv'))
  
  dnv_table$group <- rep(NA, dim(dnv_table)[1])
  dnv_table$true_labels <- true_labels$label[match(rownames(dnv_table), true_labels$X)]
  dnv_table$training_labels <- training_true_labels$label[match(rownames(dnv_table), training_true_labels$X)]
  
  dnv_table$group[dnv_table$true_labels == 1 & !is.na(dnv_table$training_labels)] = 'TP: training'
  dnv_table$group[dnv_table$true_labels == 0 & !is.na(dnv_table$training_labels)] = 'TN: training'
  dnv_table$group[dnv_table$true_labels == 1 & is.na(dnv_table$training_labels)] = 'TP'
  dnv_table$group[dnv_table$true_labels == 0 & is.na(dnv_table$training_labels)] = 'TN'
  
  all_dnv_table <- rbind(all_dnv_table, dnv_table)
}
result.folder <- paste0('figs/')
dir.create(result.folder)

dnv_table <- all_dnv_table

# check the true FDR vs FDR
dnv_table.ranked <- dnv_table[order(dnv_table$qvalue), ]
get_FDR_curve <- function(true_labels, qvalue) {
  FDR.points <- 10000
  FDR_AutoEncoder <- data.frame(FDR=rep(NA, FDR.points+1),
                                qvalue=seq(0, FDR.points)/FDR.points,
                                num_points=rep(NA, FDR.points+1))
  for (i in 1:(FDR.points+1)) {
    all_points <- sum(qvalue<=FDR_AutoEncoder$qvalue[i])
    false_points <- all_points - sum(true_labels[qvalue<=FDR_AutoEncoder$qvalue[i]])
    if (all_points==0) {
      FDR_AutoEncoder$FDR[i] <- 0
    } else {
      FDR_AutoEncoder$FDR[i] <- false_points / all_points
    }
    FDR_AutoEncoder$num_points[i] <- all_points
  }
  FDR_AutoEncoder
}
FDR_curve_1 <- get_FDR_curve(dnv_table.ranked$true_labels, 
                             dnv_table.ranked$qvalue)
FDR_curve_1$model <- 'VBASS'
dnv_table.ranked <- dnv_table[order(dnv_table$TADA_qvalue), ]
FDR_curve_2 <- get_FDR_curve(dnv_table.ranked$true_labels,
                             dnv_table.ranked$TADA_qvalue)
FDR_curve <- data.frame(FDR=c(FDR_curve_1$qvalue, FDR_curve_2$qvalue),
                        real.FDR=c(FDR_curve_1$FDR, FDR_curve_2$FDR),
                        num_points=c(FDR_curve_1$num_points, FDR_curve_2$num_points),
                        model=c(rep('VBASS', dim(FDR_curve_1)[1]),
                                rep('extTADA', dim(FDR_curve_2)[1]))
)

# plot together
p <- ggplot(FDR_curve, aes(x=FDR, y=real.FDR, col=model)) +
  # geom_point() +
  geom_abline(intercept = 0, col="red", alpha=0.5, linetype='dotdash') +
  # geom_line() +
  geom_smooth() +
  xlim(0, 0.1) +
  # ylim(0, 0.11) +
  scale_y_continuous(breaks = seq(0, 0.125, by = 0.025), limits = c(0, 0.11)) +
  theme_light()
# geom_text_repel(size=2.5, colour='black')
ggsave(plot = p, filename = paste0(result.folder, 'fig.3.A.pdf'),
       width = 5, height = 4)

# precision recall
dnv_table <- dnv_table[order(dnv_table$qvalue), ]
all.true <- sum(dnv_table$true)
VBASS.pc <- data.frame(qvalue=dnv_table$qvalue,
                       true=dnv_table$true_labels,
                       precision=dnv_table$true_labels,
                       recall=dnv_table$true_labels)
VBASS.pc <- VBASS.pc[VBASS.pc$qvalue <= 0.5, ]
dnv_table <- dnv_table[order(dnv_table$TADA_qvalue), ]
extTADA.pc <- data.frame(qvalue=dnv_table$TADA_qvalue,
                         true=dnv_table$true_labels,
                         precision=dnv_table$true_labels,
                         recall=dnv_table$true_labels)
extTADA.pc <- extTADA.pc[extTADA.pc$qvalue <= 0.5, ]
dnv_table <- dnv_table[order(dnv_table$real.qvalue), ]
real.pc <- data.frame(qvalue=dnv_table$real.qvalue,
                         true=dnv_table$true_labels,
                         precision=dnv_table$true_labels,
                         recall=dnv_table$true_labels)
real.pc <- real.pc[real.pc$qvalue <= 0.5, ]
for (i in 1:dim(VBASS.pc)[1]) {
  VBASS.pc$precision[i] <- sum(VBASS.pc$true[1:i])/i
  VBASS.pc$recall[i] <- sum(VBASS.pc$true[1:i])/all.true
}
for (i in 1:dim(extTADA.pc)[1]) {
  extTADA.pc$precision[i] <- sum(extTADA.pc$true[1:i])/i
  extTADA.pc$recall[i] <- sum(extTADA.pc$true[1:i])/all.true
}
for (i in 1:dim(real.pc)[1]) {
  real.pc$precision[i] <- sum(real.pc$true[1:i])/i
  real.pc$recall[i] <- sum(real.pc$true[1:i])/all.true
}
to.plot <- data.frame(precision=c(VBASS.pc$precision, extTADA.pc$precision, real.pc$precision),
                      recall=c(VBASS.pc$recall, extTADA.pc$recall, real.pc$recall),
                      model=c(rep('VBASS', dim(VBASS.pc)[1]),
                              rep('extTADA', dim(extTADA.pc)[1]),
                              rep('real', dim(real.pc)[1])))
p<-ggplot(to.plot, aes(x=recall,y=precision,col=model)) +
  geom_line() +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_light()
ggsave(plot = p, filename = paste0(result.folder, 'fig.3.B.pdf'),
       width = 5, height = 4)
