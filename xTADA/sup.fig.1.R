args = commandArgs(trailingOnly=TRUE)
# C <- as.numeric(args[1])
# samplesize <- as.numeric(args[1])
# seed <- args[3]
C <- 0.28242408
samplesize <- 5000
seed <- NA
result.folder = "figs/"
dir.create(paste0(result.folder, 'fig.', samplesize, '/'))

table.all <- c()
if (is.na(seed)) {
  for (i in 1:100) {
    # extTADAsim <- readRDS(paste0("RDS.files/extTADA_simulation_c.", C ,".size.", 
    #                              samplesize, ".seed.", i, ".RDS"))
    extTADAsim_post <- readRDS(paste0("RDS.files/extTADA_simulation_c.", C ,".size.", 
                                      samplesize, ".seed.", i, ".posterior.RDS"))
    # xTADAsim <- readRDS(paste0("RDS.files/xTADA_simulation_c.", C ,".size.",
    #                            samplesize, ".seed.", i, ".RDS"))
    xTADAsim_post <- readRDS(paste0("RDS.files/xTADA_simulation_c.", C ,".size.",
                                    samplesize, ".seed.", i, ".posterior.RDS"))
    # true_label <- readRDS(paste0("xTADA_simulation_c.", C ,".size.", samplesize, ".posterior.RDS"))
    
    table <- xTADAsim_post
    table$extTADA_qvalue <- extTADAsim_post$qvalue[match(table$Gene, extTADAsim_post$Gene)]
    table$extTADA_PP <- extTADAsim_post$PP[match(table$Gene, extTADAsim_post$Gene)]
    table$mutSum <- table$mut_Dmis + table$mut_LGD
    table$mut.Rank <- dplyr::percent_rank(table$mutSum)
    table$cov.Rank <- dplyr::percent_rank(table$covariates_reordered)
    table$xTADA_qvalue <- table$qvalue
    table$xTADA_PP <- table$PP
    table$True_Label <- table$true
    
    table$poisson_pvalue <- 1-ppois(table$dn_Dmis+table$dn_LGD-1, (table$mut_LGD+table$mut_Dmis)*samplesize*2)
    table$poisson_qvalue <- p.adjust(table$poisson_pvalue)
    
    table.all <- rbind(table.all, table)
  }
  table <- table.all
} 

Bayesian.FDR <- function(PP, alpha=0.05) {
  # convert BFs to PPA (posterior probability of alternative model)
  #pi <- 1-pi0
  #q <- apply(pi*BF, 1, prod)/(apply(t(1 - pi), 1, prod)+apply(pi*BF, 1, prod)) # PPA
  q0 <- 1 - PP # posterior probability of null model
  
  # the FDR at each PPA cutoff
  n <- length(PP)
  FDR <- numeric(n)
  for (i in 1:n) FDR[i] <- sum(q0[1:i]) / i 
  
  # the cutoff
  t <- 1
  while (t <= length(q0) & mean(q0[1:t]) <= alpha) { t <- t+1 }
  return (list(FDR=FDR, ND=t))
}

# table <- table[order(table$xTADA_PP), ]
# table$xTADA_qvalue <- Bayesian.FDR(table$xTADA_PP)$FDR
# table <- table[order(table$extTADA_PP), ]
# table$extTADA_qvalue <- Bayesian.FDR(table$extTADA_PP)$FDR

threshold_1 <- 0.1
threshold_2 <- 0.1

library(ggplot2)
# Covariate Rank
p3 <- ggplot(table, aes(x=extTADA_qvalue, y=xTADA_qvalue, col=cov.Rank, label = Gene)) +
  geom_point(size=0.2) +
  geom_vline(xintercept=c(threshold_1), col="red", alpha=0.5, linetype='dotdash') +
  geom_hline(yintercept=c(threshold_2), col="red", alpha=0.5, linetype='dotdash') + 
  theme_light() + ggtitle("Comparison of extTADA and xTADA q-values") +
  scale_color_gradient2(midpoint=0.5, low="blue", mid="white", high="red", space ="Lab" ) +
  ggeasy::easy_center_title()

ggsave(plot = p3, filename = paste0(result.folder, 'fig.', samplesize, '/c.', C, '.size.', samplesize, '.simulation.covRank.pdf'),
       width = 5, height = 4)

# True Label
p2 <- ggplot(table, aes(x=extTADA_qvalue, y=xTADA_qvalue, col=factor(True_Label), label = Gene)) +
  geom_point(size=0.2) +
  geom_vline(xintercept=c(threshold_1), col="red", alpha=0.5, linetype='dotdash') +
  geom_hline(yintercept=c(threshold_2), col="red", alpha=0.5, linetype='dotdash') + 
  theme_light() + ggtitle("Comparison of extTADA and xTADA q-values") +
  scale_colour_discrete(name="True Label", labels=c("Non Risk","Risk")) +
  ggeasy::easy_center_title()
  
ggsave(plot = p2, filename = paste0(result.folder, 'fig.', samplesize, '/c.', C, '.size.', samplesize, '.simulation.TrueLabel.pdf'),
       width = 5, height = 4)

# Mutation Rank
p4 <- ggplot(table, aes(x=extTADA_qvalue, y=xTADA_qvalue, col=mut.Rank, label = Gene)) +
  geom_point(size=0.2) +
  geom_vline(xintercept=c(threshold_1), col="red", alpha=0.5, linetype='dotdash') +
  geom_hline(yintercept=c(threshold_2), col="red", alpha=0.5, linetype='dotdash') + 
  theme_light() + ggtitle("Comparison of extTADA and xTADA q-values") +
  scale_color_gradient2(midpoint=0.5, low="blue", mid="white", high="red", space ="Lab" ) +
  ggeasy::easy_center_title()
ggsave(plot = p4, filename = paste0(result.folder, 'fig.', samplesize, '/c.', C, '.size.', samplesize, '.simulation.mutRank.pdf'),
       width = 5, height = 4)

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


get_local_FDR_curve <- function(true_labels, qvalue) {
  FDR.points <- 100
  FDR_AutoEncoder <- data.frame(FDR=rep(NA, FDR.points),
                                qvalue=seq(0, FDR.points-1)/FDR.points+1/2/FDR.points,
                                qvalue.upper=seq(1, FDR.points)/FDR.points,
                                qvalue.lower=seq(0, FDR.points-1)/FDR.points,
                                num_points=rep(NA, FDR.points))
  for (i in 1:(FDR.points)) {
    all_points <- sum(qvalue<=FDR_AutoEncoder$qvalue.upper[i] & qvalue>FDR_AutoEncoder$qvalue.lower[i])
    false_points <- all_points - sum(true_labels[qvalue<=FDR_AutoEncoder$qvalue.upper[i] & qvalue>FDR_AutoEncoder$qvalue.lower[i]])
    if (all_points==0) {
      FDR_AutoEncoder$FDR[i] <- 0
    } else {
      FDR_AutoEncoder$FDR[i] <- false_points / all_points
    }
    FDR_AutoEncoder$num_points[i] <- all_points
  }
  FDR_AutoEncoder
}


table <- table[order(table$xTADA_qvalue), ]
FDR_curve_1 <- get_FDR_curve(table$True_Label, 
                             table$xTADA_qvalue)
table <- table[order(table$xTADA_qvalue, decreasing = T), ]
FDR_local_curve_1 <- get_local_FDR_curve(table$True_Label, 
                                         table$xTADA_PP)

table <- table[order(table$extTADA_qvalue), ]
FDR_curve_2 <- get_FDR_curve(table$True_Label,
                             table$extTADA_qvalue)
table <- table[order(table$extTADA_qvalue, decreasing = T), ]
FDR_local_curve_2 <- get_local_FDR_curve(table$True_Label, 
                                         table$extTADA_PP)

table <- table[order(table$poisson_qvalue), ]
FDR_curve_3 <- get_FDR_curve(table$True_Label,
                             table$poisson_qvalue)
table <- table[order(table$poisson_pvalue), ]
FDR_local_curve_3 <- get_local_FDR_curve(table$True_Label, 
                                         1-table$poisson_pvalue)

FDR_curve <- data.frame(qvalue=c(FDR_curve_1$qvalue, FDR_curve_2$qvalue, FDR_curve_3$qvalue),
                        FDR=c(FDR_curve_1$FDR, FDR_curve_2$FDR, FDR_curve_3$FDR),
                        num_points=c(FDR_curve_1$num_points, FDR_curve_2$num_points, FDR_curve_3$num_points),
                        model=c(rep('xTADA', dim(FDR_curve_1)[1]),
                                rep('extTADA', dim(FDR_curve_1)[1]),
                                rep('Poisson', dim(FDR_curve_1)[1])))

p <- ggplot(FDR_curve, aes(x=qvalue, y=FDR, col=model)) +
  # geom_point() +
  geom_abline(intercept = 0, col="red", alpha=0.5, linetype='dotdash') +
  # geom_line() +
  geom_smooth() +
  xlim(0, 0.1) +
  # ylim(0, 0.11) +
  scale_y_continuous(breaks = seq(0, 0.125, by = 0.025), limits = c(0, 0.11)) +
  theme_light()
# geom_text_repel(size=2.5, colour='black')
ggsave(plot = p, filename = paste0(result.folder, 'fig.', samplesize, '/c.', C, '.size.', samplesize, '.compare.qvalue.FDR.small.pdf'),
       width = 5, height = 4)
p <- ggplot(FDR_curve, aes(x=qvalue, y=FDR, col=model)) +
  # geom_point() +
  geom_abline(intercept = 0, col="red", alpha=0.5, linetype='dotdash') +
  geom_line() +
  # geom_smooth() +
  # xlim(0, 0.1) +
  # ylim(0, 0.11) +
  # scale_y_continuous(breaks = seq(0, 0.125, by = 0.025), limits = c(0, 0.11)) +
  theme_light()
# geom_text_repel(size=2.5, colour='black')
ggsave(plot = p, filename = paste0(result.folder, 'fig.', samplesize, '/c.', C, '.size.', samplesize, '.compare.qvalue.FDR.pdf'),
       width = 5, height = 4)

# local FDR
FDR_local_curve <- data.frame(qvalue=c(FDR_local_curve_1$qvalue, FDR_local_curve_2$qvalue, FDR_local_curve_3$qvalue),
                        FDR=c(FDR_local_curve_1$FDR, FDR_local_curve_2$FDR, FDR_local_curve_3$FDR),
                        num_points=c(FDR_local_curve_1$num_points, FDR_local_curve_2$num_points, FDR_local_curve_3$num_points),
                        model=c(rep('xTADA', dim(FDR_local_curve_1)[1]),
                                rep('extTADA', dim(FDR_local_curve_1)[1]),
                                rep('Poisson', dim(FDR_local_curve_1)[1])))

p <- ggplot(FDR_local_curve, aes(x=qvalue, y=FDR, col=model)) +
  # geom_point() +
  geom_abline(intercept = 0, col="red", alpha=0.5, linetype='dotdash') +
  # geom_line() +
  geom_smooth() +
  xlim(0, 0.1) +
  # ylim(0, 0.11) +
  scale_y_continuous(breaks = seq(0, 0.125, by = 0.025), limits = c(0, 0.11)) +
  theme_light()
# geom_text_repel(size=2.5, colour='black')
ggsave(plot = p, filename = paste0(result.folder, 'fig.', samplesize, '/c.', C, '.size.', samplesize, '.compare.PP.FDR.local.small.pdf'),
       width = 5, height = 4)
p <- ggplot(FDR_local_curve[FDR_local_curve$model!="Poisson",], aes(x=qvalue, y=FDR, col=model)) +
  # geom_point() +
  # geom_abline(intercept = 0, col="red", alpha=0.5, linetype='dotdash') +
  geom_line() +
  xlab("PPA") +
  # geom_smooth() +
  # xlim(0, 0.1) +
  # ylim(0, 0.11) +
  # scale_y_continuous(breaks = seq(0, 0.125, by = 0.025), limits = c(0, 0.11)) +
  theme_light()
# geom_text_repel(size=2.5, colour='black')
ggsave(plot = p, filename = paste0(result.folder, 'fig.', samplesize, '/c.', C, '.size.', samplesize, '.compare.PP.FDR.local.pdf'),
       width = 5, height = 4)


# precision recall
table <- table[order(table$xTADA_qvalue), ]
all.true <- sum(table$true)
xTADA.pc <- data.frame(qvalue=table$xTADA_qvalue,
                       true=table$true,
                       precision=table$true,
                       recall=table$true)
xTADA.pc <- xTADA.pc[1:floor(1/6*dim(xTADA.pc)[1]), ]
table <- table[order(table$extTADA_qvalue), ]
extTADA.pc <- data.frame(qvalue=table$extTADA_qvalue,
                       true=table$true,
                       precision=table$true,
                       recall=table$true)
extTADA.pc <- extTADA.pc[1:floor(1/6*dim(extTADA.pc)[1]), ]
table <- table[order(table$poisson_pvalue), ]
poisson.pc <- data.frame(pvalue=table$poisson_pvalue,
                         true=table$true,
                         precision=table$true,
                         recall=table$true)
poisson.pc <- poisson.pc[1:floor(1/6*dim(poisson.pc)[1]), ]
for (i in 1:dim(xTADA.pc)[1]) {
  xTADA.pc$precision[i] <- sum(xTADA.pc$true[1:i])/i
  xTADA.pc$recall[i] <- sum(xTADA.pc$true[1:i])/all.true
  extTADA.pc$precision[i] <- sum(extTADA.pc$true[1:i])/i
  extTADA.pc$recall[i] <- sum(extTADA.pc$true[1:i])/all.true
  poisson.pc$precision[i] <- sum(poisson.pc$true[1:i])/i
  poisson.pc$recall[i] <- sum(poisson.pc$true[1:i])/all.true
}
# for (i in 1:dim(extTADA.pc)[1]) {
#   extTADA.pc$precision[i] <- sum(extTADA.pc$true[1:i])/i
#   extTADA.pc$recall[i] <- sum(extTADA.pc$true[1:i])/all.true
# }
# for (i in 1:dim(poisson.pc)[1]) {
#   poisson.pc$precision[i] <- sum(poisson.pc$true[1:i])/i
#   poisson.pc$recall[i] <- sum(poisson.pc$true[1:i])/all.true
# }
to.plot <- data.frame(precision=c(xTADA.pc$precision, extTADA.pc$precision, poisson.pc$precision),
                      recall=c(xTADA.pc$recall, extTADA.pc$recall, poisson.pc$recall),
                      model=c(rep('xTADA', dim(xTADA.pc)[1]), 
                              rep('extTADA', dim(extTADA.pc)[1]),
                              rep('poisson', dim(poisson.pc)[1])))
p<-ggplot(to.plot, aes(x=recall,y=precision,col=model)) +
  geom_line(alpha=0.4) +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_light()
ggsave(plot = p, filename = paste0(result.folder, 'fig.', samplesize, '/c.', C, '.size.', samplesize, '.precision.recall.pdf'),
       width = 5, height = 4)


