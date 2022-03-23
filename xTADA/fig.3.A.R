library(ggplot2)
args = commandArgs(trailingOnly=TRUE)
# C <- as.numeric(args[1])
samplesize <- as.numeric(args[1])
# seed <- args[3]
C <- 0.28242408
# samplesize <- 50000
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


# precision recall
table <- table[order(table$xTADA_qvalue), ]
all.true <- sum(table$true)
xTADA.pc <- data.frame(qvalue=table$xTADA_qvalue,
                       true=table$true,
                       precision=table$true,
                       recall=table$true)
xTADA.pc <- xTADA.pc[xTADA.pc$qvalue<=0.2, ]
table <- table[order(table$extTADA_qvalue), ]
extTADA.pc <- data.frame(qvalue=table$extTADA_qvalue,
                         true=table$true,
                         precision=table$true,
                         recall=table$true)
extTADA.pc <- extTADA.pc[extTADA.pc$qvalue<=0.2, ]
table <- table[order(table$poisson_pvalue), ]
poisson.pc <- data.frame(pvalue=table$poisson_pvalue,
                         true=table$true,
                         precision=table$true,
                         recall=table$true)
poisson.pc <- poisson.pc[poisson.pc$pvalue<=0.01, ]
for (i in 1:dim(xTADA.pc)[1]) {
  xTADA.pc$precision[i] <- sum(xTADA.pc$true[1:i])/i
  xTADA.pc$recall[i] <- sum(xTADA.pc$true[1:i])/all.true
}
for (i in 1:dim(extTADA.pc)[1]) {
  extTADA.pc$precision[i] <- sum(extTADA.pc$true[1:i])/i
  extTADA.pc$recall[i] <- sum(extTADA.pc$true[1:i])/all.true
}
for (i in 1:dim(poisson.pc)[1]) {
  poisson.pc$precision[i] <- sum(poisson.pc$true[1:i])/i
  poisson.pc$recall[i] <- sum(poisson.pc$true[1:i])/all.true
}
to.plot <- data.frame(precision=c(xTADA.pc$precision, extTADA.pc$precision, poisson.pc$precision),
                      recall=c(xTADA.pc$recall, extTADA.pc$recall, poisson.pc$recall),
                      model=c(rep('xTADA', dim(xTADA.pc)[1]), 
                              rep('extTADA', dim(extTADA.pc)[1]),
                              rep('poisson', dim(poisson.pc)[1])))
p<-ggplot(to.plot, aes(x=recall,y=precision,col=model)) +
  geom_line(alpha=0.8) +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_light()
ggsave(plot = p, filename = paste0(result.folder, 'fig.', samplesize, '/c.', C, '.size.', samplesize, '.precision.recall.pdf'),
       width = 5, height = 4)


