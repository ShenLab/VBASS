C <- 0.28242408
library(ggplot2)
result.folder <- 'figs/'
# samplesize <- 2654
seed <- NA
args = commandArgs(trailingOnly=TRUE)
# C <- as.numeric(args[1])
samplesize <- as.numeric(args[1])
# seed <- args[3]
table.all <- c()
scalarVBASSsim_post <- readRDS(paste0("RDS.files/scalarVBASS_simulation_c.", C ,".size.",
                                samplesize, ".seed.", 1, ".posterior.RDS"))
true.genes <- scalarVBASSsim_post$Gene[scalarVBASSsim_post$true]
# FDR.cutoffs <- c(0.05)
# recall.mutrate <- matrix(NA, nrow = length(true.genes), ncol = 2*length(FDR.cutoffs)+1)
scalarVBASSsim_post$mut.Rank <- dplyr::percent_rank(scalarVBASSsim_post$mutSum)
true.genes.mutRate <- scalarVBASSsim_post$mutSum[match(true.genes, scalarVBASSsim_post$Gene)]
true.genes.mutRate.Rank <- scalarVBASSsim_post$mut.Rank[match(true.genes, scalarVBASSsim_post$Gene)]

for (i in 1:100) {
  # extTADAsim <- readRDS(paste0("RDS.files/extTADA_simulation_c.", C ,".size.", 
  #                              samplesize, ".seed.", i, ".RDS"))
  extTADAsim_post <- readRDS(paste0("RDS.files/extTADA_simulation_c.", C ,".size.", 
                                    samplesize, ".seed.", i, ".posterior.RDS"))
  # scalarVBASSsim <- readRDS(paste0("RDS.files/scalarVBASS_simulation_c.", C ,".size.",
  #                            samplesize, ".seed.", i, ".RDS"))
  scalarVBASSsim_post <- readRDS(paste0("RDS.files/scalarVBASS_simulation_c.", C ,".size.",
                                  samplesize, ".seed.", i, ".posterior.RDS"))
  # true_label <- readRDS(paste0("scalarVBASS_simulation_c.", C ,".size.", samplesize, ".posterior.RDS"))
  
  table <- scalarVBASSsim_post
  table$extTADA_FDR <- extTADAsim_post$qvalue[match(table$Gene, extTADAsim_post$Gene)]
  table$extTADA_PP <- extTADAsim_post$PP[match(table$Gene, extTADAsim_post$Gene)]
  table$mut.Rank <- dplyr::percent_rank(table$mutSum)
  table$cov.Rank <- dplyr::percent_rank(table$covariates_reordered)
  table$mutSum <- table$mut_Dmis + table$mut_LGD
  table$scalarVBASS_FDR <- table$qvalue
  table$scalarVBASS_PP <- table$PP
  table$True_Label <- table$true
  
  table$poisson_pvalue <- 1-ppois(table$dn_Dmis+table$dn_LGD-1, (table$mut_LGD+table$mut_Dmis)*samplesize*2)
  table$poisson_FDR <- p.adjust(table$poisson_pvalue)
  
  table.all <- rbind(table.all, table)
}

FDR.cutoff <- c(0.05)
to.plot.all <- c()
for (j in FDR.cutoff) {
  extTADA.recall <- rep(NA, length(true.genes))
  scalarVBASS.recall <- rep(NA, length(true.genes))
  poisson.recall <- rep(NA, length(true.genes))
  for (i in 1:length(true.genes)) {
    gene.table <- table.all[table.all$Gene==true.genes[i],]
    extTADA.recall[i] <- sum(gene.table$extTADA_FDR<=j)/100
    scalarVBASS.recall[i] <- sum(gene.table$scalarVBASS_FDR<=j)/100
    poisson.recall[i] <- sum(gene.table$poisson_FDR<=j)/100
  }
  to.plot <- data.frame(recall=c(extTADA.recall, scalarVBASS.recall, poisson.recall),
                        mut.rate=rep(true.genes.mutRate, 3),
                        mut.rate.rank=rep(true.genes.mutRate.Rank, 3),
                        # cov.mean = rep(cov.mean, 2),
                        model.name=c(paste0(rep('extTADA', length(true.genes)), '.FDR.', j),
                                     paste0(rep('VBASS', length(true.genes)), '.FDR.', j),
                                     paste0(rep('poisson', length(true.genes)), '.FDR.', j)
                        )
  )
  to.plot.all <- rbind(to.plot.all, to.plot)
}

p <- ggplot(to.plot.all, aes(x=mut.rate.rank, y=recall, col=model.name)) +
  # geom_point() +
  geom_smooth() +
  # geom_line() +
  theme_light()
ggsave(plot = p, filename = paste0(result.folder, 'fig.', samplesize, '/c.', C, '.size.', samplesize, '.recall.mut.pdf'),
       width = 4.5, height = 3)

write.csv(to.plot.all, paste0(result.folder, 'fig.', samplesize, '/c.', C, '.size.', samplesize, '.recall.mut.csv'))




