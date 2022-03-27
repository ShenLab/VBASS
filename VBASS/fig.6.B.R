result.folder <- 'train.discovery/'
out.folder <- paste0('figs/')

train.input <- read.csv('AutoEncoderTADA.input.discovery/input.x.csv', row.names = 1)
train.input.var <- read.csv('AutoEncoderTADA.input.discovery/input.x_var.csv', row.names = 1)
point.estimate.train.input <- train.input[,1:80] / train.input.var[,1:80]
validation.dataset <- readRDS('data/SPARK_extTADA.WES1.validation.RDS')
validation.dataset <- validation.dataset$dataFDR

remove_duplicate_line <- function(dataFDR) {
  duplicate_lines <- as.character(dataFDR$HGNC[duplicated(dataFDR$HGNC)])
  cols.to.add <- c('dn_LGD', 'dn_Dmis', 'mut_LGD', 'mut_Dmis')
  for (gene in duplicate_lines) {
    gene.lines <- dataFDR[dataFDR$HGNC == gene,]
    dataFDR[min(which(dataFDR$HGNC==gene)), cols.to.add] <- colSums(gene.lines[,cols.to.add])
  }
  dataFDR <- dataFDR[!duplicated(dataFDR$HGNC),]
  dataFDR
}
validation.dataset <- remove_duplicate_line(validation.dataset)

input.dnv.table <- cbind(validation.dataset[,c('dn_LGD', 'dn_Dmis')],
                         validation.dataset[,c('mut_LGD', 'mut_Dmis')])
rownames(input.dnv.table) <- validation.dataset$HGNC
colnames(input.dnv.table) <- c('dn.cls1', 'dn.cls2', 'mut.cls1', 'mut.cls2')

reference.dnv.table <- input.dnv.table[,c('mut.cls1', 'mut.cls2')]

reference.dnv.table$GeneID <- validation.dataset$HGNC

train.input.label <- read.csv(paste0('AutoEncoderTADA.input.discovery/input.label.csv'),
                              row.names = 1)
train.genes <- rownames(train.input.label)[!is.na(train.input.label$label)]

reference.dnv.table <- reference.dnv.table[!reference.dnv.table$GeneID %in% train.genes,]

log_BFs <- read.csv(paste0(result.folder, 'log_BFs.csv'), row.names = 1)

log_BFs$PPA <- 1/(exp(-log_BFs$X0)+1)

genes.to.check.scores <- list()
extTADA.PP <- readRDS('data//SPARK_extTADA.WES1.discovery.RDS')
extTADA.PP <- extTADA.PP$dataFDR
dn.num <- extTADA.PP$dn_Dmis+extTADA.PP$dn_LGD
mut.num <- extTADA.PP$mut_LGD+extTADA.PP$mut_Dmis
poisson.res <- ppois(dn.num-1, mut.num*2*8308)

genes.to.check.scores.value <- data.frame(row.names = row.names(reference.dnv.table),
                                          PPA = log_BFs$PPA[match(rownames(reference.dnv.table),
                                                                  rownames(log_BFs))],
                                          extTADA.PPA = extTADA.PP$PP[match(rownames(reference.dnv.table),
                                                                            extTADA.PP$HGNC)],
                                          poisson.PP = poisson.res[match(rownames(reference.dnv.table),
                                                                         extTADA.PP$HGNC)])

genes.to.check.scores$value <- cbind(genes.to.check.scores.value, genes=row.names(reference.dnv.table))
genes.to.check.scores.rank <- apply(genes.to.check.scores.value, 2, dplyr::percent_rank)
genes.to.check.scores$rank <- cbind(genes.to.check.scores.rank, genes=row.names(reference.dnv.table))

# bin_lower=c(0, 0.2, 0.8, 0.9)
bin_lower=c(seq(0, 4)*0.2)
# bin_upper=c(0.2, 0.8, 0.9, 1)
bin_upper=c(seq(1, 5)*0.2)
scores <- genes.to.check.scores$value
bin_label=paste0(bin_lower, " < PPA < ", bin_upper)

samplenumber <- 8308
reference.dnv.table$HGNC <- reference.dnv.table$GeneID
pi.estimated <- matrix(NA, nrow = length(bin_lower), ncol = dim(scores)[2]-1)
average.score <- matrix(NA, nrow = length(bin_lower), ncol = dim(scores)[2]-1)
pars.list <- list()

library(doParallel)
library(ggplot2)
library(ggrepel)
class <- c("VBASS", "extTADA", "Poisson")
for (k in 1:2) {
  bin.by.expression <- foreach(j=1:length(bin_lower), .combine=rbind) %dopar% {
    source("geneset_burden_analysis.R")
    expression.bin.geneset <- as.character(reference.dnv.table$GeneID[scores[,k]>bin_lower[j]
                                                                      & scores[,k]<=bin_upper[j]
                                                                      & !is.na(scores[,k])])
    temp.df <- geneset_burden_analysis_dnvtable_mutrate(expression.bin.geneset, input.dnv.table, samplenumber)
    temp.df
  }
  average.score <- c()
  num.genes <- c()
  for (i in 1:length(bin_lower)) {
    average.score[i] <- mean(genes.to.check.scores$value[scores[,k]>bin_lower[i]
                                                         & scores[,k]<=bin_upper[i]
                                                         & !is.na(scores[,k]),k])
    num.genes[i] <- sum(scores[,k]>bin_lower[i]
                        & scores[,k]<=bin_upper[i]
                        & !is.na(scores[,k]))
  }
  to.plot <- data.frame(burden = c(bin.by.expression$burden_LGD,bin.by.expression$burden_Dmis),
                        burden.max = c((bin.by.expression$obs_LGD+sqrt(bin.by.expression$obs_LGD))/bin.by.expression$exp_LGD,
                                       (bin.by.expression$obs_Dmis+sqrt(bin.by.expression$obs_Dmis))/(bin.by.expression$exp_Dmis)),
                        burden.min = c((bin.by.expression$obs_LGD-sqrt(bin.by.expression$obs_LGD))/(bin.by.expression$exp_LGD),
                                       (bin.by.expression$obs_Dmis-sqrt(bin.by.expression$obs_Dmis))/(bin.by.expression$exp_Dmis)),
                        PPA.mean = rep(average.score, 2),
                        label = paste0(bin_label, "\n", num.genes, " genes"),
                        group = c(rep("LGD", dim(bin.by.expression)[1]),
                                  rep("Dmis", dim(bin.by.expression)[1])))
  to.plot$label.y = c(to.plot$burden[1:length(average.score)]*2, 
                      to.plot$burden[(length(average.score)+1):(2*length(average.score))]/2)
  
  p <- ggplot(to.plot, aes(x=PPA.mean, y=burden, label=label, col=group)) +
    geom_line(alpha=0.5) +
    geom_point() +
    geom_errorbar(aes(ymin=burden.min, ymax=burden.max), width=.02, alpha=0.5) +
    geom_text_repel(aes(y=label.y), size=2) +
    ggtitle(paste0(class[k], ' PPA')) +
    scale_y_continuous(trans='log2') +
    theme_bw() +
    ggeasy::easy_center_title()
  
  ggsave(paste0(out.folder, paste0("fig.6.B.", k, ".pdf")), p, height = 3, width = 6)
  
}
