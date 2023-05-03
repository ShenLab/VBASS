go.terms <- as.list(GSA::GSA.read.gmt(paste("~/Data/GSEA/msigdb.v7.0.symbols.gmt", sep = "")))
to.select <- c()
for (i in 1:length(go.terms$genesets)) {
  if (length(go.terms$genesets[[i]]) > 5) {
    to.select <- c(to.select, i)
  }
}

library(doParallel)
cl <- makeCluster(72)
registerDoParallel(cl)
result <- foreach (i = 1:length(to.select), .combine = rbind) %dopar% {
  chd <- readRDS(paste0('CHD.GO/', i, '.RDS'))
  asd <- readRDS(paste0('ASD.GO/', i, '.RDS'))
  res <- data.frame(geneset=go.terms$geneset.names[to.select[i]],
                    description=go.terms$geneset.descriptions[to.select[i]],
                    chd=chd$genesetInfo$pGS, asd=asd$genesetInfo$pGS)
  res
}
stopCluster(cl)

result$chd.q <- p.adjust(result$chd)
result$asd.q <- p.adjust(result$asd)

write.csv(result, file ='DECO.geneset.result.csv')

chd.genesets <- result$geneset[result$chd.q <= 0.05]
asd.genesets <- result$geneset[result$asd.q <= 0.05]

# use the most significant geneset as input.
chd.geneset <- result$geneset[result$chd == min(result$chd)]
asd.geneset <- result$geneset[result$asd == min(result$asd)]
chd.idx <- which(to.select == which(go.terms$geneset.names == chd.geneset))
asd.idx <- which(to.select == which(go.terms$geneset.names == asd.geneset))

chd.result <- readRDS(paste0('CHD.GO/', chd.idx, '.RDS'))$dataPP
asd.result <- readRDS(paste0('ASD.GO/', asd.idx, '.RDS'))$dataPP

chd.VBASS.out <- read.csv("~/ld1/VBASS_xTADA_paper/scalar_VBASS/result/table.S4.csv", row.names = 1)
asd.VBASS.out <- read.csv("~/ld1/VBASS_xTADA_paper/VBASS/figs/table.S5.csv")

chd.VBASS.out$DECO.PP <- chd.result$PP[match(chd.VBASS.out$Gene, chd.result$Gene)]
asd.VBASS.out$DECO.PP <- asd.result$PP[match(asd.VBASS.out$X, asd.result$Gene)]

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

chd.VBASS.out <- chd.VBASS.out[order(chd.VBASS.out$DECO.PP, decreasing = T),]
chd.VBASS.out$DECO.FDR <- Bayesian.FDR.fromPP(chd.VBASS.out$DECO.PP)$FDR
asd.VBASS.out <- asd.VBASS.out[order(asd.VBASS.out$DECO.PP, decreasing = T),]
asd.VBASS.out$DECO.FDR <- Bayesian.FDR.fromPP(asd.VBASS.out$DECO.PP)$FDR
asd.VBASS.out$HGNC <- asd.VBASS.out$X

plot.FDR.compare <- function(dnv_table, prefix="") {
  library(ggplot2)
  library(ggrepel)
  library(ggExtra)
  dnv_table$group <- rep(NA, dim(dnv_table)[1])
  threshold_1 <- 0.05
  threshold_2 <- 0.1
  dnv_table$group[dnv_table$FDR<=threshold_2 & dnv_table$DECO.FDR<=threshold_2] = 'both'
  dnv_table$group[dnv_table$FDR>threshold_1 & dnv_table$DECO.FDR<=threshold_1] = 'FDR_0.05_DECO_only'
  dnv_table$group[dnv_table$FDR<=threshold_1 & dnv_table$DECO.FDR>threshold_1] = 'FDR_0.05_VBASS_only'
  dnv_table$group[dnv_table$FDR>threshold_2 
                  & dnv_table$DECO.FDR<=threshold_2
                  & dnv_table$DECO.FDR>threshold_1] = 'FDR_0.1_DECO_only'
  dnv_table$group[dnv_table$FDR<=threshold_2
                  & dnv_table$FDR>threshold_1
                  & dnv_table$DECO.FDR>threshold_2] = 'FDR_0.1_VBASS_only'
  
  dnv_table$group[dnv_table$FDR>threshold_2 & dnv_table$DECO.FDR>threshold_2] = 'NA'
  dnv_table$label <- dnv_table$HGNC
  dnv_table$label[dnv_table$group %in% c('NA', 'both')] <- ''
  
  p <- ggplot(dnv_table, aes(x=FDR, y=DECO.FDR, col=group, label=label)) +
    geom_point() +
    geom_vline(xintercept=c(threshold_1), col="red", alpha=0.5, linetype='dotdash') +
    geom_hline(yintercept=c(threshold_1), col="red", alpha=0.5, linetype='dotdash') + 
    geom_vline(xintercept=c(threshold_2), col="blue", alpha=0.5, linetype='dotdash') +
    geom_hline(yintercept=c(threshold_2), col="blue", alpha=0.5, linetype='dotdash') + 
    theme_light() +
    xlim(0, 0.25) +
    ylim(0, 0.25) +
    ylab('DECO.FDR') +
    xlab('VBASS.FDR') +
    geom_text_repel(size=2.5, colour='black')
  ggsave(plot = p, filename = paste0(prefix, '.DECO.VBASS.compare.pdf'),
         width = 6, height = 4)
  dnv_table
}

asd.VBASS.out <- plot.FDR.compare(asd.VBASS.out, "ASD")
chd.VBASS.out <- plot.FDR.compare(chd.VBASS.out, "CHD")
write.csv(chd.VBASS.out, file = "chd.VBASS.DECO.csv")
write.csv(asd.VBASS.out, file = "asd.VBASS.DECO.csv")


# compare with Sifrim et al 2016
sifrim <- readxl::read_excel('Sifrim.de.novo.xlsx')
for (i in 1:dim(chd.VBASS.out)[1]) {
  chd.VBASS.out$Sifrim.LGD[i] <- sum(sifrim$symbol == chd.VBASS.out$HGNC[i] & sifrim$consequence != "missense_variant")
  chd.VBASS.out$Sifrim.mis[i] <- sum(sifrim$symbol == chd.VBASS.out$HGNC[i] & sifrim$consequence == "missense_variant")
}
write.csv(chd.VBASS.out, file = "chd.VBASS.DECO.csv")

# compare with NDD and SFARI
sfari <- read.csv('SFARI_20230413.csv')
asd.VBASS.out$sfari <- sfari$gene.score[match(asd.VBASS.out$HGNC, sfari$gene.symbol)]
asd.VBASS.out$in.gene.set <- asd.VBASS.out$HGNC %in% go.terms$genesets[[which(go.terms$geneset.names == asd.geneset)]]

kaplanis <- read.delim('kaplanis.tsv')
for (i in 1:dim(asd.VBASS.out)[1]) {
  asd.VBASS.out$Kaplanis.LGD[i] <- sum(kaplanis$symbol == asd.VBASS.out$HGNC[i] & kaplanis$consequence %in% c("frameshift_variant", "splice_acceptor_variant", "splice_donor_variant", "stop_gained", "stop_lost"))
  asd.VBASS.out$Kaplanis.mis[i] <- sum(kaplanis$symbol == asd.VBASS.out$HGNC[i] & kaplanis$consequence == "missense_variant")
}
write.csv(asd.VBASS.out, file = "asd.VBASS.DECO.csv")
