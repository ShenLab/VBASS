library(rstan)
# Compare with bulk
sc.input <- read.csv('VBASS.input/input.x.csv', row.names = 1)
GW.cols <- colnames(sc.input)[grepl('GW', colnames(sc.input))]
week.cols <- colnames(sc.input)[grepl('week', colnames(sc.input))]

# sc.input <- rowSums(sc.input[,GW.cols])
sc.number <- read.csv('VBASS.input/input.x_var.csv', row.names = 1)
# sc.number <- rowSums(sc.number[,GW.cols])
# generate pseudo bulk
pseudo.bulk.1 <- data.frame(expression=rowSums(sc.input[,GW.cols]) / rowSums(sc.number[,GW.cols]),
                          HGNC=rownames(sc.input))
pseudo.bulk.2 <- data.frame(expression=rowSums(sc.input[,week.cols]) / rowSums(sc.number[,week.cols]),
                            HGNC=rownames(sc.input))

modelfile = dir('../scalar_VBASS/scalar_VBASS/', '.R$')
modelfile = modelfile[!modelfile == 'model_src.R']
for (ii in modelfile) {
  sourcedir = paste0('../scalar_VBASS/scalar_VBASS/', ii)
  print(sourcedir)
  source(sourcedir)
}
DNscalarVBASS <- readChar(paste0("../scalar_VBASS/scalar_VBASS/", "/model.stan"),
                          file.info(paste0("../scalar_VBASS/scalar_VBASS/", "/model.stan"))$size)

cases <- read.csv('data/SPARK.WES1.csv', row.names = 1)
trios <- read.csv('data/SPARK.WES1.trios.csv', row.names = 1)
samplenumber <- dim(trios)[1]

# reference <- read.csv('VBASS.input/input.x_var.csv', row.names = 1)
# reference <- reference[,81:82]/samplenumber
# reference$GeneID <- rownames(reference)
# reference$HGNC <- rownames(reference)
reference <- read.table("data/mutrate.3mer.txt", sep = "\t", header = TRUE, na.strings = ".")
blacklist <- read.table("data/GENCODEV19_blacklist.txt", sep = "\t", header = F)
reference = reference[!reference$GeneID %in% blacklist$V1,]
reference <- reference[match(pseudo.bulk.1$HGNC, reference$HGNC),]

covariatesData.1 <- data.frame(gene=as.character(reference$GeneID),
                             expression=as.numeric(pseudo.bulk.1$expression[match(reference$HGNC, pseudo.bulk.1$HGNC)]))
covariatesData.2 <- data.frame(gene=as.character(reference$GeneID),
                             expression=as.numeric(pseudo.bulk.2$expression[match(reference$HGNC, pseudo.bulk.2$HGNC)]))

geneset <- as.character(reference$GeneID)

result.1 = gene_set_scalarVBASS(geneset, cases, samplenumber, reference, covariatesData.1)
result.2 = gene_set_scalarVBASS(geneset, cases, samplenumber, reference, covariatesData.2)
saveRDS(result.1, file = "ASD_VBASS.GW.RDS")
saveRDS(result.2, file = "ASD_VBASS.week.RDS")

# draw.xTADA.S_distribution.from.pars0(result$pars0, 'review.figs/ASD.VBASS.bulk.pdf')
table.S7 <- read.csv('figs/table.S7.csv', row.names = 1)
table.S7$bulk.dataset1.FDR <- result.1$dataFDR.full.posterior$qvalue[match(table.S7$GeneID, result.1$dataFDR.full.posterior$Gene)]
table.S7$bulk.dataset1.PPA <- result.1$dataFDR.full.posterior$PP[match(table.S7$GeneID, result.1$dataFDR.full.posterior$Gene)]

table.S7$bulk.dataset2.FDR <- result.2$dataFDR.full.posterior$qvalue[match(table.S7$GeneID, result.2$dataFDR.full.posterior$Gene)]
table.S7$bulk.dataset2.PPA <- result.2$dataFDR.full.posterior$PP[match(table.S7$GeneID, result.2$dataFDR.full.posterior$Gene)]
write.csv(table.S7, file = 'figs/table.S7.csv')


table.S7$group <- NA
# FDR of table.S7
threshold_1 <- 0.05
threshold_2 <- 0.1
table.S7$group[table.S7$FDR<=threshold_2 & table.S7$bulk.dataset1.FDR<=threshold_2] = 'both'
table.S7$group[table.S7$FDR>threshold_1 & table.S7$bulk.dataset1.FDR<=threshold_1] = 'FDR_0.05_bulk_only'
table.S7$group[table.S7$FDR<=threshold_1 & table.S7$bulk.dataset1.FDR>threshold_1] = 'FDR_0.05_sc_only'
table.S7$group[table.S7$FDR>threshold_2 
                & table.S7$bulk.dataset1.FDR<=threshold_2
                & table.S7$bulk.dataset1.FDR>threshold_1] = 'FDR_0.1_bulk_only'
table.S7$group[table.S7$FDR<=threshold_2
                & table.S7$FDR>threshold_1
                & table.S7$bulk.dataset1.FDR>threshold_2] = 'FDR_0.1_sc_only'


p <- ggplot(table.S7, aes(x=FDR, y=bulk.dataset1.FDR, col=group, label=label)) +
  geom_point() +
  geom_vline(xintercept=c(threshold_1), col="red", alpha=0.5, linetype='dotdash') +
  geom_hline(yintercept=c(threshold_1), col="red", alpha=0.5, linetype='dotdash') + 
  geom_vline(xintercept=c(threshold_2), col="blue", alpha=0.5, linetype='dotdash') +
  geom_hline(yintercept=c(threshold_2), col="blue", alpha=0.5, linetype='dotdash') + 
  theme_light() +
  xlim(0, 0.25) +
  ylim(0, 0.25) +
  ylab('bulk.FDR') +
  xlab('sc.FDR') +
  geom_text_repel(size=2.5, colour='black')
ggsave(plot = p, filename = paste0('figs/VBASS.bulk.VBASS.ASD.pdf'),
       width = 6, height = 6)


