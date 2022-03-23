TADA.output <- readRDS('SPARK_extTADA.WES1.RDS')
dir.create('AutoEncoderTADA.input/')
DNV.table <- TADA.output$dataFDR
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
DNV.table <- remove_duplicate_line(DNV.table)
genes <- DNV.table$HGNC


zhong_2018 <- read.csv('data/Zhong_2018.cell.number.csv', row.names = 1)
zhong_2018.exp <- read.csv('data/week:cell_types.fraction.num.csv', row.names = 1)

lamanno_2016 <- read.csv('data/LaManno_2016.cell.number.csv', row.names = 1)
lamanno_2016.exp <- read.csv('data/Timepoint:Cell_type.fraction.num.csv', row.names = 1)


genes.zhong_2018 <- zhong_2018.exp[match(genes, rownames(zhong_2018.exp)),]
print(sum(is.na(genes.zhong_2018))/24)
genes.zhong_2018[is.na(genes.zhong_2018)] <- 0
genes.lamanno_2016 <- lamanno_2016.exp[match(genes, rownames(lamanno_2016.exp)),]
print(sum(is.na(genes.lamanno_2016))/56)
genes.lamanno_2016[is.na(genes.lamanno_2016)] <- 0

genes.dnv <- DNV.table[, c('dn_LGD', 'dn_Dmis')]

genes.x <- cbind(genes.zhong_2018, genes.lamanno_2016, genes.dnv)
rownames(genes.x) <- genes
write.csv(genes.x, file = 'AutoEncoderTADA.input/input.x.csv')
# begin x.var
genes.cell_number <- t(matrix(rep(c(zhong_2018$week.cell_types, lamanno_2016$Timepoint.Cell_type), length(genes)),
                              ncol = length(genes)))
rownames(genes.cell_number) <- genes
colnames(genes.cell_number) <- c(colnames(zhong_2018.exp), colnames(lamanno_2016.exp))
genes.mut <- DNV.table[, c('mut_LGD', 'mut_Dmis')]*16616

genes.x_var <- cbind(genes.cell_number, genes.mut)
rownames(genes.x_var) <- genes
write.csv(genes.x_var, file = 'AutoEncoderTADA.input/input.x_var.csv')
# begin gene.label
library(dplyr)
label.file <- read.csv('data/SFARI_+_Control.csv')
# label.file.genes <- ensembldb::select(EnsDb.Hsapiens.v86, keys= as.character(label.file$HGNC),
#                                       keytype = "SYMBOL", columns = c("SYMBOL","GENEID"))
# label.file$genes <- label.file.genes$GENEID
set.seed(0)
label.file.training <- sample_n(label.file, floor(dim(label.file)[1]*0.85))
label.file.testing <- label.file[!label.file$genes%in%label.file.training$genes,]
genes.label <- data.frame(row.names = genes,
                          label = rep(NA, length(genes)))
genes.label$label[match(label.file.training$genes[label.file.training$conditions=="SPARK.risk"], genes)] = 1
genes.label$label[match(label.file.training$genes[label.file.training$conditions=="control.LGD"], genes)] = 0
write.csv(genes.label, file = 'AutoEncoderTADA.input/input.label.csv')

write.csv(label.file.testing, file = "AutoEncoderTADA.input/testing.label.csv")

