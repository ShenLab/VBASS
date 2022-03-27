seed <- 0
set.seed(as.numeric(seed))
out.dir <- paste0('AutoEncoderTADA.input.discovery/')
TADA.output <- readRDS(paste0('data/SPARK_extTADA.WES1.discovery.RDS'))
dir.create(out.dir)
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

rownames(DNV.table) <- DNV.table$HGNC
# begin x
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

genes.dnv <- DNV.table[match(genes, DNV.table$HGNC), c('dn_LGD', 'dn_Dmis')]

genes.x <- cbind(genes.zhong_2018, genes.lamanno_2016, genes.dnv)
rownames(genes.x) <- genes
write.csv(genes.x, file = paste0(out.dir, '/input.x.csv'))
# begin x.var
genes.cell_number <- t(matrix(rep(c(zhong_2018$week.cell_types, lamanno_2016$Timepoint.Cell_type), length(genes)),
                              ncol = length(genes)))
rownames(genes.cell_number) <- genes
colnames(genes.cell_number) <- c(colnames(zhong_2018.exp), colnames(lamanno_2016.exp))
genes.mut <- DNV.table[match(genes, DNV.table$HGNC), c('mut_LGD', 'mut_Dmis')]*8308

genes.x_var <- cbind(genes.cell_number, genes.mut)
write.csv(genes.x_var, file = paste0(out.dir, '/input.x_var.csv'))
# begin gene.label
library(dplyr)
label.file <- read.csv('data/SFARI_+_Control.csv')
set.seed(0)
label.file.training <- sample_n(label.file, floor(dim(label.file)[1]*0.85))
label.file.testing <- label.file[!label.file$genes%in%label.file.training$genes,]
genes.label <- data.frame(row.names = genes,
                          label = rep(NA, length(genes)))
genes.label$label[match(label.file.training$HGNC[label.file.training$conditions=="SPARK.risk"], genes)] = 1
genes.label$label[match(label.file.training$HGNC[label.file.training$conditions=="control.LGD"], genes)] = 0
write.csv(genes.label, file = paste0(out.dir, '/input.label.csv'))

write.csv(label.file.testing, file = paste0(out.dir, "/testing.label.csv"))






