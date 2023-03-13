sc.input <- read.csv('VBASS.input/input.x.csv', row.names = 1)
sc.input <- sc.input[,1:80]
sc.number <- read.csv('VBASS.input/input.x_var.csv', row.names = 1)
sc.number <- sc.number[,1:80]
sc.input <- sc.input / sc.number
gene.labels <- read.csv('VBASS.input/input.label.csv', row.names = 1)
posterior <- readxl::read_excel('figs/Supplementary Table.xlsx', sheet='S6')
colnames(posterior) <- posterior[14,]
posterior <- posterior[15:dim(posterior)[1], ]

pi <- readxl::read_excel('figs/Supplementary Table.xlsx', sheet='S5')
colnames(pi) <- pi[13,]
pi <- pi[14:dim(pi)[1], ]

posterior$FDR.0.05 <- as.numeric(posterior$FDR) <= 0.05
posterior$FDR.0.1 <- as.numeric(posterior$FDR) <= 0.1
posterior <- posterior[match(rownames(gene.labels), posterior$HGNC),]
gene.labels$posterior.label <- NA
gene.labels$posterior.label[posterior$FDR.0.1] <- "FDR <= 0.1"
# gene.labels$log_pi <- 
# gene.labels$posterior.label[posterior$FDR.0.05] <- "FDR <= 0.05"

gene.labels <- gene.labels[!is.na(gene.labels$label) | !is.na(gene.labels$posterior.label),]
# visualize training genes
training.genes.sc <- sc.input[match(rownames(gene.labels), rownames(sc.input)),]
gene.labels$label[gene.labels$label==0] <- "Non-risk"
gene.labels$label[gene.labels$label==1] <- "Risk"
colnames(gene.labels)[1] <- 'training.label'
gene.labels$posterior.label[gene.labels$training.label == "Risk"] <- "Excluded"

gene.labels[is.na(gene.labels)] <- 'NA'
gene.labels$training.label <- factor(gene.labels$training.label, c("Risk", "Non-risk", "NA"))
gene.labels$posterior.label[gene.labels$posterior.label == "NA"] <- "FDR > 0.1"
gene.labels$posterior.label <- factor(gene.labels$posterior.label, c("FDR <= 0.1", "FDR > 0.1", "Excluded"))
# gene.labels$log_pi <- log(posterior$pi)
source('~/Pipeline/plot.genes.scores.heatmap.R')

col.names <- colnames(sc.input)
week <- as.numeric(unlist(regmatches(col.names, gregexpr('[0-9]+\\.', col.names))))
col.anno <- data.frame(row.names = col.names, week=week)
split <- rep(NA, dim(gene.labels)[1])
split[gene.labels$posterior.label %in% c("FDR > 0.1")] <- "Not significant"
split[!gene.labels$posterior.label %in% c("FDR > 0.1")] <- "Significant"
split <- factor(split, levels = c('Significant', 'Not significant'))
plot.complex.heatmap(training.genes.sc, row.anno = gene.labels, 
                     numeric.row.anno = data.frame(log_pi=log(as.numeric(pi$pi[match(rownames(gene.labels), pi$HGNC)]))),
                     col.anno = col.anno, title = 'Expression', split = split,
                     width = 10, height = 15, col_font_size = 6,
                     savefilename = 'figs/sup.fig.3.pdf')
