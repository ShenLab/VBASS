# 1 is result folder
# 2 is whether to use log_BF or y_logits
# 3 is whether to use certain number of epoch
library(ggplot2)
library(ggpubr)
args = commandArgs(trailingOnly=TRUE)
result.folder <- 'train.simulation.out/'
out.folder <- paste0('figs/')
train.input <- read.csv('VBASS.simulation.with.RDS/input.x.1.csv', row.names = 1)
train.input.var <- read.csv('VBASS.simulation.with.RDS/input.x_var.1.csv', row.names = 1)
point.estimate.train.input <- train.input[,1:80] / train.input.var[,1:80]
exact.train.input <- train.input[,1:80]

input.dnv.table <- cbind(train.input[,c('dn.cls1', 'dn.cls2')],
                         train.input.var[,c('mut.cls1', 'mut.cls2')])

all.pi <- c()
for (i in 1:50) {
  log_logits <- read.csv(paste0(result.folder, i, '.out/y_logits.csv'), row.names = 1)
  log_logits$log_logits <- log_logits$X1 - log_logits$X0
  log_logits$pi <- exp(log_logits$X1)
  all.pi <- c(all.pi, log_logits$pi)
}
all.pi.matrix <- matrix(all.pi, dim(log_logits)[1], 50)
mean.pi <- rowMeans(all.pi.matrix)
real.pi <- read.csv('VBASS.simulation.with.RDS/real.pi.csv', row.names = 1)
train.x.label <- read.csv('VBASS.simulation.with.RDS/input.label.1.csv', row.names = 1)
true.x.label <- read.csv('VBASS.simulation.with.RDS/test.label.1.csv', row.names = 1)

to.plot <- data.frame(mean.pi = mean.pi, real.pi = real.pi$pi)
to.plot$group <- NA
to.plot$group[true.x.label$label==1] = "TP"
to.plot$group[true.x.label$label==0] = "TN"
to.plot$group[!is.na(train.x.label$label) & train.x.label$label==1] = "TP.train"
to.plot$group[!is.na(train.x.label$label) & train.x.label$label==0] = "TN.train"
to.plot$fake_group <- "ALL"
p <- ggplot(to.plot, aes(x=mean.pi, y=real.pi, col=group)) +
  geom_point(alpha=0.3) +
  stat_smooth(method = "lm", formula = y~x, aes(col=fake_group)) +
  stat_regline_equation(
    aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~"), col=fake_group),
    formula = y~x
  ) + ylim(-0.05, 0.6) + theme_bw()
ggsave(plot = p, filename = paste0(out.folder, "fig.3.C.pdf"), height = 5, width = 5)

correlations <- c()
cell_type <- c()
label <- c()
for (i in 1:50) {
  correlation_array <- as.data.frame(cor(point.estimate.train.input, all.pi.matrix[,i], method='spearman'))
  correlations <- c(correlations, correlation_array$V1)
  cell_type <- c(cell_type, rownames(correlation_array))
  label <- c(label, rep('Spearman Correlation', dim(correlation_array)[1]))
}
real.correlation <- as.data.frame(cor(point.estimate.train.input, real.pi$pi, method = 'spearman'))
correlations <- c(correlations, real.correlation$V1)
cell_type <- c(cell_type, rownames(real.correlation))
label <- c(label, rep('Real Spearman Correlation', dim(real.correlation)[1]))

to.plot <- data.frame(correlations,
                      cell_type,
                      label)
# correlation_array <- correlation_array[order(-correlation_array$V1), ]
library(ggplot2)
library(ggrepel)
correlation_array <- as.data.frame(cor(point.estimate.train.input, mean.pi, method='spearman'))
to.plot <- data.frame(Simulated.correlation=correlation_array$V1,
                      Real.correlation=real.correlation$V1,
                      cell_type=rownames(real.correlation))
# remove unknown cell types
to.plot <- to.plot[grep("Unk", to.plot$cell_type, invert = T),]
p <- ggplot(to.plot, aes(x=Simulated.correlation, y=Real.correlation, label=cell_type)) +
  geom_point(alpha=0.3) +
  stat_smooth(method = "lm", formula = y~x) +
  geom_text_repel(size=2.5, colour='black') +
  stat_regline_equation(
    aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),
    formula = y~x
  ) + theme_bw()
ggsave(plot = p, filename = paste0(out.folder, "fig.3.D.pdf"), height = 5, width = 5)



