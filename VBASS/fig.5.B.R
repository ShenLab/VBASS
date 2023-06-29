result.folder <- 'train.real/'
out.folder <- paste0('figs/')
train.input <- read.csv('VBASS.input/input.x.csv', row.names = 1)
train.input.var <- read.csv('VBASS.input/input.x_var.csv', row.names = 1)
point.estimate.train.input <- train.input[,1:80] / train.input.var[,1:80]

all.train.input <- rowSums(train.input[,1:80]) / rowSums(train.input.var[,1:80])
all.train.input.1 <- rowSums(train.input[,startsWith(colnames(train.input), "week")]) / 
  rowSums(train.input.var[,startsWith(colnames(train.input), "week")])
all.train.input.2 <- rowSums(train.input[,startsWith(colnames(train.input), "GW")]) / 
  rowSums(train.input.var[,startsWith(colnames(train.input), "GW")])

exact.train.input <- train.input[,1:80]

# input.dnv.table <- cbind(train.input[,c('dn.cls1', 'dn.cls2')],
#                          train.input.var[,c('mut.cls1', 'mut.cls2')])
input.dnv.table <- cbind(train.input[,c('dn_LGD', 'dn_Dmis')],
                         train.input.var[,c('mut_LGD', 'mut_Dmis')])


log_logits <- read.csv(paste0(result.folder, 'y_logits.csv'), row.names = 1)

log_logits$log_logits <- log_logits$X1 - log_logits$X0
log_logits$pi <- 1/(1+exp(-log_logits$log_logits))

correlation_array <- as.data.frame(cor(point.estimate.train.input, log_logits$pi, method='spearman'))
all.corr <- cor(all.train.input, log_logits$pi, method = 'spearman')
all.corr.1 <- cor(all.train.input.1, log_logits$pi, method = 'spearman')
all.corr.2 <- cor(all.train.input.2, log_logits$pi, method = 'spearman')
correlation_array$cell_type <- rownames(correlation_array)
# remove Unk cell types
correlation_array <- correlation_array[grep("Unk", correlation_array$cell_type, invert = T),]

cell_type_1 <- correlation_array[startsWith(correlation_array$cell_type, "week"), ]
cell_type_1_V <- correlation_array$V1[startsWith(correlation_array$cell_type, "week")]
cell_type_2 <- correlation_array[startsWith(correlation_array$cell_type, "GW"), ]
cell_type_2_V <- correlation_array$V1[startsWith(correlation_array$cell_type, "GW")]
level_1 <- cell_type_1$cell_type[order(cell_type_1_V, decreasing = T)]
level_2 <- cell_type_2$cell_type[order(cell_type_2_V, decreasing = T)]
# correlation_array <- correlation_array[order(-correlation_array$V1), ]
library(ggplot2)
level <- c(level_1, level_2)
correlation_array <- correlation_array[match(level, correlation_array$cell_type),]
correlation_array.1 <- correlation_array[1:length(level_1),]
correlation_array.2 <- correlation_array[(length(level_1)+1):(length(level_1) + length(level_2)),]
p <- ggplot(correlation_array.1, aes(x=factor(cell_type, levels = level), y=V1)) +
  geom_point() +
  ylab('Spearman Correlation') +
  theme_bw() +
  xlab('cell types') +
  ylim(0.275, 0.50) +
  geom_abline(intercept = all.corr.1, slope = 0, col = 'blue', linetype="dotted") +
  annotate(geom="text", label=paste0("bulk dataset 1 corr:", format(all.corr.1, digits = 3)), x=8, y=all.corr.1, vjust=1, hjust=-1, col = "blue") +
  # geom_abline(intercept = all.corr.2, slope = 0, col = 'blue', linetype="dotted") +
  # annotate(geom="text", label=paste0("bulk dataset 2 corr:", format(all.corr.2, digits = 3)), x=8, y=all.corr.2, vjust=1, col = "blue") +
  theme(axis.text.x = element_text(angle = 90))

ggsave(plot = p, filename = paste0(out.folder, "fig.5.B.1.pdf"), height = 4, width = 6)
write.csv(correlation_array.1, file = paste0(out.folder, 'fig.5.B.1.csv'))

p <- ggplot(correlation_array.2, aes(x=factor(cell_type, levels = level), y=V1)) +
  geom_point() +
  ylab('Spearman Correlation') +
  theme_bw() +
  xlab('cell types') +
  ylim(0.275, 0.50) +
  # geom_abline(intercept = all.corr.1, slope = 0, col = 'blue', linetype="dotted") +
  # annotate(geom="text", label=paste0("bulk dataset 1 corr:", format(all.corr.1, digits = 3)), x=8, y=all.corr.1, vjust=1, col = "blue") +
  geom_abline(intercept = all.corr.2, slope = 0, col = 'blue', linetype="dotted") +
  annotate(geom="text", label=paste0("bulk dataset 2 corr:", format(all.corr.2, digits = 3)), x=8, y=all.corr.2, vjust=1, hjust=-1, col = "blue") +
  theme(axis.text.x = element_text(angle = 90))

ggsave(plot = p, filename = paste0(out.folder, "fig.5.B.2.pdf"), height = 4, width = 6)
write.csv(correlation_array.2, file = paste0(out.folder, 'fig.5.B.2.csv'))
