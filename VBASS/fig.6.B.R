result.folder <- 'train.real/'
out.folder <- paste0('figs/')
train.input <- read.csv('AutoEncoderTADA.input/input.x.csv', row.names = 1)
train.input.var <- read.csv('AutoEncoderTADA.input/input.x_var.csv', row.names = 1)
point.estimate.train.input <- train.input[,1:80] / train.input.var[,1:80]
exact.train.input <- train.input[,1:80]

# input.dnv.table <- cbind(train.input[,c('dn.cls1', 'dn.cls2')],
#                          train.input.var[,c('mut.cls1', 'mut.cls2')])
input.dnv.table <- cbind(train.input[,c('dn_LGD', 'dn_Dmis')],
                         train.input.var[,c('mut_LGD', 'mut_Dmis')])


log_logits <- read.csv(paste0(result.folder, 'y_logits.csv'), row.names = 1)

log_logits$log_logits <- log_logits$X1 - log_logits$X0
log_logits$pi <- 1/(1+exp(-log_logits$log_logits))

correlation_array <- as.data.frame(cor(point.estimate.train.input, log_logits$pi, method='spearman'))
correlation_array$cell_type <- rownames(correlation_array)
cell_type_1 <- correlation_array[startsWith(correlation_array$cell_type, "week"), ]
cell_type_1_V <- correlation_array$V1[startsWith(correlation_array$cell_type, "week")]
cell_type_2 <- correlation_array[startsWith(correlation_array$cell_type, "GW"), ]
cell_type_2_V <- correlation_array$V1[startsWith(correlation_array$cell_type, "GW")]
level_1 <- cell_type_1$cell_type[order(cell_type_1_V, decreasing = T)]
level_2 <- cell_type_2$cell_type[order(cell_type_2_V, decreasing = T)]
# correlation_array <- correlation_array[order(-correlation_array$V1), ]
library(ggplot2)
p <- ggplot(cell_type_1, aes(x=factor(cell_type, levels = level_1), y=V1)) +
  geom_point() +
  ylab('Spearman Correlation') +
  theme_bw() +
  xlab('cell types') +
  theme(axis.text.x = element_text(angle = 90))

ggsave(plot = p, filename = paste0(out.folder, "fig.6.B.1.pdf"), height = 4, width = 8)
p <- ggplot(cell_type_2, aes(x=factor(cell_type, levels = level_2), y=V1)) +
  geom_point() +
  ylab('Spearman Correlation') +
  theme_bw() +
  xlab('cell types') +
  theme(axis.text.x = element_text(angle = 90))

ggsave(plot = p, filename = paste0(out.folder, "fig.6.B.2.pdf"), height = 4, width = 8)

