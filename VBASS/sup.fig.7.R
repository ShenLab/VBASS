library(ggplot2)
dnv_table <- read.csv(file = 'figs/table.S6.csv', row.names = 1)
# plot PPA distributions
dnv_table$trainining.label <- "NA"
dnv_table$trainining.label[dnv_table$train.labels==0] <- "negative"
dnv_table$trainining.label[dnv_table$train.labels==1] <- "positive"
p1 <- ggplot(dnv_table, aes(x=PPA, col=trainining.label)) +
  geom_density() + ggtitle('VBASS') + theme_bw() + ggeasy::easy_center_title()
# ggsave(p, filename = 'figs/VBASS.PPA.density.pdf', width = 6, height = 3)
p2 <- ggplot(dnv_table, aes(x=TADA.PPA, col=trainining.label)) +
  geom_density() + ggtitle('extTADA') + theme_bw() + ggeasy::easy_center_title() + xlab('PPA')
# ggsave(p, filename = 'figs/TADA.PPA.density.pdf', width = 6, height = 3)
library(patchwork)
p <- p1 / p2
ggsave(p, filename = 'figs/sup.fig.7.pdf', width = 8, height = 6)
write.csv(dnv_table, file = 'figs/sup.fig.7.csv')