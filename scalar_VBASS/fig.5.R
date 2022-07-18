VBASS_result <- readRDS('result/CHD_VBASS_expression.RDS')
source('VBASS/draw.VBASS.S_distribution.R')
draw.VBASS.S_distribution.from.pars0(VBASS_result$pars0, savefilename = "figs/fig.5.A.pdf")

dnv_table <- read.csv('result/table.S4.csv')
threshold_1 <- 0.1
threshold_2 <- 0.1
dnv_table$group[dnv_table$FDR<=threshold_1 & dnv_table$extTADA.FDR<=threshold_2] = 'both'
dnv_table$group[dnv_table$FDR>threshold_1 & dnv_table$extTADA.FDR<=threshold_2] = 'extTADA_only'
dnv_table$group[dnv_table$FDR<=threshold_1 & dnv_table$extTADA.FDR>threshold_2] = 'VBASS_only'
dnv_table$group[dnv_table$FDR>threshold_1 & dnv_table$extTADA.FDR>threshold_2] = 'NA'
dnv_table$label <- dnv_table$HGNC
dnv_table$label[dnv_table$group %in% c('NA', 'both')] <- ''
library(ggplot2)
library(ggrepel)
p <- ggplot(dnv_table, aes(x=FDR, y=extTADA.FDR, col=group, label=label)) +
  geom_point() +
  geom_vline(xintercept=c(threshold_1), col="red", alpha=0.5, linetype='dotdash') +
  geom_hline(yintercept=c(threshold_2), col="red", alpha=0.5, linetype='dotdash') + 
  theme_light() +
  xlab('VBASS.FDR') +
  geom_text_repel(size=2, colour='black', max.overlaps=50)
ggsave(plot = p, filename = paste0('figs/fig.5.C.pdf'), width = 6, height = 5)
