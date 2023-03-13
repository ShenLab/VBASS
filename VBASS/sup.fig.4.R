library(patchwork)
# Using all genes PCA, use K-means to cluster genes
sc.input <- read.csv('VBASS.input/input.x.csv', row.names = 1)
sc.input <- sc.input[,1:80]
sc.number <- read.csv('VBASS.input/input.x_var.csv', row.names = 1)
sc.number <- sc.number[,1:80]
sc.input <- sc.input / sc.number

gene.labels <- read.csv('VBASS.input/input.label.csv', row.names = 1)

pi <- read.csv('train.real/y_logits.csv', row.names = 1)
# training.genes.sc <- sc.input[!is.na(gene.labels$label),]
# 
# training.genes.sc.pca <- summary(prcomp(training.genes.sc, scale. = F))

# all.genes.pca <- (as.matrix(sc.input) - matrix(1, nrow = dim(sc.input)[1], ncol = 1) 
#                   %*% training.genes.sc.pca$center) %*% training.genes.sc.pca$rotation
# sum(all.genes.pca[!is.na(gene.labels$label),] - training.genes.sc.pca$x)

all.genes.pca <- summary(prcomp(sc.input, scale. = F))$x

library(Seurat)
fake_sc <- CreateSeuratObject(t(sc.input))
fake_sc <- FindVariableFeatures(fake_sc)
fake_sc <- ScaleData(fake_sc)
fake_sc <- RunPCA(fake_sc)
fake_sc <- FindNeighbors(fake_sc)
fake_sc <- FindClusters(fake_sc, resolution = 0.1)

i = 2

all.genes.pca.df <- data.frame(all.genes.pca[,i:i+1])
all.genes.pca.df$louvain_clust <- fake_sc@meta.data$seurat_clusters


# do a k-means clustering
wss <- function(k) {
  kmeans(all.genes.pca[,i:i+1], k, nstart = 50, iter.max = 100)$tot.withinss
}
k.values <- 1:15
wss_values <- purrr::map_dbl(k.values, wss)
plot(k.values, wss_values,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")

result <- kmeans(all.genes.pca[,i:i+1], 3, nstart = 50, iter.max = 100)
all.genes.pca.df$k_means_clust <- as.factor(result$cluster)


pca.df <- data.frame(PC0 = all.genes.pca[,i-1],
                     PC1=all.genes.pca[,i],
                     PC2=all.genes.pca[,i+1],
                     training.label=as.character(gene.labels$label),
                     louvain.clust=all.genes.pca.df$louvain_clust,
                     k.means.clust=all.genes.pca.df$k_means_clust,
                     log_pi=pi$X1)
pca.df$training.label[pca.df$training.label=="0"] <- "Non-Risk"
pca.df$training.label[pca.df$training.label=="1"] <- "ASD-Risk"

posterior <- readxl::read_excel('figs/Supplementary Table.xlsx', sheet='S6')
colnames(posterior) <- posterior[14,]
posterior <- posterior[15:dim(posterior)[1], ]
pca.df$posterior.label <- posterior$group[match(rownames(pca.df), posterior$HGNC)]
pca.df$posterior.label[pca.df$posterior.label=="NA"] <- NA

p1 <- ggplot(pca.df, aes(x=PC1, y=PC2, col=k.means.clust)) +
  geom_point(alpha=0.7, size=0.1) + xlab(paste0("PC ", i)) + ylab(paste0("PC ", i+1)) +
  theme_bw()
# p1 <- ggExtra::ggMarginal(p1, groupColour = T, type="density")
# ggsave(paste0("review.figs/review.fig.5.k.means.clust.1.pdf"), plot = p1, width = 6, height = 5)

pca.result <- prcomp(sc.input, scale. = F)
summary.pca <- summary(pca.result)
pca.loadings <- data.frame(Variables = row.names(pca.result$rotation),
                           PC0 = pca.result$rotation[,i-1],
                           PC1 = pca.result$rotation[,i],
                           PC2 = pca.result$rotation[,i+1])
arrow.length.1 <- max(abs(c(pca.df$PC0, pca.df$PC1))) / max(abs(c(pca.loadings$PC0, pca.loadings$PC1)))/2
mid.point <- log(mean(exp(pca.df$log_pi)))
p0.5 <- ggplot(data = pca.df, aes(x=PC0, y=PC1)) +
  geom_point(alpha=0.5, size=1, stroke=0, col='grey') +
  # xlim(-7.5, 2) + ylim(-3.5, 5) +
  theme_bw() + 
  xlab(paste("PC", i-1, ":", summary.pca$importance[2,i-1], sep = "")) + 
  ylab(paste("PC", i, ":", summary.pca$importance[2,i], sep = "")) + 
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = mid.point) +
  geom_segment(data = pca.loadings, aes(x = 0, y = 0, xend = (PC0*arrow.length.1),
                                        yend = (PC1*arrow.length.1)), arrow = arrow(length = unit(1/2, "picas")),
               color = "red", alpha = 0.3) +
  geom_text_repel(data = pca.loadings, 
                  aes(x = (PC0*arrow.length.1), y = (PC1*arrow.length.1), label = Variables),
                  max.overlaps = 70,
                  size = 0.8)
arrow.length <- max(abs(c(pca.df$PC1, pca.df$PC2))) / max(abs(c(pca.loadings$PC1, pca.loadings$PC2)))/2
mid.point <- log(mean(exp(pca.df$log_pi)))
p1 <- ggplot(data = pca.df, aes(x=PC1, y=PC2)) +
  geom_point(alpha=0.5, size=1, stroke=0, col='grey') +
  # xlim(-7.5, 2) + ylim(-3.5, 5) +
  theme_bw() + 
  xlab(paste("PC", i, ":", summary.pca$importance[2,i], sep = "")) + 
  ylab(paste("PC", i+1, ":", summary.pca$importance[2,i+1], sep = "")) + 
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = mid.point) +
  geom_segment(data = pca.loadings, aes(x = 0, y = 0, xend = (PC1*arrow.length),
                                        yend = (PC2*arrow.length)), arrow = arrow(length = unit(1/2, "picas")),
               color = "red", alpha = 0.3) +
  geom_text_repel(data = pca.loadings, 
                  aes(x = (PC1*arrow.length), y = (PC2*arrow.length), label = Variables),
                  size = 0.8)
ggsave(paste0("figs/sup.fig.4.B.pdf"), plot = p0.5 + p1, width = 6, height = 3)

p1.5 <- ggplot(pca.df, aes(x=PC0, y=log_pi)) +
  geom_point(alpha=0.7, size=0.1, col='light blue') + xlab(paste0("PC", i-1)) + ylab(paste0("log(pi)")) +
  geom_smooth(fill='blue') +
  theme_bw()
p2 <- ggplot(pca.df, aes(x=PC1, y=log_pi)) +
  geom_point(alpha=0.7, size=0.1, col='light blue') + xlab(paste0("PC", i)) + ylab(paste0("log(pi)")) +
  geom_smooth() +
  theme_bw()
p3 <- ggplot(pca.df, aes(x=PC2, y=log_pi)) +
  geom_point(alpha=0.7, size=0.1, col='light blue') + xlab(paste0("PC", i+1)) + ylab(paste0("log(pi)")) +
  geom_smooth() +
  theme_bw()
# ggsave(paste0("review.figs/review.fig.5.pi.k.means.1.pdf"), plot = p2, width = 6, height = 3)
p <- p1.5 + p2 + p3
ggsave("figs/sup.fig.4.A.pdf", plot = p, width = 10, height = 3)




