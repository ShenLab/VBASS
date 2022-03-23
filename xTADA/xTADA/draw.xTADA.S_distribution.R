draw.xTADA.S_distribution <- function(A, B, C, pi0, savefilename, title=NULL, ylim=NULL, width=5, height=3) {
  library(ggplot2)
  x = seq(0, 1, 0.001)
  L = (1-C) * A / (log(exp(A)+exp(B*A)) - log(exp(B*A)+1))
  y1 = 1 + pi0/(1-pi0)*(1-C-L/(1+exp(-A*(x-B))))
  y2 = C + L / (1+exp(-A*(x-B)))
  
  to.plot <- data.frame(x=c(x, x), y=c(y1, y2),
                        Hypothesis=c(rep("Null", length(x)), rep("Alternative", length(x))))
  if (is.null(title)) {
    title = "xTADA S distribution"
  }
  p <- ggplot(to.plot, aes(x=x, y=y, col=Hypothesis)) +
    geom_line() + 
    xlim(0, 1) +
    ylim +
    theme_bw() + theme(panel.border = element_blank(),
                       panel.grid.major = element_blank(),
                       axis.line = element_line(colour = "black"),
                       text = element_text(size=10)) +
    xlab('Covariates Rank Percentile') + ylab('Probability Density Function') +
    ggtitle(title) + ggeasy::easy_center_title()
  ggsave(savefilename, plot = p, width=width, height=height)
}


draw.xTADA.S_distribution.from.pars0 <- function(pars0, savefilename, title=NULL, ylim=NULL, width=5, height=3) {
  library(ggplot2)
  x = seq(0, 1, 0.001)
  A <- pars0[2,1]
  B <- pars0[3,1]
  C <- pars0[4,1]
  pi0 <- pars0[1,1]
  L = (1-C) * A / (log(exp(A)+exp(B*A)) - log(exp(B*A)+1))
  y1 = 1 + pi0/(1-pi0)*(1-C-L/(1+exp(-A*(x-B))))
  y2 = C + L / (1+exp(-A*(x-B)))
  
  to.plot <- data.frame(x=c(x, x), y=c(y1, y2),
                        Hypothesis=c(rep("Null", length(x)), rep("Alternative", length(x))))
  if (is.null(title)) {
    title = "xTADA S distribution"
  }
  p <- ggplot(to.plot, aes(x=x, y=y, col=Hypothesis)) +
    geom_line() + 
    xlim(0, 1) +
    ylim +
    theme_bw() + theme(panel.border = element_blank(),
                       panel.grid.major = element_blank(),
                       axis.line = element_line(colour = "black"),
                       text = element_text(size=10)) +
    xlab('Covariates Rank Percentile') + ylab('Probability Density Function') +
    ggtitle(title) + ggeasy::easy_center_title()
  ggsave(savefilename, plot = p, width=width, height=height)
}
