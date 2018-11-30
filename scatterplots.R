library(latex2exp)
library(ggplot2)
library(isotone)
library(gridExtra)
library(grid)
source('utilities.R')

read.functions <- c(read.GSE46980, 
                    read.GSE70580, 
                    read.GSE95601,
                    read.GSE86977,
                    read.EMTAB2805,
                    read.EMTAB3929,
                    read.EMTAB5522)

titles = c('GSE46980', 
           'GSE70580', 
           'GSE95601',
           'GSE86977',
           'E-MTAB-2805',
           'E-MTAB-3929',
           'E-MTAB-5522')

plots <- c()
plots_corrected <- c()

for (i in 1:7) {
  dat <- read.functions[[i]]()
  endo.counts <- dat$endo.counts
  ercc.counts <- dat$ercc.counts
  endo.means <- rowMeans(endo.counts)
  endo.cvs <- (rowVars(endo.counts)-endo.means)/endo.means^2+1
  ercc.means <- rowMeans(ercc.counts)
  ercc.cvs <- (rowVars(ercc.counts)-ercc.means)/ercc.means^2+1
  df1 <- data.frame(x=log10(endo.means), y=log10(endo.cvs), origin='endogenous')
  df2 <- data.frame(x=log10(ercc.means), y=log10(ercc.cvs), origin='spike-in')
  df <- rbind(df1, df2)
  df <- df[df$x > 1, ]
  df$rank <- 1
  df$rank[order(df$x)] = 1:length(df$x)
  plots[[i]] <- ggplot() + 
                theme(text = element_text(size = 20), 
                      legend.text = element_text(size=20), 
                      axis.text = element_text(size=20),
                      legend.position = c(0.5,0.88))+
                geom_point(data = df, mapping=aes(x=rank, y=y, color=origin)) +
                xlab(TeX("Gene index")) + 
                ylab(TeX("$ln\\left[\\frac{Var(X_g)-E(X_g)}{E(X_g)^2}+1\\right]$")) + 
                scale_color_manual(name=NULL, values=c('black', 'red'))+
                ggtitle(titles[[i]])
  
  res <- gpava(-df2$x/2, df2$y, ties = "primary", solver = weighted.fractile, p = 0.5)
  fit.spline<-smooth.spline(df2$x/2, (res$x), df = 11)
  df$y_corrected <- df$y - predict(fit.spline, df$x/2)$y
  plots_corrected[[i]] <- ggplot() + 
    theme(text = element_text(size = 20), 
          legend.text = element_text(size=20), 
          axis.text = element_text(size=20),
          legend.position = c(0.5,0.88))+
    geom_point(data = df, mapping=aes(x=rank, y=y_corrected, color=origin)) +
    xlab(TeX("Gene index")) + 
    ylab(TeX("$ln\\left[\\frac{Var(X_g)-E(X_g)}{E(X_g)^2}+1\\right]$")) + 
    scale_color_manual(name=NULL, values=c('black', 'red'))+
    ggtitle(paste0(titles[[i]], ' (corrected)'))
}

for (i in 1:7) {
  grid.arrange(plots[[i]], plots_corrected[[i]], nrow=1)
}
