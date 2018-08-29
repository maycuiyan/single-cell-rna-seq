library(latex2exp)
library(ggplot2)
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
           'E-MTAB5522')
plots <- c()

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
  
  plots[[i]] <- ggplot() + 
                theme(text = element_text(size = 20), 
                      legend.text = element_text(size=20), 
                      axis.text = element_text(size=20),
                      legend.position = c(0.8,0.88))+
                geom_point(data = df, mapping=aes(x=x, y=y, color=origin)) +
                xlim(0, 5) +
                xlab(TeX("$ln\\[E(x_g)\\]$")) + 
                ylab(TeX("$ln\\left[\\frac{Var(X_g)-E(X_g)}{E(X_g)^2}+1\\right]$")) + 
                scale_color_manual(name=NULL, values=c('black', 'red'))+
                ggtitle(titles[[i]])
}


plots[[1]]
plots[[2]]
plots[[3]]
plots[[4]]
plots[[5]]
plots[[6]]
plots[[7]]
