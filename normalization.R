source('utilities.R')
library(isotone)
library(ggplot2)

read.functions <- c(read.GSE46980, read.GSE70580, read.GSE95601)
titles = c('GSE46980', 'GSE70580', 'GSE95601')
plots <- c()

for (i in 1:3) {
  dat <- read.functions[[i]]()
  endo.counts <- dat$endo.counts
  ercc.counts <- dat$ercc.counts
  
  # obtain lowly variable genes (lvg)
  lvg <- detect.lvg(endo.counts)
  
  # normalization size factors 
  s.lvg <- norm.libsize(endo.counts[lvg, ])$s
  s.libsize <- norm.libsize(endo.counts)$s
  s.scran <- norm.scran(endo.counts)$s
  
  # plot the result
  spike.list <- split(ercc.counts, 1:nrow(ercc.counts)) # boxplot requires a list as input
  names(spike.list) <- rownames(ercc.counts)
  spike.mean <- sapply(spike.list, mean)
  spike.list <- spike.list[order(spike.mean, decreasing = T)] # order the list according to average count
  num.genes <- 10
  tmp <- lapply(spike.list[1:num.genes], function(x) log(x/s.lvg))
  df1 <- data.frame(norm = 'Proposed', y = do.call(c, tmp), x=rep(names(spike.list[1:num.genes]), each=length(s.lvg)))
  tmp <- lapply(spike.list[1:num.genes], function(x) log(x/s.libsize))
  df2 <- data.frame(norm = 'library size', y = do.call(c, tmp), x=rep(names(spike.list[1:num.genes]), each=length(s.libsize)))
  tmp <- lapply(spike.list[1:num.genes], function(x) log(x/s.scran))
  df3 <- data.frame(norm = 'deconvolution', y = do.call(c, tmp), x=rep(names(spike.list[1:num.genes]), each=length(s.scran)))
  df <- rbind(df1, df2, df3)
  df$x <- factor(df$x, levels = names(sort(spike.mean, decreasing = T))[1:10])
  plots[[i]] <- ggplot(df, aes(x=x, y=y, fill=norm)) +
                geom_boxplot() + 
                theme(axis.text.x = element_text(angle = 70, hjust = 1))+
                xlab('spike-ins')+
                ylab('normalized counts')+
                ggtitle(titles[[i]])
}

plots[[1]]
plots[[2]]
plots[[3]]
