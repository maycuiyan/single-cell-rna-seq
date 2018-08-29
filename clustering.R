library(GEOquery)
library(matrixStats)
library(NMF)
library(Seurat)
library(Matrix)
library(tsne)
require(mclust)
library(SingleCellExperiment)
library(scran)
library(isotone)
library(readr)
library(statmod)
library(org.Hs.eg.db)
library(VennDiagram)
library(gridExtra)
source('utilities.R')

counts.all <- list()
group.all <- list()

##--------Li et al. dataset---------
dat <- readRDS("./data/li.rds") ## good
annot <- dat@colData
counts.all[[1]] <- counts(dat)
annot$cell_type1 <- gsub('_B[12]', '', as.character(annot$cell_type1))
group.all[[1]] <- as.factor(annot$cell_type1)

##--------Klein et al. dataset---------
dat <- readRDS("./data/klein.rds") ## good
annot <- dat@colData
counts.all[[2]] <- counts(dat)
group.all[[2]] <- as.factor(annot$cell_type1)

##--------Bjorklund et al. dataset---------
filenames <- dir('./data/GSE70580_RAW/')
counts <- list()
for (i in 1:length(filenames)) {
  count <- read_delim(paste0('./data/GSE70580_RAW/', filenames[i]), "\t", escape_double = FALSE, trim_ws = TRUE)
  counts[[i]] <- count$Reads
}
counts <- do.call(cbind, counts)
rownames(counts) <- count$`#Gene symbol`
colnames(counts) <- gsub('_.*', '', filenames)
rownames(counts) <- toupper(rownames(counts))
symbol2entrez <- as.list(org.Hs.egALIAS2EG)
symbol2entrez <- symbol2entrez[!is.na(symbol2entrez)]
endo.counts <- counts[rownames(counts) %in% names(symbol2entrez), ] # keep genes with identified symbols
gene.freq <- table(rownames(endo.counts))
unique.genes <- names(gene.freq)[gene.freq==1]
endo.counts <- endo.counts[unique.genes,]
counts.all[[3]] <- endo.counts
gset <- getGEO('GSE70580')
annot <- pData(gset[[1]])
group.all[[3]] <- as.factor(annot$`facs gating:ch1`)


num.genes <- 1000
titles = c('Li', 'Klein', 'Bjorklund')
silhouettes = c()

for (i in c(1:3)) {
  
  counts <- counts.all[[i]]
  group <- group.all[[i]]
  genes.keep <- which(rowMeans(counts) >= 1)
  counts <- as.matrix(counts[genes.keep, ])
  normcounts <- norm.libsize(counts)$normcounts
  
  colors = rainbow(length(unique(group)))
  names(colors) = unique(group)
  ecb = function(x,y){ plot(x,t='n'); text(x,labels=as.numeric(group), col=colors[group]) }
  
  # method 1: proposed
  hvg <- order(hvg.detection(counts), decreasing = T)[1:num.genes]
  X <- log10(t(normcounts[hvg, ])+1)
  tsne = tsne(X)
  si <- silhouette(as.numeric(group), dist(tsne))
  silhouettes = c(silhouettes, mean(si[,3]))
  plot(tsne[,1], tsne[,2], col=colors[group], pch=as.numeric(group)+14,
       xlab='TSNE[:,1]', 
       ylab='TSNE[:,2]', 
       main = paste('Proposed (silhouette = ', round(mean(si[,3]), 3), ')'), 
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, 
       xlim = c(-60, 60))
  legend('bottomright', legend=unique(group), col=colors[unique(group)], pch=as.numeric(unique(group))+14, 
         cex = 1)
  
  # method 2: seurat
  hvg.seurat <- order(hvg.detection.seurat(counts), decreasing = T)[1:num.genes]
  X <- log10(t(normcounts[hvg.seurat, ])+1)
  tsne_seurat = tsne(X)
  si <- silhouette(as.numeric(group), dist(tsne_seurat))
  silhouettes = c(silhouettes, mean(si[,3]))
  plot(tsne_seurat[,1], tsne_seurat[,2], col=colors[group], pch=as.numeric(group)+14,
       xlab='TSNE[:,1]', 
       ylab='TSNE[:,2]', 
       main = paste('Seurat (silhouette = ', round(mean(si[,3]), 3), ')'), 
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
  legend('bottomright', legend=unique(group), col=colors[unique(group)], pch=as.numeric(unique(group))+14, 
         cex = 1)
  
  # method 3: scran
  hvg.scran <- order(hvg.detection.scran(counts))[1:num.genes]
  X <- log10(t(normcounts[hvg.scran, ])+1)
  tsne_scran = tsne(X)
  si <- silhouette(as.numeric(group), dist(tsne_scran))
  silhouettes = c(silhouettes, mean(si[,3]))
  plot(tsne_scran[,1], tsne_scran[,2], col=colors[group], pch=as.numeric(group)+14,
       xlab='TSNE[:,1]', 
       ylab='TSNE[:,2]', 
       main = paste('Scran (silhouette = ', round(mean(si[,3]), 3), ')'), 
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
  legend('bottomright', legend=unique(group), col=colors[unique(group)], pch=as.numeric(unique(group))+14, 
         cex = 1)
  
  # method 4: brennecke
  hvg.brennecke <- order(hvg.detection.brennecke(counts),decreasing=T)[1:num.genes]
  X <- log10(t(normcounts[hvg.brennecke, ])+1)
  tsne_brennecke = tsne(X)
  si <- silhouette(as.numeric(group), dist(tsne_brennecke))
  silhouettes = c(silhouettes, mean(si[,3]))
  plot(tsne_brennecke[,1], tsne_brennecke[,2], col=colors[group], pch=as.numeric(group)+14,
       xlab='TSNE[:,1]', 
       ylab='TSNE[:,2]', 
       main = paste('Brennecke (silhouette = ', round(mean(si[,3]), 3), ')'), 
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
  legend('bottomright', legend=unique(group), col=colors[unique(group)], pch=as.numeric(unique(group))+14, 
         cex = 1)
  
  # venn diagram
  venn.diagram(list(hvg, hvg.seurat, hvg.scran, hvg.brennecke),
               filename = paste0(titles[[i]], '.tiff'),
               fill = c("orange", "red", "green", "blue"),
               category=c('Proposed', 'Seurat', 'Scran', 'Brennecke'),
               main = paste(titles[[i]], 'Dataset'))
}

# barplot for all datasets
data.m <- data.frame(Datasets=c('Li et al.', 
                                'Li et al.', 
                                'Li et al.', 
                                'Li et al.', 
                                'Klein et al.',
                                'Klein et al.',
                                'Klein et al.',
                                'Klein et al.',
                                'Bjorklund et al.',
                                'Bjorklund et al.',
                                'Bjorklund et al.',
                                'Bjorklund et al.'), 
                     Methods=c('Proposed', 'Seurat', 'Scran', 'Brennecke',
                               'Proposed', 'Seurat', 'Scran', 'Brennecke',
                               'Proposed', 'Seurat', 'Scran', 'Brennecke'), 
                     Silhouette=silhouettes)


data.m$Datasets <- factor(data.m$Datasets, levels = c('Li et al.', 'Klein et al.', 'Bjorklund et al.'))
data.m$Methods <- factor(data.m$Methods, levels = c('Proposed', 'Seurat', 'Scran', 'Brennecke'))

ggplot(data.m, aes(Datasets, Silhouette)) +   
  theme(axis.text = element_text(size=20), text = element_text(size = 20))+
  geom_bar(aes(fill = Methods), position = "dodge", stat="identity")
