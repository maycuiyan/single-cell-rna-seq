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
library(dplyr)
library(org.Hs.eg.db)
library(VennDiagram)
library(gridExtra)
source('utilities.R')


counts.all <- list()
group.all <- list()

##--------Li et al. dataset---------
dat <- readRDS("./data/li.rds")
annot <- dat@colData
counts.all[[1]] <- counts(dat)
annot$cell_type1 <- gsub('_B[12]', '', as.character(annot$cell_type1))
group.all[[1]] <- as.factor(annot$cell_type1)

##--------Klein et al. dataset---------
dat <- readRDS("./data/klein.rds")
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
endo.counts <- counts[rownames(counts) %in% names(symbol2entrez), ]
endo.counts <- endo.counts %>%
  as.data.frame() %>%
  mutate(genename = row.names(.)) %>%
  group_by(genename) %>%
  summarise_all(funs(sum))
endo.counts <- as.data.frame(endo.counts)
rownames(endo.counts) <- endo.counts[,1]
counts.all[[3]] <- as.matrix(endo.counts[,-1])
gset <- getGEO('GSE70580')
annot <- pData(gset[[1]])
group.all[[3]] <- as.factor(annot$`facs gating:ch1`)


num.genes <- 1000
methods <- c(hvg.detection,
             hvg.detection.seurat,
             hvg.detection.scran,
             hvg.detection.brennecke)
titles = c('Li', 'Klein', 'Bjorklund')
ARI = c()

for (i in c(1:3)) {
  counts <- counts.all[[i]]
  group <- group.all[[i]]
  genes.keep <- which(rowMeans(counts) > 1)
  counts <- as.matrix(counts[genes.keep, ])
  normcounts <- norm.libsize(counts)$normcounts
  
  hvgs <- list()
  for (j in 1:4) {
    method <- methods[[j]]
    hvg <- order(method(counts), decreasing = T)[1:num.genes]
    hvgs[[j]] <- hvg
    X <- log10(t(normcounts[hvg, ])+1)
    X <- prcomp(X, scale. = F)$x[, 1:10]
    sil <- c()
    for (k in 2:10){
      clust <- Mclust(X, G = k:k, modelNames = c('VVV'))
      si <- silhouette(as.numeric(clust$classification), dist(X))
      sil <- c(sil, mean(si[,3]))
    }
    k <- which.max(sil)+1
    clust <- Mclust(X, G = k:k, modelNames = c('VVV'))
    ARI <- c(ARI, adjustedRandIndex(group, clust$classification))
  }
  
  venn.diagram(hvgs,
               filename = paste0(titles[[i]], '.tiff'),
               fill = c("orange", "red", "green", "blue"),
               category=c('Proposed', 'Seurat', 'Scran', 'Brennecke'),
               main = paste(titles[[i]], 'Dataset'))
}



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
                     ARI=ARI)


data.m$Datasets <- factor(data.m$Datasets, levels = c('Li et al.', 'Klein et al.', 'Bjorklund et al.'))
data.m$Methods <- factor(data.m$Methods, levels = c('Proposed', 'Seurat', 'Scran', 'Brennecke'))

ggplot(data.m, aes(Datasets, ARI)) +   
  coord_cartesian(ylim=c(0.5,0.9)) + 
  theme(axis.text = element_text(size=20), text = element_text(size = 20))+
  geom_bar(aes(fill = Methods), position = "dodge", stat="identity")
