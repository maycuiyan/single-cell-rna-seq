library(isotone)
library(readr)
library(SingleCellExperiment)
library(Matrix)
library(ggplot2)
require(DESeq)
require(statmod)
library(Seurat)
library(scran)
library(ggplot2)
library(gridExtra)
library(grid)
library(GGally)
source('utilities.R')

all.counts <- read_delim("./data/GSE75790_ziegenhain_complete_data.txt",
                         "\t", escape_double = FALSE, trim_ws = TRUE, skip=1, col_names = F)
sample.names <- read_delim("./data/GSE75790_ziegenhain_complete_data.txt",
                           "\t", escape_double = FALSE, trim_ws = TRUE, col_names = F, n_max = 1)
all.counts <- as.data.frame(all.counts)
rownames(all.counts) <- all.counts[,1]
all.counts <- all.counts[,-1]
colnames(all.counts) <- sample.names
counts.A <- list(as.matrix(all.counts[, grep('CELseq2A', colnames(all.counts))]),
                 as.matrix(all.counts[, grep('DropSeqA', colnames(all.counts))]),
                 as.matrix(all.counts[, grep('MARSseqA', colnames(all.counts))]),
                 as.matrix(all.counts[, grep('SCRBseqA', colnames(all.counts))]),
                 as.matrix(all.counts[, grep('SmartSeqA', colnames(all.counts))]),
                 as.matrix(all.counts[, grep('SmartSeq2A', colnames(all.counts))]))
counts.B <- list(as.matrix(all.counts[, grep('CELseq2B', colnames(all.counts))]),
                 as.matrix(all.counts[, grep('DropSeqB', colnames(all.counts))]),
                 as.matrix(all.counts[, grep('MARSseqB', colnames(all.counts))]),
                 as.matrix(all.counts[, grep('SCRBseqB', colnames(all.counts))]),
                 as.matrix(all.counts[, grep('SmartSeqB', colnames(all.counts))]),
                 as.matrix(all.counts[, grep('SmartSeq2B', colnames(all.counts))]))

num.genes <- 1000
A <- c()
cnt <- 1
for (i in 1:6) {
  for (j in i:6) {
    if (i==j) {
      counts1 <- counts.A[[i]]
      counts2 <- counts.B[[j]]
    } else {
      counts1 <- counts.A[[i]]
      counts2 <- counts.A[[j]]
    }
    genes.keep <- which(rowMeans(counts1) >= 1 & rowMeans(counts2) >= 1)
    counts1 <- counts1[genes.keep, ]
    counts2 <- counts2[genes.keep, ]
    # method 1: proposed
    hvg1 <- order(hvg.detection(counts1), decreasing = T)[1:num.genes]
    hvg2 <- order(hvg.detection(counts2), decreasing = T)[1:num.genes]
    A[[cnt]] <- length(intersect(hvg1, hvg2))
    cnt <- cnt + 1
    # method 2: seurat
    hvg1.seurat <- order(hvg.detection.seurat(counts1), decreasing = T)[1:num.genes]
    hvg2.seurat <- order(hvg.detection.seurat(counts2), decreasing = T)[1:num.genes]
    A[[cnt]] <- length(intersect(hvg1.seurat, hvg2.seurat))
    cnt <- cnt + 1
    # method 3: scran
    hvg1.scran <- order(hvg.detection.scran(counts1))[1:num.genes]
    hvg2.scran <- order(hvg.detection.scran(counts2))[1:num.genes]
    A[[cnt]] <- length(intersect(hvg1.scran, hvg2.scran))
    cnt <- cnt + 1
    # method 4: brenneke
    hvg1.brennecke <- order(hvg.detection.brennecke(counts1), decreasing = T)[1:num.genes]
    hvg2.brennecke <- order(hvg.detection.brennecke(counts2), decreasing = T)[1:num.genes]
    A[[cnt]] <- length(intersect(hvg1.brennecke, hvg2.brennecke))
    cnt <- cnt + 1
  }
}

A <- matrix(A, nrow=4)

plotList <- list()
cnt <- 1
for (i in 1:6) {
  for (j in 1:6) {
    if (j >= i) {
      df <- data.frame(method=c("Proposed", "Seurat", "Scran", "Brennecke"),
                       overlap=A[,cnt])
      plotList[[(i-1)*6+j]] <- ggplot(df, aes(x=method, y=overlap, fill=method, label=overlap))+
        geom_bar(position = 'dodge', stat="identity", color="black") +
        scale_x_discrete(limits=c("Proposed", "Seurat", "Scran", "Brennecke")) +
        scale_fill_discrete(breaks=c("Proposed", "Seurat", "Scran", "Brennecke"))+
        coord_cartesian(ylim=c(min(A[,cnt])-50, max(A[,cnt])+50)) +
        ylab(NULL)+
        xlab(NULL)+
        theme_gray()+
        theme(legend.position="none",
              axis.text.x = element_blank(),
              axis.line=element_blank(),
              axis.text.y=element_blank(),
              axis.title.x=element_blank(),
              axis.title.y=element_blank(),
              axis.ticks=element_blank())+
        geom_text(position=position_dodge(width=0.9), vjust=-0.25)
      cnt = cnt+1
    }
    else {
      plotList[[(i-1)*6+j]] <- ggplot()
    }
  }
}

pm <- ggmatrix(
  plotList,
  6, 6,
  c("CEL-seq2", "Drop-seq", "MARS-seqs", 'SCRB-seq', 'Smart-seq', 'Smart-seq2'),
  c("CEL-seq2", "Drop-seq", "MARS-seqs", 'SCRB-seq', 'Smart-seq', 'Smart-seq2'),
  byrow = T,
  showXAxisPlotLabels = T
)

pm
