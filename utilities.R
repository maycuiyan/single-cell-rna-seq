library(Seurat)
library(org.Hs.eg.db)
library(SummarizedExperiment)
library(scone)
library(readr)
library(scran)


hvg.detection <- function(counts) {
  xx <- rowMeans(counts)^2
  yy <- apply(counts, 1, var)-rowMeans(counts)
  res <- gpava(-log(xx)/2, log(yy/xx+1), ties = "primary", solver = weighted.fractile, p = 0.5)
  fit.spline<-smooth.spline(log(xx)/2, (res$x), df = 11)
  fitted.spline <- predict(fit.spline, log(xx)/2)
  score <- log(yy/xx+1)-fitted.spline$y
  score
}


hvg.detection.seurat <- function(counts) {
  seurat <- new("seurat", raw.data = Matrix(counts))
  seurat <- Setup(seurat, min.cells = 0, min.genes = 0, do.scale = F, project = "sim", do.center = F)
  seurat <- MeanVarPlot(seurat)
  score <- seurat@mean.var[,'data.norm.y']
  score
}


hvg.detection.scran <- function(counts) {
  sce <- SingleCellExperiment(list(counts=counts))
  sce <- computeSumFactors(sce)
  sce <- normalize(sce)
  var.fit <- trendVar(sce, use.spikes=F)
  var.out <- decomposeVar(sce, var.fit)
  score <- -var.out$p.value
  score
}


hvg.detection.brennecke <- function(counts) {
  lib.size <- norm.libsize(counts)$s
  ed <- t(t(counts)/lib.size)
  means <- rowMeans(ed)
  vars <- apply(ed,1,var)
  cv2 <- vars/means^2
  minMeanForFit <- unname(quantile(means[which(cv2 > 0.3)], .95))
  useForFit <- means >= minMeanForFit
  fit <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/means[useForFit] ),cv2[useForFit])
  a0 <- unname( fit$coefficients["a0"] )
  a1 <- unname( fit$coefficients["a1tilde"])
  afit <- a1/means+a0
  score <- vars/(afit*means^2)
  score
}


detect.lvg <- function(counts) {
  xx <- rowMeans(counts)^2
  yy <- apply(counts, 1, var)-rowMeans(counts)
  res <- gpava(-log(xx)/2, log(yy/xx+1), ties = "primary", solver = weighted.fractile, p = 0.1)
  fit.spline<-smooth.spline(log(xx)/2, (res$x), df = 11)
  fitted.spline <- predict(fit.spline, log(xx)/2)
  lvg <- which(log(yy/xx+1) <= fitted.spline$y)
  lvg
}


norm.scran <- function(counts) {
  clust<-quickCluster(counts, min.size=20, method='igraph')
  s <- computeSumFactors(counts, clusters=clust)
  s <- s/mean(s)
  normcounts <- t(t(counts)/s)
  list(s=s, normcounts=normcounts)
}


norm.libsize <- function(counts) {
  s <- colMeans(counts)
  s <- s/mean(s)
  normcounts <- t(t(counts)/s)
  list(s=s, normcounts=normcounts)
}


cell.metrics <- function(endo.counts, mito.counts = NULL, ercc.counts = NULL) {
  if (is.null(mito.counts)) {
    mito.counts <- matrix(0, nrow=1, ncol=ncol(endo.counts))
  }
  if (is.null(ercc.counts)) {
    ercc.counts <- matrix(0, nrow=1, ncol=ncol(endo.counts))
  }
  libsize <- colSums(endo.counts) + colSums(mito.counts) + colSums(ercc.counts)
  totfeat <- colSums(endo.counts>0) + colSums(mito.counts>0) + colSums(ercc.counts>0)
  spikeprop <- colSums(ercc.counts) / libsize
  mitoprop <- colSums(mito.counts) / libsize
  list(libsize=libsize, totfeat=totfeat, spikeprop=spikeprop, mitoprop=mitoprop)
}


gene.metrics <- function(endo.counts) {
  mean <- rowMeans(endo.counts)
  exprprop <- rowMeans(endo.counts>0)
  list(mean=mean, exprprop=exprprop)
}


read.GSE70580 <- function() {
  filenames <- dir('./data/GSE70580_RAW/')
  counts <- list()
  for (i in 1:length(filenames)) {
    count <- read_delim(paste0('./data/GSE70580_RAW/', filenames[i]), "\t", escape_double = FALSE, trim_ws = TRUE)
    counts[[i]] <- count$Reads
  }
  counts <- do.call(cbind, counts)
  rownames(counts) <- count$`#Gene symbol`
  colnames(counts) <- gsub('_.*', '', filenames)
  ercc.counts <- counts[tail(rownames(counts), 92), ]
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
  endo.counts <- as.matrix(endo.counts[,-1])
  endo.counts <- endo.counts[rowMeans(endo.counts)>1, ]
  cm <- cell.metrics(endo.counts = endo.counts, ercc.counts = ercc.counts)
  sample.keep <- cm$libsize < 2e06 & cm$totfeat > 1100 & cm$spikeprop < 0.33 & cm$spikeprop > 0.12
  endo.counts <- endo.counts[, sample.keep]
  ercc.counts <- ercc.counts[, sample.keep]
  gm <- gene.metrics(endo.counts)
  gene.keep <- gm$mean > 1 & gm$exprprop > 0.1
  endo.counts <- endo.counts[gene.keep, ]
  ercc.counts <- ercc.counts[rowMeans(ercc.counts)>1,]
  list(endo.counts=endo.counts, ercc.counts=ercc.counts)
}


read.GSE86977 <- function() {
  counts <- read_csv("./data/GSE86977_UMI.2684.csv", col_names = F, skip = 1)
  counts <- as.data.frame(counts)
  gene.names <- counts[,1]
  index.keep <- which(!duplicated(gene.names))
  counts <- counts[index.keep,c(-1,-2)]
  rownames(counts) <- gene.names[index.keep]
  colnames(counts) <- read_csv("./data/GSE86977_UMI.2684.csv", col_names = F, n_max = 1)[-1]
  endo.counts <- as.matrix(counts[1:23709,])
  ercc.counts <- as.matrix(counts[23710:23801,])
  cm <- cell.metrics(endo.counts = endo.counts, ercc.counts = ercc.counts)
  sample.keep <- cm$libsize < 2e05 & cm$totfeat > 5000 & cm$spikeprop < 0.9 & cm$spikeprop > 0.01
  endo.counts <- endo.counts[, sample.keep]
  ercc.counts <- ercc.counts[, sample.keep]
  gm <- gene.metrics(endo.counts)
  gene.keep <- gm$mean > 1 & gm$exprprop > 0.1
  endo.counts <- endo.counts[gene.keep, ]
  ercc.counts <- ercc.counts[rowMeans(ercc.counts)>1,] 
  list(endo.counts=endo.counts, ercc.counts=ercc.counts)
}


read.GSE46980 <- function() {
  counts <- read_delim("./data/GSE46980_CombinedMoleculeCounts.tab", 
                           "\t", escape_double = FALSE, trim_ws = TRUE, skip=7, col_names = F)
  counts <- as.data.frame(counts)
  rownames(counts) <- counts$X1
  counts <- counts[,-c(1:7)]
  sample.names <- read_delim("./data/GSE46980_CombinedMoleculeCounts.tab", "\t", 
                             escape_double = FALSE, 
                             trim_ws = TRUE, 
                             skip=5, col_names = F, n_max = 1)
  colnames(counts) <- sample.names[,8:103]
  ercc.counts <- as.matrix(counts[1:96,])
  endo.counts <- as.matrix(counts[97:25914,])
  cm <- cell.metrics(endo.counts = endo.counts, ercc.counts = ercc.counts)
  sample.keep <- cm$libsize > 2e04 & cm$totfeat > 400 & cm$spikeprop < 0.9 & cm$spikeprop > 0.01
  endo.counts <- endo.counts[, sample.keep]
  ercc.counts <- ercc.counts[, sample.keep]
  gm <- gene.metrics(endo.counts)
  gene.keep <- gm$mean > 1 & gm$exprprop > 0.1
  endo.counts <- endo.counts[gene.keep, ]
  ercc.counts <- ercc.counts[rowMeans(ercc.counts)>1,] 
  list(endo.counts=endo.counts, ercc.counts=ercc.counts)
}


read.EMTAB5522 <- function() {
  counts <- read_delim("./data/counts_Calero_20160113.tsv", 
                           "\t", escape_double = FALSE, trim_ws = TRUE)
  counts <- as.data.frame(counts)
  rownames(counts) <- counts$GeneID
  counts <- counts[,-c(1,2)]
  ercc.counts <- as.matrix(counts[46604:46702,])
  endo.counts <- as.matrix(counts[1:46603,])
  cm <- cell.metrics(endo.counts = endo.counts, ercc.counts = ercc.counts)
  sample.keep <- cm$libsize > 5e05 & cm$totfeat > 5000 & cm$spikeprop < 0.9 & cm$spikeprop > 0.01
  endo.counts <- endo.counts[, sample.keep]
  ercc.counts <- ercc.counts[, sample.keep]
  gm <- gene.metrics(endo.counts)
  gene.keep <- gm$mean > 1 & gm$exprprop > 0.1
  endo.counts <- endo.counts[gene.keep, ]
  ercc.counts <- ercc.counts[rowMeans(ercc.counts)>1,] 
  endo.counts <- endo.counts[rowMeans(endo.counts)>1, ]
  list(endo.counts=endo.counts, ercc.counts=ercc.counts)
}


read.EMTAB3929 <- function() {
  endo.counts <- read_delim("./data/counts.txt",
                       "\t", escape_double = FALSE, trim_ws = TRUE)
  endo.counts <- as.data.frame(endo.counts)
  rownames(endo.counts) <- endo.counts$X1
  endo.counts <- as.matrix(endo.counts[,-1])
  ercc.counts <- read_delim("./data/ercc.counts.txt",
                     "\t", escape_double = FALSE, trim_ws = TRUE)
  ercc.counts <- as.data.frame(ercc.counts)
  rownames(ercc.counts) <- ercc.counts$X1
  ercc.counts <- as.matrix(ercc.counts[,-1])
  cm <- cell.metrics(endo.counts = endo.counts, ercc.counts = ercc.counts)
  sample.keep <- cm$libsize > 5e04 & cm$totfeat > 400 & cm$spikeprop < 0.9 & cm$spikeprop > 0.02
  endo.counts <- endo.counts[, sample.keep]
  ercc.counts <- ercc.counts[, sample.keep]
  gm <- gene.metrics(endo.counts)
  gene.keep <- gm$mean > 1 & gm$exprprop > 0.1
  endo.counts <- endo.counts[gene.keep, ]
  ercc.counts <- ercc.counts[rowMeans(ercc.counts)>1,]
  list(endo.counts=endo.counts, ercc.counts=ercc.counts)
}


read.EMTAB2805 <- function() {
  read.data <- function(filepath) {
    counts <- read_delim(filepath,
                             "\t", escape_double = FALSE, trim_ws = TRUE)
    counts <- counts[1:38385,]
    counts <- as.data.frame(counts)
    rownames(counts) <- counts$EnsemblGeneID
    counts <- counts[,-c(1,2,3,4)]
    endo.counts <- as.matrix(counts[1:38293,])
    ercc.counts <- as.matrix(counts[38294:38385,])
    list(endo.counts=endo.counts, ercc.counts=ercc.counts)
  }
  dat.G1 <- read.data("./data/G1_singlecells_counts.txt")
  dat.G2 <- read.data("./data/G2M_singlecells_counts.txt")
  dat.S <- read.data("./data/S_singlecells_counts.txt")
  endo.counts <- do.call(cbind,list(dat.G1$endo.counts, dat.G2$endo.counts, dat.S$endo.counts))
  ercc.counts <- do.call(cbind,list(dat.G1$ercc.counts, dat.G2$ercc.counts, dat.S$ercc.counts))
  cm <- cell.metrics(endo.counts = endo.counts, ercc.counts = ercc.counts)
  sample.keep <- cm$libsize > 1e06 & cm$totfeat > 2000 & cm$spikeprop < 0.8
  endo.counts <- endo.counts[, sample.keep]
  ercc.counts <- ercc.counts[, sample.keep]
  gm <- gene.metrics(endo.counts)
  gene.keep <- gm$mean > 1 & gm$exprprop > 0.1
  endo.counts <- endo.counts[gene.keep, ]
  ercc.counts <- ercc.counts[rowMeans(ercc.counts)>1,] 
  list(endo.counts=endo.counts, ercc.counts=ercc.counts)
}


read.GSE95601 <- function() {
  load("./data/GSE95601_oeHBCdiff_Cufflinks_eSet.Rda")
  E <- assayData(Cufflinks_eSet)$counts_table
  E <- na.omit(E)
  E <- E[rowSums(E)>0,]
  cre <- E["CreER",]
  ercc.counts <- E[grep("^ERCC-", rownames(E)),]
  E <- E[grep("^ERCC-", rownames(E), invert = TRUE), ]
  E <- E[-which(rownames(E)=="CreER"), ]
  qc <- as.matrix(protocolData(Cufflinks_eSet)@data)[,c(1:5, 10:18)]
  qc <- cbind(qc, CreER = cre, ERCC_reads = colSums(ercc.counts))
  batch <- droplevels(pData(Cufflinks_eSet)$MD_c1_run_id)
  bio <- droplevels(pData(Cufflinks_eSet)$MD_expt_condition)
  clusterLabels <- read.table("./data/oeHBCdiff_clusterLabels.txt",
                              sep = "\t", stringsAsFactors = FALSE)
  m <- match(colnames(E), clusterLabels[, 1])
  metadata <- data.frame("Experiment" = bio,
                         "Batch" = batch,
                         "publishedClusters" = clusterLabels[m,2],
                         qc)
  metadata$publishedClusters[is.na(metadata$publishedClusters)] <- -2
  se <- SummarizedExperiment(assays = list(counts = E),
                             colData = metadata)
  data("housekeeping")
  hk = rownames(se)[toupper(rownames(se)) %in% housekeeping$V1]
  mfilt <- metric_sample_filter(assay(se),
                                nreads = colData(se)$NREADS,
                                ralign = colData(se)$RALIGN,
                                pos_controls = rownames(se) %in% hk,
                                zcut = 3, mixture = FALSE,
                                plot = F)
  mfilt <- !apply(simplify2array(mfilt[!is.na(mfilt)]), 1, any)
  se <- se[, mfilt]
  endo.counts <- assay(se)
  ercc.counts <- ercc.counts[,colnames(endo.counts)]
  batch <- colData(se)$Batch
  endo.counts <- endo.counts[, order(batch)]
  ercc.counts <- ercc.counts[, order(batch)]
  batch <- batch[order(batch)]
  thresh <- 10
  endo.counts <- endo.counts[, colMeans(ercc.counts)>thresh]
  batch <- batch[colMeans(ercc.counts)>thresh]
  ercc.counts <- ercc.counts[, colMeans(ercc.counts)>thresh]
  batch <- droplevels(batch)
  bats <- c(3,4,5,6,8)
  endo.counts <- endo.counts[, as.numeric(batch) %in% bats]
  ercc.counts <- ercc.counts[, as.numeric(batch) %in% bats]
  cm <- cell.metrics(endo.counts = endo.counts, ercc.counts = ercc.counts)
  sample.keep <- cm$libsize > 5000 & cm$totfeat > 1000 & cm$spikeprop < 0.8 & cm$spikeprop > 0.05
  endo.counts <- endo.counts[, sample.keep]
  ercc.counts <- ercc.counts[, sample.keep]
  gm <- gene.metrics(endo.counts)
  gene.keep <- gm$mean > 1 & gm$exprprop > 0.1
  endo.counts <- endo.counts[gene.keep, ]
  ercc.counts <- ercc.counts[rowMeans(ercc.counts)>1,] 
  list(endo.counts=endo.counts, ercc.counts=ercc.counts)
}
