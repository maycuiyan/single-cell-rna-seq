library(DESeq2)
library(copula)
library(isotone)
library(pROC)
library(plotrix)
library(Seurat)
library(Matrix)
library(scran)
require(statmod)
source('utilities.R')

#---------------------simulate scRNAseq data---------------------
set.seed(123456) # for reproducibility
n1 <- 1000
n2 <- 1000
n3 <- 1000
n <- n1+n2+n3 # number of samples per group
d <- 5000 # number of genes

## shape parameter and mean of each biological gene expr
alpha <- 10*runif(d)
mu <- exp(rnorm(d, 1, 2))

## outliers
outlier_rate = 0.02
ind <- sample(1:d, size=d*outlier_rate, replace = F)
mu[ind] <- median(mu) * 5

## parameters for technical noise
a <- 1
b <- 1
c <- 1
alpha_s <- a*mu^b+c # rate parameter of technical noise
mu_s <- 1
beta_s <- alpha_s/mu_s

## Three groups with different DE proportion
mu1 <- mu * sample(c(0.5, 1, 2), d, replace=T, prob= c(0.6, 0.2, 0.2))
mu2 <- mu * sample(c(0.5, 1, 2), d, replace=T, prob= c(0.2, 0.6, 0.2))
mu3 <- mu * sample(c(0.5, 1, 2), d, replace=T, prob= c(0.2, 0.2, 0.6))

## from generative model
X1 <- matrix(rgamma(d*n1, shape=alpha, rate=alpha/mu1), nrow=d)
X2 <- matrix(rgamma(d*n2, shape=alpha, rate=alpha/mu2), nrow=d)
X3 <- matrix(rgamma(d*n3, shape=alpha, rate=alpha/mu3), nrow=d)
X <- cbind(X1, X2, X3)

#-------S[j,] and S[k,] are fully correlated-----
C <- matrix(rep(runif(d), each = n), nrow=d)
S <- qgamma(C, shape=alpha_s, rate=beta_s)

## Read counts
Y <- X * S
counts <- matrix(rpois(d*n, Y), nrow=d)
rownames(counts) <- paste('Gene', 1:nrow(counts), sep='_')
colnames(counts) <- paste('Cell', 1:ncol(counts), sep='_')

## Drop-out
Drop <- matrix(rbinom(d*n, 1, 0.8), nrow = d)
counts <- counts * Drop


## keep expressed genes
genes.keep <- which(rowMeans(counts)>0)
counts <- counts[genes.keep, ]
X <- X[genes.keep,]
Y <- Y[genes.keep,]
S <- S[genes.keep,]
alpha <- alpha[genes.keep]
mu <- mu[genes.keep]
alpha_s <- alpha_s[genes.keep]
mu_s <- mu_s[genes.keep]
beta_s <- beta_s[genes.keep]

## True CV^2
cv_j <- (rowSds(X) / rowMeans(X))^2
cv_s <- 1/alpha_s

## proposed method
score <- hvg.detection(counts)

## seurat
score.seurat <- hvg.detection.seurat(counts)

# scran
score.scran <- hvg.detection.scran(counts)

## brennecke
score.brenneke <- hvg.detection.brennecke(counts)

## plot ROC cuves
plot(roc(cv_j > 0.5, score), col='red', lty=1, lwd=3, xlim=c(1, 0.0), ylim=c(0.0, 1), identity=F)
plot(roc(cv_j > 0.5, score.seurat), add=T, col='green', lty=2, lwd=3)
plot(roc(cv_j > 0.5, score.scran), add=T, col='blue', lty=3, lwd=3)
plot(roc(cv_j > 0.5, score.brenneke), add=T, col='purple', lty=4, lwd=3)

legend("bottomright", 
       legend = c(paste('Proposed (AUC = ', round(auc(roc(cv_j > 0.5, score)), 3), ')'), 
                  paste('Seurat (AUC = ', round(auc(roc(cv_j > 0.5, score.seurat)), 3), ')'),
                  paste('Scran (AUC = ', round(auc(roc(cv_j > 0.5, score.scran)), 3), ')'),
                  paste('Brennecke (AUC = ', round(auc(roc(cv_j > 0.5, score.brenneke)), 3), ')')), 
       col = c("red", "green", "blue", "purple"),
       lty = c(1,2,3,4),
       lwd = 2)

