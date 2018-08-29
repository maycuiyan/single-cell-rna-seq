# single-cell-rna-seq

This repo contains the R source codes for reproducing the results presented in the following paper submitted for review:

*Yi Cui, Bailiang Li, and Ruijiang Li. "Dissecting biological and technical variability in single-cell RNA-seq data".*

The R source codes here contains everything that is required to reproduce the results presented in our paper. To use it, follow the steps below.

+ Run download_datasets.R. This will download all the necessary scRNAseq data used in our paper, (create if not exists) and populate /data folder under the current directory.

+ Run scatterplots.R. This will reproduce Figure 1 and Supplementary Figure 1 in our paper.

+ Run simulation.R. This will reproduce Figure 3 and Supplmntary Figure 2 (by changing the number of cells to simulate) in our paper.

+ Run clustering.R. This will reproduce Figure 4, Figure 5, and Supplementary Figure 3 in our paper.

+ Run reproducibility.R. This will reproduce Figure 6 in our paper.

+ Run normalization.R. This will reproduce Figuree 7 and Supplementary Figure 4 in our paper. 
