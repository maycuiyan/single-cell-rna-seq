library(readr)

##---------------------directory to store the downloaded datasets---------------------
dir.create('./data')

#---------------------------------------GSE70580--------------------------------------
if (!'GSE70580_RAW.tar' %in% list.files('./data')) {
  download.file(url='https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE70580&format=file', 
                destfile='./data/GSE70580_RAW.tar',
                mode='wb')
  untar(tarfile='./data/GSE70580_RAW.tar', exdir='./data/GSE70580_RAW')
}

#---------------------------------------GSE95601--------------------------------------
if (!'GSE95601_oeHBCdiff_Cufflinks_eSet.Rda' %in% list.files('./data')) {
  download.file(url='https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE95601&format=file&file=GSE95601%5FoeHBCdiff%5FCufflinks%5FeSet%2ERda%2Egz', 
                destfile='./data/GSE95601_oeHBCdiff_Cufflinks_eSet.Rda.gz',
                mode='wb')
  download.file(url='https://raw.githubusercontent.com/rufletch/p63-HBC-diff/master/ref/oeHBCdiff_clusterLabels.txt',
                destfile='./data/oeHBCdiff_clusterLabels.txt',
                mode='wb')
  gunzip(filename='./data/GSE95601_oeHBCdiff_Cufflinks_eSet.Rda.gz', 
         destname='./data/GSE95601_oeHBCdiff_Cufflinks_eSet.Rda')
}

#---------------------------------------E-MTAB-5522--------------------------------------
if (!'E-MTAB-5522.processed.1.zip' %in% list.files('./data')) {
  download.file(url='https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-5522/E-MTAB-5522.processed.1.zip', 
                destfile='./data/E-MTAB-5522.processed.1.zip',
                mode='wb') 
  unzip(zipfile='./data/E-MTAB-5522.processed.1.zip', exdir='./data')
}

#---------------------------------------GSE86977--------------------------------------
if (!'GSE86977_UMI.2684.csv' %in% list.files('./data')) {
  download.file(url='https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE86977&format=file&file=GSE86977%5FUMI%2E2684%2Ecsv%2Egz', 
                destfile='./data/GSE86977_UMI.2684.csv.gz',
                mode='wb')
  gunzip(filename='./data/GSE86977_UMI.2684.csv.gz', 
         destname='./data/GSE86977_UMI.2684.csv')
}

#---------------------------------------GSE46980--------------------------------------
if (!'GSE46980_CombinedMoleculeCounts.tab.gz' %in% list.files('./data')) {
  download.file(url='https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE46980&format=file&file=GSE46980%5FCombinedMoleculeCounts%2Etab%2Egz', 
                destfile='./data/GSE46980_CombinedMoleculeCounts.tab.gz',
                mode='wb') 
  gunzip(filename='./data/GSE46980_CombinedMoleculeCounts.tab.gz', 
         destname='./data/GSE46980_CombinedMoleculeCounts.tab')
}

#---------------------------------------E-MTAB-3929--------------------------------------
if (!'E-MTAB-3929.processed.1.zip' %in% list.files('./data')) {
  download.file(url='https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-3929/E-MTAB-3929.processed.1.zip', 
                destfile='./data/E-MTAB-3929.processed.1.zip',
                mode='wb') 
  download.file(url='https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-3929/E-MTAB-3929.processed.3.zip', 
                destfile='./data/E-MTAB-3929.processed.3.zip',
                mode='wb') 
  unzip(zipfile='./data/E-MTAB-3929.processed.1.zip', exdir='./data')
  unzip(zipfile='./data/E-MTAB-3929.processed.3.zip', exdir='./data')
}

#---------------------------------------E-MTAB-2805--------------------------------------
if (!'E-MTAB-2805.processed.1.zip' %in% list.files('./data')) {
  download.file(url='https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-2805/E-MTAB-2805.processed.1.zip', 
                destfile='./data/E-MTAB-2805.processed.1.zip',
                mode='wb') 
  unzip(zipfile='./data/E-MTAB-2805.processed.1.zip', exdir='./data')
}

#---------------------------------------Li et al dataset--------------------------------------
if (!'li.rds' %in% list.files('./data')) {
  download.file(url='https://scrnaseq-public-datasets.s3.amazonaws.com/scater-objects/li.rds', 
                destfile='./data/li.rds',
                mode='wb')
}

#---------------------------------------Klein et al dataset--------------------------------------
if (!'klein.rds' %in% list.files('./data')) {
  download.file(url='https://scrnaseq-public-datasets.s3.amazonaws.com/scater-objects/klein.rds', 
                destfile='./data/klein.rds',
                mode='wb')
}

#---------------------------------------GSE75790--------------------------------------
if (!'GSE75790_ziegenhain_complete_data.txt' %in% list.files('./data')) {
  download.file(url='https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE75790&format=file&file=GSE75790%5Fziegenhain%5Fcomplete%5Fdata%2Etxt%2Egz', 
                destfile='./data/GSE75790_ziegenhain_complete_data.txt.gz',
                mode='wb') 
  gunzip(filename='./data/GSE75790_ziegenhain_complete_data.txt.gz', 
         destname='./data/GSE75790_ziegenhain_complete_data.txt')
}
