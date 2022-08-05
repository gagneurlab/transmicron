set.seed(123)

#load packages
library(ggplot2)
library(data.table)
library(stringr)
library(seqLogo)
library(GenomicRanges)
library(Biostrings)
library(BSgenome.Mmusculus.UCSC.mm10)
library(BSgenome.Mmusculus.UCSC.mm9)
library(BSgenome)
#BiocManager::install("ggseqlogo")
library(ggseqlogo)
#library(TFBSTools)
require(ggplot2)
library(SRAdb)
library(ShortRead)
library(pROC)
library(reshape2)
library(regioneR)
library(data.table)
#library(SimRAD)
library(caret)
library(dplyr)
library(e1071)
library(readr)
library(PRROC)
#BiocManager::install("motifRG")
#install.packages("Bioconductor")
library(motifRG)

library(precrec)
library(plotROC)
library(stringr)
#library(strip)
library(randomForest)
#BiocManager::install("regioneR")
#BiocManager::install("plyranges")
library(plyranges)
library(tidyverse)
require(tidyquant)
library(knitr)
library(plyr)
#insSB.packages("caret")
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(plyranges)
#insSB.packages("pheatmap")
library(pheatmap)
#insSB.packages("e1071")
library(e1071)
library(readxl)
#insSB.packages("glmnet")
library(glmnet)
#insSB.packages("mvtnorm")
library(mvtnorm) 
#insSB.packages("HDCI")
#library(HDCI)
#insSB.packages("ROCR")
library(ROCR)
#BiocManager::insSB("org.Mm.eg.db")
#library(biglm)

library(org.Mm.eg.db)
#insSB.packages("randomForest")
library(randomForest)
#insSB.packages("bigglm")
#library(biglm)
#insSB.packages("Boruta")
#insSB.packages("randomForestSRC")

library(mlbench)
library(caret)
library(ggplot2)
library(ranger)
library(performance)
library(Boruta)
library(randomForestSRC)
library(vita)
library(MLmetrics)
library(doParallel)

#input sleeping beauty vs piggy bac
PiggyBac_insertions <-fread("/s/project/transposon/Output/Snakemake/Input/20210212TcellPB.BED")


#specify genome
genome=BSgenome.Mmusculus.UCSC.mm10

######get sequence
#convert to Granges
insertion_table_to_Granges<- function(x){
  y<-data.table()
  y<-y[,chr:= x$Chromosome]
  y<-y[,start:= as.numeric(as.character(x$TransposonIntegrationSite))]
  y<-y[,end:= as.numeric(as.character(x$TransposonIntegrationSite))]
  y<-y[,orientation:= x$orientation]
  y<-y[,gene:= x$Gene]
  y<-y[,ReadCoverage:= x$ReadCoverage]
  y<-y[,normCov:= x$normCov]
  y<-y[,TumorID:= x$TumorID]
  y<-y[,organ:=x$organ]
  y<-y[,diagnosis:=x$diagnosis]
  y<-y[,transposon_type:=x$transposon_type]
  y<-regioneR::toGRanges(y, genome=genome)
  #strand(y)<-y$orientation
  return(y)
}

PiggyBac_Granges<-insertion_table_to_Granges(PiggyBac_insertions)
names(PiggyBac_Granges) <- 1:length(PiggyBac_Granges)

#liftover to mm9
  chain<-import.chain(con="/s/project/transposon/Input/mm10ToMm9.over.chain")
  PiggyBac_Granges<- unlist(liftOver(PiggyBac_Granges,chain))


genome=BSgenome.Mmusculus.UCSC.mm9
seqlevels(PiggyBac_Granges) <- seqlevels(genome)
seqinfo(PiggyBac_Granges) <- seqinfo(genome)
#PiggyBac_Granges <- trim(PiggyBac_Granges)
#PiggyBac_Granges <- PiggyBac_Granges[-which(width(PiggyBac_Granges)!=1)]


#check whether insertions are within TTAA (sequencing error)
check_TTAA<-function(x){
  
  start(x) <- start(x) - 1
  end(x) <- end(x) + 2
  return(x)
}


PiggyBac_Granges<-check_TTAA(PiggyBac_Granges)

sequences_PiggyBac<-(DNAStringSet(getSeq(BSgenome.Mmusculus.UCSC.mm9, PiggyBac_Granges)))
basematrix_PiggyBac<-as.matrix(sequences_PiggyBac)
# consensusMatrix(basematrix_PiggyBac)
  #check which insertions are shifted by 1
PiggyBac_shifted_Granges<-PiggyBac_Granges

shift_to_TTAA<-function(x){
  original<-x
  y<-x
  end(y) <- end(y) -1
  to_shift<-which(DNAStringSet(getSeq(BSgenome.Mmusculus.UCSC.mm9, y))=="TAA")
  z<-x[to_shift]
  
  start(z)<-start(z)-1
  end(z)<-end(z)-1
  x[to_shift]<-z
  start(original) <- start(original) +1
  
  to_shift2<-which(DNAStringSet(getSeq(BSgenome.Mmusculus.UCSC.mm9, original))=="TTA")
  z<-x[to_shift2]
  
  start(z)<-start(z)+1
  end(z)<-end(z)+1
  
  x[to_shift2]<-z
  return(x)
}

PiggyBac_Granges_correct<-shift_to_TTAA(PiggyBac_shifted_Granges)
names(PiggyBac_Granges_correct) <- 1:length(PiggyBac_Granges_correct)
strand(PiggyBac_Granges_correct) <- PiggyBac_Granges_correct$orientation


saveRDS(PiggyBac_Granges_correct, "/s/project/transposon/Output/Snakemake/TALL_mm9/correct_ranges.RData")


sequences_PiggyBac<-(DNAStringSet(getSeq(BSgenome.Mmusculus.UCSC.mm9, PiggyBac_Granges_correct)))
basematrix_PiggyBac<-as.matrix(sequences_PiggyBac)
test<- sample(1:nrow(basematrix_PiggyBac), 100000)
ggseqlogo(consensusMatrix(basematrix_PiggyBac[test,]))
