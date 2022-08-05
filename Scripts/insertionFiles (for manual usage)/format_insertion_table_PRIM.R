setwd("/s/project/transposon/Output/Snakemake/PRIM_dataset")
set.seed(123)
print("worked")
#install packages

#load packages
library(ggplot2)
library(data.table)
library(stringr)
library(seqLogo)
library(GenomicRanges)
library(Biostrings)
library(BSgenome)
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
library(motifRG)
library(precrec)
library(plotROC)
library(stringr)
#library(strip)
library(randomForest)
#BiocManager::install("regioneR")
library(plyranges)
#install.packages("tidyverse")
library(tidyverse)
require(tidyquant)
library(knitr)
library(plyr)
#install.packages("caret")
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(plyranges)
#install.packages("pheatmap")
library(pheatmap)
#install.packages("e1071")
library(e1071)
library(readxl)
#install.packages("glmnet")
library(glmnet)
#install.packages("mvtnorm")
library(mvtnorm) 
#install.packages("HDCI")
#library(HDCI)
#install.packages("ROCR")
library(ROCR)
#BiocManager::install("org.Mm.eg.db")
#library(biglm)

library(org.Mm.eg.db)
#install.packages("randomForest")
library(randomForest)
#install.packages("bigglm")
#library(biglm)
#install.packages("Boruta")
#install.packages("randomForestSRC")

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
library(BSgenome.Mmusculus.UCSC.mm9)
#library(BSgenome.Mmusculus.UCSC.mm10)
#SleepingBeauty_insertions<-fread("/s/project/transposon/Output/Snakemake/Input/TS_colon_combined_annot.txt")
sleep <- read.table("/s/project/transposon/Output/Snakemake/Input/TS_colon_combined_annot.txt", 
                    header=T,sep='\t')

#Pulling out the CIS insertions that were in the tail data
#This one doesn't have GI insertions
tails <- ((sleep$MouseChrom==2 & sleep$StartInsertionbp > 98502000 & sleep$StartInsertionbp < 98557231) | (sleep$MouseChrom==4 & sleep$StartInsertionbp > 134752971 & sleep$StartInsertionbp < 134809202) | (sleep$MouseChrom==7 & sleep$StartInsertionbp > 37314639 & sleep$StartInsertionbp < 37330249) | (sleep$MouseChrom==16 & sleep$StartInsertionbp > 15656308 & sleep$StartInsertionbp < 15728591) | (sleep$MouseChrom==16 & sleep$StartInsertionbp > 79341661 & sleep$StartInsertionbp < 79474573))
sleep <- sleep[!tails,]

iss <- paste(sleep$MouseChrom,sleep$StartInsertionbp,sep="/")
poo <- table(sleep$MouseID,iss)
tempiss <- dimnames(poo)$iss
poo <- (table(sleep$MouseID,iss) > 0)
fakeReps <- apply(poo,2,sum)
uniqueIS <- tempiss[fakeReps==1]

sleepTemp <- sleep[is.element(iss,uniqueIS),]
iss <- paste(sleepTemp$MouseChrom,sleepTemp$StartInsertionbp,
             sep="/")
sleepInsert <- sleepTemp[!duplicated(iss),]
rm(fakeReps,iss,poo,sleep,sleepTemp,tempiss,uniqueIS)
#BiocManager::install("BSgenome.Mmusculus.UCSC.mm10", force=TRUE) 

######get sequence
#convert to Granges
Mmusculus=BSgenome.Mmusculus.UCSC.mm9

insertion_table_to_Granges<- function(x){
  y<-data.table()
  y<-y[,chr:= x$MouseChrom]
  y<-y[,start:= as.numeric(as.character(x$StartInsertionbp))]
  y<-y[,end:= as.numeric(as.character(x$StartInsertionbp))]
  y<-y[,orientation:= as.character(x$Strand)]
  y<-y[,gene:= x$GeneName]
  y<-y[,ReadCoverage:= x$freq]
  y<-y[,normCov:= x$normCov]
  y<-y[,TumorID:= x$Library]
  y<-y[,organ:=x$Tissue]
  y<-y[,diagnosis:=x$diagnosis]
  y<-y[,transposon_type:=x$transposon_type]
  y<-regioneR::toGRanges(y, genome=Mmusculus)
  strand(y)<-y$orientation
  names(y)<-1:length(y)
  return(y)
}

SB_Granges <- insertion_table_to_Granges(sleepInsert)
chain<-import.chain(con="/s/project/transposon/Input/mm9ToMm10.over.chain")
SB_Granges <- unlist(liftOver(SB_Granges,chain))
seqlevels(SB_Granges) <- seqlevels(BSgenome.Mmusculus.UCSC.mm10)
seqinfo(SB_Granges) <- seqinfo(BSgenome.Mmusculus.UCSC.mm10)


retrieve_sequence<-function(x){
  
  index_minus<-which(as.vector(strand(x))=="-")
  index_plus<-which(as.vector(strand(x))=="+")
  
  end(x[index_plus]) <- end(x[index_plus]) + lengtharound
  
  
  start(x[index_minus]) <- start(x[index_minus]) -lengtharound
  return(x)
}
# 


SB_Granges<- retrieve_sequence(SB_Granges)
sequencesallSB<-(DNAStringSet(getSeq(BSgenome.Mmusculus.UCSC.mm10, SB_Granges)))

basematrix_positives<-as.matrix(sequencesallSB)
ggseqlogo(consensusMatrix(basematrix_positives), plot=TRUE)

saveRDS(SB_Granges, "/s/project/transposon/Output/Snakemake/PRIM_dataset/correct_ranges.RData")

