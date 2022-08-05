set.seed(123)

#load packages
library(ggplot2)
library(data.table)
library(stringr)
library(seqLogo)
library(GenomicRanges)
library(Biostrings)
library(BSgenome.Mmusculus.UCSC.mm10.UCSC.mm10)
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
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

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
library(BSgenome.Mmusculus.UCSC.mm10)
library(BSgenome.Mmusculus.UCSC.mm9)
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
pancreas13<-fread("/s/project/transposon/Output/Snakemake/Input/CRIUK_PPM_PANCREAS_20110408_SB13_INSERTIONS.bed")
uniqueN(pancreas13$V4)
t <- as.data.frame(tstrsplit(pancreas13$V4, "_"))

pancreas13$TumorID <-t$c..6.273....2.118....6.273....3.166....2.121....4.215....6.309...

#subset for sb
SB_insertions<-pancreas13
#tail_samples


#specify genome
genome=BSgenome.Mmusculus.UCSC.mm9

######get sequence
#convert to Granges
insertion_table_to_Granges<- function(x){
  y<-data.table()
  y<-y[,chr:= x$V1]
  y<-y[,start:= as.numeric(as.character(x$V2))]
  y<-y[,end:= as.numeric(as.character(x$V3))]
  y<-y[,orientation:= as.character(x$V6)]
  y<-y[,ReadCoverage:= x$V5]
  y<-y[,TumorID:= x$TumorID]
  y<- toGRanges(y, keep.extra.columns = TRUE, genome=genome)
  strand(y)<-y$orientation
  return(y)
}

SB_Granges<-insertion_table_to_Granges(SB_insertions)
chain<-import.chain(con="/s/project/transposon/Input/mm9ToMm10.over.chain")
SB_Granges<-unlist(liftOver(SB_Granges,chain))
seqlevels(SB_Granges) <- seqlevels(BSgenome.Mmusculus.UCSC.mm10)
seqinfo(SB_Granges) <- seqinfo(BSgenome.Mmusculus.UCSC.mm10)


sequencesallSB<-(DNAStringSet(getSeq(BSgenome.Mmusculus.UCSC.mm10, SB_Granges)))
basematrix_positives<-as.matrix(sequencesallSB)
ggseqlogo(consensusMatrix(basematrix_positives), plot=TRUE)

SB_Granges_correct <- SB_Granges

#remove cancer genes
TCGA_pancancer_all <- fread("/s/project/transposon/Output/Snakemake/Input/oncovar/pan_cancer/TCGA.PanCancer.all.genes.OncoVar.tsv.gz")
TCGA_pancancer_consensusscore <- TCGA_pancancer_all[-which(TCGA_pancancer_all$Consensus_Score==0),]$Gene_symbol
cancer_gene_catalogue <- TCGA_pancancer_consensusscore

genes <-genes(TxDb.Mmusculus.UCSC.mm10.knownGene)
genes$symbol <- (mapIds(org.Mm.eg.db,
                                keys = genes$gene_id,
                                column = "SYMBOL",
                                keytype = "ENTREZID",
                                multiVals = "first"))

cancerGenes <- genes[which(genes$symbol%in%(str_to_title(cancer_gene_catalogue)))]
insertionInCancerGenes <- which(overlapsAny(SB_Granges_correct,cancerGenes)==TRUE)
length(cancerGenes)/length(genes)
length(insertionInCancerGenes)/length(SB_Granges_correct)

SB_Granges_correct<- SB_Granges_correct[-insertionInCancerGenes]
saveRDS(SB_Granges_correct, "/s/project/transposon/Output/Snakemake/PancreasNoCancer/correct_ranges.RData")

