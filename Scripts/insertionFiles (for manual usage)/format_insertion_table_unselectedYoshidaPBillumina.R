set.seed(123)

#load packages
library(data.table)
library(GenomicRanges)
library(Biostrings)
library(BSgenome.Mmusculus.UCSC.mm10)
library(readxl)
#input sleeping beauty vs piggy bac
insertions <-read.table("/s/project/transposon/Output/Snakemake/Input/P-17b.rmdup.IS.cov1.unique.bed")

colnames(insertions) <- c("chr", "start", "end", "ID", "ReadCoverage", "orientation")
insertions$TumorID <-tstrsplit(insertions$ID, "_")[[2]]
#subset for unselected samples
#specify genome
genome=BSgenome.Mmusculus.UCSC.mm10

######get sequence
#convert to Granges
insertion_table_to_Granges<- function(x){
  y<-data.table()
  y<-y[,chr:= x$chr]
  y<-y[,start:= as.numeric(as.character(x$start))]
  y<-y[,end:= as.numeric(as.character(x$end))]
  y<-y[,orientation:= x$orientation]
  y<-y[,gene:= x$Gene]
  y<-y[,ReadCoverage:= x$ReadCoverage]
  y<-y[,normCov:= x$normCov]
  y<-y[,TumorID:= x$TumorID]
  y<-y[,organ:=x$organ]
  y<-y[,time:=x$time]
  y<-y[,transposon_type:=x$transposon_type]
  y<-regioneR::toGRanges(y)
  strand(y)<-y$orientation
  return(y)
}
unselected_Granges<-insertion_table_to_Granges(insertions)
names(unselected_Granges) <- 1:length(unselected_Granges)

chain<-import.chain(con="//s/project/transposon/Output/Snakemake/Input/mm8ToMm10.over.chain")
unselected_Granges<-unlist(liftOver(unselected_Granges,chain))
seqlevels(unselected_Granges) <- seqlevels(BSgenome.Mmusculus.UCSC.mm10)
seqinfo(unselected_Granges) <- seqinfo(BSgenome.Mmusculus.UCSC.mm10)
#end(unselected_Granges) <- end(unselected_Granges)-1


shiftRanges <- function(inputRanges){
  minusStrand <-which(inputRanges$orientation=="-")
  end(inputRanges)[minusStrand] <- end(inputRanges)[minusStrand]+3
  start(inputRanges)[minusStrand] <- start(inputRanges)[minusStrand]+3
  
  minusStrand <-which(inputRanges$orientation=="+")
  end(inputRanges)[minusStrand] <- end(inputRanges)[minusStrand]-2
  start(inputRanges)[minusStrand] <- start(inputRanges)[minusStrand]-2
  
  
  return(inputRanges)
}

#check whether insertions are within TTAA (sequencing error)
check_TTAA<-function(x){
  
  start(x) <- start(x)-1
  end(x) <- end(x)  +1
  return(x)
}
#PiggyBac_Granges<- unselected_Granges
PiggyBac_Granges<- shiftRanges(unselected_Granges)
PiggyBac_Granges<-check_TTAA(PiggyBac_Granges)


#unselected_Grangesold <- unselected_Granges
#PiggyBac_Granges <- shiftRanges(PiggyBac_Granges)

sequences_PiggyBac<-(DNAStringSet(getSeq(BSgenome.Mmusculus.UCSC.mm10, PiggyBac_Granges)))
basematrix_PiggyBac<-as.matrix(sequences_PiggyBac)
ggseqlogo(consensusMatrix(basematrix_PiggyBac), plot=TRUE)

PiggyBac_Granges_correct <-PiggyBac_Granges
names(PiggyBac_Granges_correct) <- 1:length(PiggyBac_Granges_correct)
strand(PiggyBac_Granges_correct) <- PiggyBac_Granges_correct$orientation

saveRDS(PiggyBac_Granges_correct, "/s/project/transposon/Output/Snakemake/unselectedYohsibaPBIllumina/correct_ranges.RData")

sequencesallunselected<-(DNAStringSet(getSeq(BSgenome.Mmusculus.UCSC.mm10, PiggyBac_Granges_correct)))
basematrix_positives<-as.matrix(sequencesallunselected)
ggseqlogo(consensusMatrix(basematrix_positives), plot=TRUE)

