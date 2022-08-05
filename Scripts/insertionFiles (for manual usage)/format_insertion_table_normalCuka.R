library(data.table)
library(readxl)
library(BSgenome.Mmusculus.UCSC.mm9)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ggseqlogo)

#input sleeping beauty vs piggy bac
Cuka<-as.data.frame(read_xlsx("/s/project/transposon/Output/Snakemake/Input/cuka.xlsx"))

#specify genome
genome=BSgenome.Mmusculus.UCSC.mm9

######get sequence
#convert to Granges
insertion_table_to_Granges<- function(x){
  x <- Cuka
  y<-data.table()
  y<-y[,chr:= x$Chr]
  y<-y[,start:= as.numeric(as.character(x$start))]
  y<-y[,end:= as.numeric(as.character(x$end))]
  y<-y[,orientation:= as.character(x$strand)]
  y<-y[,ReadCoverage:= x$reads]
  y<-y[,TumorID:= x$TumorID]
  y<- toGRanges(y, keep.extra.columns = TRUE, genome=genome)
  strand(y)<-y$orientation
  return(y)
}

SB_Granges<-insertion_table_to_Granges(Cuka)
chain<-import.chain(con="/s/project/transposon/Input/mm9ToMm10.over.chain")
SB_Granges<-unlist(liftOver(SB_Granges,chain))
seqlevels(SB_Granges) <- seqlevels(BSgenome.Mmusculus.UCSC.mm10)
seqinfo(SB_Granges) <- seqinfo(BSgenome.Mmusculus.UCSC.mm10)


sequencesallSB<-(DNAStringSet(getSeq(BSgenome.Mmusculus.UCSC.mm10, SB_Granges)))
basematrix_positives<-as.matrix(sequencesallSB)
ggseqlogo(consensusMatrix(basematrix_positives), plot=TRUE)


#remove donor 
SB_Granges <- SB_Granges[-which(seqnames(SB_Granges)%in%c("chr9"))]
SB_Granges_correct <- SB_Granges
saveRDS(SB_Granges_correct, "/s/project/transposon/Output/Snakemake/Cuka/correct_ranges.RData")

