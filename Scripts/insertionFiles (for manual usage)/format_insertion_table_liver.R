library(data.table)
library(readxl)
library(BSgenome.Mmusculus.UCSC.mm9)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ggseqlogo)

#input sleeping beauty vs piggy bac
Liver<-fread("/s/project/transposon/Output/Snakemake/Input/Liver.txt")

#remove donors

t <- as.data.frame(tstrsplit(Liver$V4, "_"))
uniqueN(t$c..2.bard99....2.bard99....2.bard99....2.bard99....2.bard99...)
Liver$tumorID <- t$c..2.bard99....2.bard99....2.bard99....2.bard99....2.bard99...

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
  y<-y[,TumorID:= x$tumorID]
  y<- toGRanges(y, keep.extra.columns = TRUE, genome=genome)
  strand(y)<-y$orientation
  return(y)
}

SB_Granges<-insertion_table_to_Granges(Liver)
chain<-import.chain(con="/s/project/transposon/Input/mm9ToMm10.over.chain")
SB_Granges<-unlist(liftOver(SB_Granges,chain))
seqlevels(SB_Granges) <- seqlevels(BSgenome.Mmusculus.UCSC.mm10)
seqinfo(SB_Granges) <- seqinfo(BSgenome.Mmusculus.UCSC.mm10)


sequencesallSB<-(DNAStringSet(getSeq(BSgenome.Mmusculus.UCSC.mm10, SB_Granges)))
basematrix_positives<-as.matrix(sequencesallSB)
ggseqlogo(consensusMatrix(basematrix_positives), plot=TRUE)

#remove donor 
SB_Granges <- SB_Granges[-which(seqnames(SB_Granges)=="chr4")]
SB_Granges_correct <- SB_Granges
SB_Granges_correct <- SB_Granges_correct[-which(seqnames(SB_Granges_correct)=="chrX")]
SB_Granges_correct <- SB_Granges_correct[-which(seqnames(SB_Granges_correct)=="chrY")]

SBDT <- as.data.table(SB_Granges_correct)
SBDT <- SBDT[,c("seqnames","start", "end", "orientation", "TumorID", "ReadCoverage")]

write.table(SBDT, "/s/project/transposon/Output/PublishedPipeline/Input/BEDInsertionTesting/Liver.BED")
write.table(SBDT, "/s/project/transposon/Output/PaperPipeline/Input/BEDInsertionTesting/Liver.BED")
