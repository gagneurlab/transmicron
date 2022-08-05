set.seed(123)

#load packages
library(data.table)
library(GenomicRanges)
library(Biostrings)
library(BSgenome.Mmusculus.UCSC.mm10)

#input sleeping beauty vs piggy bac
SBvsPB<-fread("/s/project/transposon/rawdata/SBvsPB/SBvsPB.BED")
#subset for sb
SB_insertions<-SBvsPB[which(SBvsPB$transposon_type=="SB"),]
SB_insertions_edited<-SB_insertions
#analyze
length(which(SB_insertions$diagnosis=="TAIL"))
uniqueN(SB_insertions$TumorID)

#tail_samples
SB_insertions<-SB_insertions[-which(SB_insertions$diagnosis=="TAIL"),]

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
  strand(y)<-y$orientation
  return(y)
}
SB_Granges<-insertion_table_to_Granges(SB_insertions)

#check whether insertions are within TA (sequencing error)
check_TA<-function(x){
  start(x) <- start(x)
  end(x) <- end(x) + 1
  return(x)
}

SB_Granges<-check_TA(SB_Granges)

sequencesallSB<-(DNAStringSet(getSeq(BSgenome.Mmusculus.UCSC.mm10, SB_Granges)))
basematrix_positives<-as.matrix(sequencesallSB)
ggseqlogo(consensusMatrix(basematrix_positives), plot=TRUE)

#check which insertions are shifted by 1
SB_shifted_Granges<-SB_Granges

#shift the wrong insertions
shift_to_TA<-function(x){
  original<-x
  y<-original
  end(y) <- end(y) -1
  start(y) <- start(y) -1
  to_shift<-which(DNAStringSet(getSeq(BSgenome.Mmusculus.UCSC.mm10, y))=="TA")
  start(x[to_shift])<-start(x[to_shift])-1
  end(x[to_shift])<-end(x[to_shift])-1
  
  return(x)
}


SB_Granges_correct<-shift_to_TA(SB_shifted_Granges)
names(SB_Granges_correct) <- 1:length(SB_Granges_correct)

sequencesallSB<-(DNAStringSet(getSeq(BSgenome.Mmusculus.UCSC.mm10, SB_Granges_correct)))
basematrix_positives<-as.matrix(sequencesallSB)
ggseqlogo(consensusMatrix(basematrix_positives), plot=TRUE)

#remove donor 
SB_Granges_correct <- SB_Granges_correct[-which(seqnames(SB_Granges_correct)=="chr14")]
SB_Granges_correct <- SB_Granges_correct[-which(seqnames(SB_Granges_correct)=="chrX")]
SB_Granges_correct <- SB_Granges_correct[-which(seqnames(SB_Granges_correct)=="chrY")]


SBDT <- as.data.table(SB_Granges_correct)
SBDT <- SBDT[,c("seqnames","start", "end", "orientation", "TumorID", "ReadCoverage")]


write.table(SBDT, "/s/project/transposon/Output/PublishedPipeline/Input/BEDInsertionTesting/DLBCLSB.BED")
write.table(SBDT, "/s/project/transposon/Output/PaperPipeline/Input/BEDInsertionTesting/DLBCLSB.BED")
