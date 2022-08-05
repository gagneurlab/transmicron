set.seed(123)

#load packages
library(data.table)
library(GenomicRanges)
library(Biostrings)
library(BSgenome.Mmusculus.UCSC.mm10)
library(readxl)
#input sleeping beauty vs piggy bac
insertions <- fread("/s/project/transposon/Output/Snakemake/Input/SB+G418_mm8.spacing.bed")
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
  #strand(y)<-y$orientation
  return(y)
}
unselected_Granges<-insertion_table_to_Granges(insertions)
names(unselected_Granges) <- 1:length(unselected_Granges)

chain<-import.chain(con="//s/project/transposon/Output/Snakemake/Input/mm8ToMm10.over.chain")
unselected_Granges<-unlist(liftOver(unselected_Granges,chain))
seqlevels(unselected_Granges) <- seqlevels(BSgenome.Mmusculus.UCSC.mm10)
seqinfo(unselected_Granges) <- seqinfo(BSgenome.Mmusculus.UCSC.mm10)
end(unselected_Granges) <- end(unselected_Granges)+1
start(unselected_Granges) <- start(unselected_Granges)+1

sequences_SleeepinBeauty<-(DNAStringSet(getSeq(BSgenome.Mmusculus.UCSC.mm10, unselected_Granges)))
basematrix_SleeepinBeauty<-as.matrix(sequences_SleeepinBeauty)
ggseqlogo(consensusMatrix(basematrix_SleeepinBeauty), plot=TRUE)

SleeepinBeauty_Granges_correct <-unselected_Granges
names(SleeepinBeauty_Granges_correct) <- 1:length(SleeepinBeauty_Granges_correct)
strand(SleeepinBeauty_Granges_correct) <- SleeepinBeauty_Granges_correct$orientation

SleeepinBeauty_Granges_correct <- SleeepinBeauty_Granges_correct[-which(seqnames(SleeepinBeauty_Granges_correct)=="chrX")]
SleeepinBeauty_Granges_correct <- SleeepinBeauty_Granges_correct[-which(seqnames(SleeepinBeauty_Granges_correct)=="chrY")]

SBDT <- as.data.table(SleeepinBeauty_Granges_correct)
SBDT <- SBDT[,c("seqnames","start", "end", "orientation", "TumorID", "ReadCoverage")]


write.table(SBDT, "/s/project/transposon/Output/PublishedPipeline/Input/BEDInsertionTesting/SBmESC.BED")
write.table(SBDT, "/s/project/transposon/Output/PaperPipeline/Input/BEDInsertionTesting/SBmESC.BED")
