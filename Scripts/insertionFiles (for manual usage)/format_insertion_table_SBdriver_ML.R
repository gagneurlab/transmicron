library(data.table)
library(GenomicRanges)
library(regioneR)
library(Biostrings)
library(BSgenome.Mmusculus.UCSC.mm9)
Mmusculus <- BSgenome.Mmusculus.UCSC.mm9
#read input
ML_input <- read.table("/data/nasif12/home_if12/bredthau/workspace/transposon/Benchmarks/SBdriver/sbdriver-kit-20180531/results/ML/ML.BENCHMARK.bed")

#convert to Granges
head(ML_input)
insertion_table_to_Granges<- function(x){
  y<-data.table()
  y<-y[,chr:= x$V1]
  y<-y[,start:= as.numeric(as.character(x$V2))]
  y<-y[,end:= as.numeric(as.character(x$V3))]
  y<-y[,orientation:= x$V6]
  y<-y[,ReadCoverage:= x$V5]
  y<-y[,TumorID:= x$V4]
  y<-toGRanges(y, genome=Mmusculus)
  strand(y)<-y$orientation
  return(y)
}

SB_Granges<-insertion_table_to_Granges(ML_input)
names(SB_Granges) <- 1:length(SB_Granges)




chain<-import.chain(con="/s/project/transposon/Input/mm9ToMm10.over.chain")
SB_Granges<-unlist(liftOver(SB_Granges,chain))
seqlevels(SB_Granges) <- seqlevels(BSgenome.Mmusculus.UCSC.mm10)
seqinfo(SB_Granges) <- seqinfo(BSgenome.Mmusculus.UCSC.mm10)


sequencesallSB<-(DNAStringSet(getSeq(BSgenome.Mmusculus.UCSC.mm10, SB_Granges)))
basematrix_positives<-as.matrix(sequencesallSB)
ggseqlogo(consensusMatrix(basematrix_positives), plot=TRUE)

SB_Granges_correct <- SB_Granges

donor_chromosome ="chr18"
on_donor <- which(as.vector(seqnames(SB_Granges_correct))==donor_chromosome)
if(length(on_donor)!=0){
  SB_Granges_correct <- SB_Granges_correct[-on_donor]
  print("excluded insertions:")
  print(length(on_donor))}

saveRDS(SB_Granges_correct, "/s/project/transposon/Output/Snakemake/SBdriver_ML/correct_ranges.RData")

