library(data.table)
library(readxl)
library(BSgenome.Mmusculus.UCSC.mm9)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ggseqlogo)

#input sleeping beauty vs piggy bac
melanoma<-as.data.frame(read_xls("/s/project/transposon/Output/Snakemake/Input/melanoma.xls"))
uniqueN(melanoma$tumor)

#remove donors
t <- as.data.frame(tstrsplit(melanoma$tumor, "_"))
uniqueN(t$c..2.MBM0064....5.MBM0140....5.MBM0147....1.MBM0006....5.MBM0151...) #-> 77 melanomas https://www-nature-com.eaccess.ub.tum.de/articles/ng.3275#MOESM100
melanoma$tumorID <- t$c..2.MBM0064....5.MBM0140....5.MBM0147....1.MBM0006....5.MBM0151...

#specify genome
genome=BSgenome.Mmusculus.UCSC.mm9

######get sequence
#convert to Granges
insertion_table_to_Granges<- function(x){
  y<-data.table()
  y<-y[,chr:= x$chromosome]
  y<-y[,start:= as.numeric(as.character(x$`nucleotide start`))]
  y<-y[,end:= as.numeric(as.character(x$`nucleotide end`))]
  y<-y[,orientation:= as.character(x$strand)]
  y<-y[,ReadCoverage:= x$frequency]
  y<-y[,TumorID:= x$tumorID]
  y<- toGRanges(y, keep.extra.columns = TRUE, genome=genome)
  strand(y)<-y$orientation
  return(y)
}

SB_Granges<-insertion_table_to_Granges(melanoma)
chain<-import.chain(con="/s/project/transposon/Input/mm9ToMm10.over.chain")
SB_Granges<-unlist(liftOver(SB_Granges,chain))
seqlevels(SB_Granges) <- seqlevels(BSgenome.Mmusculus.UCSC.mm10)
seqinfo(SB_Granges) <- seqinfo(BSgenome.Mmusculus.UCSC.mm10)


sequencesallSB<-(DNAStringSet(getSeq(BSgenome.Mmusculus.UCSC.mm10, SB_Granges)))
basematrix_positives<-as.matrix(sequencesallSB)
ggseqlogo(consensusMatrix(basematrix_positives), plot=TRUE)


#remove donor 
SB_Granges <- SB_Granges[-which(seqnames(SB_Granges)%in%c("chr1","chr4"))]


SB_Granges_correct <- SB_Granges
saveRDS(SB_Granges_correct, "/s/project/transposon/Output/Snakemake/Melanoma/correct_ranges.RData")

