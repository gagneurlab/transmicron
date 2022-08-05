set.seed(123)

#load packages
library(data.table)
library(GenomicRanges)
library(Biostrings)
library(BSgenome.Mmusculus.UCSC.mm10)
library(readxl)
#input sleeping beauty vs piggy bac
insertions <- read_excel("/s/project/transposon/Output/Snakemake/Input/unselected.xls")
insertions$start <- insertions$pos
insertions$end <- insertions$pos
insertions$pos<-NULL
insertions$chr <- as.character(insertions$chr)
insertions[which(is.na(insertions$chr)),1] <- "X"
col_order <- c("chr", "start", "end",
               "strand", "sample", "time")
insertions <- insertions[, col_order]

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
  y<-y[,orientation:= x$strand]
  y<-y[,gene:= x$Gene]
  y<-y[,ReadCoverage:= x$ReadCoverage]
  y<-y[,normCov:= x$normCov]
  y<-y[,TumorID:= x$sample]
  y<-y[,organ:=x$organ]
  y<-y[,time:=x$time]
  y<-y[,transposon_type:=x$transposon_type]
  y<-regioneR::toGRanges(y, genome=genome)
  strand(y)<-y$orientation
  return(y)
}
unselected_Granges<-insertion_table_to_Granges(insertions)
names(unselected_Granges) <- 1:length(unselected_Granges)

chain<-import.chain(con="/s/project/transposon/Input/mm9ToMm10.over.chain")
unselected_Granges<-unlist(liftOver(unselected_Granges,chain))
seqlevels(unselected_Granges) <- seqlevels(BSgenome.Mmusculus.UCSC.mm10)
seqinfo(unselected_Granges) <- seqinfo(BSgenome.Mmusculus.UCSC.mm10)


shiftRanges <- function(inputRanges){
  minusStrand <-which(inputRanges$orientation=="-")
  
  start(inputRanges)[minusStrand] <- start(inputRanges)[minusStrand]-3
  end(inputRanges)[minusStrand] <- end(inputRanges)[minusStrand]-3

  return(inputRanges)
}
#unselected_Grangesold <- unselected_Granges
unselected_Granges <- shiftRanges(unselected_Granges)

#check whether insertions are within TTAA (sequencing error)
check_TTAA<-function(x){
  
  start(x) <- start(x) 
  end(x) <- end(x)  +3
  return(x)
}


PiggyBac_Granges<-check_TTAA(unselected_Granges)

sequences_PiggyBac<-(DNAStringSet(getSeq(BSgenome.Mmusculus.UCSC.mm10, PiggyBac_Granges)))
basematrix_PiggyBac<-as.matrix(sequences_PiggyBac)
head(basematrix_PiggyBac)
ggseqlogo(consensusMatrix(basematrix_PiggyBac), plot=TRUE)

PiggyBac_Granges_correct <-PiggyBac_Granges
names(PiggyBac_Granges_correct) <- 1:length(PiggyBac_Granges_correct)
strand(PiggyBac_Granges_correct) <- PiggyBac_Granges_correct$orientation


#remove donor chromosome

donorLookup <- matrix(ncol=2, nrow=2)
donorLookup[1,1] <- "H"
donorLookup[2,1] <- "G"
donorLookup[1,2] <- "chrX"
donorLookup[2,2] <- "chr11"
colnames(donorLookup) <- c("TumorID", "donor")
donorLookup <- as.data.frame(donorLookup)
donorDF <- as.data.frame(PiggyBac_Granges_correct)
donorDF <- left_join(donorDF, donorLookup)
onDonor <- which(donorDF$seqnames==donorDF$donor)
PiggyBac_Granges_correct <- PiggyBac_Granges_correct[-onDonor]
names(PiggyBac_Granges_correct) <- 1:length(PiggyBac_Granges_correct)

saveRDS(PiggyBac_Granges_correct, "/s/project/transposon/Output/Snakemake/unselected/correct_ranges.RData")

sequencesallunselected<-(DNAStringSet(getSeq(BSgenome.Mmusculus.UCSC.mm10, PiggyBac_Granges_correct)))
basematrix_positives<-as.matrix(sequencesallunselected)
ggseqlogo(consensusMatrix(basematrix_positives), plot=TRUE)

