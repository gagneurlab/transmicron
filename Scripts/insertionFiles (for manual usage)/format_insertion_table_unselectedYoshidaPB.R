set.seed(123)

#load packages
library(data.table)
library(GenomicRanges)
library(Biostrings)
library(BSgenome.Mmusculus.UCSC.mm10)
library(readxl)
#input sleeping beauty vs piggy bac
insertions <- fread("/s/project/transposon/Output/Snakemake/Input/PB+G418_mm8.spacing.bed")
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
end(unselected_Granges) <- end(unselected_Granges)-1



#check whether insertions are within TTAA (sequencing error)
check_TTAA<-function(x){
  
  start(x) <- start(x)-1
  end(x) <- end(x)  +2
  return(x)
}


PiggyBac_Granges<-check_TTAA(unselected_Granges)

shiftRanges <- function(inputRanges){
  minusStrand <-which(inputRanges$orientation=="+")
  start(inputRanges)[minusStrand] <- start(inputRanges)[minusStrand]+2
  end(inputRanges)[minusStrand] <- end(inputRanges)[minusStrand]+2
  
  return(inputRanges)
}

#unselected_Grangesold <- unselected_Granges
PiggyBac_Granges <- shiftRanges(PiggyBac_Granges)
sequences_PiggyBac<-(DNAStringSet(getSeq(BSgenome.Mmusculus.UCSC.mm10, PiggyBac_Granges)))
basematrix_PiggyBac<-as.matrix(sequences_PiggyBac)
ggseqlogo(consensusMatrix(basematrix_PiggyBac), plot=TRUE)

PiggyBac_Granges_correct <-PiggyBac_Granges
names(PiggyBac_Granges_correct) <- 1:length(PiggyBac_Granges_correct)
strand(PiggyBac_Granges_correct) <- PiggyBac_Granges_correct$orientation

#this function is used to convert a granges object containing TTAA / TA sites to a granges object with single bp insertion loci
get_insertion_site_fromPattern <-function(granges_pattern, matching_pattern){
  
  
  strand<-as.vector(strand(granges_pattern))
  index_minus<-which(strand=="-")
  index_plus<-which(strand=="+")
  
  if(matching_pattern == "TTAA"){
    
    end(granges_pattern[index_plus]) <- end(granges_pattern[index_plus]) -1
    start(granges_pattern[index_minus]) <- start(granges_pattern[index_minus]) +1
    
    start(granges_pattern[index_plus]) <- start(granges_pattern[index_plus]) +1
    end(granges_pattern[index_minus]) <- end(granges_pattern[index_minus]) -1
  }
  
  else if(matching_pattern == "TA"){
    end(granges_pattern[index_plus]) <- end(granges_pattern[index_plus]) -1
    start(granges_pattern[index_minus]) <- start(granges_pattern[index_minus]) +1
  }
  
  return(granges_pattern)
}

PiggyBac_Granges_correct <- get_insertion_site_fromPattern(PiggyBac_Granges_correct, "TTAA")

PiggyBacDT <- as.data.table(PiggyBac_Granges_correct)
PiggyBacDT <- PiggyBacDT[,c("seqnames","start", "end", "orientation", "TumorID", "ReadCoverage")]
PiggyBac_Granges_correct <- PiggyBac_Granges_correct[-which(seqnames(PiggyBac_Granges_correct)=="chrX")]
PiggyBac_Granges_correct <- PiggyBac_Granges_correct[-which(seqnames(PiggyBac_Granges_correct)=="chrY")]
write.table(PiggyBacDT, "/s/project/transposon/Output/PublishedPipeline/Input/BEDInsertionTesting/PBmESC.BED")
write.table(PiggyBacDT, "/s/project/transposon/Output/PaperPipeline/Input/BEDInsertionTesting/PBmESC.BED")

sequencesallunselected<-(DNAStringSet(getSeq(BSgenome.Mmusculus.UCSC.mm10, PiggyBac_Granges_correct)))
basematrix_positives<-as.matrix(sequencesallunselected)
ggseqlogo(consensusMatrix(basematrix_positives), plot=TRUE)

