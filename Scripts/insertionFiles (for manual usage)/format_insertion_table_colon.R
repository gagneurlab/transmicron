set.seed(123)

#load packages
library(data.table)
library(GenomicRanges)
library(Biostrings)
library(BSgenome.Mmusculus.UCSC.mm10)

#input sleeping beauty vs piggy bac
colon<-fread("/s/project/transposon/Output/Snakemake/Input/insertion_inputs/Anja_colon/anja_download/ColonPBGene_allcohorts-ext.BED")

#subset for colon samples
colon <- colon[which(grepl(colon$tissue, pattern="TU")),]

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
  #strand(y)<-y$orientation
  return(y)
}
colon_Granges<-insertion_table_to_Granges(colon)
check_TTAA<-function(x){
  
  start(x) <- start(x) - 1
  end(x) <- end(x) + 2
  return(x)
}


PiggyBac_Granges<-check_TTAA(colon_Granges)

sequences_PiggyBac<-(DNAStringSet(getSeq(BSgenome.Mmusculus.UCSC.mm10, PiggyBac_Granges)))
basematrix_PiggyBac<-as.matrix(sequences_PiggyBac)

#check which insertions are shifted by 1
PiggyBac_shifted_Granges<-PiggyBac_Granges

shift_to_TTAA<-function(x){
  original<-x
  y<-x
  end(y) <- end(y) -1
  to_shift<-which(DNAStringSet(getSeq(BSgenome.Mmusculus.UCSC.mm10,y))=="TAA")
  z<-x[to_shift]
  
  start(z)<-start(z)-1
  end(z)<-end(z)-1
  x[to_shift]<-z
  start(original) <- start(original) +1
  
  to_shift2<-which(DNAStringSet(getSeq(BSgenome.Mmusculus.UCSC.mm10,original))=="TTA")
  z<-x[to_shift2]
  
  start(z)<-start(z)+1
  end(z)<-end(z)+1
  
  x[to_shift2]<-z
  return(x)
}

PiggyBac_Granges_correct<-shift_to_TTAA(PiggyBac_shifted_Granges)
names(PiggyBac_Granges_correct) <- 1:length(PiggyBac_Granges_correct)
strand(PiggyBac_Granges_correct) <- PiggyBac_Granges_correct$orientation

saveRDS(PiggyBac_Granges_correct, "/s/project/transposon/Output/Snakemake/colonRAD/correct_ranges.RData")

sequencesallcolon<-(DNAStringSet(getSeq(BSgenome.Mmusculus.UCSC.mm10, PiggyBac_Granges_correct)))
basematrix_positives<-as.matrix(sequencesallcolon)
ggseqlogo(consensusMatrix(basematrix_positives), plot=TRUE)
