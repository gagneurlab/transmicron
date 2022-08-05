library(data.table)
library(BSgenome.Mmusculus.UCSC.mm10)

#input sleeping beauty vs piggy bac
SBvsPB<-fread("/s/project/transposon/rawdata/SBvsPB/SBvsPB.BED")
#subset for PB
PiggyBac_insertions<-SBvsPB[which(SBvsPB$transposon_type=="PB"),]

#analyze
length(which(PiggyBac_insertions$diagnosis=="TAIL"))
uniqueN(PiggyBac_insertions$TumorID)

#remove tail_samples
PiggyBac_insertions<-PiggyBac_insertions[-which(PiggyBac_insertions$diagnosis=="TAIL"),]

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
  y<- makeGRangesFromDataFrame(y, keep.extra.columns = TRUE)
  #strand(y)<-y$orientation
  return(y)
}

PiggyBac_Granges<-insertion_table_to_Granges(PiggyBac_insertions)

#check whether insertions are within TTAA (sequencing error)
check_TTAA<-function(x){
  
  start(x) <- start(x) - 1
  end(x) <- end(x) + 2
  return(x)
}


PiggyBac_Granges<-check_TTAA(PiggyBac_Granges)

sequences_PiggyBac<-(DNAStringSet(getSeq(BSgenome.Mmusculus.UCSC.mm10,PiggyBac_Granges)))
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
  
  to_shift2<-which(DNAStringSet(getSeq(BSgenome.Mmusculus.UCSC.mm10, original))=="TTA")
  z<-x[to_shift2]
  
  start(z)<-start(z)+1
  end(z)<-end(z)+1
  
  x[to_shift2]<-z
  return(x)
}

PiggyBac_Granges_correct<-shift_to_TTAA(PiggyBac_shifted_Granges)
strand(PiggyBac_Granges_correct) <- PiggyBac_Granges_correct$orientation

names(PiggyBac_Granges_correct) <-  1:length(PiggyBac_Granges_correct)

#this function is used to convert a granges object containing TTAA / TA sites to a granges object with single bp insertion loci (width = 2)
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

#remove donor 
PiggyBac_Granges_correct <- PiggyBac_Granges_correct[-which(seqnames(PiggyBac_Granges_correct)=="chr14")]
PiggyBac_Granges_correct <- PiggyBac_Granges_correct[-which(seqnames(PiggyBac_Granges_correct)=="chrX")]
PiggyBac_Granges_correct <- PiggyBac_Granges_correct[-which(seqnames(PiggyBac_Granges_correct)=="chrY")]

PiggyBacDT <- as.data.table(PiggyBac_Granges_correct)
PiggyBacDT <- PiggyBacDT[,c("seqnames","start", "end", "orientation", "TumorID", "ReadCoverage")]

write.table(PiggyBacDT, "/s/project/transposon/Output/PublishedPipeline/Input/BEDInsertionTesting/DLBCLPB.BED")
write.table(PiggyBacDT, "/s/project/transposon/Output/PaperPipeline/Input/BEDInsertionTesting/DLBCLPB.BED")

#precentage remaining
nrow(PiggyBacDT)/nrow(PiggyBac_insertions)
