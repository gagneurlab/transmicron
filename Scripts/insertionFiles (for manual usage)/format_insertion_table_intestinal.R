set.seed(123)

#load packages
library(ggplot2)
library(data.table)
library(stringr)
library(seqLogo)
library(GenomicRanges)
library(Biostrings)
library(BSgenome.Mmusculus.UCSC.mm10.UCSC.mm10)
library(BSgenome)
library(ggseqlogo)


#input sleeping beauty vs piggy bac
intestinal<-fread("/s/project/transposon/Output/Snakemake/Input/insertion_inputs/intestinal_cancer/intestinal_cancer.txt")
intestinal
t <- as.data.frame(tstrsplit(intestinal$V4, "_"))
intestinal$TumorID <- t$c..7.115....7.115....7.115....7.115....7.115....7.115....7.115...

#subset for sb
SB_insertions<-intestinal


#specify genome
genome=BSgenome.Mmusculus.UCSC.mm9

######get sequence
#convert to Granges
insertion_table_to_Granges<- function(x){
  y<-data.table()
  y<-y[,chr:= x$`track type=bed name=Final_SB_Insertions_Sites description=Final_SB_Insertion_Sites visibility=2 color=0,0,255`]
  y<-y[,start:= as.numeric(as.character(x$V2))]
  y<-y[,end:= as.numeric(as.character(x$V3))]
  y<-y[,orientation:= as.character(x$V6)]
  y<-y[,ReadCoverage:= x$V5]
  y<-y[,TumorID:= x$TumorID]
  y<- toGRanges(y, keep.extra.columns = TRUE, genome=genome)
  strand(y)<-y$orientation
  return(y)
}

SB_Granges<-insertion_table_to_Granges(SB_insertions)
chain<-import.chain(con="/s/project/transposon/Input/mm9ToMm10.over.chain")
SB_Granges<-unlist(liftOver(SB_Granges,chain))
seqlevels(SB_Granges) <- seqlevels(BSgenome.Mmusculus.UCSC.mm10)
seqinfo(SB_Granges) <- seqinfo(BSgenome.Mmusculus.UCSC.mm10)

#out_of_bounds <- GenomicRanges:::get_out_of_bound_index(SB_Granges)
#SB_Granges <- SB_Granges[-out_of_bounds]
sequencesallSB<-(DNAStringSet(getSeq(BSgenome.Mmusculus.UCSC.mm10, SB_Granges)))
basematrix_positives<-as.matrix(sequencesallSB)
ggseqlogo(consensusMatrix(basematrix_positives), plot=TRUE)


############
# donor was already removed by the authors
############
SB_Granges_correct <- SB_Granges

plyr::count(as.vector(seqnames(SB_Granges_correct)))
#remove sex chr
SB_Granges_correct <- SB_Granges_correct[-which(seqnames(SB_Granges_correct)=="chrX")]
SB_Granges_correct <- SB_Granges_correct[-which(seqnames(SB_Granges_correct)=="chrY")]

SBDT <- as.data.table(SB_Granges_correct)
SBDT <- SBDT[,c("seqnames","start", "end", "orientation", "TumorID", "ReadCoverage")]

fwrite(SBDT,  "/s/project/transposon/Output/PublishedPipeline/Input/BEDInsertionTesting/Intestinal.BED")
write.table(SBDT,  "/s/project/transposon/Output/PaperPipeline/Input/BEDInsertionTesting/Intestinal.BED")




